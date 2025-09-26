!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module antitracer_mod

!BOP
! !MODULE: antitracer_mod
!
! Module for antitracer
!
! The units of concentration for these tracers are 1/cm^3.
!
! The units of surface fluxes for these tracers are
!        1/cm^3 * cm/s == 1/cm^2/s.
!
! !DESCRIPTION:
! This module implements a generic "deficit" tracer or set of tracers
! for ocean carbon interventions. It allows for multiple independent
! tracers, each with its own forcing data and air-sea gas exchange properties.
!
! !REVISION HISTORY:
! SVN:$Id: $
! 2024-07-30: Updated for multi-tracer implementation with dynamic allocation
!             and refined shr_strdata integration.

! !USES:

    use POP_KindsMod
    use POP_ErrorMod

    use kinds_mod
    use blocks,      only: nx_block, ny_block, block
    use domain_size, only: max_blocks_clinic, km
    use domain,      only: nblocks_clinic, distrb_clinic
    use exit_mod,    only: sigAbort, exit_POP
    use communicate, only: my_task, master_task
    use constants,   only: c0, c1
    use io_types,    only: stdout
    use io_tools,    only: document
    use tavg,        only: define_tavg_field, accumulate_tavg_field
    use time_management, only : iyear, iday_of_year, frac_day, days_in_year, thour00

    use passive_tracer_tools, only: forcing_monthly_every_ts, ind_name_pair
    use passive_tracer_tools, only : read_field, tracer_read
    use broadcast
    use strdata_interface_mod, only : strdata_input_type, POP_strdata_type_set, &
                                      POP_strdata_type_match, POP_strdata_type_append_field, &
                                      POP_strdata_type_field_count, POP_strdata_create, POP_strdata_advance, &
                                      POP_strdata_type_cp
    use POP_HaloMod, only : POP_HaloUpdate
    use POP_GridHorzMod, only : POP_gridHorzLocCenter
    use POP_FieldMod, only : POP_fieldKindScalar
    use POP_CommMod, only : POP_communicator
    use POP_ErrorMod, only : POP_Success
    use domain, only : POP_haloClinic
    use blocks, only : get_block

    implicit none
    save

!-----------------------------------------------------------------------
! public/private declarations
!-----------------------------------------------------------------------

    private

! !PUBLIC MEMBER FUNCTIONS:

    public :: &
        antitracer_tracer_cnt, &
        antitracer_init, &
        antitracer_set_sflux,  &
        antitracer_tavg_forcing, &
        antitracer_column_integral_tavg, &
        tavg_ANTITRACER_COLUMN_INTEGRAL



!EOP
!BOC

    integer (int_kind), dimension(:), allocatable :: tavg_ANTITRACER_FORCING
    integer (int_kind), dimension(:), allocatable :: tavg_ANTITRACER_COLUMN_INTEGRAL

!-----------------------------------------------------------------------
! module variables required by passive_tracers
!-----------------------------------------------------------------------

    integer(int_kind), parameter :: antitracer_tracer_cnt = ANTITRACER_TRACER_CNT

!-----------------------------------------------------------------------
! relative tracer indices
!-----------------------------------------------------------------------
    integer (int_kind), parameter :: antitracer_global_offset_dummy = 1

!-----------------------------------------------------------------------
! derived type & parameter for tracer index lookup
!-----------------------------------------------------------------------

    type(ind_name_pair), dimension(:), allocatable :: ind_name_table

    ! New derived type to store forcing information for each antitracer
    type antitracer_forcing_info_type
        character(char_len) :: name
        integer(int_kind)   :: tracer_local_idx
        character(char_len) :: filename
        character(char_len) :: file_varname
        integer(int_kind)   :: year_first
        integer(int_kind)   :: year_last
        integer(int_kind)   :: year_align
        real(r8)            :: scale_factor
        integer(int_kind)   :: surface_strdata_inputlist_ind
        integer(int_kind)   :: surface_strdata_var_ind
    end type antitracer_forcing_info_type

    type(antitracer_forcing_info_type), dimension(:), allocatable :: all_antitracer_forcing_info

!-----------------------------------------------------------------------
! mask that eases avoidance of computation over land
!-----------------------------------------------------------------------

    logical (log_kind), dimension(:,:,:), allocatable :: LAND_MASK

!-----------------------------------------------------------------------
! Forcing data streams
!-----------------------------------------------------------------------

    type (strdata_input_type), pointer :: surface_strdata_inputlist_ptr(:)
    logical(log_kind), save :: antitracer_io_initialized = .false.

!-----------------------------------------------------------------------
! define tavg id for 2d fields related to surface fluxes
!-----------------------------------------------------------------------

    real (r8), dimension(:,:,:,:), allocatable :: ANTITRACER_SFLUX_TAVG

    integer (int_kind) ::   &
      tavg_ANTITRACER_IFRAC,      &
      tavg_ANTITRACER_XKW,        &
      tavg_ANTITRACER_SCHMIDT,    &
      tavg_ANTITRACER_PV

!-----------------------------------------------------------------------
! timers
!-----------------------------------------------------------------------

    integer (int_kind) :: antitracer_sflux_timer

    real(r8), parameter :: BETA = 1.0_r8

!EOC
!***********************************************************************

contains

!***********************************************************************
!BOP
! !IROUTINE: antitracer_init
! !INTERFACE:

  subroutine antitracer_init(antitracer_ind_begin, init_ts_file_fmt, read_restart_filename, &
                             tracer_d_module, TRACER_MODULE, errorCode)

! !DESCRIPTION:
! Initialize antitracer module.

! !USES:

    use constants,  only: char_blank, delim_fmt
    use prognostic, only: curtime, oldtime, tracer_field
    use grid,       only: KMT, n_topo_smooth, fill_points
    use io_types,   only: nml_in, nml_filename
    use timers,     only: get_timer
    use time_management, only: int_to_char
    use passive_tracer_tools, only: rest_read_tracer_block, file_read_tracer_block

! !INPUT PARAMETERS:

    integer (int_kind), intent(in) :: antitracer_ind_begin
    character (*), intent(in)      :: init_ts_file_fmt, read_restart_filename

! !INPUT/OUTPUT PARAMETERS:

    type (tracer_field), dimension(:), intent(inout) :: tracer_d_module
    real (r8), dimension(nx_block,ny_block,km,antitracer_tracer_cnt,3,max_blocks_clinic), &
      intent(inout) :: TRACER_MODULE

! !OUTPUT PARAMETERS:

    integer (POP_i4), intent(out) :: errorCode

!EOP
!BOC
!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------
    character(*), parameter :: subname = 'antitracer_mod:antitracer_init'

    character(char_len) :: init_antitracer_option, init_antitracer_init_file
    character(char_len) :: init_antitracer_init_file_fmt, n_char
    character(char_len) :: antitracer_restart_filename
    integer (int_kind)  :: n, k, iblock, nml_error
    type(tracer_read), dimension(:), allocatable :: tracer_init_ext

    type antitracer_forcing_nml_type
        character(char_len) :: name, file, varname
        integer(int_kind)   :: year_first, year_last, year_align
        real(r8)            :: scale_factor
    end type antitracer_forcing_nml_type

    type(antitracer_forcing_nml_type), dimension(:), allocatable :: antitracer_forcing_nml_array

    namelist /antitracer_nml/ &
      init_antitracer_option, init_antitracer_init_file, init_antitracer_init_file_fmt, &
      tracer_init_ext, antitracer_forcing_nml_array

!-----------------------------------------------------------------------
! default namelist settings
!-----------------------------------------------------------------------
    init_antitracer_option      = 'unknown'
    init_antitracer_init_file     = 'unknown'
    init_antitracer_init_file_fmt = 'bin'

    allocate(tracer_init_ext(antitracer_tracer_cnt))
    allocate(antitracer_forcing_nml_array(antitracer_tracer_cnt))
    allocate(ind_name_table(antitracer_tracer_cnt))
    allocate(all_antitracer_forcing_info(antitracer_tracer_cnt))

    do n = 1, antitracer_tracer_cnt
      tracer_init_ext(n)%mod_varname  = 'unknown'
      tracer_init_ext(n)%filename     = 'unknown'
      tracer_init_ext(n)%file_varname = 'unknown'
      tracer_init_ext(n)%scale_factor = c1
      tracer_init_ext(n)%default_val  = c0
      tracer_init_ext(n)%file_fmt     = 'bin'

      call int_to_char(3, n, n_char)
      antitracer_forcing_nml_array(n)%name         = 'ANTITRACER' // n_char
      antitracer_forcing_nml_array(n)%file         = 'unknown'
      antitracer_forcing_nml_array(n)%varname      = 'antitracer_forcing' // n_char
      antitracer_forcing_nml_array(n)%year_first   = 1999
      antitracer_forcing_nml_array(n)%year_last    = 2019
      antitracer_forcing_nml_array(n)%year_align   = 347
      antitracer_forcing_nml_array(n)%scale_factor = 1.0e4_r8
    end do

    if (my_task == master_task) then
       open (nml_in, file=nml_filename, status='old', iostat=nml_error)
       if (nml_error /= 0) then
         nml_error = -1
       else
         nml_error =  1
       endif
       !*** keep reading until find right namelist
       do while (nml_error > 0)
         read(nml_in, nml=antitracer_nml,iostat=nml_error)
       end do
       if (nml_error == 0) close(nml_in)
    end if

    call broadcast_scalar(nml_error, master_task)
    if (nml_error /= 0) then
       call exit_POP(sigAbort,'ERROR reading antitracer namelist')
    endif

    call broadcast_scalar(init_antitracer_option, master_task)
    call broadcast_scalar(init_antitracer_init_file, master_task)
    call broadcast_scalar(init_antitracer_init_file_fmt, master_task)

    do n = 1, antitracer_tracer_cnt
      call broadcast_scalar(tracer_init_ext(n)%mod_varname, master_task)
      call broadcast_scalar(tracer_init_ext(n)%filename, master_task)
      call broadcast_scalar(tracer_init_ext(n)%file_varname, master_task)
      call broadcast_scalar(tracer_init_ext(n)%scale_factor, master_task)
      call broadcast_scalar(tracer_init_ext(n)%default_val, master_task)
      call broadcast_scalar(tracer_init_ext(n)%file_fmt, master_task)
      call broadcast_scalar(antitracer_forcing_nml_array(n)%name, master_task)
      call broadcast_scalar(antitracer_forcing_nml_array(n)%file, master_task)
      call broadcast_scalar(antitracer_forcing_nml_array(n)%varname, master_task)
      call broadcast_scalar(antitracer_forcing_nml_array(n)%year_first, master_task)
      call broadcast_scalar(antitracer_forcing_nml_array(n)%year_last, master_task)
      call broadcast_scalar(antitracer_forcing_nml_array(n)%year_align, master_task)
      call broadcast_scalar(antitracer_forcing_nml_array(n)%scale_factor, master_task)
    end do

    do n = 1, antitracer_tracer_cnt
      all_antitracer_forcing_info(n)%name           = antitracer_forcing_nml_array(n)%name
      all_antitracer_forcing_info(n)%tracer_local_idx = n
      all_antitracer_forcing_info(n)%filename       = antitracer_forcing_nml_array(n)%file
      all_antitracer_forcing_info(n)%file_varname   = antitracer_forcing_nml_array(n)%varname
      all_antitracer_forcing_info(n)%year_first     = antitracer_forcing_nml_array(n)%year_first
      all_antitracer_forcing_info(n)%year_last      = antitracer_forcing_nml_array(n)%year_last
      all_antitracer_forcing_info(n)%year_align     = antitracer_forcing_nml_array(n)%year_align
      all_antitracer_forcing_info(n)%scale_factor   = antitracer_forcing_nml_array(n)%scale_factor
      ind_name_table(n) = ind_name_pair(n, antitracer_forcing_nml_array(n)%name)
    end do

    if (size(tracer_d_module) < antitracer_tracer_cnt) then
      call exit_POP(sigAbort, 'TRACER_MODULE allocation error in ' // subname)
    endif

    do n = 1, antitracer_tracer_cnt
      tracer_d_module(n)%short_name = ind_name_table(n)%name
      tracer_d_module(n)%long_name  = ind_name_table(n)%name
      tracer_d_module(n)%units      = '1/cm^3'
      tracer_d_module(n)%tend_units = '1/cm^3/s'
      tracer_d_module(n)%flux_units = '1/cm^2/s'
    end do

    select case (trim(init_antitracer_option))
    case ('zero')
      TRACER_MODULE(:,:,:,1:antitracer_tracer_cnt,curtime,:) = c0
      TRACER_MODULE(:,:,:,1:antitracer_tracer_cnt,oldtime,:) = c0

    case ('restart')
      if (trim(init_antitracer_init_file) == 'same_as_TS') then
        antitracer_restart_filename = read_restart_filename
        init_antitracer_init_file_fmt = init_ts_file_fmt
      else
        antitracer_restart_filename = trim(init_antitracer_init_file)
      endif
      call rest_read_tracer_block(antitracer_ind_begin, &
            init_antitracer_init_file_fmt, antitracer_restart_filename, &
            tracer_d_module(1:antitracer_tracer_cnt), &
            TRACER_MODULE(:,:,:,1:antitracer_tracer_cnt,:,:))

    case ('file')
      call file_read_tracer_block(init_antitracer_init_file_fmt, &
            init_antitracer_init_file, tracer_d_module(1:antitracer_tracer_cnt), &
            ind_name_table(1:antitracer_tracer_cnt), tracer_init_ext(1:antitracer_tracer_cnt), &
            TRACER_MODULE(:,:,:,1:antitracer_tracer_cnt,:,:))
      if (n_topo_smooth > 0) then
        do n = 1, antitracer_tracer_cnt
          do k = 1, km
            call fill_points(k,TRACER_MODULE(:,:,k,n,curtime,:), errorCode)
            if (errorCode /= POP_Success) return
          end do
        end do
      endif

    case default
      call exit_POP(sigAbort, 'unknown init_antitracer_option')
    end select

    do iblock = 1, nblocks_clinic
      do n = 1, antitracer_tracer_cnt
        do k = 1, km
          where (k > KMT(:,:,iblock))
            TRACER_MODULE(:,:,k,n,curtime,iblock) = c0
            TRACER_MODULE(:,:,k,n,oldtime,iblock) = c0
          end where
        end do
      end do
    end do

    allocate( LAND_MASK(nx_block,ny_block,max_blocks_clinic) )
    LAND_MASK = (KMT .gt. 0)

    call get_timer(antitracer_sflux_timer, 'ANTITRACER_SFLUX', 1, distrb_clinic%nprocs)

    call antitracer_init_tavg
    call antitracer_init_sflux

!EOC
  end subroutine antitracer_init

!***********************************************************************
!BOP
! !IROUTINE: antitracer_init_tavg
! !INTERFACE:

  subroutine antitracer_init_tavg

! !DESCRIPTION:
! Define tavg fields not automatically handled by the base model.
!EOP
!BOC
!-----------------------------------------------------------------------
    integer (int_kind) :: var_cnt, n
    character(char_len) :: sname, lname, units, coordinates
!-----------------------------------------------------------------------

    var_cnt = 0
    call define_tavg_field(tavg_ANTITRACER_IFRAC,'ANTITRACER_IFRAC',2, &
          long_name='Ice Fraction for ANTITRACER fluxes', &
          units='fraction', grid_loc='2110', coordinates='TLONG TLAT time')
    var_cnt = var_cnt+1

    call define_tavg_field(tavg_ANTITRACER_XKW,'ANTITRACER_XKW',2, &
          long_name='XKW for ANTITRACER fluxes', &
          units='cm/s', grid_loc='2110', coordinates='TLONG TLAT time')
    var_cnt = var_cnt+1

    call define_tavg_field(tavg_ANTITRACER_SCHMIDT,'ANTITRACER_SCHMIDT',2, &
          long_name='ANTITRACER Schmidt Number', &
          units='none', grid_loc='2110', coordinates='TLONG TLAT time')
    var_cnt = var_cnt+1

    call define_tavg_field(tavg_ANTITRACER_PV,'ANTITRACER_PV',2, &
          long_name='ANTITRACER piston velocity', &
          units='cm/s', grid_loc='2110', coordinates='TLONG TLAT time')
    var_cnt = var_cnt+1

    allocate(ANTITRACER_SFLUX_TAVG(nx_block,ny_block,var_cnt,max_blocks_clinic))
    ANTITRACER_SFLUX_TAVG = c0

    allocate(tavg_ANTITRACER_FORCING(antitracer_tracer_cnt))
    do n = 1, antitracer_tracer_cnt
        sname = trim(ind_name_table(n)%name) // '_FORCING'
        lname = 'Forcing for ' // trim(ind_name_table(n)%name)
        units = '1/cm^2/s'
        coordinates = 'TLONG TLAT time'
        call define_tavg_field(tavg_ANTITRACER_FORCING(n), sname, 2, &
              long_name=lname, units=units, grid_loc='2110', coordinates=coordinates)
    end do

    ! === Define TAVG fields for column-integrated tracers ===
    allocate(tavg_ANTITRACER_COLUMN_INTEGRAL(antitracer_tracer_cnt))

    do n = 1, antitracer_tracer_cnt
        sname = trim(ind_name_table(n)%name) // '_COL_INT'
        lname = 'Column Integrated ' // trim(ind_name_table(n)%name)
        units = '1/cm^2'  ! units are (1/cm^3) * cm
        coordinates = 'TLONG TLAT time'

        call define_tavg_field(tavg_ANTITRACER_COLUMN_INTEGRAL(n), sname, 2, &
                               long_name=lname, units=units, grid_loc='2110', &
                               coordinates=coordinates)
    end do
!EOC
  end subroutine antitracer_init_tavg

!***********************************************************************
!BOP
! !IROUTINE: antitracer_init_sflux
! !INTERFACE:

  subroutine antitracer_init_sflux()

! !DESCRIPTION:
! Initialize surface flux computations.
!EOP
!BOC
!-----------------------------------------------------------------------
    integer (int_kind) :: n_tracer, m, n_strdata_entries
    type (strdata_input_type) :: surface_strdata_input_var
    type (strdata_input_type), pointer :: surface_strdata_inputlist_tmp_ptr(:)
!-----------------------------------------------------------------------

    if (.not. associated(surface_strdata_inputlist_ptr)) then
      allocate(surface_strdata_inputlist_ptr(0))
    end if

    do n_tracer = 1, antitracer_tracer_cnt
      associate(forcing_info => all_antitracer_forcing_info(n_tracer))
        call POP_strdata_type_set(surface_strdata_input_var, &
          file_name   = forcing_info%filename, &
          field       = forcing_info%file_varname, &
          timer_label = 'antitracer_file_' // trim(forcing_info%name), &
          year_first  = forcing_info%year_first, &
          year_last   = forcing_info%year_last, &
          year_align  = forcing_info%year_align, &
          depth_flag  = .false., &
          tintalgo    = 'linear', &
          taxMode     = 'cycle')

        n_strdata_entries = size(surface_strdata_inputlist_ptr)
        forcing_info%surface_strdata_inputlist_ind = 0

        do m = 1, n_strdata_entries
          if (POP_strdata_type_match(surface_strdata_input_var, surface_strdata_inputlist_ptr(m))) then
            call POP_strdata_type_append_field(forcing_info%file_varname, surface_strdata_inputlist_ptr(m))
            forcing_info%surface_strdata_inputlist_ind = m
            exit
          end if
        end do

        if (forcing_info%surface_strdata_inputlist_ind == 0) then
          n_strdata_entries = n_strdata_entries + 1
          allocate(surface_strdata_inputlist_tmp_ptr(n_strdata_entries))
          if (associated(surface_strdata_inputlist_ptr)) then
             if (size(surface_strdata_inputlist_ptr) > 0) then
                do m = 1, n_strdata_entries - 1
                   call POP_strdata_type_cp(surface_strdata_inputlist_ptr(m), &
                                            surface_strdata_inputlist_tmp_ptr(m))
                end do
             endif
             deallocate(surface_strdata_inputlist_ptr)
          endif
          surface_strdata_inputlist_ptr => surface_strdata_inputlist_tmp_ptr
          call POP_strdata_type_cp(surface_strdata_input_var, surface_strdata_inputlist_ptr(n_strdata_entries))
          forcing_info%surface_strdata_inputlist_ind = n_strdata_entries
        endif

        forcing_info%surface_strdata_var_ind = POP_strdata_type_field_count( &
          surface_strdata_inputlist_ptr(forcing_info%surface_strdata_inputlist_ind))
      end associate
    end do

!EOC
  end subroutine antitracer_init_sflux

!***********************************************************************
!BOP
! !IROUTINE: antitracer_set_sflux
! !INTERFACE:

  subroutine antitracer_set_sflux(U10_SQR,IFRAC,SST, &
                                  SURF_VALS,STF_MODULE)

! !DESCRIPTION:
! Compute ANTITRACER surface flux.

! !USES:
    use constants, only: xkw_coeff
    use timers, only: timer_start, timer_stop
    use domain,    only: blocks_clinic


! !INPUT PARAMETERS:
    ! FIX: Changed from assumed-size (:) to explicit-shape (max_blocks_clinic)
    real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: U10_SQR
    real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: IFRAC
    real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: SST
    real (r8), dimension(nx_block,ny_block,antitracer_tracer_cnt,max_blocks_clinic), intent(in) :: SURF_VALS

! !OUTPUT PARAMETERS:
    ! FIX: Changed from assumed-size (:) to explicit-shape (max_blocks_clinic)
    real (r8), dimension(nx_block,ny_block,antitracer_tracer_cnt,max_blocks_clinic), intent(inout) :: STF_MODULE

!EOP
!BOC
!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------
    character(*), parameter :: subname = 'antitracer_mod:antitracer_set_sflux'
    integer (int_kind)       :: iblock, n_tracer, m, i, j, n_idx
    integer(POP_i4)          :: errorCode
    type(block)              :: this_block
    logical(log_kind), save   :: first_call_this_timestep = .true.

    real (r8), dimension(nx_block,ny_block) :: &
      IFRAC_USED, XKW_USED, ANTITRACER_SCHMIDT, XKW_ICE, PV

    real(r8), dimension(nx_block,ny_block,max_blocks_clinic) :: tracer_forcing_data
    
    ! Storage for the pre-computed, tracer-independent piston velocity
    real(r8), dimension(nx_block,ny_block,max_blocks_clinic) :: PV_field

!-----------------------------------------------------------------------

    call timer_start(antitracer_sflux_timer)

    ! Create shr_strdata objects on the first call
    if (.not. antitracer_io_initialized) then
      do m = 1, size(surface_strdata_inputlist_ptr)
          call POP_strdata_create(surface_strdata_inputlist_ptr(m))
      end do
      antitracer_io_initialized = .true.
    endif

    ! Advance all unique shr_strdata streams once per timestep
    if (first_call_this_timestep) then
        do m = 1, size(surface_strdata_inputlist_ptr)
            call POP_strdata_advance(surface_strdata_inputlist_ptr(m))
        end do
        first_call_this_timestep = .false.
    end if

    !=======================================================================
    ! STEP 1: Pre-compute tracer-independent fields for all blocks first.
    !=======================================================================
    do iblock = 1, nblocks_clinic
        ! Pre-compute common gas exchange parameters for this block
        where (LAND_MASK(:,:,iblock))
            IFRAC_USED = IFRAC(:,:,iblock)
            XKW_USED   = xkw_coeff * U10_SQR(:,:,iblock)
        elsewhere
            IFRAC_USED = c0
            XKW_USED   = c0
        endwhere
        ! Clamp ice fraction values
        where (LAND_MASK(:,:,iblock) .and. IFRAC_USED < c0) IFRAC_USED = c0
        where (LAND_MASK(:,:,iblock) .and. IFRAC_USED > c1) IFRAC_USED = c1
    
        ! Compute Schmidt number from Sea Surface Temperature (SST)
        call comp_antitracer_schmidt(LAND_MASK(:,:,iblock), SST(:,:,iblock), ANTITRACER_SCHMIDT)
    
        ! Compute final piston velocity (PV)
        where (LAND_MASK(:,:,iblock))
            XKW_ICE = (c1 - IFRAC_USED) * XKW_USED
            PV      = XKW_ICE * sqrt(660.0_r8 / ANTITRACER_SCHMIDT)
        elsewhere
            XKW_ICE = c0
            PV      = c0
        endwhere
    
        ! Store the computed piston velocity for this block for later use
        PV_field(:,:,iblock) = PV(:,:)
    
        ! Optional: Store intermediate fields for time averaging/diagnostics
        ANTITRACER_SFLUX_TAVG(:,:,1,iblock) = IFRAC_USED(:,:)
        ANTITRACER_SFLUX_TAVG(:,:,2,iblock) = XKW_USED(:,:)
        ANTITRACER_SFLUX_TAVG(:,:,3,iblock) = ANTITRACER_SCHMIDT(:,:)
        ANTITRACER_SFLUX_TAVG(:,:,4,iblock) = PV(:,:)
    end do

    !=======================================================================
    ! STEP 2: Loop over each tracer individually, using the pre-computed fields.
    !=======================================================================
    do n_tracer = 1, antitracer_tracer_cnt
        associate(forcing_info => all_antitracer_forcing_info(n_tracer))
    
            ! a) Get forcing data for THIS tracer into the reusable 3D array
            do iblock = 1, nblocks_clinic
                this_block = get_block(blocks_clinic(iblock), iblock)
                n_idx = 0
                do j = this_block%jb, this_block%je
                    do i = this_block%ib, this_block%ie
                        n_idx = n_idx + 1
                        tracer_forcing_data(i,j,iblock) = &
                            surface_strdata_inputlist_ptr(forcing_info%surface_strdata_inputlist_ind)%sdat%avs(forcing_info%surface_strdata_var_ind)%rAttr(forcing_info%surface_strdata_var_ind, n_idx)
                    enddo
                enddo
                ! Accumulate time average for this tracer's forcing
                call accumulate_tavg_field(tracer_forcing_data(:,:,iblock), tavg_ANTITRACER_FORCING(n_tracer), iblock, 1)
            end do
    
            ! b) Apply halo update to THIS tracer's forcing data
            call POP_HaloUpdate(tracer_forcing_data, POP_haloClinic, &
                                POP_gridHorzLocCenter, POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
            if (errorCode /= POP_Success) then
                call document(subname, 'error updating halo for antitracer forcing field')
                call exit_POP(sigAbort, 'Stopping in ' // subname)
            endif
    
            ! c) Compute the final flux using the pre-computed PV_field
            do iblock = 1, nblocks_clinic
                where (LAND_MASK(:,:,iblock))
                    STF_MODULE(:,:,n_tracer,iblock) = forcing_info%scale_factor * tracer_forcing_data(:,:,iblock) - &
                                                      (PV_field(:,:,iblock) / BETA) * SURF_VALS(:,:,n_tracer,iblock)
                elsewhere
                    STF_MODULE(:,:,n_tracer,iblock) = c0
                endwhere
            end do
    
        end associate
    end do

    call timer_stop(antitracer_sflux_timer)

!-----------------------------------------------------------------------
!EOC

  end subroutine antitracer_set_sflux
!***********************************************************************
!BOP
! !IROUTINE: comp_antitracer_schmidt
! !INTERFACE:

  subroutine comp_antitracer_schmidt(LAND_MASK, SST_IN, ANTITRACER_SCHMIDT)

! !DESCRIPTION:
! Compute Schmidt numbers of ANTITRACERs.

! !USES:
    use shr_infnan_mod, only : shr_infnan_isnan
    use exit_mod,       only : exit_POP, sigAbort
    use communicate,    only : my_task, master_task
    use io_types,       only : stdout

! !INPUT PARAMETERS:

    logical (log_kind), intent(in)  :: LAND_MASK(nx_block,ny_block)
    real (r8)         , intent(in)  :: SST_IN(nx_block,ny_block)

! !OUTPUT PARAMETERS:

    real (r8)         , intent(out) :: ANTITRACER_SCHMIDT(nx_block,ny_block)

!EOP
!BOC
!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------
    integer(int_kind) :: i, j
    real (r8)         :: SST(nx_block,ny_block)

    real (r8), parameter :: a =  2116.8_r8
    real (r8), parameter :: b = -136.25_r8
    real (r8), parameter :: c =   4.7353_r8
    real (r8), parameter :: d =  -0.092307_r8
    real (r8), parameter :: e =   0.0007555_r8

!-----------------------------------------------------------------------
    do j = 1, ny_block
      do i = 1, nx_block
        if (shr_infnan_isnan(SST_IN(i,j))) then
           if (my_task == master_task) write(stdout,*) 'NaN in SST_IN in comp_antitracer_schmidt'
           call exit_POP(sigAbort, 'NaN in comp_antitracer_schmidt')
        endif
        if (LAND_MASK(i,j)) then
            SST(i,j) = max(-2.0_r8, min(40.0_r8, SST_IN(i,j)))
            ANTITRACER_SCHMIDT(i,j) = a + SST(i,j) * (b + SST(i,j) * (c + SST(i,j) * (d + SST(i,j) * e)))
        else
            ANTITRACER_SCHMIDT(i,j) = c0
        endif
      end do
    end do

!EOC
  end subroutine comp_antitracer_schmidt

!***********************************************************************
!BOP
! !IROUTINE: antitracer_tavg_forcing
! !INTERFACE:

  subroutine antitracer_tavg_forcing

! !DESCRIPTION:
! Make accumulation calls for forcing related tavg fields.
!EOP
!BOC
!-----------------------------------------------------------------------
    integer (int_kind) :: iblock
!-----------------------------------------------------------------------
    do iblock = 1, nblocks_clinic
      call accumulate_tavg_field(ANTITRACER_SFLUX_TAVG(:,:,1,iblock),tavg_ANTITRACER_IFRAC,iblock,1)
      call accumulate_tavg_field(ANTITRACER_SFLUX_TAVG(:,:,2,iblock),tavg_ANTITRACER_XKW,iblock,1)
      call accumulate_tavg_field(ANTITRACER_SFLUX_TAVG(:,:,3,iblock),tavg_ANTITRACER_SCHMIDT,iblock,1)
      call accumulate_tavg_field(ANTITRACER_SFLUX_TAVG(:,:,4,iblock),tavg_ANTITRACER_PV,iblock,1)
    end do
!EOC
  end subroutine antitracer_tavg_forcing

!***********************************************************************
!BOP
! !IROUTINE: antitracer_column_integral_tavg
! !INTERFACE:

  subroutine antitracer_column_integral_tavg(TRACER_MODULE, DZ, KMT)

! !DESCRIPTION:
! Compute and accumulate the time-average for the column-integrated
! amount of each antitracer. The result has units of 1/cm^2.

! !USES:
    ! DZ and KMT are now passed as arguments, so they are removed from here.
    use prognostic, only: curtime

! !INPUT PARAMETERS:
    ! This declaration matches the slice being passed from passive_tracers
    real (r8), dimension(nx_block,ny_block,km,antitracer_tracer_cnt,3,max_blocks_clinic), &
      intent(in) :: TRACER_MODULE

    ! << ADDED >>: Explicit declarations for the dummy arguments
    real (r8), dimension(km), intent(in) :: DZ
    integer (int_kind), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: KMT

!EOP
!BOC
!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------
    integer (int_kind) :: iblock, n_tracer, k
    real (r8), dimension(nx_block,ny_block) :: col_int_field
!-----------------------------------------------------------------------

    do iblock = 1, nblocks_clinic
      do n_tracer = 1, antitracer_tracer_cnt

        ! Zero the temporary field for this block and tracer
        col_int_field = c0

        ! Sum (concentration * layer thickness) over the water column
        do k = 1, km
          where (k <= KMT(:,:,iblock))
            col_int_field(:,:) = col_int_field(:,:) + &
                                 TRACER_MODULE(:,:,k,n_tracer,curtime,iblock) * DZ(k)
          endwhere
        enddo

        ! Accumulate the result into the tavg field for this tracer
        call accumulate_tavg_field(col_int_field, &
                                   tavg_ANTITRACER_COLUMN_INTEGRAL(n_tracer), &
                                   iblock, 1)
      enddo ! end n_tracer loop
    enddo ! end iblock loop

!EOC
  end subroutine antitracer_column_integral_tavg

!***********************************************************************

end module antitracer_mod
