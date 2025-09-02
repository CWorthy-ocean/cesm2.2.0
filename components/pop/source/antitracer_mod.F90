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
!       1/cm^3 * cm/s == 1/cm^2/s.
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
                                       POP_strdata_type_field_count, POP_strdata_create, POP_strdata_advance
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
        antitracer_tavg_forcing

!EOP
!BOC

!-----------------------------------------------------------------------
! module variables required by passive_tracers
!-----------------------------------------------------------------------

    integer(int_kind), parameter :: antitracer_tracer_cnt = ANTITRACER_TRACER_CNT

!-----------------------------------------------------------------------
! relative tracer indices
!-----------------------------------------------------------------------
    ! This parameter indicates the first antitracer index in the global tracer array.
    ! Individual antitracer indices will be 1, 2, ..., antitracer_tracer_cnt relative to this start.
    ! Note: `antitracer_ind` would typically be `n` in loops, representing the local index.
    ! If `antitracer_ind` was used for a single fixed tracer, it's now replaced by `n_tracer`.
    integer (int_kind), parameter :: antitracer_global_offset_dummy = 1 ! A placeholder; actual offset passed to `antitracer_init`

!-----------------------------------------------------------------------
! derived type & parameter for tracer index lookup
!-----------------------------------------------------------------------

    type(ind_name_pair), dimension(:), allocatable :: ind_name_table

    ! New derived type to store forcing information for each antitracer
    type antitracer_forcing_info_type
        character(char_len) :: name            ! Name of the antitracer (e.g., 'ANTITRACER1')
        integer(int_kind)   :: tracer_local_idx ! Local index of this antitracer (1, 2, ...)
        character(char_len) :: filename        ! Forcing file name
        character(char_len) :: file_varname    ! Variable name in forcing file
        integer(int_kind)   :: year_first
        integer(int_kind)   :: year_last
        integer(int_kind)   :: year_align
        real(r8)            :: scale_factor
        integer(int_kind)   :: surface_strdata_inputlist_ind ! Index in module's surface_strdata_inputlist_ptr
        integer(int_kind)   :: surface_strdata_var_ind       ! Index for the variable within the strdata entry
    end type antitracer_forcing_info_type

    type(antitracer_forcing_info_type), dimension(:), allocatable :: all_antitracer_forcing_info

!-----------------------------------------------------------------------
! mask that eases avoidance of computation over land
!-----------------------------------------------------------------------

    logical (log_kind), dimension(:,:,:), allocatable :: LAND_MASK

!-----------------------------------------------------------------------
! Forcing data streams. This pointer array will hold unique shr_strdata objects
! representing input files/streams for various tracers.
!-----------------------------------------------------------------------

    type (strdata_input_type), pointer :: surface_strdata_inputlist_ptr(:)
    logical(log_kind), save :: antitracer_io_initialized = .false.

!-----------------------------------------------------------------------
! define tavg id for 2d fields related to surface fluxes
!-----------------------------------------------------------------------

    real (r8), dimension(:,:,:,:), allocatable :: ANTITRACER_SFLUX_TAVG

    integer (int_kind) ::    &
      tavg_ANTITRACER_IFRAC,      & ! tavg id for ice fraction
      tavg_ANTITRACER_XKW,        & ! tavg id for xkw
      tavg_ANTITRACER_SCHMIDT,    & ! tavg id for antitracer Schmidt number
      tavg_ANTITRACER_PV          ! tavg id for antitracer piston velocity

!-----------------------------------------------------------------------
! timers
!-----------------------------------------------------------------------

    integer (int_kind) :: antitracer_sflux_timer

    ! Define BETA here, assuming it's a general constant for gas exchange
    real(r8), parameter :: BETA = 1.0_r8 ! Placeholder value

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
! Initialize antitracer module. This involves setting metadata, reading
! the modules namelist and setting initial conditions.

! !REVISION HISTORY:
! same as module

! !USES:

    use constants,  only: char_blank, delim_fmt
    use prognostic, only: curtime, oldtime, tracer_field
    use grid,       only: KMT, n_topo_smooth, fill_points
    use io_types,   only: nml_in, nml_filename
    use timers,     only: get_timer
    use time_management, only: int_to_char

    use passive_tracer_tools, only: init_forcing_monthly_every_ts, &
        rest_read_tracer_block, file_read_tracer_block

    use io_read_fallback_mod, only: io_read_fallback_register_tracer

! !INPUT PARAMETERS:

    integer (int_kind), intent(in) :: &
      antitracer_ind_begin            ! starting index of antitracers in global tracer array
                                      ! passed through to rest_read_tracer_block

    character (*), intent(in) ::  &
      init_ts_file_fmt,    & ! format (bin or nc) for input file
      read_restart_filename  ! file name for restart file

! !INPUT/OUTPUT PARAMETERS:

    ! tracer_d_module and TRACER_MODULE now expect variable size based on antitracer_tracer_cnt
    type (tracer_field), dimension(:), intent(inout) :: &
      tracer_d_module    ! descriptors for each tracer

    ! TRACER_MODULE is typically allocated in the main model driver,
    ! so its last dimension is passed in as `:`, and its actual size
    ! for the tracer dimension must be >= antitracer_tracer_cnt.
    real (r8), dimension(nx_block,ny_block,km,antitracer_tracer_cnt,3,max_blocks_clinic), &
      intent(inout) :: TRACER_MODULE

! !OUTPUT PARAMETERS:

    integer (POP_i4), intent(out) :: &
      errorCode          ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------

    character(*), parameter :: subname = 'antitracer_mod:antitracer_init'

    character(char_len) :: &
      init_antitracer_option,      &   ! option for initialization of bgc
      init_antitracer_init_file,    &  ! filename for option 'file'
      init_antitracer_init_file_fmt, & ! file format for option 'file'
      n_char 

    integer (int_kind) :: &
      n,                       & ! index for looping over tracers
      k,                       & ! index for looping over depth levels
      iblock,                  & ! index for looping over blocks
      nml_error                  ! namelist i/o error flag

    type(tracer_read), dimension(:), allocatable :: & ! Now allocatable
      tracer_init_ext            ! namelist variable for initializing tracers

    ! New namelist array for antitracer forcing definitions
    type antitracer_forcing_nml_type
        character(char_len) :: name
        character(char_len) :: file
        character(char_len) :: varname
        integer(int_kind)   :: year_first
        integer(int_kind)   :: year_last
        integer(int_kind)   :: year_align
        real(r8)            :: scale_factor
    end type antitracer_forcing_nml_type

    type(antitracer_forcing_nml_type), dimension(:), allocatable :: antitracer_forcing_nml_array

    namelist /antitracer_nml/ &
      init_antitracer_option, init_antitracer_init_file, init_antitracer_init_file_fmt, &
      tracer_init_ext, &
      antitracer_forcing_nml_array ! New namelist array for individual tracer forcing

    character (char_len) :: &
      antitracer_restart_filename      ! modified file name for restart file

!-----------------------------------------------------------------------
! default namelist settings - for antitracer_tracer_cnt and forcing_nml_array
!   these are only defaults on master, will be read/broadcast
!-----------------------------------------------------------------------
    init_antitracer_option      = 'unknown'
    init_antitracer_init_file      = 'unknown'
    init_antitracer_init_file_fmt = 'bin'
    
    ! Allocate module-level arrays based on antitracer_tracer_cnt
    allocate(tracer_init_ext(antitracer_tracer_cnt))
    allocate(antitracer_forcing_nml_array(antitracer_tracer_cnt))
    allocate(ind_name_table(antitracer_tracer_cnt))
    allocate(all_antitracer_forcing_info(antitracer_tracer_cnt))

    ! Set defaults for the allocated namelist arrays
    do n = 1, antitracer_tracer_cnt
      tracer_init_ext(n)%mod_varname  = 'unknown'
      tracer_init_ext(n)%filename     = 'unknown'
      tracer_init_ext(n)%file_varname = 'unknown'
      tracer_init_ext(n)%scale_factor = c1
      tracer_init_ext(n)%default_val  = c0
      tracer_init_ext(n)%file_fmt     = 'bin'

      ! Default forcing info for each antitracer.
      ! Names are generated automatically if not provided in namelist.
      call int_to_char(3, n, n_char)
      antitracer_forcing_nml_array(n)%name         = 'ANTITRACER' // n_char
      antitracer_forcing_nml_array(n)%file         = 'unknown'
      antitracer_forcing_nml_array(n)%varname      = 'antitracer_forcing' // n_char
      antitracer_forcing_nml_array(n)%year_first   = 1999
      antitracer_forcing_nml_array(n)%year_last    = 2019
      antitracer_forcing_nml_array(n)%year_align   = 347
      antitracer_forcing_nml_array(n)%scale_factor = 1.0e4_r8 ! convert from 1/m^2/s to 1/cm^2/s

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


!-----------------------------------------------------------------------
! broadcast all namelist variables (including the newly allocated arrays)
!-----------------------------------------------------------------------

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

!-----------------------------------------------------------------------
! Initialize internal forcing info structure and ind_name_table
!-----------------------------------------------------------------------

    do n = 1, antitracer_tracer_cnt
        all_antitracer_forcing_info(n)%name = antitracer_forcing_nml_array(n)%name
        all_antitracer_forcing_info(n)%tracer_local_idx = n ! Store local index
        all_antitracer_forcing_info(n)%filename = antitracer_forcing_nml_array(n)%file
        all_antitracer_forcing_info(n)%file_varname = antitracer_forcing_nml_array(n)%varname
        all_antitracer_forcing_info(n)%year_first = antitracer_forcing_nml_array(n)%year_first
        all_antitracer_forcing_info(n)%year_last = antitracer_forcing_nml_array(n)%year_last
        all_antitracer_forcing_info(n)%year_align = antitracer_forcing_nml_array(n)%year_align
        all_antitracer_forcing_info(n)%scale_factor = antitracer_forcing_nml_array(n)%scale_factor

        ! Populate ind_name_table for tracer descriptor creation
        ind_name_table(n) = ind_name_pair(n, antitracer_forcing_nml_array(n)%name)
    end do


!-----------------------------------------------------------------------
! initialize tracer_d values (descriptors)
!-----------------------------------------------------------------------

    ! tracer_d_module is passed in from the calling routine. It should be
    ! allocated to a size sufficient to hold all antitracer_tracer_cnt descriptors.
    ! Check that the passed in array is large enough.
    if (size(tracer_d_module) < antitracer_tracer_cnt) then
        call document(subname, 'tracer_d_module size mismatch. Expected at least', antitracer_tracer_cnt)
        call document(subname, 'Actual size of tracer_d_module', size(tracer_d_module))
        call exit_POP(sigAbort, 'TRACER_MODULE allocation error in ' // subname)
    endif

    do n = 1, antitracer_tracer_cnt
      tracer_d_module(n)%short_name = ind_name_table(n)%name
      tracer_d_module(n)%long_name  = ind_name_table(n)%name
      tracer_d_module(n)%units      = '1/cm^3'
      tracer_d_module(n)%tend_units = '1/cm^3/s'
      tracer_d_module(n)%flux_units = '1/cm^2/s'
    end do

!-----------------------------------------------------------------------
!   initialize 3D tracer fields (TRACER_MODULE)
!-----------------------------------------------------------------------

    ! TRACER_MODULE is passed in from the calling routine. Its tracer dimension
    ! must be correctly sized to accommodate all antitracers.
    ! The compiler will handle the `:` dimension mapping based on the calling context.

    select case (init_antitracer_option)

    case ('zero')
      ! Initialize all antitracer fields to zero
      TRACER_MODULE(:,:,:,1:antitracer_tracer_cnt,curtime,:) = c0
      TRACER_MODULE(:,:,:,1:antitracer_tracer_cnt,oldtime,:) = c0
      if (my_task == master_task) then
           write(stdout,delim_fmt)
           write(stdout,*) ' Initial 3-d antitracers set to all zeros'
           write(stdout,delim_fmt)
      endif

    case ('restart')

      antitracer_restart_filename = char_blank

      if (init_antitracer_init_file == 'same_as_TS') then
          if (read_restart_filename == 'undefined') then
             call document(subname, 'no restart file to read ANTITRACERs from')
             call exit_POP(sigAbort, 'stopping in ' // subname)
          endif
          antitracer_restart_filename = read_restart_filename
          init_antitracer_init_file_fmt = init_ts_file_fmt

      else  ! do not read from TS restart file

          antitracer_restart_filename = trim(init_antitracer_init_file)

      endif

      ! rest_read_tracer_block should be able to handle multiple tracers based on tracer_d_module size
      call rest_read_tracer_block(antitracer_ind_begin,        &
                                  init_antitracer_init_file_fmt, &
                                  antitracer_restart_filename,    &
                                  tracer_d_module(1:antitracer_tracer_cnt), & ! Pass relevant part
                                  TRACER_MODULE(:,:,:,1:antitracer_tracer_cnt,:,:)) ! Pass relevant part

    case ('file')

      call document(subname, 'ANTITRACERs being read from separate file')

      ! file_read_tracer_block should handle multiple tracers based on tracer_d_module size
      call file_read_tracer_block(init_antitracer_init_file_fmt, &
                                  init_antitracer_init_file,      &
                                  tracer_d_module(1:antitracer_tracer_cnt), & ! Pass relevant part
                                  ind_name_table(1:antitracer_tracer_cnt),  & ! Pass relevant part
                                  tracer_init_ext(1:antitracer_tracer_cnt), & ! Pass relevant part
                                  TRACER_MODULE(:,:,:,1:antitracer_tracer_cnt,:,:)) ! Pass relevant part

      if (n_topo_smooth > 0) then
          do n = 1, antitracer_tracer_cnt
            do k = 1, km
               call fill_points(k,TRACER_MODULE(:,:,k,n,curtime,:), &
                                 errorCode)

               if (errorCode /= POP_Success) then
                  call POP_ErrorSet(errorCode, &
                     'antitracer_init: error in fill_points')
                  return
               endif
            end do
          end do
      endif

    case default
      call document(subname, 'init_antitracer_option', init_antitracer_option)
      call exit_POP(sigAbort, 'unknown init_antitracer_option')

    end select

!-----------------------------------------------------------------------
! apply land mask to tracers
!-----------------------------------------------------------------------

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

!-----------------------------------------------------------------------
! allocate and initialize LAND_MASK (true for ocean points)
!-----------------------------------------------------------------------

    allocate( LAND_MASK(nx_block,ny_block,max_blocks_clinic) )
    LAND_MASK = (KMT.gt.0)

    call get_timer(antitracer_sflux_timer, 'ANTITRACER_SFLUX', 1, distrb_clinic%nprocs)

!-----------------------------------------------------------------------
! call other initialization subroutines
!-----------------------------------------------------------------------

    call antitracer_init_tavg
    call antitracer_init_sflux

!-----------------------------------------------------------------------
!EOC

  end subroutine antitracer_init

!***********************************************************************
!BOP
! !IROUTINE: antitracer_init_tavg
! !INTERFACE:

  subroutine antitracer_init_tavg

! !DESCRIPTION:
! Define tavg fields not automatically handled by the base model.

! !REVISION HISTORY:
! same as module

!EOP
!BOC
!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------

    integer (int_kind) :: &
      var_cnt            ! how many tavg variables are defined

!-----------------------------------------------------------------------

    var_cnt = 0

    call define_tavg_field(tavg_ANTITRACER_IFRAC,'ANTITRACER_IFRAC',2,        &
                             long_name='Ice Fraction for ANTITRACER fluxes',&
                             units='fraction', grid_loc='2110',      &
                             coordinates='TLONG TLAT time')
    var_cnt = var_cnt+1

    call define_tavg_field(tavg_ANTITRACER_XKW,'ANTITRACER_XKW',2,            &
                             long_name='XKW for ANTITRACER fluxes',          &
                             units='cm/s', grid_loc='2110',          &
                             coordinates='TLONG TLAT time')
    var_cnt = var_cnt+1

    call define_tavg_field(tavg_ANTITRACER_SCHMIDT,'ANTITRACER_SCHMIDT',2,    &
                             long_name='ANTITRACER Schmidt Number',          &
                             units='none', grid_loc='2110',          &
                             coordinates='TLONG TLAT time')
    var_cnt = var_cnt+1

    call define_tavg_field(tavg_ANTITRACER_PV,'ANTITRACER_PV',2,              &
                             long_name='ANTITRACER piston velocity',          &
                             units='cm/s', grid_loc='2110',          &
                             coordinates='TLONG TLAT time')
    var_cnt = var_cnt+1

!-----------------------------------------------------------------------

    allocate(ANTITRACER_SFLUX_TAVG(nx_block,ny_block,var_cnt,max_blocks_clinic))
    ANTITRACER_SFLUX_TAVG = c0
           
    write(stdout,*) ' Done with antitracer_init_tavg'

!-----------------------------------------------------------------------
!EOC

  end subroutine antitracer_init_tavg

!***********************************************************************
!BOP
! !IROUTINE: antitracer_init_sflux
! !INTERFACE:

  subroutine antitracer_init_sflux()

! !DESCRIPTION:
! Initialize surface flux computations for all antitracer tracer modules.
! Sets up strdata_input_type for each antitracer's forcing.
! This routine is called once during initialization.

! !USES:

    use strdata_interface_mod, only : POP_strdata_type_set
    use strdata_interface_mod, only : POP_strdata_type_match
    use strdata_interface_mod, only : POP_strdata_type_append_field
    use strdata_interface_mod, only : POP_strdata_type_cp
    use strdata_interface_mod, only : POP_strdata_type_field_count
    ! Add these for the debug statements
    use communicate, only : my_task, master_task
    use io_types,    only : stdout

!EOP
!BOC
!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------

    character(*), parameter :: subname = 'antitracer_mod:antitracer_init_sflux'
    integer (int_kind) :: &
      n_tracer, m, n_strdata_entries ! Loop indices and size tracking
    type (strdata_input_type)           :: surface_strdata_input_var ! temporary for setting
    type (strdata_input_type), pointer  :: surface_strdata_inputlist_tmp_ptr(:) ! temporary for reallocation

!-----------------------------------------------------------------------
! Loop through each antitracer and set up its forcing data source.
! This will create/append to surface_strdata_inputlist_ptr.
!-----------------------------------------------------------------------
    if (my_task == master_task) then
        write(stdout,'(A)') '==> DEBUG: Entering antitracer_init_sflux'
    end if

    ! Initialize surface_strdata_inputlist_ptr if not already done.
    ! It should be allocated to size 0 initially, then grow as needed.
    if (.not. associated(surface_strdata_inputlist_ptr)) then
        if (my_task == master_task) write(stdout,'(A)') '... DEBUG: Allocating surface_strdata_inputlist_ptr to size 0'
        allocate(surface_strdata_inputlist_ptr(0))
    end if

    do n_tracer = 1, antitracer_tracer_cnt
        if (my_task == master_task) then
            write(stdout,'(A,I0,A)') '---------------------------------------------------'
            write(stdout,'(A,I0)') '... DEBUG: Processing tracer #', n_tracer
        end if
        ! Get info for the current antitracer from the module-level array
        associate(forcing_info => all_antitracer_forcing_info(n_tracer))
            if (my_task == master_task) then
                write(stdout,'(A,A)')   '... DEBUG: Forcing Filename = ', trim(forcing_info%filename)
                write(stdout,'(A,A)')   '... DEBUG: Forcing Varname  = ', trim(forcing_info%file_varname)
                write(stdout,'(A,I0)')  '... DEBUG: Forcing Year First= ', forcing_info%year_first
                write(stdout,'(A,I0)')  '... DEBUG: Forcing Year Last = ', forcing_info%year_last
            end if

            ! Set up a temporary strdata_input_type for the current antitracer forcing.
            ! This defines the file, variable, and time parameters for shr_strdata.

            call POP_strdata_type_set(surface_strdata_input_var, &
                                      file_name = forcing_info%filename,      &
                                      field = forcing_info%file_varname,      &
                                      timer_label = 'antitracer_file_' // trim(forcing_info%name), &
                                      year_first = forcing_info%year_first,   &
                                      year_last = forcing_info%year_last,     &
                                      year_align = forcing_info%year_align,   &
                                      depth_flag = .false.,     & ! Forcing is 2D surface flux
                                      tintalgo = 'linear',      &
                                      taxMode = 'cycle')        ! Cycle or extend time series

            n_strdata_entries = size(surface_strdata_inputlist_ptr)
            if (my_task == master_task) write(stdout,'(A,I0)') '... DEBUG: Current number of unique strdata entries = ', n_strdata_entries
            forcing_info%surface_strdata_inputlist_ind = 0 ! Initialize to 'not found' state

            ! Check if a matching strdata entry (same file, same years, etc.) already exists
            ! in the global `surface_strdata_inputlist_ptr` array.
            do m = 1, n_strdata_entries
                if (my_task == master_task) write(stdout,'(A,I0)') '... DEBUG:   Checking for match with existing entry #', m
                if (POP_strdata_type_match(surface_strdata_input_var, surface_strdata_inputlist_ptr(m))) then
                    ! If a match is found, append the current tracer's variable name to the existing entry's field list.
                    if (my_task == master_task) write(stdout,'(A,I0,A)') '... DEBUG:   MATCH FOUND with entry #', m, '. Appending field.'
                    call POP_strdata_type_append_field(forcing_info%file_varname, surface_strdata_inputlist_ptr(m))
                    forcing_info%surface_strdata_inputlist_ind = m ! Store the index of this shared entry
                    exit ! Found a match, stop searching for this tracer
                end if
            end do

            ! If no match was found, create a new entry in `surface_strdata_inputlist_ptr`.
            if (forcing_info%surface_strdata_inputlist_ind == 0) then
                if (my_task == master_task) write(stdout,'(A)') '... DEBUG:   NO match found. Creating a new strdata entry.'
                n_strdata_entries = n_strdata_entries + 1
                if (my_task == master_task) write(stdout,'(A,I0)') '... DEBUG:   New total entries will be = ', n_strdata_entries
                ! Reallocate the pointer array to accommodate the new entry
                allocate(surface_strdata_inputlist_tmp_ptr(n_strdata_entries))
                ! Copy existing entries to the new larger array
                do m = 1, n_strdata_entries - 1
                    call POP_strdata_type_cp(surface_strdata_inputlist_ptr(m), surface_strdata_inputlist_tmp_ptr(m))
                end do
                deallocate(surface_strdata_inputlist_ptr)
                surface_strdata_inputlist_ptr => surface_strdata_inputlist_tmp_ptr

                ! Copy the new `strdata_input_var` (current tracer's details) to the new slot
                call POP_strdata_type_cp(surface_strdata_input_var, surface_strdata_inputlist_ptr(n_strdata_entries))
                forcing_info%surface_strdata_inputlist_ind = n_strdata_entries ! Store the index of the new entry
            endif

            ! Store the variable's index within the `shr_strdata` object's field list.
            ! This tells `sdat%avs` which variable to retrieve.
            forcing_info%surface_strdata_var_ind = POP_strdata_type_field_count( &
                surface_strdata_inputlist_ptr(forcing_info%surface_strdata_inputlist_ind))

            if (my_task == master_task) then
                write(stdout,'(A,I0)') '... DEBUG: This tracer will use strdata entry index: ', forcing_info%surface_strdata_inputlist_ind
                write(stdout,'(A,I0)') '... DEBUG: Its variable index within that entry is: ', forcing_info%surface_strdata_var_ind
            end if

        end associate ! forcing_info
    end do ! n_tracer

    if (my_task == master_task) then
        write(stdout,'(A,I0,A)') '---------------------------------------------------'
        write(stdout,'(A,I0,A)') '... DEBUG: Finished processing all tracers. Final number of unique strdata entries is ', &
                                 size(surface_strdata_inputlist_ptr), '.'
    end if

!-----------------------------------------------------------------------
!EOC

end subroutine antitracer_init_sflux

!***********************************************************************
!BOP
! !IROUTINE: antitracer_set_sflux
! !INTERFACE:

  subroutine antitracer_set_sflux(U10_SQR,IFRAC,SST, &
                                  SURF_VALS,STF_MODULE)

! !DESCRIPTION:
! Compute ANTITRACER surface flux and store related tavg fields for
! subsequent accumulating. This routine iterates over all defined
! antitracers, reads their forcing data, and computes air-sea gas exchange.

! !REVISION HISTORY:
! same as module

! !USES:

    use constants, only: xkw_coeff !, p5
    use timers, only: timer_start, timer_stop
    use domain,                 only : blocks_clinic

! !INPUT PARAMETERS:

    real (r8), dimension(nx_block,ny_block,:), intent(in) :: &
      U10_SQR,    & ! 10m wind speed squared (cm/s)**2
      IFRAC,      & ! sea ice fraction (non-dimensional)
      SST           ! sea surface temperature (C)

    ! SURF_VALS contains the current concentration of all antitracers at the surface.
    real (r8), dimension(nx_block,ny_block,antitracer_tracer_cnt,:), &
              intent(in) :: SURF_VALS

! !OUTPUT PARAMETERS:

    ! STF_MODULE will store the computed surface flux for all antitracers.
    real (r8), dimension(nx_block,ny_block,antitracer_tracer_cnt,:), &
              intent(inout) :: STF_MODULE

!EOP
!BOC
!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------
    character(*), parameter :: subname = 'antitracer_mod:antitracer_set_sflux'

    integer (int_kind) :: &
      iblock, n_tracer, m ! block and tracer indices

    ! These variables are computed per-block and are common for all antitracers
    ! for gas exchange calculations.
    real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      IFRAC_USED,      & ! used ice fraction (non-dimensional)
      XKW_USED           ! part of piston velocity (cm/s)

    real (r8), dimension(nx_block,ny_block) :: &
      ANTITRACER_SCHMIDT,    & ! ANTITRACER Schmidt number
      XKW_ICE,          & ! common portion of piston vel., (1-fice)*xkw (cm/s)
      PV                   ! piston velocity (cm/s)

    type(block)           :: this_block ! block info for the current block
    integer(int_kind)     :: i, j       ! loop indices for spatial iteration
    integer(int_kind)     :: n_idx      ! local linear index for rAttr array access
    integer(POP_i4)       :: errorCode ! error code for POP_HaloUpdate

    ! Temporary array to hold the forcing data from shr_strdata for a single tracer
    real(r8), dimension(nx_block,ny_block,max_blocks_clinic) :: current_tracer_forcing_data

    ! Static flag for first call within a timestep to prevent redundant `strdata_advance` calls.
    logical(log_kind), save :: first_call_this_timestep = .true.

!-----------------------------------------------------------------------

    call timer_start(antitracer_sflux_timer)

    !-----------------------------------------------------------------------
    ! Create the shr_strdata stream objects on the first timestep.
    ! This must be done here rather than in the init sequence to ensure all
    ! parallel components (PIO, MCT maps) are fully initialized.
    !-----------------------------------------------------------------------

    if (.not. antitracer_io_initialized) then

      do m = 1, size(surface_strdata_inputlist_ptr)
          if (my_task == master_task) then
              write(stdout,'(A,I0,A,A)') '... DEBUG: Calling POP_strdata_create for entry #', m, ' (File: ', &
                                        trim(surface_strdata_inputlist_ptr(m)%file_name), ')'
          end if
          call POP_strdata_create(surface_strdata_inputlist_ptr(m))
      end do
      antitracer_io_initialized = .true.
    endif

    if (my_task == master_task) then
        write(stdout,'(A,I0,A)') '---------------------------------------------------'
        write(stdout,'(A)') '==> DEBUG: Exiting antitracer_init_sflux'
    end if

!-----------------------------------------------------------------------
! Advance all unique shr_strdata streams.
! This is done only once per model timestep to ensure data is updated.
! `surface_strdata_inputlist_ptr` contains all unique shr_strdata objects.
!-----------------------------------------------------------------------

    if (first_call_this_timestep) then
        do m = 1, size(surface_strdata_inputlist_ptr)
            call POP_strdata_advance(surface_strdata_inputlist_ptr(m))
        end do
        first_call_this_timestep = .false.
    end if

!-----------------------------------------------------------------------
! Pre-compute common gas exchange parameters (not tracer-dependent)
! These only need to be computed once per timestep, per block.
! Moved calculation outside the `n_tracer` loop for efficiency.
!-----------------------------------------------------------------------

    do iblock = 1, nblocks_clinic
      where (LAND_MASK(:,:,iblock))
          IFRAC_USED(:,:,iblock) = IFRAC(:,:,iblock)
          XKW_USED(:,:,iblock) = xkw_coeff * U10_SQR(:,:,iblock)
      endwhere
      ! Clamp IFRAC_USED to valid range [0, 1]
      where (LAND_MASK(:,:,iblock) .and. IFRAC_USED(:,:,iblock) < c0) &
          IFRAC_USED(:,:,iblock) = c0
      where (LAND_MASK(:,:,iblock) .and. IFRAC_USED(:,:,iblock) > c1) &
          IFRAC_USED(:,:,iblock) = c1

      call comp_antitracer_schmidt(LAND_MASK(:,:,iblock), SST(:,:,iblock), &
                                    ANTITRACER_SCHMIDT)

      ! Calculate XKW_ICE and PV based on common parameters
      where (LAND_MASK(:,:,iblock))
          XKW_ICE = (c1 - IFRAC_USED(:,:,iblock)) * XKW_USED(:,:,iblock)
          PV = XKW_ICE * sqrt(660.0_r8 / ANTITRACER_SCHMIDT)
      elsewhere
          XKW_ICE(:,:) = c0
          PV(:,:) = c0
      endwhere
    end do

    ! Accumulate tavg fields related to general surface fluxes.
    ! These are accumulated once per timestep, using the common calculated values.
    ! They don't need to be inside the n_tracer loop.
    ! Note: A dedicated call `antitracer_tavg_forcing` handles this after all fluxes are set.
    ! However, if these were supposed to be accumulated here, it would be:
    ! do iblock = 1, nblocks_clinic
    !     call accumulate_tavg_field(IFRAC_USED(:,:,iblock),tavg_ANTITRACER_IFRAC,iblock,1)
    !     call accumulate_tavg_field(XKW_USED(:,:,iblock),tavg_ANTITRACER_XKW,iblock,1)
    !     call accumulate_tavg_field(ANTITRACER_SCHMIDT,tavg_ANTITRACER_SCHMIDT,iblock,1) ! Schmidt is 2D
    !     call accumulate_tavg_field(PV(:,:,iblock),tavg_ANTITRACER_PV,iblock,1)
    ! end do
    ! But since `antitracer_tavg_forcing` exists, we'll let it handle the module-level ANTITRACER_SFLUX_TAVG.

!-----------------------------------------------------------------------
! Loop over each antitracer to compute its specific surface flux
!-----------------------------------------------------------------------

    !$OMP PARALLEL DO PRIVATE(n_tracer, iblock, this_block, i, j, n_idx, errorCode)
    do n_tracer = 1, antitracer_tracer_cnt
        associate(forcing_info => all_antitracer_forcing_info(n_tracer))

            do iblock = 1, nblocks_clinic
                this_block = get_block(blocks_clinic(iblock), iblock)
                n_idx = 0 ! Local linear index within the rAttr array for this block
                do j = this_block%jb, this_block%je
                    do i = this_block%ib, this_block%ie
                        n_idx = n_idx + 1
                        current_tracer_forcing_data(i,j,iblock) = &
                            surface_strdata_inputlist_ptr(forcing_info%surface_strdata_inputlist_ind)%sdat%avs(forcing_info%surface_strdata_var_ind)%rAttr(forcing_info%surface_strdata_var_ind, n_idx)
                    enddo
                enddo
            enddo

            ! Apply halo update for this tracer's forcing field
            call POP_HaloUpdate(current_tracer_forcing_data, POP_haloClinic, &
                                POP_gridHorzLocCenter, POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
            if (errorCode /= POP_Success) then
                call document(subname, 'error updating halo for antitracer forcing field for '//trim(forcing_info%name))
                call exit_POP(sigAbort, 'Stopping in ' // subname)
            endif

            ! Initialize/Set STF_MODULE for this tracer with its forcing term
            do iblock = 1, nblocks_clinic
                where (LAND_MASK(:,:,iblock))
                    STF_MODULE(:,:,n_tracer,iblock) = forcing_info%scale_factor * current_tracer_forcing_data(:,:,iblock)
                elsewhere
                    STF_MODULE(:,:,n_tracer,iblock) = c0 ! No forcing over land
                endwhere
            enddo

            ! Apply air-sea gas exchange term (loss to atmosphere) for this tracer
            do iblock = 1, nblocks_clinic
                where (LAND_MASK(:,:,iblock))
                    ! This subtracts the loss due to gas exchange from the initial forcing.
                    ! STF_MODULE represents the NET surface flux into the ocean.
                    STF_MODULE(:,:,n_tracer,iblock) = &
                        STF_MODULE(:,:,n_tracer,iblock) - &
                        (PV / BETA) * SURF_VALS(:,:,n_tracer,iblock)
                elsewhere
                    STF_MODULE(:,:,n_tracer,iblock) = c0 ! Ensure zero flux over land
                endwhere
            end do ! iblock for gas exchange

        end associate ! forcing_info
    end do ! n_tracer
    !$OMP END PARALLEL DO

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
! This is a common utility function, not tracer-specific.
!
! range of validity of fit is -2:40
!
! Ref : Wanninkhof 2014, Relationship between wind speed
!         and gas exchange over the ocean revisited,
!         Limnol. Oceanogr.: Methods, 12,
!         doi:10.4319/lom.2014.12.351
!
! !REVISION HISTORY:
! same as module

! !USES:
    ! Add these for the NaN check and a clean exit
    use shr_infnan_mod, only : shr_infnan_isnan
    use exit_mod,       only : exit_POP, sigAbort
    use communicate,    only : my_task, master_task
    use io_types,       only : stdout

! !INPUT PARAMETERS:

    logical (log_kind), intent(in)  :: LAND_MASK(nx_block,ny_block)    ! land mask for this block
    real (r8)          , intent(in)  :: SST_IN(nx_block,ny_block)      ! sea surface temperature (C)

! !OUTPUT PARAMETERS:

    real (r8)          , intent(out) :: ANTITRACER_SCHMIDT(nx_block,ny_block)  ! Schmidt number of ANTITRACER (non-dimensional)

!EOP
!BOC
!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------
    integer(int_kind)    :: i, j
    real (r8)            :: SST(nx_block,ny_block)

    real (r8), parameter :: a = 2116.8_r8
    real (r8), parameter :: b = -136.25_r8
    real (r8), parameter :: c =    4.7353_r8
    real (r8), parameter :: d =   -0.092307_r8
    real (r8), parameter :: e =    0.000000697_r8

!-----------------------------------------------------------------------
    do j = 1, ny_block
      do i = 1, nx_block
          if (LAND_MASK(i,j)) then
             SST(i,j) = max(-2.0_r8, min(40.0_r8, SST_IN(i,j)))
             ANTITRACER_SCHMIDT(i,j) = a + SST(i,j) * (b + SST(i,j) * (c + SST(i,j) * (d + SST(i,j) * e)))
          else
             ANTITRACER_SCHMIDT(i,j) = c0
          endif
      end do
    end do

!-----------------------------------------------------------------------
!EOC

  end subroutine comp_antitracer_schmidt

!***********************************************************************

!BOP
! !IROUTINE: antitracer_tavg_forcing
! !INTERFACE:

  subroutine antitracer_tavg_forcing

! !DESCRIPTION:
! Make accumulation calls for forcing related tavg fields. This is
! necessary because the forcing routines are called before tavg flags
! are set. These tavg fields are generally for the overall "antitracer"
! system's physical parameters (like piston velocity), not for each
! individual antitracer specifically.

! !REVISION HISTORY:
! same as module

!EOP
!BOC
!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------

    integer (int_kind) :: &
      iblock             ! block loop index

!-----------------------------------------------------------------------

    !$OMP PARALLEL DO PRIVATE(iblock)

    do iblock = 1, nblocks_clinic
        ! The ANTITRACER_SFLUX_TAVG array holds quantities (IFRAC, XKW, SCHMIDT, PV)
        ! that are common to all antitracers' gas exchange calculations.
        ! These are calculated in `antitracer_set_sflux` and are independent of `n_tracer`.
        ! We accumulate them here once per timestep.
        call accumulate_tavg_field(ANTITRACER_SFLUX_TAVG(:,:,1,iblock),tavg_ANTITRACER_IFRAC,iblock,1)
        call accumulate_tavg_field(ANTITRACER_SFLUX_TAVG(:,:,2,iblock),tavg_ANTITRACER_XKW,iblock,1)
        call accumulate_tavg_field(ANTITRACER_SFLUX_TAVG(:,:,3,iblock),tavg_ANTITRACER_SCHMIDT,iblock,1)
        call accumulate_tavg_field(ANTITRACER_SFLUX_TAVG(:,:,4,iblock),tavg_ANTITRACER_PV,iblock,1)
    end do

    !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!EOC

  end subroutine antitracer_tavg_forcing

!***********************************************************************

end module antitracer_mod
