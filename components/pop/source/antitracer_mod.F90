!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module antitracer_mod

!BOP
! !MODULE: antitracer_mod
!
!  Module for antitracer 
!
!  The units of concentration for these tracers are 1/cm^3.
!
!  The units of surface fluxes for these tracers are
!     1/cm^3 * cm/s == 1/cm^2/s.
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  SVN:$Id: $

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
   use time_management, only : iyear, iday_of_year, frac_day, days_in_year

   use passive_tracer_tools, only: forcing_monthly_every_ts, ind_name_pair
   use passive_tracer_tools, only : read_field, tracer_read
   use broadcast

   implicit none
   save

!-----------------------------------------------------------------------
!  public/private declarations
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
!  module variables required by passive_tracers
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
       antitracer_tracer_cnt = 1

!-----------------------------------------------------------------------
!  relative tracer indices
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      antitracer_ind =  1

!-----------------------------------------------------------------------
!  derived type & parameter for tracer index lookup
!-----------------------------------------------------------------------
   
   type(ind_name_pair), dimension(antitracer_tracer_cnt) :: &
        ind_name_table = (/ &
        ind_name_pair(antitracer_ind, 'ANTITRACER')/)
   
!-----------------------------------------------------------------------
!  mask that eases avoidance of computation over land
!-----------------------------------------------------------------------

   logical (log_kind), dimension(:,:,:), allocatable :: &
      LAND_MASK

!-----------------------------------------------------------------------
!  forcing related variables
!-----------------------------------------------------------------------

   integer(int_kind)   :: antitracer_forcing_shr_stream_year_first   ! first year in stream to use
   integer(int_kind)   :: antitracer_forcing_shr_stream_year_last    ! last year in stream to use
   integer(int_kind)   :: antitracer_forcing_shr_stream_year_align   ! align ndep_shr_stream_year_first with this model year
   character(char_len) :: antitracer_forcing_shr_stream_file         ! file containing domain and input data
   real(r8)            :: antitracer_forcing_shr_stream_scale_factor ! unit conversion factor

   type (strdata_input_type), pointer :: surface_strdata_inputlist_ptr(:)

!-----------------------------------------------------------------------
!  define tavg id for 2d fields related to surface fluxes
!-----------------------------------------------------------------------

   real (r8), dimension(:,:,:,:), allocatable ::   &
      ANTITRACER_SFLUX_TAVG

   integer (int_kind) ::  &
      tavg_ANTITRACER_IFRAC,     & ! tavg id for ice fraction
      tavg_ANTITRACER_XKW,       & ! tavg id for xkw
      tavg_ANTITRACER_SCHMIDT,   & ! tavg id for antitracer Schmidt number
      tavg_ANTITRACER_PV,        & ! tavg id for antitracer piston velocity

!-----------------------------------------------------------------------
!  timers
!-----------------------------------------------------------------------

   integer (int_kind) :: antitracer_sflux_timer

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
!  Initialize antitracer module. This involves setting metadata, reading
!  the modules namelist and setting initial conditions.

! !REVISION HISTORY:
!  same as module

! !USES:

   use constants,  only: char_blank, delim_fmt
   use prognostic, only: curtime, oldtime, tracer_field
   use grid,       only: KMT, n_topo_smooth, fill_points
   use io_types,   only: nml_in, nml_filename
   use timers,     only: get_timer

   use passive_tracer_tools, only: init_forcing_monthly_every_ts, &
       rest_read_tracer_block, file_read_tracer_block

   use io_read_fallback_mod, only: io_read_fallback_register_tracer

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      antitracer_ind_begin          ! starting index of antitracers in global tracer array
                                    ! passed through to rest_read_tracer_block

   character (*), intent(in) ::  &
      init_ts_file_fmt,    & ! format (bin or nc) for input file
      read_restart_filename  ! file name for restart file

! !INPUT/OUTPUT PARAMETERS:

   type (tracer_field), dimension(antitracer_tracer_cnt), intent(inout) :: &
      tracer_d_module   ! descriptors for each tracer

   real (r8), dimension(nx_block,ny_block,km,antitracer_tracer_cnt,3,max_blocks_clinic), &
      intent(inout) :: TRACER_MODULE

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode         ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: subname = 'antitracer_mod:antitracer_init'

   character(char_len) :: &
      init_antitracer_option,        & ! option for initialization of bgc
      init_antitracer_init_file,     & ! filename for option 'file'
      init_antitracer_init_file_fmt    ! file format for option 'file'

   integer (int_kind) :: &
      n,                      & ! index for looping over tracers
      k,                      & ! index for looping over depth levels
      iblock,                 & ! index for looping over blocks
      nml_error                 ! namelist i/o error flag

   type(tracer_read), dimension(antitracer_tracer_cnt) :: &
      tracer_init_ext           ! namelist variable for initializing tracers

   namelist /antitracer_nml/ &
      init_antitracer_option, init_antitracer_init_file, init_antitracer_init_file_fmt, &
      tracer_init_ext, &
      antitracer_forcing_shr_stream_year_first,     &
      antitracer_forcing_shr_stream_year_last, antitracer_forcing_shr_stream_year_align, &
      antitracer_forcing_shr_stream_file, antitracer_forcing_shr_stream_scale_factor,    &


   character (char_len) ::  &
      antitracer_restart_filename      ! modified file name for restart file

!-----------------------------------------------------------------------
!  initialize tracer_d values
!-----------------------------------------------------------------------

   do n = 1, antitracer_tracer_cnt
      tracer_d_module(n)%short_name = ind_name_table(n)%name
      tracer_d_module(n)%long_name  = ind_name_table(n)%name
      tracer_d_module(n)%units      = '1/cm^3'
      tracer_d_module(n)%tend_units = '1/cm^3/s'
      tracer_d_module(n)%flux_units = '1/cm^2/s'
   end do

!-----------------------------------------------------------------------
!  default namelist settings
!-----------------------------------------------------------------------

   init_antitracer_option        = 'unknown'
   init_antitracer_init_file     = 'unknown'
   init_antitracer_init_file_fmt = 'bin'

   do n = 1, antitracer_tracer_cnt
      tracer_init_ext(n)%mod_varname  = 'unknown'
      tracer_init_ext(n)%filename     = 'unknown'
      tracer_init_ext(n)%file_varname = 'unknown'
      tracer_init_ext(n)%scale_factor = c1
      tracer_init_ext(n)%default_val  = c0
      tracer_init_ext(n)%file_fmt     = 'bin'
   end do

    antitracer_forcing_shr_stream_year_first = 1999
    antitracer_forcing_shr_stream_year_last = 2019
    antitracer_forcing_shr_stream_year_align = 347
    antitracer_forcing_shr_stream_file = '/glade/work/mclong/o-nets/data/forcing/alk-forcing.001.nc'
    antitracer_forcing_shr_stream_scale_factor = 1.0e4_r8  ! convert from 1/m^2/s to 1/cm^2/s

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=antitracer_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call document(subname, 'antitracer_nml not found')
      call exit_POP(sigAbort, 'stopping in ' /&
                           &/ subname)
   endif

!-----------------------------------------------------------------------
!  broadcast all namelist variables
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
   end do

!-----------------------------------------------------------------------
!   initialize tracers
!-----------------------------------------------------------------------

   select case (init_antitracer_option)

   case ('zero')
      TRACER_MODULE = c0
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
            call exit_POP(sigAbort, 'stopping in ' /&
                                 &/ subname)
         endif
         antitracer_restart_filename = read_restart_filename
         init_antitracer_init_file_fmt = init_ts_file_fmt

      else  ! do not read from TS restart file

         antitracer_restart_filename = trim(init_antitracer_init_file)

      endif

      call rest_read_tracer_block(antitracer_ind_begin,          &
                                  init_antitracer_init_file_fmt, &
                                  antitracer_restart_filename,   &
                                  tracer_d_module,        &
                                  TRACER_MODULE)

   case ('file')

      call document(subname, 'ANTITRACERs being read from separate file')

      call file_read_tracer_block(init_antitracer_init_file_fmt, &
                                  init_antitracer_init_file,     &
                                  tracer_d_module,        &
                                  ind_name_table,         &
                                  tracer_init_ext,        &
                                  TRACER_MODULE)

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
!  apply land mask to tracers
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
!  allocate and initialize LAND_MASK (true for ocean points)
!-----------------------------------------------------------------------

   allocate( LAND_MASK(nx_block,ny_block,max_blocks_clinic) )
   LAND_MASK = (KMT.gt.0)

   call get_timer(antitracer_sflux_timer, 'ANTITRACER_SFLUX', 1, distrb_clinic%nprocs)

!-----------------------------------------------------------------------
!  call other initialization subroutines
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
!  Define tavg fields not automatically handled by the base model.

! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      var_cnt             ! how many tavg variables are defined

!-----------------------------------------------------------------------

   var_cnt = 0

   call define_tavg_field(tavg_ANTITRACER_IFRAC,'ANTITRACER_IFRAC',2,           &
                          long_name='Ice Fraction for ANTITRACER fluxes',&
                          units='fraction', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_ANTITRACER_XKW,'ANTITRACER_XKW',2,               &
                          long_name='XKW for ANTITRACER fluxes',         &
                          units='cm/s', grid_loc='2110',          &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_ANTITRACER_SCHMIDT,'ANTITRACER_SCHMIDT',2,   &
                          long_name='ANTITRACER Schmidt Number',       &
                          units='none', grid_loc='2110',          &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_ANTITRACER_PV,'ANTITRACER_PV',2,             &
                          long_name='ANTITRACER piston velocity',      &
                          units='cm/s', grid_loc='2110',          &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

!-----------------------------------------------------------------------

   allocate(ANTITRACER_SFLUX_TAVG(nx_block,ny_block,var_cnt,max_blocks_clinic))
   ANTITRACER_SFLUX_TAVG = c0

!-----------------------------------------------------------------------
!EOC

 end subroutine antitracer_init_tavg

!***********************************************************************
!BOP
! !IROUTINE: antitracer_init_sflux
! !INTERFACE:

 subroutine antitracer_init_sflux(filename, file_varname, rank, &
      year_first, year_last, year_align, tintalgo, taxMode, strdata_inputlist_ptr)

! !DESCRIPTION:
!  Initialize surface flux computations for antitracer tracer module.

! !USES:

    use strdata_interface_mod, only : POP_strdata_type_set
    use strdata_interface_mod, only : POP_strdata_type_match
    use strdata_interface_mod, only : POP_strdata_type_append_field
    use strdata_interface_mod, only : POP_strdata_type_cp
    use strdata_interface_mod, only : POP_strdata_type_field_count

    character(len=*),                            intent(in)    :: filename
    character(len=*),                            intent(in)    :: file_varname
    integer (kind=int_kind),                     intent(in)    :: rank
    integer(kind=int_kind),            optional, intent(in)    :: year_first
    integer(kind=int_kind),            optional, intent(in)    :: year_last
    integer(kind=int_kind),            optional, intent(in)    :: year_align
    character(len=*),                  optional, intent(in)    :: tintalgo
    character(len=*),                  optional, intent(in)    :: taxMode
    type(strdata_input_type), pointer, optional, intent(inout) :: strdata_inputlist_ptr(:)

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: subname = 'antitracer_mod:antitracer_init_sflux'

   integer (int_kind) :: &
      n,                 & ! index for looping over tracers
      iblock               ! index for looping over blocks

   real (r8), dimension (nx_block,ny_block,12,max_blocks_clinic), target :: &
      WORK_READ            ! temporary space to read in fields

!-----------------------------------------------------------------------
   do n = 1, size(interior_strdata_inputlist_ptr)
     call POP_strdata_create(surface_strdata_inputlist_ptr(n))
   end do

   call POP_strdata_type_set(strdata_input_var, &
     file_name   = this%filename,     &
     field       = this%file_varname, &
     timer_label = 'marbl_file',      &
     year_first  = this%year_first,   &
     year_last   = this%year_last,    &
     year_align  = this%year_align,   &
     depth_flag  = (rank == 3),       &
     tintalgo    = tintalgo,          &
     taxMode     = taxMode)

!-----------------------------------------------------------------------
!EOC

 end subroutine antitracer_init_sflux

!***********************************************************************
!BOP
! !IROUTINE: antitracer_set_sflux
! !INTERFACE:

 subroutine antitracer_set_sflux(U10_SQR,IFRAC,SST, &
                          SURF_VALS,STF_MODULE)
! subroutine antitracer_set_sflux(U10_SQR,IFRAC,SST, &
!                          SURF_VALS_OLD,SURF_VALS_CUR,STF_MODULE)

! !DESCRIPTION:
!  Compute ANTITRACER surface flux and store related tavg fields for
!  subsequent accumulating.

! !REVISION HISTORY:
!  same as module

! !USES:

   use constants, only: xkw_coeff !, p5
   use timers, only: timer_start, timer_stop

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: &
      U10_SQR,   & ! 10m wind speed squared (cm/s)**2
      IFRAC,     & ! sea ice fraction (non-dimensional)
      SST,       & ! sea surface temperature (C)

   real (r8), dimension(nx_block,ny_block,antitracer_tracer_cnt,max_blocks_clinic), &
         intent(in) :: SURF_VALS ! module tracers
!         intent(in) :: SURF_VALS_OLD, SURF_VALS_CUR ! module tracers

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,antitracer_tracer_cnt,max_blocks_clinic), &
         intent(inout) :: STF_MODULE

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      iblock             ! block index

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      IFRAC_USED,      & ! used ice fraction (non-dimensional)
      XKW_USED,        & ! part of piston velocity (cm/s)

   real (r8), dimension(nx_block,ny_block) :: &
      !SURF_VALS,       & ! filtered surface tracer values
      ANTITRACER_SCHMIDT,   & ! ANTITRACER Schmidt number
      XKW_ICE,         & ! common portion of piston vel., (1-fice)*xkw (cm/s)
      PV,              & ! piston velocity (cm/s)

!-----------------------------------------------------------------------

   call timer_start(antitracer_sflux_timer)

!-----------------------------------------------------------------------
!   read antitracer forcing data
!-----------------------------------------------------------------------

   call POP_strdata_advance(surface_strdata_inputlist_ptr(:))

   stream_index = metadata%field_file_info%strdata_inputlist_ind
   var_ind      = metadata%field_file_info%strdata_var_ind

   n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock), iblock)
      do j = this_block%jb, this_block%je
         do i = this_block%ib, this_block%ie
            n = n + 1
            shr_stream(i,j,iblock) = surface_strdata_inputlist_ptr(stream_index)%sdat%avs(1)%rAttr(var_ind,n)
         enddo
      enddo
   enddo

   call POP_HaloUpdate(shr_stream, POP_haloClinic, &
        POP_gridHorzLocCenter, POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
   if (errorCode /= POP_Success) then
      call document(subname, 'error updating halo for shr_stream field')
      call exit_POP(sigAbort, 'Stopping in ' // subname)
   endif

   do iblock = 1, nblocks_clinic
      where (land_mask(:,:,iblock))
         forcing_field%field_0d(:,:,iblock) = shr_stream(:,:,iblock)
      endwhere
   enddo

   if (metadata%ltime_varying) then
      do iblock = 1, nblocks_clinic
         call apply_unit_conv_factor(land_mask(:,:,iblock), forcing_field, iblock)
      enddo
   end if

!-----------------------------------------------------------------------
!   compute air-sea gas exchange
!-----------------------------------------------------------------------

   do iblock = 1, nblocks_clinic
      IFRAC_USED(:,:,iblock) = c0
      XKW_USED(:,:,iblock) = c0
   end do

   !$OMP PARALLEL DO PRIVATE(iblock,SURF_VALS,ANTITRACER_SCHMIDT, &
   !$OMP                     XKW_ICE,PV)
   do iblock = 1, nblocks_clinic

      where (LAND_MASK(:,:,iblock))
         IFRAC_USED(:,:,iblock) = IFRAC(:,:,iblock)
         XKW_USED(:,:,iblock) = xkw_coeff * U10_SQR(:,:,iblock)
      endwhere
      where (LAND_MASK(:,:,iblock) .and. IFRAC_USED(:,:,iblock) < c0) &
         IFRAC_USED(:,:,iblock) = c0
      where (LAND_MASK(:,:,iblock) .and. IFRAC_USED(:,:,iblock) > c1) &
         IFRAC_USED(:,:,iblock) = c1

      call comp_antitracer_schmidt(LAND_MASK(:,:,iblock), SST(:,:,iblock), &
                            ANTITRACER_SCHMIDT)

      where (LAND_MASK(:,:,iblock))
         ANTITRACER_SFLUX_TAVG(:,:,1,iblock) = IFRAC_USED(:,:,iblock)
         ANTITRACER_SFLUX_TAVG(:,:,2,iblock) = XKW_USED(:,:,iblock)
         ANTITRACER_SFLUX_TAVG(:,:,3,iblock) = ANTITRACER_SCHMIDT

         XKW_ICE = (c1 - IFRAC_USED(:,:,iblock)) * XKW_USED(:,:,iblock)
         PV = XKW_ICE * sqrt(660.0_r8 / ANTITRACER_SCHMIDT)
         ANTITRACER_SFLUX_TAVG(:,:,4,iblock) = PV
                  
         !SURF_VALS = p5*(SURF_VALS_OLD(:,:,sf6_ind,iblock) + &
         !                SURF_VALS_CUR(:,:,sf6_ind,iblock))

         STF_MODULE(:,:,antitracer_ind,iblock) = &
            (PV / BETA) * SURF_VALS

      elsewhere
         STF_MODULE(:,:,antitracer_ind,iblock) = c0
      endwhere

   end do
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
!  Compute Schmidt numbers of ANTITRACERs.
!
!  range of validity of fit is -2:40
!
!  Ref : Wanninkhof 2014, Relationship between wind speed 
!        and gas exchange over the ocean revisited,
!        Limnol. Oceanogr.: Methods, 12, 
!        doi:10.4319/lom.2014.12.351
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   logical (log_kind), intent(in)  :: LAND_MASK(nx_block,ny_block)    ! land mask for this block
   real (r8)         , intent(in)  :: SST_IN(nx_block,ny_block)       ! sea surface temperature (C)

! !OUTPUT PARAMETERS:

   real (r8)         , intent(out) :: ANTITRACER_SCHMIDT(nx_block,ny_block)  ! Schmidt number of ANTITRACER (non-dimensional)

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
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
!  Make accumulation calls for forcing related tavg fields. This is
!  necessary because the forcing routines are called before tavg flags
!  are set.

! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      iblock              ! block loop index

!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock)

   do iblock = 1, nblocks_clinic
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

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
