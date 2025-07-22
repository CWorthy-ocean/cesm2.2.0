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
   use forcing_timeseries_mod, only: forcing_timeseries_dataset
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

   character(char_len) :: &
      antitracer_formulation,     & ! how to calculate flux (ocmip or model)

   integer (int_kind) ::  &
      model_year,             & ! arbitrary model year
      data_year,              & ! year in data that corresponds to model_year

   type (forcing_timeseries_dataset) :: &
      pantitracer_atm_forcing_dataset  ! data structure for atm pantitracer timeseries

   real (r8), dimension(:,:,:,:), allocatable :: &
      INTERP_WORK            ! temp array for interpolate_forcing output

   type(forcing_monthly_every_ts) :: &
      fice_file,           & ! ice fraction, if read from file
      xkw_file,            & ! a * wind-speed ** 2, if read from file
      ap_file                ! atmoshperic pressure, if read from file

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

   type(tracer_read) :: &
      gas_flux_fice,          & ! ice fraction for gas fluxes
      gas_flux_ws,            & ! wind speed for gas fluxes

   namelist /antitracer_nml/ &
      init_antitracer_option, init_antitracer_init_file, init_antitracer_init_file_fmt, &
      tracer_init_ext, model_year, data_year, &
      antitracer_formulation, gas_flux_fice, gas_flux_ws

   real (r8) :: &
      mapped_date               ! date of current model timestep mapped to data timeline

   character (char_len) ::  &
      antitracer_restart_filename      ! modified file name for restart file

!-----------------------------------------------------------------------
!  initialize forcing_monthly_every_ts variables
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call init_forcing_monthly_every_ts(fice_file)
   call init_forcing_monthly_every_ts(xkw_file)
   call init_forcing_monthly_every_ts(ap_file)

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

   gas_flux_fice%filename     = 'unknown'
   gas_flux_fice%file_varname = 'FICE'
   gas_flux_fice%scale_factor = c1
   gas_flux_fice%default_val  = c0
   gas_flux_fice%file_fmt     = 'bin'

   gas_flux_ws%filename     = 'unknown'
   gas_flux_ws%file_varname = 'XKW'
   gas_flux_ws%scale_factor = c1
   gas_flux_ws%default_val  = c0
   gas_flux_ws%file_fmt     = 'bin'

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

   call broadcast_scalar(model_year, master_task)
   call broadcast_scalar(data_year, master_task)
   call broadcast_scalar(antitracer_formulation, master_task)

   call broadcast_scalar(gas_flux_fice%filename, master_task)
   call broadcast_scalar(gas_flux_fice%file_varname, master_task)
   call broadcast_scalar(gas_flux_fice%scale_factor, master_task)
   call broadcast_scalar(gas_flux_fice%default_val, master_task)
   call broadcast_scalar(gas_flux_fice%file_fmt, master_task)

   fice_file%input = gas_flux_fice

   call broadcast_scalar(gas_flux_ws%filename, master_task)
   call broadcast_scalar(gas_flux_ws%file_varname, master_task)
   call broadcast_scalar(gas_flux_ws%scale_factor, master_task)
   call broadcast_scalar(gas_flux_ws%default_val, master_task)
   call broadcast_scalar(gas_flux_ws%file_fmt, master_task)

   xkw_file%input = gas_flux_ws

   call broadcast_scalar(gas_flux_ap%filename, master_task)
   call broadcast_scalar(gas_flux_ap%file_varname, master_task)
   call broadcast_scalar(gas_flux_ap%scale_factor, master_task)
   call broadcast_scalar(gas_flux_ap%default_val, master_task)
   call broadcast_scalar(gas_flux_ap%file_fmt, master_task)

   ap_file%input = gas_flux_ap

!-----------------------------------------------------------------------
!   initialize tracers
!-----------------------------------------------------------------------

   select case (init_antitracer_option)

   case ('ccsm_startup', 'zero', 'ccsm_startup_spunup')
      TRACER_MODULE = c0
      if (my_task == master_task) then
          write(stdout,delim_fmt)
          write(stdout,*) ' Initial 3-d antitracers set to all zeros'
          write(stdout,delim_fmt)
      endif

   case ('restart', 'ccsm_continue', 'ccsm_branch', 'ccsm_hybrid' )

      ! if mapped_date is less than pantitracer_first_nonzero_year then
      ! register c0 as an io_read_fallback option
      mapped_date = iyear + (iday_of_year-1+frac_day)/days_in_year &
                    - model_year + data_year
      if (mapped_date < pantitracer_first_nonzero_year) then
         call io_read_fallback_register_tracer(tracername='ANTITRACER', &
            fallback_opt='const', const_val=c0)
      endif

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

   call define_tavg_field(tavg_ANTITRACER_ATM_PRESS,'ANTITRACER_ATM_PRESS',2,   &
                          long_name='Atmospheric Pressure for ANTITRACER fluxes',&
                          units='atmospheres', grid_loc='2110',   &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_pANTITRACER,'pANTITRACER',2,                 &
                          long_name='ANTITRACER atmospheric partial pressure',&
                          units='pmol/mol', grid_loc='2110',      &
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

   call define_tavg_field(tavg_ANTITRACER_surf_sat,'ANTITRACER_surf_sat',2, &
                          long_name='ANTITRACER Saturation',           &
                          units='fmol/cm^3', grid_loc='2110',     &
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

 subroutine antitracer_init_sflux

! !USES:

   use forcing_tools, only: find_forcing_times
   use forcing_timeseries_mod, only: forcing_timeseries_init_dataset

! !DESCRIPTION:
!  Initialize surface flux computations for antitracer tracer module.
! !REVISION HISTORY:
!  same as module

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

   call forcing_timeseries_init_dataset(pantitracer_file, &
      varnames      = (/ 'ANTITRACERNH', 'ANTITRACERSH' /), &
      model_year    = model_year, &
      data_year     = data_year, &
      taxmode_start = 'endpoint', &
      taxmode_end   = 'extrapolate', &
      dataset       = pantitracer_atm_forcing_dataset)

!-----------------------------------------------------------------------
!  read gas flux forcing (if required)
!  otherwise, use values passed in
!-----------------------------------------------------------------------

   select case (antitracer_formulation)

   case ('ocmip')

!-----------------------------------------------------------------------
!  allocate space for interpolate_forcing
!-----------------------------------------------------------------------

      allocate(INTERP_WORK(nx_block,ny_block,max_blocks_clinic,1))

!-----------------------------------------------------------------------
!  first, read ice file
!-----------------------------------------------------------------------

      allocate(fice_file%DATA(nx_block,ny_block,max_blocks_clinic,1,12))

      call read_field(fice_file%input%file_fmt, &
                      fice_file%input%filename, &
                      fice_file%input%file_varname, &
                      WORK_READ)
      !$OMP PARALLEL DO PRIVATE(iblock, n)
      do iblock=1,nblocks_clinic
      do n=1,12
         fice_file%DATA(:,:,iblock,1,n) = WORK_READ(:,:,n,iblock)
         where (.not. LAND_MASK(:,:,iblock)) &
            fice_file%DATA(:,:,iblock,1,n) = c0
         fice_file%DATA(:,:,iblock,1,n) = &
            fice_file%DATA(:,:,iblock,1,n) * fice_file%input%scale_factor
      end do
      end do
      !$OMP END PARALLEL DO

      call find_forcing_times(fice_file%data_time, &
                              fice_file%data_inc, fice_file%interp_type, &
                              fice_file%data_next, fice_file%data_time_min_loc, &
                              fice_file%data_update, fice_file%data_type)

!-----------------------------------------------------------------------
!  next, read piston velocity file
!-----------------------------------------------------------------------

      allocate(xkw_file%DATA(nx_block,ny_block,max_blocks_clinic,1,12))

      call read_field(xkw_file%input%file_fmt, &
                      xkw_file%input%filename, &
                      xkw_file%input%file_varname, &
                      WORK_READ)

      !$OMP PARALLEL DO PRIVATE(iblock, n)
      do iblock=1,nblocks_clinic
      do n=1,12
         xkw_file%DATA(:,:,iblock,1,n) = WORK_READ(:,:,n,iblock)
         where (.not. LAND_MASK(:,:,iblock)) &
            xkw_file%DATA(:,:,iblock,1,n) = c0
         xkw_file%DATA(:,:,iblock,1,n) = &
            xkw_file%DATA(:,:,iblock,1,n) * xkw_file%input%scale_factor
      end do
      end do
      !$OMP END PARALLEL DO

      call find_forcing_times(xkw_file%data_time, &
                              xkw_file%data_inc, xkw_file%interp_type, &
                              xkw_file%data_next, xkw_file%data_time_min_loc, &
                              xkw_file%data_update, xkw_file%data_type)

!-----------------------------------------------------------------------
!  last, read atmospheric pressure file
!-----------------------------------------------------------------------

      allocate(ap_file%DATA(nx_block,ny_block,max_blocks_clinic,1,12))

      call read_field(ap_file%input%file_fmt, &
                      ap_file%input%filename, &
                      ap_file%input%file_varname, &
                      WORK_READ)

      !$OMP PARALLEL DO PRIVATE(iblock, n)
      do iblock=1,nblocks_clinic
      do n=1,12
         ap_file%DATA(:,:,iblock,1,n) = WORK_READ(:,:,n,iblock)
         where (.not. LAND_MASK(:,:,iblock)) &
            ap_file%DATA(:,:,iblock,1,n) = c0
         ap_file%DATA(:,:,iblock,1,n) = &
            ap_file%DATA(:,:,iblock,1,n) * ap_file%input%scale_factor
      end do
      end do
      !$OMP END PARALLEL DO

      call find_forcing_times(ap_file%data_time, &
                              ap_file%data_inc, ap_file%interp_type, &
                              ap_file%data_next, ap_file%data_time_min_loc, &
                              ap_file%data_update, ap_file%data_type)

   case ('model')

      if (my_task == master_task) then
         write(stdout,*)  &
            ' Using fields from model forcing for calculating ANTITRACER flux'
      endif

   case default
      call document(subname, 'antitracer_formulation', antitracer_formulation)

      call exit_POP(sigAbort, &
                    'antitracer_init_sflux: Unknown value for antitracer_formulation')

   end select

!-----------------------------------------------------------------------
!EOC

 end subroutine antitracer_init_sflux

!***********************************************************************
!BOP
! !IROUTINE: antitracer_set_sflux
! !INTERFACE:

 subroutine antitracer_set_sflux(U10_SQR,IFRAC,PRESS,SST,SSS, &
                          SURF_VALS_OLD,SURF_VALS_CUR,STF_MODULE)

! !DESCRIPTION:
!  Compute ANTITRACER surface flux and store related tavg fields for
!  subsequent accumulating.

! !REVISION HISTORY:
!  same as module

! !USES:

   use constants, only: field_loc_center, field_type_scalar, p5, xkw_coeff
   use time_management, only: thour00
   use forcing_tools, only: update_forcing_data, interpolate_forcing
   use timers, only: timer_start, timer_stop
   use forcing_timeseries_mod, only: forcing_timeseries_dataset_update_data

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: &
      U10_SQR,   & ! 10m wind speed squared (cm/s)**2
      IFRAC,     & ! sea ice fraction (non-dimensional)
      PRESS,     & ! sea level atmospheric pressure (dyne/cm**2)
      SST,       & ! sea surface temperature (C)
      SSS          ! sea surface salinity (psu)

   real (r8), dimension(nx_block,ny_block,antitracer_tracer_cnt,max_blocks_clinic), &
         intent(in) :: SURF_VALS_OLD, SURF_VALS_CUR ! module tracers

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
      AP_USED            ! used atm pressure (converted from dyne/cm**2 to atm)

   real (r8), dimension(nx_block,ny_block) :: &
      SURF_VALS,       & ! filtered surface tracer values
      pANTITRACER,          & ! atmospheric ANTITRACER mole fraction (pmol/mol)
      ANTITRACER_SCHMIDT,   & ! ANTITRACER Schmidt number
      ANTITRACER_SOL_0,     & ! solubility of ANTITRACER at 1 atm (mol/l/atm)
      XKW_ICE,         & ! common portion of piston vel., (1-fice)*xkw (cm/s)
      PV,              & ! piston velocity (cm/s)
      ANTITRACER_surf_sat       ! ANTITRACER surface saturation (fmol/cm^3)

   character (char_len) :: &
      tracer_data_label          ! label for what is being updated

   character (char_len), dimension(1) :: &
      tracer_data_names          ! short names for input data fields

   integer (int_kind), dimension(1) :: &
      tracer_bndy_loc,          &! location and field type for ghost
      tracer_bndy_type           !    cell updates

!-----------------------------------------------------------------------

   call timer_start(antitracer_sflux_timer)

   do iblock = 1, nblocks_clinic
      IFRAC_USED(:,:,iblock) = c0
      XKW_USED(:,:,iblock) = c0
      AP_USED(:,:,iblock) = c0
   end do

!-----------------------------------------------------------------------
!  Interpolate gas flux forcing data if necessary
!-----------------------------------------------------------------------

   call forcing_timeseries_dataset_update_data(pantitracer_atm_forcing_dataset)

   if (antitracer_formulation == 'ocmip') then
       if (thour00 >= fice_file%data_update) then
          tracer_data_names = fice_file%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'Ice Fraction'
          call update_forcing_data(          fice_file%data_time,   &
               fice_file%data_time_min_loc,  fice_file%interp_type, &
               fice_file%data_next,          fice_file%data_update, &
               fice_file%data_type,          fice_file%data_inc,    &
               fice_file%DATA(:,:,:,:,1:12), fice_file%data_renorm, &
               tracer_data_label,            tracer_data_names,     &
               tracer_bndy_loc,              tracer_bndy_type,      &
               fice_file%filename,           fice_file%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK, &
            fice_file%DATA(:,:,:,:,1:12), &
            fice_file%data_time,         fice_file%interp_type, &
            fice_file%data_time_min_loc, fice_file%interp_freq, &
            fice_file%interp_inc,        fice_file%interp_next, &
            fice_file%interp_last,       0)
       IFRAC_USED = INTERP_WORK(:,:,:,1)

       if (thour00 >= xkw_file%data_update) then
          tracer_data_names = xkw_file%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'Piston Velocity'
          call update_forcing_data(         xkw_file%data_time,   &
               xkw_file%data_time_min_loc,  xkw_file%interp_type, &
               xkw_file%data_next,          xkw_file%data_update, &
               xkw_file%data_type,          xkw_file%data_inc,    &
               xkw_file%DATA(:,:,:,:,1:12), xkw_file%data_renorm, &
               tracer_data_label,           tracer_data_names,    &
               tracer_bndy_loc,             tracer_bndy_type,     &
               xkw_file%filename,           xkw_file%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK, &
            xkw_file%DATA(:,:,:,:,1:12), &
            xkw_file%data_time,         xkw_file%interp_type, &
            xkw_file%data_time_min_loc, xkw_file%interp_freq, &
            xkw_file%interp_inc,        xkw_file%interp_next, &
            xkw_file%interp_last,       0)
       XKW_USED = INTERP_WORK(:,:,:,1)

       if (thour00 >= ap_file%data_update) then
          tracer_data_names = ap_file%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'Atmospheric Pressure'
          call update_forcing_data(        ap_file%data_time,   &
               ap_file%data_time_min_loc,  ap_file%interp_type, &
               ap_file%data_next,          ap_file%data_update, &
               ap_file%data_type,          ap_file%data_inc,    &
               ap_file%DATA(:,:,:,:,1:12), ap_file%data_renorm, &
               tracer_data_label,          tracer_data_names,   &
               tracer_bndy_loc,            tracer_bndy_type,    &
               ap_file%filename,           ap_file%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK, &
            ap_file%DATA(:,:,:,:,1:12), &
            ap_file%data_time,         ap_file%interp_type, &
            ap_file%data_time_min_loc, ap_file%interp_freq, &
            ap_file%interp_inc,        ap_file%interp_next, &
            ap_file%interp_last,       0)
       AP_USED = INTERP_WORK(:,:,:,1)
   endif

   !$OMP PARALLEL DO PRIVATE(iblock,SURF_VALS,pANTITRACER,ANTITRACER_SCHMIDT, &
   !$OMP                     ANTITRACER_SOL_0,XKW_ICE,&
   !$OMP                     PV,ANTITRACER_surf_sat)
   do iblock = 1, nblocks_clinic

      if (antitracer_formulation == 'ocmip') then
         where (LAND_MASK(:,:,iblock) .and. IFRAC_USED(:,:,iblock) < 0.2000_r8) &
            IFRAC_USED(:,:,iblock) = 0.2000_r8
         where (LAND_MASK(:,:,iblock) .and. IFRAC_USED(:,:,iblock) > 0.9999_r8) &
            IFRAC_USED(:,:,iblock) = 0.9999_r8
      endif

      if (antitracer_formulation == 'model') then
         where (LAND_MASK(:,:,iblock))
            IFRAC_USED(:,:,iblock) = IFRAC(:,:,iblock)

            XKW_USED(:,:,iblock) = xkw_coeff * U10_SQR(:,:,iblock)

            AP_USED(:,:,iblock) = PRESS(:,:,iblock)
         endwhere
         where (LAND_MASK(:,:,iblock) .and. IFRAC_USED(:,:,iblock) < c0) &
            IFRAC_USED(:,:,iblock) = c0
         where (LAND_MASK(:,:,iblock) .and. IFRAC_USED(:,:,iblock) > c1) &
            IFRAC_USED(:,:,iblock) = c1
      endif

!-----------------------------------------------------------------------
!  assume PRESS is in cgs units (dyne/cm**2) since that is what is
!    required for pressure forcing in barotropic
!  want units to be atmospheres
!  convertion from dyne/cm**2 to Pascals is P(mks) = P(cgs)/10.
!  convertion from Pascals to atm is P(atm) = P(Pa)/101.325e+3_r8
!-----------------------------------------------------------------------

      AP_USED(:,:,iblock) = AP_USED(:,:,iblock) * (c1 / 1013.25e+3_r8)

      call comp_pantitracer(iblock, LAND_MASK(:,:,iblock), pANTITRACER)

      call comp_antitracer_schmidt(LAND_MASK(:,:,iblock), SST(:,:,iblock), &
                            ANTITRACER_SCHMIDT)

      call comp_antitracer_sol_0(LAND_MASK(:,:,iblock), SST(:,:,iblock), SSS(:,:,iblock), &
                          ANTITRACER_SOL_0)

      where (LAND_MASK(:,:,iblock))
         ANTITRACER_SFLUX_TAVG(:,:,1,iblock) = IFRAC_USED(:,:,iblock)
         ANTITRACER_SFLUX_TAVG(:,:,2,iblock) = XKW_USED(:,:,iblock)
         ANTITRACER_SFLUX_TAVG(:,:,3,iblock) = AP_USED(:,:,iblock)
         ANTITRACER_SFLUX_TAVG(:,:,4,iblock) = pANTITRACER
         ANTITRACER_SFLUX_TAVG(:,:,5,iblock) = ANTITRACER_SCHMIDT

         XKW_ICE = (c1 - IFRAC_USED(:,:,iblock)) * XKW_USED(:,:,iblock)
         PV = XKW_ICE * sqrt(660.0_r8 / ANTITRACER_SCHMIDT)
         ANTITRACER_surf_sat = AP_USED(:,:,iblock) * ANTITRACER_SOL_0 * pANTITRACER
         SURF_VALS = p5*(SURF_VALS_OLD(:,:,antitracer_ind,iblock) + &
                         SURF_VALS_CUR(:,:,antitracer_ind,iblock))
         STF_MODULE(:,:,antitracer_ind,iblock) = &
            PV * (ANTITRACER_surf_sat - SURF_VALS)

         ANTITRACER_SFLUX_TAVG(:,:,6,iblock) = PV
         ANTITRACER_SFLUX_TAVG(:,:,7,iblock) = ANTITRACER_surf_sat

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
! !IROUTINE: comp_pantitracer
! !INTERFACE:

 subroutine comp_pantitracer(iblock, LAND_MASK, pANTITRACER)

! !DESCRIPTION:
!  Compute atmospheric mole fractions of ANTITRACERs
!  Linearly interpolate hemispheric values to current time step
!  Spatial pattern is determined by :
!     Northern Hemisphere value is used North of 10N
!     Southern Hemisphere value is used North of 10S
!     Linear Interpolation (in latitude) is used between 10N & 10S

! !REVISION HISTORY:
!  same as module

! !USES:

   use grid, only : TLATD
   use constants, only : c10
   use forcing_timeseries_mod, only: forcing_timeseries_dataset_get_var

! !INPUT PARAMETERS:

   logical (log_kind), dimension(nx_block,ny_block), intent(in) :: &
      LAND_MASK          ! land mask for this block

   integer (int_kind) :: &
      iblock          ! block index

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(out) :: &
      pANTITRACER  ! atmospheric ANTITRACER mole fraction (pmol/mol)

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j              ! loop indices

   real (r8) :: &
      pantitracer_nh_curr,   & ! pantitracer_nh for current time step (pmol/mol)
      pantitracer_sh_curr      ! pantitracer_sh for current time step (pmol/mol)

!-----------------------------------------------------------------------
!  Generate hemisphere values for current time step.
!
!  varind in the following calls must match varname ordering in 
!  call to forcing_timeseries_init_dataset in subroutine antitracer_init_sflux
!-----------------------------------------------------------------------

   call forcing_timeseries_dataset_get_var(pantitracer_atm_forcing_dataset, varind=1, data_1d=pantitracer_nh_curr)
   call forcing_timeseries_dataset_get_var(pantitracer_atm_forcing_dataset, varind=2, data_1d=pantitracer_sh_curr)

!-----------------------------------------------------------------------
!     Merge hemisphere values.
!-----------------------------------------------------------------------

   do j = 1, ny_block
      do i = 1, nx_block
         if (LAND_MASK(i,j)) then
            if (TLATD(i,j,iblock) < -c10) then
               pANTITRACER(i,j) = pantitracer_sh_curr
            else if (TLATD(i,j,iblock) > c10) then
               pANTITRACER(i,j) = pantitracer_nh_curr
            else
               pANTITRACER(i,j) = pantitracer_sh_curr + (TLATD(i,j,iblock)+c10) &
                  * 0.05_r8 * (pantitracer_nh_curr - pantitracer_sh_curr)
            endif
         endif
      end do
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine comp_pantitracer

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

   real (r8), parameter :: a = 3177.5_r8
   real (r8), parameter :: b = -200.57_r8
   real (r8), parameter :: c =    6.8865_r8
   real (r8), parameter :: d =   -0.13335_r8
   real (r8), parameter :: e =    0.0010877_r8

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
! !IROUTINE: comp_antitracer_sol_0
! !INTERFACE:

 subroutine comp_antitracer_sol_0(LAND_MASK, SST, SSS, ANTITRACER_SOL_0)

! !DESCRIPTION:
!  Compute solubilities of ANTITRACERs at 1 atm.
!  Ref: Bullister et al., 2002: The solubility of sulfur 
!       hexafluoride in water and seawater, DSR, 49(1),
!       doi:10.1016/S0967-0637(01)00051-6.
!
! !REVISION HISTORY:
!  same as module

! !USES:

   use constants, only: T0_Kelvin

! !INPUT PARAMETERS:

   logical (log_kind), dimension(nx_block,ny_block) :: &
      LAND_MASK          ! land mask for this block

   real (r8), dimension(nx_block,ny_block) :: &
      SST,             & ! sea surface temperature (C)
      SSS                ! sea surface salinity (psu)

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(out) :: &
      ANTITRACER_SOL_0  ! solubility of ANTITRACER at 1 atm (mol/l/atm)

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (r8), parameter :: &
      a1 = -96.5975_r8,    &
      a2 = 139.883_r8,     &
      a3 =  37.8193_r8,    &
      a4 =   0.00000_r8,   &
      b1 =   0.0310693_r8, &
      b2 =  -0.0356385_r8, &
      b3 =   0.00743254_r8

   real (r8), dimension(nx_block,ny_block) :: &
      SSTKp01  ! .01 * sea surface temperature (in Kelvin)

!-----------------------------------------------------------------------

   SSTKp01 = merge( ((SST + T0_Kelvin)* 0.01_r8), c1, LAND_MASK)

   where (LAND_MASK)
      ANTITRACER_SOL_0 = EXP(a1 + a2 / SSTKp01 &
                        + a3 * LOG(SSTKp01) + a4 * SSTKp01 ** 2 &
                        + SSS * (b1 + SSTKp01 * (b2 + b3 * SSTKp01)))
   elsewhere
      ANTITRACER_SOL_0 = c0
   endwhere

!-----------------------------------------------------------------------
!EOC

 end subroutine comp_antitracer_sol_0

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
         call accumulate_tavg_field(ANTITRACER_SFLUX_TAVG(:,:,3,iblock),tavg_ANTITRACER_ATM_PRESS,iblock,1)
         call accumulate_tavg_field(ANTITRACER_SFLUX_TAVG(:,:,4,iblock),tavg_pANTITRACER,iblock,1)
         call accumulate_tavg_field(ANTITRACER_SFLUX_TAVG(:,:,6,iblock),tavg_ANTITRACER_SCHMIDT,iblock,1)
         call accumulate_tavg_field(ANTITRACER_SFLUX_TAVG(:,:,6,iblock),tavg_ANTITRACER_PV,iblock,1)
         call accumulate_tavg_field(ANTITRACER_SFLUX_TAVG(:,:,7,iblock),tavg_ANTITRACER_surf_sat,iblock,1)
   end do

   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!EOC

 end subroutine antitracer_tavg_forcing

!***********************************************************************

end module antitracer_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
