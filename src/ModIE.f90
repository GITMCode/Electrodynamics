
MODULE ModIE

  use ModCharSize
  use ModW05_read_data
  use w05sc, only: setmodel, epotval
  use EIE_ModWeimer, only: get_tilt
  use ModFtaModel
  use ModIHP
  use ModTimeConvert
  use ModKind
  
  implicit none

  private

  integer, parameter, public :: iZero_ = 0
  integer, parameter, public :: iWeimer05_ = 1
  integer, parameter, public :: iMillstone_ = 2
  integer, parameter, public :: iHepMay_ = 3
  integer, parameter, public :: iFTA_ = 1
  integer, parameter, public :: iFRE_ = 2
  integer, parameter, public :: iPEM_ = 3
  integer, parameter, public :: iOvationPrime_ = 3
  integer, parameter, public :: iOvationSme_ = 4
  real, parameter, public :: rBadValue = -1.0e32

  integer, external :: efield_interpret_name
  integer, external :: aurora_interpret_name

  type, public :: ieModel

     logical :: isOk = .true.
     
     integer :: iDebugLevel = 0
     integer :: iEfield_ = -1
     integer :: iAurora_ = -1

     character (len = iCharLenIE_) :: modelDir = "data/ext/"
     character (len = iCharLenIE_) :: modelDirFta = "FTA/"     
     character (len = iCharLenIE_) :: northFile = "none"
     character (len = iCharLenIE_) :: southFile = "none"
     
     ! ----------------------------------------------------------------
     ! These are the states that the models has, if we either read in
     ! a file or we are coupling to some other model.
     ! They are 3D, because there is latitude, local time, and block,
     ! where block can be north/south or whatever.
     ! ----------------------------------------------------------------
     integer :: havenLats = 0
     integer :: havenMLTs = 0
     integer :: havenBLKs = 0
     real, allocatable, dimension(:,:,:) :: haveLats
     real, allocatable, dimension(:,:,:) :: haveMLTs 
     real, allocatable, dimension(:,:,:) :: havePotential
     real, allocatable, dimension(:,:,:) :: haveDiffuseEeFlux
     real, allocatable, dimension(:,:,:) :: haveDiffuseEAveE
     real, allocatable, dimension(:,:,:) :: haveDiffuseIeFlux
     real, allocatable, dimension(:,:,:) :: haveDiffuseIAveE

     ! ----------------------------------------------------------------
     ! These are what the code that is calling this library needs
     ! ----------------------------------------------------------------

     integer :: neednLats = 0
     integer :: neednMLTs = 0
     real, allocatable, dimension(:,:) :: needLats
     real, allocatable, dimension(:,:) :: needMLTs 
     
     integer :: iProc = 0

     integer :: havenTimes = 0
     real (kind = Real8_) :: currentTime = rBadValue

     real :: needImfBz = rBadValue
     real :: needImfBy = rBadValue
     real :: needSwV = rBadValue
     real :: needSwN = rBadValue
     real :: needHp = rBadValue
     real :: needHpN = rBadValue
     real :: needHpS = rBadValue
     real :: needAu = rBadValue
     real :: needAl = rBadValue
     real :: needAe = rBadValue
     real :: needKp = rBadValue
     logical :: useAeForHp = .false.

     real :: weimerTilt = 0.0

   contains

     ! Set verbose level:
     procedure :: verbose => set_verbose

     ! Set model types:
     procedure :: efield_model => set_efield_model
     procedure :: aurora_model => set_aurora_model
     procedure :: filename_north => set_filename_north
     procedure :: filename_south => set_filename_south
     
     ! Where to find the data files for empirical models:
     procedure :: model_dir => set_model_dir
     
     ! Initialize the library:
     procedure :: init => initialize

     ! set indices to run empirical models:
     procedure :: imfBz => set_bz
     procedure :: imfBy => set_by
     procedure :: swV => set_swv
     procedure :: swN => set_swn
     procedure :: hp => set_hp
     procedure :: hpN => set_hpn
     procedure :: hpS => set_hps
     procedure :: au => set_au
     procedure :: al => set_al
     procedure :: ae => set_ae
     procedure :: kp => set_kp
     procedure :: useAeHp => set_useAeForHp
     procedure :: dontAeHp => unset_useAeForHp
     procedure :: aehp => set_hp_from_ae
     procedure :: check_indices => run_check_indices

     ! Grid information for the calling code:
     procedure :: nMlts => set_nMlts
     procedure :: nLats => set_nLats
     procedure :: grid => set_grid
     procedure :: mlts => set_mlts
     procedure :: lats => set_lats

     procedure :: time_real => set_time_real
     procedure :: time_ymdhms => set_time_ymdhms
     procedure :: check_time => run_check_time

     ! Get model results:
     procedure :: get_potential => run_potential_model
     procedure :: get_aurora => run_aurora_model
     procedure :: weimer05 => run_weimer05_model
     procedure :: hepmay => run_heppner_maynard_model
     procedure :: fta => run_fta_model
     procedure :: hpi_pem => run_hpi_pem_model

  end type ieModel

contains

  ! ------------------------------------------------------------
  ! Initialize all of the different models
  subroutine initialize(this)
    class(ieModel) :: this
    character (len = iCharLenIE_) :: modelDirTotal
    character (len = iCharLenIE_) :: inFileNameTotal
    character (len = iCharLenIE_) :: name
    integer :: UnitTmp_ = 76
    integer :: iError = 0

    if (this%iDebugLevel > 1) &
         write(*,*) "==> Model data directory : ", trim(this%modelDir)

    !\
    ! --------------------------------------------------------------------
    ! Electric Field Models
    ! --------------------------------------------------------------------
    !/
    if (this % iEfield_ == iWeimer05_) &
         call read_all_files(this%modelDir)

    if (this % iEfield_ == iHepMay_) then
       inFileNameTotal = 'hmr89.cofcnts'
       call merge_str(this%modelDir, inFileNameTotal)
       open(UnitTmp_, file = inFileNameTotal, status='old', iostat = iError)
       if (iError /= 0) then
          call set_error('Error opening heppner maynard file :')
          call set_error(inFileNameTotal)
       else
          call gethmr(UnitTmp_)
          close(UnitTmp_)
       endif
    endif

    !\
    ! --------------------------------------------------------------------
    ! Aurora Models
    ! --------------------------------------------------------------------
    !/
    if (this % iAurora_ == iFta_) then
       modelDirTotal = this%modelDirFta
       call merge_str(this%modelDir, modelDirTotal)
       call initialize_fta(modelDirTotal)
    endif

    if (this % iAurora_ == iFRE_) then
       name = 'ihp'
       call read_conductance_model(name, this%modelDir, this%iDebugLevel)       
    endif
    if (this % iAurora_ == iPEM_) then
       name = 'pem'
       call read_conductance_model(name, this%modelDir, this%iDebugLevel)       
    endif
    
  end subroutine initialize

  ! ------------------------------------------------------------
  ! Set the verbose level for the library:
  subroutine set_verbose(this, level)
    class(ieModel) :: this
    integer, intent(in) :: level
    this%iDebugLevel = level
  end subroutine set_verbose
  
  ! ------------------------------------------------------------
  ! Check indices to see if the model has what it needs:
  subroutine run_check_time(ie)
    class(ieModel) :: ie
    if (ie%currentTime == rBadValue) then
       call set_error("current time not set!")
       call set_error("   --> need to call one of the set_time routines!")
    endif
    return
  end subroutine run_check_time

  ! ------------------------------------------------------------
  ! Set the current time to run the empirical model
  subroutine set_time_real(this, ut)
    class(ieModel) :: this
    real (kind = Real8_), intent(in) :: ut
    integer :: iYear
    integer :: iMonth
    integer :: iDay
    integer, dimension(7) :: itime
    real :: rHour
    
    this%currentTime = ut

    if (this % iEfield_ == iWeimer05_) then
       call time_real_to_int(ut, itime)
       iYear = itime(1)
       iMonth = itime(2)
       iDay = itime(3)
       rHour = float(iTime(4)) + float(iTime(5))/60.0
       this % weimerTilt = get_tilt(iYear, iMonth, iDay, rHour)
    endif
  end subroutine set_time_real
  
  ! ------------------------------------------------------------
  ! Set the current time to run the empirical model
  subroutine set_time_ymdhms(this, iYear, iMonth, iDay, iHour, iMin, iSec)
    class(ieModel) :: this
    integer, intent(in) :: iYear
    integer, intent(in) :: iMonth
    integer, intent(in) :: iDay
    integer, intent(in) :: iHour
    integer, intent(in) :: iMin
    integer, intent(in) :: iSec
    integer :: iTime(1:7)
    real (kind = Real8_) :: ut

    iTime(1) = iYear
    iTime(2) = iMonth
    iTime(3) = iDay
    iTime(4) = iHour
    iTime(5) = iMin
    iTime(6) = iSec
    iTime(7) = 0
    call time_int_to_real(iTime, ut)
    call this % time_real(ut)
  end subroutine set_time_ymdhms
  
  ! ------------------------------------------------------------
  ! Set indices routines
  ! ------------------------------------------------------------
  INCLUDE "model_subroutines.f90"
  
  ! ------------------------------------------------------------
  ! Set indices routines
  ! ------------------------------------------------------------
  INCLUDE "indices_subroutines.f90"
  
  ! ------------------------------------------------------------
  ! Set grid routines
  ! ------------------------------------------------------------
  INCLUDE "grid_subroutines.f90"
  
end MODULE ModIE
