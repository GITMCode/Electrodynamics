
MODULE ModIE

  use ModCharSize
  use ModW05_read_data
  use ModFtaModel

  implicit none

  private

  integer, parameter, public :: iZero_ = 0
  integer, parameter, public :: iWeimer05_ = 1
  integer, parameter, public :: iMillstone_ = 2
  integer, parameter, public :: iFTA_ = 1
  integer, parameter, public :: iFRE_ = 2
  integer, parameter, public :: iOvationPrime_ = 3
  integer, parameter, public :: iOvationSme_ = 4

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

     integer :: iProc = 0

     integer :: havenTimes = 0
     real*8 :: currentTime = 0.0

     real :: needImfBz = -1e32
     real :: needImfBy = -1e32
     real :: needSwV = -1e32
     real :: needSwN = -1e32
     real :: needHp = -1e32
     real :: needHpN = -1e32
     real :: needHpS = -1e32
     real :: needAu = -1e32
     real :: needAl = -1e32
     real :: needAe = -1e32
     logical :: useAeForHp = .false.

   contains

     procedure :: verbose => set_verbose
     procedure :: efield_model => set_efield_model
     procedure :: aurora_model => set_aurora_model
     procedure :: filename_north => set_filename_north
     procedure :: filename_south => set_filename_south
     procedure :: model_dir => set_model_dir
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
     procedure :: useAeHp => set_useAeForHp
     procedure :: dontAeHp => unset_useAeForHp
     procedure :: aehp => set_hp_from_ae

  end type ieModel

contains

  ! ------------------------------------------------------------
  ! Initialize all of the different models
  subroutine initialize(this)
    class(ieModel) :: this
    character (len = iCharLenIE_) :: modelDirTotal

    if (this%iDebugLevel > 1) &
         write(*,*) "==> Model data directory : ", trim(this%modelDir)

    !\
    ! --------------------------------------------------------------------
    ! Electric Field Models
    ! --------------------------------------------------------------------
    !/
    if (this % iEfield_ == iWeimer05_) &
         call read_all_files(this%modelDir)

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
    
  end subroutine initialize

  
  ! ------------------------------------------------------------
  ! Set the verbose level for the library:
  subroutine set_verbose(this, level)
    class(ieModel) :: this
    integer, intent(in) :: level
    this%iDebugLevel = level
  end subroutine set_verbose
  
  ! ------------------------------------------------------------
  ! set efield model
  subroutine set_efield_model(this, efield_model)
    class(ieModel) :: this
    character (len = *), intent(in) :: efield_model
    if (this%iDebugLevel > 0) &
         write(*,*) "=> Setting efield model to : ", efield_model
    this%iEfield_ = efield_interpret_name(efield_model)
    if (this%iDebugLevel > 0) &
         write(*,*) "=> That is model : ", this%iEfield_
  end subroutine set_efield_model
  
  ! ------------------------------------------------------------
  ! set aurora model
  subroutine set_aurora_model(this, aurora_model)
    class(ieModel) :: this
    character (len = *), intent(in) :: aurora_model
    if (this%iDebugLevel > 0) &
         write(*,*) "=> Setting aurora model to : ", aurora_model
    this%iAurora_ = aurora_interpret_name(aurora_model)
    if (this%iDebugLevel > 0) &
         write(*,*) "=> That is model : ", this%iAurora_
  end subroutine set_aurora_model
  
  ! ------------------------------------------------------------
  ! filename for north
  subroutine set_filename_north(this, filename)
    class(ieModel) :: this
    character (len = *), intent(in) :: filename
    if (this%iDebugLevel > 0) &
         write(*,*) "=> Setting north file name to : ", filename
    this%northFile = filename
  end subroutine set_filename_north
  
  ! ------------------------------------------------------------
  ! filename for south
  subroutine set_filename_south(this, filename)
    class(ieModel) :: this
    character (len = *), intent(in) :: filename
    if (this%iDebugLevel > 0) &
         write(*,*) "=> Setting south file name to : ", filename
    this%southFile = filename
  end subroutine set_filename_south
  
  ! ------------------------------------------------------------
  ! set the directory for finding all of the coef files
  subroutine set_model_dir(this, dir)
    class(ieModel) :: this
    character (len = *), intent(in) :: dir
    if (this%iDebugLevel > 0) &
         write(*,*) "=> Setting model directory to : ", dir
    this%modelDir = dir
  end subroutine set_model_dir

  ! ------------------------------------------------------------
  ! Set indices routines
  ! ------------------------------------------------------------

  ! ------------------------------------------------------------
  ! set IMF Bz
  subroutine set_bz(this, value)
    class(ieModel) :: this
    real, intent(in) :: value
    if (this%iDebugLevel > 2) &
         write(*,*) "=> Setting imf bz : ", value
    this%needImfBz = value
  end subroutine set_bz

  ! ------------------------------------------------------------
  ! set IMF By
  subroutine set_by(this, value)
    class(ieModel) :: this
    real, intent(in) :: value
    if (this%iDebugLevel > 2) &
         write(*,*) "=> Setting imf by : ", value
    this%needImfBy = value
  end subroutine set_by

  ! ------------------------------------------------------------
  ! set solar wind velocity
  subroutine set_swv(this, value)
    class(ieModel) :: this
    real, intent(in) :: value
    if (this%iDebugLevel > 2) &
         write(*,*) "=> Setting Solar Wind Velocity: ", value
    ! Make sure that this velocity is a positive value:
    this%needSwV = abs(value)
  end subroutine set_swv

  ! ------------------------------------------------------------
  ! set solar wind density
  subroutine set_swn(this, value)
    class(ieModel) :: this
    real, intent(in) :: value
    if (this%iDebugLevel > 2) &
         write(*,*) "=> Setting Solar Wind Density: ", value
    this%needSwN = value
  end subroutine set_swn

  ! ------------------------------------------------------------
  ! set the hemispheric power (in gigawatts)
  subroutine set_hp(this, value)
    class(ieModel) :: this
    real, intent(in) :: value
    if (this%iDebugLevel > 2) &
         write(*,*) "=> Setting Hemispheric Power (in both hems) : ", value
    this%needHp = value
    ! If the user calls this routine, then it sets the north and south HPs
    this%needHpN = value
    this%needHpS = value
  end subroutine set_hp

  ! ------------------------------------------------------------
  ! set north hemispheric power (in gigawatts)
  subroutine set_hpn(this, value)
    class(ieModel) :: this
    real, intent(in) :: value
    if (this%iDebugLevel > 2) &
         write(*,*) "=> Setting Hemispheric Power (in north) : ", value
    this%needHpN = value
  end subroutine set_hpn

  ! ------------------------------------------------------------
  ! set south hemispheric power (in gigawatts)
  subroutine set_hps(this, value)
    class(ieModel) :: this
    real, intent(in) :: value
    if (this%iDebugLevel > 2) &
         write(*,*) "=> Setting Hemispheric Power (in south) : ", value
    this%needHpS = value
  end subroutine set_hps

  ! ------------------------------------------------------------
  ! set Hemispheric Power from AE
  subroutine set_hp_from_ae(this, value)
    class(ieModel) :: this
    real, intent(in) :: value
    if (this%useAeForHp) then
       ! Set the hemispheric power based on the AE index:
       this%needHp = 0.102 * value + 8.953
       this%needHpN = this%needHp
       this%needHpS = this%needHp
       if (this%iDebugLevel > 2) then
          write(*,*) "=> Setting HP from AE. AE = ", value
          write(*,*) "                       HP = ", this%needHp
       endif
    endif
  end subroutine set_hp_from_ae
         
  ! ------------------------------------------------------------
  ! set AU (and derive AE)
  subroutine set_au(this, value)
    class(ieModel) :: this
    real, intent(in) :: value
    if (this%iDebugLevel > 2) &
         write(*,*) "=> Setting AU (and deriving AE): ", value
    this%needAu = value
    ! derive AE, done in both functions so that the last one called works
    this%needAe = this%needAu - this%needAl
    call this%aehp(this%needAe)
  end subroutine set_au

  ! ------------------------------------------------------------
  ! set AL (and derive AE)
  subroutine set_al(this, value)
    class(ieModel) :: this
    real, intent(in) :: value
    if (this%iDebugLevel > 2) &
         write(*,*) "=> Setting AL (and deriving AE): ", value
    this%needAl = value
    ! derive AE, done in both functions so that the last one called works
    this%needAe = this%needAu - this%needAl
    call this%aehp(this%needAe)
  end subroutine set_al

  ! ------------------------------------------------------------
  ! set AE
  subroutine set_ae(this, value) 
    class(ieModel) :: this
    real, intent(in) :: value
    if (this%iDebugLevel > 2) &
         write(*,*) "=> Setting AE : ", value
    this%needAe = value
    call this%aehp(this%needAe)
  end subroutine set_ae

  ! ------------------------------------------------------------
  ! set Use AE to determine HP to true
  subroutine set_useAeForHp(this)
    class(ieModel) :: this
    if (this%iDebugLevel > 2) &
         write(*,*) "=> Using AE to specify the HP!"
    this%useAeForHp = .true.
  end subroutine set_useAeForHp

  ! ------------------------------------------------------------
  ! set Use AE to determine HP to false
  subroutine unset_useAeForHp(this)
    class(ieModel) :: this
    if (this%iDebugLevel > 2) &
         write(*,*) "=> NOT Using AE to specify the HP!"
    this%useAeForHp = .false.
  end subroutine unset_useAeForHp

  
end MODULE ModIE
