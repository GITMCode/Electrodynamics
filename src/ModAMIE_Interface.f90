!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
Module ModAMIE_Interface

  use ModCharSize
  character (len=iCharLenIE_) :: AMIE_FileName
  integer :: AMIE_nLats, AMIE_nMlts, AMIE_nTimes

  ! For a single file
  real*4, allocatable,dimension(:) :: AMIE_Lats, AMIE_MLTs
  real*8, allocatable,dimension(:,:) :: AMIE_Time

  real*4, allocatable,dimension(:,:,:,:) :: AMIE_Potential
  real*4, allocatable,dimension(:,:,:,:) :: AMIE_PotentialY

  real*4, allocatable,dimension(:,:,:,:) :: AMIE_Value

  ! Electron Diffuse Aurora
  real*4, allocatable,dimension(:,:,:,:) :: AMIE_eleDiffEFlux
  real*4, allocatable,dimension(:,:,:,:) :: AMIE_eleDiffAveE

  ! Electron Discrete / Monoenergetic Aurora
  real*4, allocatable,dimension(:,:,:,:) :: AMIE_eleMonoEFlux
  real*4, allocatable,dimension(:,:,:,:) :: AMIE_eleMonoAveE

  ! Electron Broadband / Wave-driven Aurora
  real*4, allocatable,dimension(:,:,:,:) :: AMIE_eleWaveEFlux
  real*4, allocatable,dimension(:,:,:,:) :: AMIE_eleWaveAveE

  ! Ion Diffuse Aurora
  real*4, allocatable,dimension(:,:,:,:) :: AMIE_IonEFlux
  real*4, allocatable,dimension(:,:,:,:) :: AMIE_IonAveE

  integer, parameter :: AMIE_Closest_     = 1
  integer, parameter :: AMIE_After_       = 2
  integer, parameter :: AMIE_Interpolate_ = 3

  integer :: AMIE_iDebugLevel = 0

  integer :: AMIE_South_ = 1
  integer :: AMIE_North_ = 2

  integer, parameter :: potential_ = 1
  integer, parameter :: potentialy_ = 2
  integer, parameter :: ele_diff_eflux_ = 3
  integer, parameter :: ele_diff_avee_ = 4
  integer, parameter :: ion_diff_eflux_ = 5
  integer, parameter :: ion_diff_avee_ = 6
  integer, parameter :: ele_mono_eflux_ = 7
  integer, parameter :: ele_mono_avee_ = 8
  integer, parameter :: ele_wave_eflux_ = 9
  integer, parameter :: ele_wave_avee_ = 10

end Module ModAMIE_Interface
