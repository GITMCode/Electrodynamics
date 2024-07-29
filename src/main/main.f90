
program test_01

  use ModCharSize
  use ModErrors
  use ModIE
  implicit none

  type(ieModel), pointer :: model
  character (len = iCharLenIE_) :: filename

  integer, parameter :: nMlts = 4
  integer, parameter :: nLats = 14
  real :: lats(nMlts, nLats)
  real :: mlts(nMlts, nLats)
  real :: potential(nMlts, nLats)
  real :: avee(nMlts, nLats)
  real :: eflux(nMlts, nLats)
  integer :: iMlt, iLat

  do iMlt = 1, nMlts
     do iLat = 1, nLats
        lats(iMlt, iLat) = 60.0 + iLat*2
        mlts(iMlt, iLat) = 0.0 + iMlt
     enddo
  enddo
  
  allocate(model)
  model = ieModel()

  call model % verbose(10)
  call model % efield_model("hepmay")
  call model % aurora_model("pem")
  call model % model_dir("data/ext/")
  call model % init()

  call model % imfBz(-5.0)
  call model % imfBy(0.0)
  call model % swv(450.0)
  call model % swn(5.0)
  call model % kp(4.0)
  call model % useAeHp()
  call model % au(100.0)
  call model % al(-500.0)

  call model % time_ymdhms(2001, 10, 29, 12, 0, 0)

  call model % nMlts(nMlts)
  call model % nLats(nLats)
  
  call model % grid(mlts, lats)
  call model % get_potential(potential)
  call model % get_aurora(eflux, avee)

  do iMlt = 1, nMlts
     do iLat = 1, nLats
        write(*,*) lats(iMlt, iLat), mlts(iMlt,iLat), &
             potential(iMlt, iLat), eflux(iMlt, iLat), avee(iMlt, iLat)
     enddo
  enddo
  
  call report_errors
  
end program test_01
