
program test_01

  use ModMain
  use ModCharSize
  use ModErrors
  use ModIE
  use ModTimeConvert
  implicit none

  type(ieModel), pointer :: model

  integer :: iMlt, iLat, iLon

  real :: dLat, dMlt, dLon
  real (kind = Real8_) :: dt, jan1Time
  real :: date
  real :: apexHeight, apexLat, apexLon, Bmag, Bx, By, Bz, V
  real :: sbsllat, sbsllon, magMlt, sec
  real :: coLatPole, LonPole, vPol
  
  integer :: iTime_I(7), iDay

  ! This is to figure out the grid.
  ! The dMlt is the spacing in magnetic local time.  If you want geographic,
  ! then dMLT is the spacing in local time.
  isGeographic = .true.
  minLat = 50.
  dLat = 0.5
  dMlt = 0.25
  dLon = 15.0 * dMlt
  
  ! Define on edges:
  nLats = (90.0 - minLat)/dLat + 1
  nLons = 24.0/dMlt + 1

  call allocate_all_variables
  
  allocate(model)
  model = ieModel()

  call model % verbose(10)
  !call model % efield_model("hepmay")
  call model % efield_model("weimer05")

  ! set aurora model:
  !call model % aurora_model("hpi")
  !call model % aurora_model("pem")
  call model % aurora_model("fta")

  ! Set where the code can find the model files:
  call model % model_dir("data/ext/")
  call model % init()
  
  ! Handle Start Time:
  iTime_I(1) = 2011
  iTime_I(2) = 8
  iTime_I(3) = 11
  iTime_I(4) = 0
  iTime_I(5) = 0
  iTime_I(6) = 0
  iTime_I(7) = 0
  call time_int_to_real(iTime_I, startTime)

  ! Handle End Time:
  iTime_I(1) = 2011
  iTime_I(2) = 8
  iTime_I(3) = 12
  call time_int_to_real(iTime_I, endTime)

  if (isGeographic) then
     iTime_I(2) = 1
     iTime_I(3) = 1
     call time_int_to_real(iTime_I, jan1Time)
     date = iTime_I(1) + (startTime - jan1Time)/(24.0 * 3600.0)/365.25
     write(*,*) 'Date for Apex Calcs : ', date
     do iLon = 1, nLons
        do iLat = 1, nLats
           geoLats(iLon, iLat) = minLat + (iLat - 1) * dLat
           geoLons(iLon, iLat) = 0.0 + (iLon - 1) * dLon
           call apex(date, geoLats(iLon, iLat), geoLons(iLon, iLat), 120.0, &
                apexHeight, apexLat, apexLon, Bmag, Bx, By, Bz, V)
           magLats(iLon, iLat) = apexLat
           magLons(iLon, iLat) = apexLon
        enddo
     enddo
     call dypol (coLatPole, LonPole, vPol)
  else
     do iMlt = 1, nLons
        do iLat = 1, nLats
           magLats(iMlt, iLat) = minLat + (iLat-1) * dLat
           magMlts(iMlt, iLat) = 0.0 + (iMlt-1) * dMlt
        enddo
     enddo
  endif
  
  call model % nMlts(nLons)
  call model % nLats(nLats)
  if (.not. isGeographic) &
       call model % grid(magMlts, magLats)
  
  dt = (endTime - startTime) / nTimes
  
  do iT = 1, nTimes

     currentTime = startTime + dt * (iT - 1)
     call time_real_to_int(currentTime, iTime_I)
     write(*,*) " --> iTime : ", iTime_I
     ! Set the time in the IE library:
     call model % time_ymdhms( &
          iTime_I(1), iTime_I(2), iTime_I(3), &
          iTime_I(4), iTime_I(5), iTime_I(6))

     write(*,*) " --> currentTime : ", currentTime

     if (isGeographic) then
        ! Calculate the subsolar point:
        iDay = floor((startTime - jan1Time)/(24.0 * 3600.0)) + 1
        sec = iTime_I(6)
        call SUBSOLR (iTime_I(1), iDay, iTime_I(4), iTime_I(5), sec, &
             sbsllat, sbsllon)
        do iLon = 1, nLons
           do iLat = 1, nLats              
              apexLon = magLons(iLon, iLat)
              call magloctm(apexLon, sbsllat, sbsllon, coLatPole, LonPole, &
                   magMlt)
              magMlts(iLon, iLat) = magMlt
           enddo
        enddo
        call model % grid(magMlts, magLats)
     endif
     
     imfBx = 0.0
     imfBy = -2.5
     imfBz = -5.0
     swV = 400.0
     swN = 5.0
     au = 100.0
     al = -500.0
     ae = au - al
     hp = 0.102 * ae + 8.953

     ! Set indices in the IE library:
     call model % imfBz(imfBz)
     call model % imfBy(imfBy)
     call model % swv(swV)
     call model % swn(swN)
     call model % kp(4.0)
     call model % useAeHp()
     call model % au(au)
     call model % al(al)
          
     ! Get potential from the IE library:
     call model % get_potential(potential)
     write(*,*) '  --> CPCP : ', maxval(potential) - minval(potential)
     ! Get electron diffuse aurora from the IE library:
     call model % get_aurora(eDiffuseEflux, eDiffuseAvee)
     write(*,*) '  --> max eflux : ', maxval(eDiffuseEflux)
     write(*,*) '  --> max ave-e : ', maxval(eDiffuseAveE)
     call output_ie

  enddo

  write(LunMainOutput) 0.1
  close(LunMainOutput)

  call report_errors
  
end program test_01
