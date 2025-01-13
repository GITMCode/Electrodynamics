
  ! ------------------------------------------------------------
  ! set efield model
  subroutine set_efield_model(this, efield_model)
    class(ieModel) :: this
    character (len = *), intent(in) :: efield_model
    if (this%iDebugLevel > 0) &
         write(*,*) "=> Setting efield model to : ", trim(efield_model)
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
         write(*,*) "=> Setting aurora model to : ", trim(aurora_model)
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
         write(*,*) "=> Setting north file name to : ", trim(filename)
    this%northFile = filename
  end subroutine set_filename_north
  
  ! ------------------------------------------------------------
  ! filename for south
  subroutine set_filename_south(this, filename)
    class(ieModel) :: this
    character (len = *), intent(in) :: filename
    if (this%iDebugLevel > 0) &
         write(*,*) "=> Setting south file name to : ", trim(filename)
    this%southFile = filename
  end subroutine set_filename_south
  
  ! ------------------------------------------------------------
  ! set the directory for finding all of the coef files
  subroutine set_model_dir(this, dir)
    class(ieModel) :: this
    character (len = *), intent(in) :: dir
    if (this%iDebugLevel > 0) &
         write(*,*) "=> Setting model directory to : ", trim(dir)
    this%modelDir = dir
  end subroutine set_model_dir

  ! ------------------------------------------------------------
  ! run efield model

  subroutine run_potential_model(ie, potential)

    use ModAMIE_Interface, only: get_amie_potential
    class(ieModel) :: ie
    real, dimension(ie%neednMlts, &
                    ie%neednLats), intent(out) :: potential
    real :: currentTilt = rBadValue, lastTilt = rBadValue
    real :: potVal

    integer :: iMLT, iLat
    potential = 0.0

    call ie%check_time()
    call ie%check_indices()

    if (.not. isOk) then
      call set_error("Not ok after checking time & indices in run_potential_model.")
      call report_errors
      ! return ! debated whether to return or not. Probably let the code try and run...
    endif

    if (ie%iEfield_ == iWeimer05_) call ie%weimer05(potential)
    if (ie%iEfield_ == iHepMay_) call ie%hepmay(potential)

    if (ie%iEfield_ == iAmiePot_) call get_amie_potential(potential)

    return
  end subroutine run_potential_model


  ! ------------------------------------------------------------
  ! run Weimer05
  subroutine run_weimer05_model(ie, potential)
    class(ieModel) :: ie
    real, dimension(ie%neednMlts, &
                    ie%neednLats), intent(inout) :: potential
    real :: currentTilt = rBadValue, lastTilt = rBadValue
    real :: potVal

    integer :: iMLT, iLat

    do iMLT = 1, ie%neednMLTs
       do iLat = 1, ie%neednLats
          if (abs(ie%needLats(iMlt, iLat)) > 45.0) then
             ! this is to check if we have changed hemispheres:
             currentTilt = sign(ie%weimerTilt, ie%needLats(iMlt, iLat))
             if (currentTilt .ne. lastTilt) then
                ! Only need to set up the model once, when everything
                ! stays the same (including the hemisphere!):
                call setmodel( &
                     ie%needIMFBy, &
                     ie%needIMFBz, &
                     currentTilt, &
                     ie%needSWV, &
                     ie%needSWN, 'epot')
                lastTilt = currentTilt
             endif
             ! Run Weimer for specific lat and mlt:
             call epotval( &
                  abs(ie%needLats(iMlt, iLat)), &
                  ie%needMlts(iMlt, iLat), &
                  0.0, &
                  potVal)
             ! Store potential and convert to V:
             potential(iMlt, iLat) = potVal * 1000.0
          endif
       enddo
    enddo

    return

  end subroutine run_weimer05_model
  
  ! ------------------------------------------------------------
  ! run Heppner Maynard model
  subroutine run_heppner_maynard_model(ie, potential)
    class(ieModel) :: ie
    real, dimension(ie%neednMlts, &
                    ie%neednLats), intent(inout) :: potential
    real :: currentTilt = rBadValue, lastTilt = rBadValue
    real :: potVal, eTheta, ePhi

    integer :: iMLT, iLat, iFirst

    iFirst = 1
    do iMLT = 1, ie%neednMLTs
       do iLat = 1, ie%neednLats
          if (abs(ie%needLats(iMlt, iLat)) > 50.0) then
             call hmrepot( &
                  ie%needLats(iMlt, iLat), &
                  ie%needMlts(iMlt, iLat), &
                  ie%needIMFBy, &
                  ie%needIMFBz, &
                  ie%needKp, &
                  iFirst, &
                  eTheta, &
                  ePhi, &
                  potVal)
             potential(iMlt, iLat) = potVal * 1000.0
             iFirst = iFirst + 1
          endif
       enddo
    enddo

    return
  end subroutine run_heppner_maynard_model

  ! ------------------------------------------------------------
  ! run aurora model

  subroutine run_aurora_model(ie, eflux, avee)
    use ModAMIE_Interface, only: &
         get_amie_electron_diffuse_eflux, &
         get_amie_electron_diffuse_avee
    class(ieModel) :: ie
    real, dimension(ie%neednMlts, &
                    ie%neednLats), intent(out) :: eFlux
    real, dimension(ie%neednMlts, &
                    ie%neednLats), intent(out) :: AveE

    eFlux = 0.001 ! ergs/cm2/s
    AveE = 3.0 ! keV
    
    call ie%check_time()
    call ie%check_indices()

    if (.not. isOk) then ! I think this is worthwhile to have, but masks errors in potential.
      call set_error("Not ok after checking time & indices in run_aurora_model.")
      call report_errors
      ! return
    endif

    if (ie%iAurora_ == iFTA_) call ie%fta(eFlux, AveE)
    ! These two models are the same, because they use the same
    if (ie%iAurora_ == iFRE_) call ie%hpi_pem(eFlux, AveE)
    if (ie%iAurora_ == iPEM_) call ie%hpi_pem(eFlux, AveE)

    if (ie%iAurora_ == iAmieAur_) then
       call get_amie_electron_diffuse_eflux(eFlux)
       call get_amie_electron_diffuse_avee(AveE)
    endif
    
    return
  end subroutine run_aurora_model


  ! ------------------------------------------------------------
  ! run FTA Model, providing aveE and E-Flux
  subroutine run_fta_model(ie, eFlux, AveE)
    class(ieModel) :: ie
    real, dimension(ie%neednMlts, &
                    ie%neednLats), intent(inout) :: eFlux
    real, dimension(ie%neednMlts, &
                    ie%neednLats), intent(inout) :: AveE
    real :: eFluxVal, AveEVal
    integer :: iError = 0, iMlt, iLat
    
    call update_fta_model(ie%needAu, ie%needAl, iError)
    if (iError /= 0) then
       call set_error('FTA Model update has an error!')
       return
    endif

    do iMLT = 1, ie%neednMLTs -1
       do iLat = 1, ie%neednLats -1
          call get_fta_model_result( &
               ie%needMlts(iMlt, iLat), &
               ie%needLats(iMlt, iLat), &
               eFluxVal, AveEVal)
          eFlux(iMlt, iLat) = eFluxVal
          AveE(iMlt, iLat) = AveEVal
       enddo
    enddo
    return
  end subroutine run_fta_model


  ! ------------------------------------------------------------
  ! run HPI Model, providing aveE and E-Flux
  subroutine run_hpi_pem_model(ie, eFlux, AveE)
    class(ieModel) :: ie
    real, dimension(ie%neednMlts, &
                    ie%neednLats), intent(inout) :: eFlux
    real, dimension(ie%neednMlts, &
                    ie%neednLats), intent(inout) :: AveE
    real :: eFluxVal, AveEVal, hp
    integer :: iError = 0, iMlt, iLat
    
    do iMLT = 1, ie%neednMLTs
       do iLat = 1, ie%neednLats
          if (ie%needLats(iMlt, iLat) > 0) then
             hp = ie%needHpN
          else
             hp = ie%needHpS
          endif
          call get_auroral_conductance( &
               ie%needMlts(iMlt, iLat), &
               ie%needLats(iMlt, iLat), &
               hp, &
               eFluxVal, AveEVal)
          eFlux(iMlt, iLat) = eFluxVal
          AveE(iMlt, iLat) = AveEVal
       enddo
    enddo
    return
  end subroutine run_hpi_pem_model

  
  
