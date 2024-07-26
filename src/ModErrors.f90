!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModErrors

  use ModCharSize

  implicit none
  
  integer, parameter :: nErrorsMax = 1000
  character (len=iCharLenIE_), dimension(nErrorsMax) :: cErrorCodes

  integer :: nErrors = 0
  logical :: isOk = .true.

contains

  subroutine set_error(cError)
    character (len = *), intent(in) :: cError
    nErrors = nErrors + 1
    cErrorCodes(nErrors) = cError
    isOk = .false.
  end subroutine set_error

  subroutine report_errors()
    integer :: iError
    if (nErrors == 0) write(*, *) "No errors to report!"
    do iError = 1, nErrors
       write(*, *) "--> Error : ", trim(cErrorCodes(iError))
    enddo
  end subroutine report_errors

  
end module ModErrors
