
program test_01

  use ModCharSize
  use ModErrors
  use ModIE
  implicit none

  type(ieModel), pointer :: model
  character (len = iCharLenIE_) :: filename

  allocate(model)
  model = ieModel()

  call model % verbose(2)
  call model % efield_model("weimer")
  call model % aurora_model("fta")
  call model % model_dir("data/ext/")
  call model % init()
  
  call report_errors
  
end program test_01
