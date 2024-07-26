
program test_01

  use ModCharSize
  use ModErrors
  use ModIE
  implicit none

  type(ieModel), pointer :: model
  character (len = iCharLenIE_) :: filename

  allocate(model)
  model = ieModel()

  call model % verbose(10)
  call model % efield_model("weimer")
  call model % aurora_model("fta")
  call model % model_dir("data/ext/")
  call model % init()

  call model % imfBz(-5.0)
  call model % useAeHp()
  call model % au(100.0)
  call model % al(-500.0)
  
  call report_errors
  
end program test_01
