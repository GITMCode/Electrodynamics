
default : TEST

include ../../build/Makefile.def

SEARCH = -I..

EXE = ../ie_test.exe

OBJECTS = \
	main.o

${EXE}: ${OBJECTS}
	${LINK.f90} -o ${EXE} ${OBJECTS} -L../../lib -lIE 

TEST:
	@make ${EXE}

