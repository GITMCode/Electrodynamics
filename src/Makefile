
default : LIB

include ../build/Makefile.def

MYLIB = ../lib/libIE.a

OBJECTS = \
	ModCharSize.o\
	ModErrors.o\
	ModWeimer05.o\
	ModFtaModel.o\
	ModIE.o\
	tools.o\
	efield_interpret_name.o

LIB:
	make ${MYLIB}

${MYLIB}: ${OBJECTS}
	rm -rf ${MYLIB}
	${AR} ${MYLIB} ${OBJECTS}


