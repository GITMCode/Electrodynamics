
default : LIB

include ../build/Makefile.def

MYLIB = ../lib/libIE.a

OBJECTS = \
	ModKind.o\
	ModNumConst.o\
	ModConst.o\
	ModCharSize.o\
	ModErrors.o\
	ModTimeConvert.o\
	ModWeimer.o\
	ModWeimer05.o\
	hmr89.o\
	ModFtaModel.o\
	ModIHP.o\
	ModIE.o\
	tools.o\
	efield_interpret_name.o

LIB:
	make ${MYLIB}

${MYLIB}: ${OBJECTS}
	rm -rf ${MYLIB}
	${AR} ${MYLIB} ${OBJECTS}


