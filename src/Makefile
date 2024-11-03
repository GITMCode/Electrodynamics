
default : LIB

include ../build/Makefile.def

-include ../../../Makefile.def
-include Makefile.DEPEND

MYLIB = ../lib/libIE.a

OBJECTS = \
	ModKind.o\
	ModNumConst.o\
	ModConst.o\
	ModCharSize.o\
	ModErrors.o\
	ModTimeAmie.o\
	ModWeimer.o\
	ModWeimer05.o\
	hmr89.o\
	ModFtaModel.o\
	ModIHP.o\
	ModAMIE_Interface.o\
	ModIE.o\
	tools.o\
	interpret_names.o

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl -I${SHAREDIR} ${OBJECTS}

LIB:
	make ${MYLIB}

${MYLIB}: ${OBJECTS}
	rm -rf ${MYLIB}
	${AR} ${MYLIB} ${OBJECTS}

SHARELIB: DEPEND
	make ${MYLIB}
	@ln -s ../lib/libIE.a ${LIBDIR}/libIE.a