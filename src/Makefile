
default : LIB

include ../build/Makefile.def

-include ../../../Makefile.def
-include ../../../Makefile.conf
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
	@perl ${SCRIPTDIR}/depend.pl  ${SEARCH} ${LIBDIR} ${OBJECTS}

LIB:
	make ${MYLIB}

${MYLIB}: ${MF} ${OBJECTS}
	@echo
	@echo 'creating LibIE'
	@echo
	${AR} ${MYLIB} ${OBJECTS} ${MODULES}

SHARELIB: DEPEND
	make ${MYLIB}
	@cp ../lib/libIE.a ${LIBDIR}/libIE.a
