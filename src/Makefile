
default : LIB

#include ../build/Makefile.def

include ../../../Makefile.def
include ../../../Makefile.conf

include Makefile.DEPEND

MYLIB = ../lib/libIE.a

OBJECTS = \
	ModKind.o\
	ModNumConst.o\
	ModConst.o\
	ModCharSize.o\
	ModErrors.o\
	ModNewell.o\
	ModTimeAmie.o\
	ModWeimer.o\
	ModWeimer05.o\
	hmr89.o\
	ModFtaModel.o\
	ModIHP.o\
	ModIoUnit.o\
	ModAMIE_Interface.o\
	ModIE.o\
	tools.o\
	interpret_names.o

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl  ${SEARCH} ${LIBDIR} ${OBJECTS}

LIB:
	@make ${MYLIB}
	@echo

${MYLIB}: ${OBJECTS}
	@echo
	@echo 'creating LibIE'
	@echo
	@rm -f ${MYLIB}
	${AR} ${MYLIB} ${OBJECTS}

SHARELIB: DEPEND
	@make ${MYLIB}
