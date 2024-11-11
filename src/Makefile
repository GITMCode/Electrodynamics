
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
        ModIoUnit.o\
	ModAMIE_Interface.o\
	ModIE.o\
	tools.o\
	ED_ModInterface.o\
	ED_Interface.o\
        ED_r_elec_ed_lup_subs.o\
	interpret_names.o

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl  ${SEARCH} ${LIBDIR} ${OBJECTS}

LIB:
	make ${MYLIB}

${MYLIB}: ${OBJECTS}
	@echo
	@echo 'creating LibIE'
	@echo
	${AR} ${MYLIB} ${OBJECTS}

SHARELIB: DEPEND
	make ${MYLIB}
	@cp ../lib/libIE.a ${LIBDIR}/libIE.a
