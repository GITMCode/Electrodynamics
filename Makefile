
default : ALL

include build/Makefile.def

ALL:
	@cd src; make LIB
	@cd src/main; make TEST

cleanall:
	@echo "--> Cleaning Electrodynamics"
	rm -f lib/*.a
	cd src; make clean
	cd src/main; make clean
	@echo "--> Done Cleaning Electrodynamics"

rundir:
	rm -rf run
	mkdir run
	cd run ; ln -s ../src/ie_test.exe ; ln -s ../data .

LIB:
	@cd src; make SHARELIB
