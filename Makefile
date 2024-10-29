
default : ALL

include build/Makefile.def

ALL:
	@cd src; make LIB
	@cd src/main; make TEST

allclean:
	rm -f lib/*.a
	@cd src; make clean
	@cd src/main; make clean

rundir:
	rm -rf run
	mkdir run
	cd run ; ln -s ../src/ie_test.exe ; ln -s ../data .


