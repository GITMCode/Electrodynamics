
.PHONY: default ALL LIB cleanall rundir

default : ALL

# This defines the directory structures:
include build/Makefile.dirs

ALL:
	@cd src; make LIB
	@cd src/main; make TEST

cleanall:
	@echo "--> Cleaning Electrodynamics"
	@echo "  --> Current directory : `pwd`"
	@echo "  --> Removing library file in lib:"
	@rm -f lib/*.a
	@echo "  --> Removing files in src:"
	@cd src; make --no-print-directory DIRSFILE=${DIRSFILE} clean
	@echo "  --> Removing files in src/main:"
	@cd src/main; make --no-print-directory DIRSFILE=${DIRSFILE} clean
	@echo "--> Done Cleaning Electrodynamics"

rundir:
	rm -rf run
	mkdir run
	cd run ; ln -s ../src/ie_test.exe ; ln -s ../data .

LIB:
	@cd src; make --no-print-directory DIRSFILE=${DIRSFILE} SHARELIB

clean:	cleanall

