
default : ALL

include build/Makefile.def

ALL:
	@cd src; make LIB
	@cd src/main; make TEST

clean:
	rm -f lib/*.a
	@cd src; make clean
	@cd src/main; make clean
