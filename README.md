
This is a library of codes that allows running different models of the high latitude ionospheric electrodynamics.


Available auroral models:
- Feature Tracking of the Auroral (FTA): Requires AU and AL indices
- Fuller Rowell and Evans [1986]: Requires hemispheric power, which can be derived using the AE index.
- OVATION-Prime: Requires IMF and solar wind.
- OVATION-SME: Requires AE.

Available electric potential models:
- Weimer: Requires IMF and solar wind.
- Heppner-Maynard: Requires IMF and Kp.

Available file formats:
- Binary file from AMIE files.


Products that can be provided:
- High-latitude electric potential (in kV)
- Diffuse auroral electron total energy flux (in mW/m2 or ergs/cm2/s)
- Diffuse auroral electron average energy (in keV)
- Diffuse auroral ion total energy flux (in mW/m2 or ergs/cm2/s)
- Diffuse auroral ion average energy (in keV)
- Monoenergetic auroral electron total energy flux (in mW/m2 or ergs/cm2/s)
- Monoenergetic auroral electron average energy (in keV)
- Wave-driven auroral electron total energy flux (in mW/m2 or ergs/cm2/s)
- Wave-driven auroral electron average energy (in keV)
- Future: high-latitude field-aligned currents


-----

## Configuring & Compiling

Electrodynamics can be run on two modes, either standalone (mostly for reference and
debugging) or coupled to other models.

### Standalone

Run the config script after cloning:

    ./config.sh --compiler=gfortran   # or nagfor

This copies the appropriate compiler rules into `build/Makefile.conf` and
writes `build/Makefile.local`, which the build system requires.  Then build:

    make

To rebuild from scratch:

    make clean
    make

### Coupled to a host model (GITM, TIEGCM, etc.)

The host model's configure step writes `build/Makefile.local` automatically —
no separate step inside this repo is needed.  Just run the host model's
configure script (e.g. `./Config.pl -install` in GITM) and then build the
host model as usual.

### Coupling to a new host model

Three steps in your host model's configure script:

1. Write `build/Makefile.local` pointing to your build config

```make
    BUILDDIR  := /path/to/your/build/dir       # must contain Makefile.conf
    DIRSFILE := /path/to/your/Makefile.dirs    # path/directory variables

    # Or using the SWMF-style build system:
    DIRSFILE := /path/to/your/Makefile.def    # path/directory variables
```

2. Touch `src/Makefile.DEPEND` so make can include it:

        touch /path/to/ext/Electrodynamics/src/Makefile.DEPEND

3. Build the library (from your configure script or top-level Makefile):

        cd /path/to/ext/Electrodynamics && make LIB

### How it works

The build system requires `build/Makefile.local` to exist before compiling.
It contains two variables:

    DIRSFILE := /path/to/Makefile.dirs    # directory/path variables
    BUILDDIR  := /path/to/build/dir       # directory also searched for Makefile.conf

`Makefile.conf` holds compiler flags and suffix rules.

For standalone builds all three files live in this repo's `build/` directory.
For coupled builds, `DIRSFILE` and `BUILDDIR` point into the host model's tree
so that the host model's compiler settings are used.
