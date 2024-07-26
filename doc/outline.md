
The general idea of this library is to make a variety of ionospheric electrodynamics codes available to the user.  There are three types of methods that will be available:

- Empirical models, including a variety of electric potential and auroral precipitation models.
- Files, including custom binary, netCDF, and HDF5 formats.
- Coupling to codes, including to the Space Weather Modeling Framework.


The library works in the following order:
- Initialize: set what type of model is desired.  Provide filenames or empirical models to use. Provide the number of latitudes and local times to use.
- Set Indices: for empirical models, the library needs the indices to use.
- Set Time: tell the library the current time.
- Set grid: tell the library the magnetic latitudes and magnetic local times to get results on.
- Get result: get the specific result (potential, energy flux, etc.) on the specific grid at the specific time.


