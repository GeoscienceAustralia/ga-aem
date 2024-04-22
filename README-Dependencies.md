# GA-AEM third-party library dependencies

## Description
For full functionality and to build all programs the following packages are required: FFTW, MPI, NetCDF, GDAL and PETSc.

### Installing dependencies on Linux
Depending on your flavour of Linux, you are likely to need to install the packages (or later versions) for Ubuntu as follows:

```bash
> sudo apt install cmake
> sudo apt install gcc-9
> sudo apt install libfftw3-dev
> sudo apt install libopenmpi-dev
> sudo apt install libnetcdf-dev
> sudo apt install libnetcdf-c++4-dev
> sudo apt install libgdal-dev
> sudo apt install libpetsc-real-dev
```

### Installing dependencies on Windows
On Windows the packages can be downloaded and installed as outlined below.
- Note that 64 bit versions are required with the include and library files.
- Note that the directory path `%LocalAppData%` is the Windows environment variable where non-administrator users can install applications for their own use and is typically located at `C:\Users\<your-username>\AppData\Local\`.
- Note that `%LocalAppData%` can be substituted for any other suitable installation path, for example `C:\Program Files\` or `C:\Win10Dev\`.

#### FFTW
- Download the 64 bit precompiled FFTW 3.3.5 Windows DLLs - https://fftw.org/pub/fftw/fftw-3.3.5-dll64.zip.
- Unzip these to `%LocalAppData%\fftw-3.3.5-dll64`.

#### MPI
- Go to the Microsoft MPI v10.0 webpage - https://www.microsoft.com/en-us/download/details.aspx?id=57467.
- Follow the install instructions.
- For building the programs you will need to install the SDK - 'msmpisdk.msi'
- For running the built programs you will need to install the runtime `MSMpiSetup.exe`

#### NETCDF
- Download the netCDF-4 64-bit prebuilt Windows binaries installer from - https://downloads.unidata.ucar.edu/netcdf-c/4.9.2/netCDF4.9.2-NC4-64.exe.
- This can be installed to the default path `C:\Program Files\netCDF 4.9.2`.
- Earlier version of netCDF 4 will also work.
- The netCDF C++ bindings library source code is included in this repository as a submodule and does NOT need to be downloaded separately.

#### GDAL
- From the GIS Internals GDAL v3.7.1 MapServer v8-0-1 (MSVC 2019 64bit) pre-built packages webpage - https://gisinternals.com/query2.html?content=filelist&file=release-1928-x64-gdal-3-7-1-mapserver-8-0-1.zip.
- Download the - Compiled binaries - https://build2.gisinternals.com/sdk/downloads/release-1928-x64-gdal-3-7-1-mapserver-8-0-1.zip.
- Download the - Compiled libraries and headers - https://build2.gisinternals.com/sdk/downloads/release-1911-x64-gdal-3-0-4-mapserver-7-4-3-libs.zip.
- Unzip both these into a single directory e.g. `gdal-3.7.1-mapserver-8-0-1` that should then have bin, doc, include and lib sub-directories.
- Place that folder into a suitable install path, for example `%LocalAppData%\gdal-3.7.1-mapserver-8-0-1`.
- Note that earlier versions of GDAL can also be similarly used, for example, v3.0.4 - https://gisinternals.com/query2.html?content=filelist&file=release-1911-x64-gdal-3-0-4-mapserver-7-4-3.zip.

#### PETSc
- Download the prebuilt PETSc 3.9.4 Windows 64 bit libraries - https://github.com/rcb547/petsc-3.9.4-vs2017-win64-libraries/releases/download/v1.0.0/petsc-3.9.4-vs2017-win64-libraries.zip
- Unzip the zip file to, for example, `%LocalAppData%` such that you then have a directory named `%LocalAppData%\petsc\3.9.4\vs2017`.
