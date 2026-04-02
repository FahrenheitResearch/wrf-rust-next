# Build And Run

This file is the stable WRF/WPS workflow.

## Recommended philosophy

- Build CPU WRF first.
- Build WPS second.
- Run a tiny sanity case.
- Only then move to real-data prep and larger domains.

## Official anchor

The closest official workflow is the WRF compile tutorial:

- compile WRF first
- compile WPS second
- for WPS `v4.4+`, use `./configure --build-grib2-libs`

## Community fast path on Ubuntu or WSL

Install the usual build stack first:

```bash
sudo apt update
sudo apt install -y \
  build-essential gcc g++ gfortran make m4 perl csh tcsh flex bison \
  curl file git cmake \
  libnetcdf-dev libnetcdff-dev \
  libopenmpi-dev openmpi-bin \
  zlib1g-dev libpng-dev
```

If that path becomes weird or inconsistent, fall back to the stricter official-library build workflow from the WRF compile tutorial.

## Build order

### WRF

```bash
git clone https://github.com/wrf-model/WRF.git
cd WRF
export NETCDF=/usr
export WRFIO_NCD_LARGE_FILE_SUPPORT=1
./configure
./compile em_real -j "$(nproc)"
```

Guide the user to choose:

- GNU serial for the smallest and simplest starting point
- or GNU dmpar if they already intend to use MPI

Do not rely on fixed menu numbers. They vary.

The success condition is that `main/wrf.exe` and `main/real.exe` exist.

### WPS

```bash
git clone https://github.com/wrf-model/WPS.git
cd WPS
./configure --build-grib2-libs
export WRF_DIR=/path/to/WRF
./compile
```

The success condition is:

- `geogrid.exe`
- `ungrib.exe`
- `metgrid.exe`

## Static geography is part of "WPS works"

Do not leave out the WPS geographical data.

- Download the official WPS static geography package.
- The highest-resolution mandatory package is the normal default and is about `29 GB` uncompressed.
- Point `geog_data_path` in `namelist.wps` at that directory before running `geogrid.exe`.

## Stable run sequence

1. Set up geog data and `namelist.wps`
2. Run `geogrid.exe`
3. Link meteorological GRIB input
4. Choose the correct Vtable
5. Run `ungrib.exe`
6. Run `metgrid.exe`
7. Copy `met_em*` into the WRF run directory
8. Prepare `namelist.input`
9. Run `real.exe`
10. Confirm `wrfinput_d01` and `wrfbdy_d01` exist
11. Run `wrf.exe`

## Sanity checks

After `real.exe`:

- `wrfinput_d01` exists
- `wrfbdy_d01` exists
- start and end dates match the forcing

After `wrf.exe`:

- `rsl.out.0000` ends with success
- CFL is not exploding
- output files are actually being written

## Use wrf-rust to validate output fast

```bash
pip install wrf-rust
python -m wrf info wrfout_d01_YYYY-MM-DD_HH:MM:SS
python -m wrf stats wrfout_d01_YYYY-MM-DD_HH:MM:SS sbcape slp temp
```

If the file opens and core fields look sane, the run is at least structurally alive.
