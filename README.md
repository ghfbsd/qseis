# QSEIS

This is a slightly modified version of R. Wang's QSEIS (2006 version).  See
[this url](https://www.gfz-potsdam.de/en/section/physics-of-earthquakes-and-volcanoes/infrastructure/tool-development-lab/) 
for information and the most up-to-date version.
The programs are all written in Fortran 77.

The modifications include:

- Added Unix Makefile to build programs.

- Changed input to read input file name from the command line.

- Changed input so that synthetics for other planets may be made (e.g. Mars).

- Merged some bug fixes with QSEIS source code from Nov. 2006 version.

- Output receiver positions in their units of input (km or deg).

- Added program that will write SAC output files.

First make the programs:
```
make qsmain qseissac
```

See these annotated test input files included with the source code for details
on preparing input.

- **Simple test** (regional crustal propagation synthetics).
  There will be 60 seismograms written at 5 km offsets out to 300 km distance.
  ```

   qsmain qs6testinput.dat          ## runs synthetics
   qseissac qs6testinput.dat        ## writes SAC files

  # After this, examine one of the files with SAC:
   sac
   r 0001.seis.z    ;* reads data file
   lh               ;* lists header information
   p1               ;* plots seismogram
   quit
  ```

- **Complex test** (teleseismic propagation, synthetics take longer).
  There will be 31 seismograms written between 60 and 90 degrees distance.
  ```
   qsmain qs6inp.dat                ## runs synthetics
   qseissac qs6inp.dat              ## writes SAC files
   sac
   r 0001.seis-2.z  ;* reads data file
   m ttsac          ;* plots seismogram, marks phase arrivals
   r 0011.seis-2.z  ;* reads data file
   m ttsac          ;* plots seismogram, marks phase arrivals
   r 0021.seis-2.z  ;* reads data file
   m ttsac          ;* plots seismogram, marks phase arrivals
   r 0031.seis-2.z  ;* reads data file
   m ttsac          ;* plots seismogram, marks phase arrivals
   quit
  ```
  See [this url](https://members.elsi.jp/~george/sac-bugs.html#ttimes) and
  [this url](https://github.com/ghfbsd/sacbook/tree/master/methods/split) and
  for making the ttsac macro available to your version of SAC.

- **Mars test** (vertical incidence, ambient noise autocorrelation synthetic).
  There will be 8 seismograms written at offsets from 0-5 km, 10 and 50 km.
  Explosive source at 60 km depth.
  ```
   qsmain qs6testmars.dat                ## runs synthetics
   qseissac qs6testmars.dat              ## writes SAC files
   sac
   r D*K.mex.z      ;* reads data files
   xlim 250 400     ;* window arrival of 410 km disc. analog and PcP
   p1               ;* plots seismograms
   quit
  ```
