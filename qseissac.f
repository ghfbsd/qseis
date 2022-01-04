C     QSEISSAC -- Convert QSEIS reflectivity program output to SAC format, and
C                write SAC seismograms.
C
C     Usage:
C        qseissac [-dir xxxx] [-debug] [-out {ex|ss|ds|clvd|fz|fh|seis}] infile
C            where `infile' is the input file name for qsmain (the QSEIS
C            reflectivity program)
C
C     Use -d[ir] option to write output files in another directory rather than
C        current one.
C
C     Use -in option to read input files in another directory rather than
C        current one.
C
C     Use -out {ex|ss|ds|clvd|fz|fh|seis} to select output type:
C        ex - explosion source-time function
C        ss - strike-slip source-time function
C        ds - dip-slip source-time function
C        clvd - clvd source-time function
C        fz - vertical force source-time function
C        fh - horizontal force source-time function
C        seis - resulting seismogram [default]
C
C     Use -debug to see random debug output.
C
C     G. Helffrich, ELSI / Tokyo Tech
C        last update 31 Dec. 2021

      parameter (irmax=100,nsmax=2**15)
      parameter (mtok=20, mwav=256)
      real rst(irmax),azst(irmax),raz(irmax),tt(irmax)
      real wav(mwav), mtc(6)
      character datfil*64,seifil*32,stnm(irmax)*5,arg*80
      character outdir*32,indir*32,tok(mtok)*16,fmt*16,outnm(7)*8
      character outcod(7)*8,file_name*80, bigl*8192
      real datt(nsmax)
      real dataz(nsmax,irmax),datar(nsmax,irmax),datat(nsmax,irmax)
      integer otyp,outsw(7)
      logical ofake,ospec,odeb
      data outsw/7*0/
      data outcod/'ex','ss','ds','clvd','fz','fh','seis'/
c----------------------------------------------------------------------
      pi=4.*atan(1.)
      datfil = ' '
      outdir = '.'
      indir = '.'
      ofake = .false.
      ospec = .false.
      odeb = .false.
      otyp = 7
      iskip = 0
      do 5 i=1,iargc()
         if (i .le. iskip) go to 5
	 call getarg(i,arg)
	 if (arg .eq. '-spike') then
	    ofake = .true.
	 else if (arg .eq. '-spec') then
	    ospec = .true.
	 else if (arg .eq. '-debug') then
	    odeb = .true.
	 else if (arg .eq. '-out') then
c           seis, ex 
	    call getarg(i+1,arg)
	    iskip = i+1
            do j=1,7
               if (arg.eq.outcod(j)) exit
            enddo
            if (j.gt.7) then
               write(*,*) '**Bad output type: ',arg(1:nblen(arg)),
     &            ', using seis'
               j = 7
            endif
            otyp = j
	 else if (arg(1:2) .eq. '-d') then
	    call getarg(i+1,outdir)
	    iskip = i+1
	 else if (arg .eq. '-in') then
	    call getarg(i+1,indir)
	    iskip = i+1
	 else if (arg(1:1) .eq. '-') then
	    write(0,*) '**Unrecognized option "',
     &         arg(1:index(arg,' ')-1),'", skipping.'
         else
	    datfil = arg
	 endif
5     continue
      if (datfil .eq. ' ')
     &   stop '**No input file given on command line.'
      open(8,file=datfil,status='old',iostat=ios)
      if (ios .ne. 0)
     &   stop '**Invalid input file given on command line.'

c     ----------------------------------    process input file
      lnum = 0
      call in(lnum,arg)
c10.0                    |dble: source_depth;
      read(arg,*,iostat=ios) smin
      call in(lnum,arg)
c0.000                 |dble: receiver_depth;
      read(arg,*,iostat=ios) rmin
      call in(lnum,arg)
c1  1                  |int: sw_equidistant, sw_d_unit;  0 = deg, 1 = km
      read(arg,*,iostat=ios) ieq, iun
      call in(lnum,arg)
c60                    |int: no_distances;
      read(arg,*,iostat=ios) ndist
      call in(lnum,arg)
c5.0 300.0             |dble: d_1,d_n; or d_1,d_2, ...(no comments in between!);
      if (ieq.eq.1) then
         read(arg,*,iostat=ios) rmin,rmax
         dr = (rmax-rmin)/max(1,ndist-1)
         do i=1,ndist
            rst(i) = rmin + (i-1)*dr
         enddo
      else
         i = 0
         do 
            ix = nclen(arg)
            do j=1,ix
               if (arg(j:j).eq.',') arg(j:j) = ' '
            enddo
            call tokens(arg(1:ix),mtok,n,tok)
            do j=1,n
               i = i + 1
               read(tok(j),*,iostat=ios) rst(i)
            enddo
            if (i .ge. ndist) exit
            call in(lnum,arg)
         enddo
      endif
      call in(lnum,arg)
c0.0  8.0  1024        |dble: t_start,t_window; int: no_t_samples;
      read(arg,*,iostat=ios) t0,tl,lt
      dt = tl/(lt-1)
      call in(lnum,arg)
c1  6.5                |int: sw_t_reduce; dble: t_reduce;
      read(arg,*,iostat=ios) ired, tred
      do i=1,ndist
         if (tred.gt.0.0) then
            if (iun.ne.ired) then
               print*,'**Velocity reduction unit mismatch to distance'
               tt(i) = t0
            else if (iun.eq.0) then
               tt(i) = t0 + rst(i)*tred
            else
               tt(i) = t0 + rst(i)/tred
            endif
         else
            tt(i) = t0
         endif
      enddo
      call in(lnum,arg)
c2                                   |int: sw_algorithm;
      read(arg,*,iostat=ios) ialg
      call in(lnum,arg)
c0.000  0.000  0.000  0.000          |dble: slw(1-4);
      read(arg,*,iostat=ios) cl,cn,cm,cs
      call in(lnum,arg)
c1.00                                |dble: sample_rate;
      read(arg,*,iostat=ios) srat
      call in(lnum,arg)
c0.01                                |dble: supp_factor;
      read(arg,*,iostat=ios) sfac
      call in(lnum,arg)
c0                            |int: isur
      call in(lnum,arg)
c0  560.0                     |int: sw_path_filter; dble:shallow_depth_limit;
      call in(lnum,arg)
c0                            |int: no_of_depth_ranges;
      read(arg,*,iostat=ios) nskip
      do i=1,nskip
         call in(lnum,arg)
      enddo
      call in(lnum,arg)
c4.0   1                         |int:dble: wavelet_duration; sw_wavelet;
      read(arg,*,iostat=ios) wdur,iwtp
      if (iwtp.eq.0) then
         read(arg,*,iostat=ios) nwp
         i = 0
         do 
            ix = nclen(arg)
            do j=1,ix
               if (arg(j:j).eq.',') arg(j:j) = ' '
            enddo
            call tokens(arg(1:ix),mtok,n,tok)
            do j=1,n
               i = min(i + 1,mwav)
               read(tok(j),*,iostat=ios) wav(i)
            enddo
            if (i .ge. nwp) exit
            call in(lnum,arg)
         enddo
      endif
c     Skip filter norm factor, zeroes and poles
      call in(lnum,arg)
      do j=1,2
         call in(lnum,arg)
         read(arg,*,iostat=ios) nskip
         do i=1,nskip
            call in(lnum,arg)
         enddo
      enddo
      call in(lnum,arg)
c     Skip output of source components
c  1           1           1          1          1           1         |int
c  'ex'        'ss'        'ds'       'cl'       'fz'        'fh'      |char
      read(arg,*,iostat=ios) (outsw(i),i=1,6)
      call in(lnum,arg)
      read(arg,*,iostat=ios) (outnm(i),i=1,6)
      call in(lnum,arg)
c 1  -0.36e+019 -5.12e+019  5.48e+019 -6.21e+019  2.40e+019 -3.84e+019 'seis'  
      read(arg,*,iostat=ios) imty,(mtc(i),i=1,6),outnm(7)
      outsw(7) = imty
      call in(lnum,arg)
      read(arg,*,iostat=ios) iaz
      call in(lnum,arg)
      if (iaz.eq.0) then
         iaz = 1
         read(arg,*,iostat=ios) azst(1)
         do i=2,ndist
            azst(i) = azst(1)
         enddo
      else
         i = 0
         do 
            ix = nclen(arg)
            do j=1,ix
               if (arg(j:j).eq.',') arg(j:j) = ' '
            enddo
            call tokens(arg(1:ix),mtok,n,tok)
            do j=1,n
               i = i + 1
               read(tok(j),*,iostat=ios) azst(i)
            enddo
            if (i .ge. ndist) exit
            call in(lnum,arg)
         enddo
      endif
      call in(lnum,arg)
c0                               |int: sw_flat_earth_transform;
      call tokens(arg(1:nclen(arg)),mtok,n,tok)
      prad = 6371.
      if (n.gt.1) then
         read(arg,*,iostat=ios) iefa,prad
         if (ios.ne.0) prad = 6371.
      endif
      read(arg,*,iostat=ios) iefa
      degkm = pi*prad/180
      call in(lnum,arg)
c0.25  0.25  5.0                 |dble: vp_res, vs_res, ro_res;  ## skip
      do j=1,2
c7                               |int: no_model_lines;
         call in(lnum,arg)
         read(arg,*,iostat=ios) nskip
c                                skip reading model
         do i=1,nskip
            call in(lnum,arg)
         enddo
      enddo
c     -------------------------- end of input

c     Calculate receiver positions, azimuths and times
      lst = ndist
      write(*,4003) lst
4003  format(' number of stations',i5)
      if (lst.le.10 .or. odeb) then
         write(*,'('' dist--azi (epi to sta)--travel time'')')
      endif
      do i=1,lst
         if (lst.le.10 .or. odeb) write(*,4004) rst(i),azst(i),tt(i)
4004     format(f10.1,2(1x,f9.3))
      enddo

c     ---------------------------------

c     Write out source-time function time series if desired
      if (odeb .and. iwtp.eq.0) then
         if (nwp.gt.mwav) write(0,*) '**Wavelet too long, truncated.'
	 call newhdr
	 call wsac1('qs6wav.z',wav,min(nwp,mwav),0,dt,nerr)
	 if (nerr .ne. 0)
     &      write(0,*) '**Warning: Error writing source wavelet.'
      endif

c     Extract data from output file and process

      if (lt.gt.nsmax) stop '**Too many samples, recompile.'
      if (lst.gt.irmax) stop '**Too many receivers, recompile.'
      if (outsw(otyp).eq.0) then
         write(*,*) '**QSEIS didn''t generate ' // outnm(otyp)
         stop
      endif

      seifil = indir(1:nblen(indir)) // '/' // outnm(otyp)
      datfil = seifil(1:nblen(seifil)) // '.tz'
      open(9,file=datfil,status='old')
      datfil = seifil(1:nblen(seifil)) // '.tr'
      open(10,file=datfil,status='old')
      datfil = seifil(1:nblen(seifil)) // '.tt'
      open(11,file=datfil,status='old',iostat=ios)
      if (ios.ne.0) then
         ic = 1
      else
         ic = 2
      endif
      do i=0,ic
         read(9+i,'(a)',iostat=ios) bigl
         call tokens(bigl,mtok,n,tok)
         if (n .ne. lst+1) then
            write(0,'(a,i4,a,i4,a)') '**Data-input mismatch:',
     &         lst,' traces expected, ',n,' found.'
            stop
         endif
      enddo
      write(fmt,'(a,i4,a)') '(',1+lst,'f12.0)'
      if (fmt(2:4).eq.' ') fmt(2:) = fmt(5:)
      if (fmt(2:3).eq.' ') fmt(2:) = fmt(4:)
      if (fmt(2:2).eq.' ') fmt(2:) = fmt(3:)
      if (odeb) write(*,*) '..output format ',fmt(1:nblen(fmt))
      do i=1,lt
         read(9, fmt,iostat=ios) datt(i),(dataz(i,j),j=1,lst)
         read(10,fmt,iostat=ios) datt(i),(datar(i,j),j=1,lst)
         read(11,fmt,iostat=ios) datt(i),(datat(i,j),j=1,lst)
      enddo
      if (abs(1 - (datt(2)-datt(1))/dt) .gt. 1e-4)
     &   write(0,*) '**FFT dT versus sample dT mismatch -- check!',
     &   datt(2)-datt(1),dt, abs(1 - (datt(2)-datt(1))/dt)
      close(9)
      close(10)
      if(ic.gt.1)close(11)

c     Prepare to write seismogram
      dmin=smin
      do i=1,lst
         evla=0.001
         evlo=0.001
         write(stnm(i),'(i4.4)') i
         write(file_name,'(a,1h/,a,1h.,a)')
     *      outdir(1:nblen(outdir)),stnm(i)(1:4),outnm(otyp)
         kif = index(file_name,' ')
         if (iun.eq.0) then
            gcarc = rst(i)
            distkm = acos(cos(gcarc*pi/180))*180/pi*degkm
         else
            distkm = rst(i)
            gcarc = rst(i)/degkm
         endif
         if (gcarc.eq.0) then
            gcnew = 0.0
         else
C     Use Newton-Raphson to locate exact elliptical earth distance to yield
C        specified distance in kilometers.  Only two iterations are usually
C        required to find distance within 0.1%.
            do j=1,5
               gclo=0.95*gcarc
               gchi=1.05*gcarc
               call dazll(evla,evlo,gcarc,azst(i),stla,stlo)
               call gcd(evla,evlo,stla,stlo,delta,delkm,saz,baz)
               call dazll(evla,evlo,gclo,azst(i),stla,stlo)
               call gcd(evla,evlo,stla,stlo,delta,dello,saz,baz)
               call dazll(evla,evlo,gchi,azst(i),stla,stlo)
               call gcd(evla,evlo,stla,stlo,delta,delhi,saz,baz)
               ddddel = (delhi-dello)/(gchi-gclo)
               gcnew = gcarc - (delkm-distkm)/ddddel
               if (abs(gcnew-gcarc) .lt. 1e-3) exit
               gcarc = gcnew
            enddo
            if (j.gt.5) write(0,*) '**Odd - could not find GC distance.'
         endif

         call dazll(evla,evlo,gcnew,azst(i),stla,stlo)
         call gcd(evla,evlo,stla,stlo,delta,delkm,saz,baz)
         call newhdr
         call setnhv('npts',lt,nerr)
         call setfhv('b',tt(i),nerr)
         call setfhv('o',0.0,nerr)
         call setfhv('delta',dt,nerr)
         call setfhv('gcarc',delta,nerr)
         call setfhv('evla',evla,nerr)
         call setfhv('evlo',evlo,nerr)
         call setfhv('evdp',dmin,nerr)
         call setfhv('az',saz,nerr)
         call setfhv('baz',baz,nerr)
         if (baz .gt. 180.0) then
            raz = baz - 180.0
         else
            raz = baz + 180.0
         endif
         call setfhv('cmpaz',raz(i),nerr)
         call setkhv('kstnm',stnm(i),nerr)
         call setfhv('stla',stla,nerr)
         call setfhv('stlo',stlo,nerr)
         call setnhv('nzyear',1999,nerr)
         call setnhv('nzjday',365,nerr)
         call setnhv('nzhour',23,nerr)
         call setnhv('nzmin',59,nerr)
         call setnhv('nzsec',59,nerr)
         call setnhv('nzmsec',999,nerr)
C        P-SV seismograms
         call setfhv('b',tt(i),nerr)
         call setfhv('cmpinc',0.0,nerr)
         call setkhv('kcmpnm','VERT',nerr)
         file_name(kif:) = '.z'
         call wsac0(file_name,dataz(1,i),dataz(1,i),nerr)
         if (nerr .ne. 0) write(0,*) '**Trouble writing ',
     &         file_name(1:nblen(file_name)),'.'
         call setfhv('b',tt(i),nerr)
         call setfhv('cmpinc',90.0,nerr)
         call setkhv('kcmpnm','RAD',nerr)
         file_name(kif:) = '.r'
         call wsac0(file_name,datar(1,i),datar(1,i),nerr)
C        SH seismograms
         call setfhv('b',tt(i),nerr)
         if (raz(i) .gt. 270.0) then
            raz(i) = raz(i) - 270.0
         else
            raz(i) = raz(i) + 90.0
         endif
         call setfhv('cmpaz',raz(i),nerr)
         call setfhv('cmpinc',90.0,nerr)
         call setkhv('kcmpnm','TAN',nerr)
         file_name(kif:) = '.t'
         call wsac0(file_name,datat(1,i),datat(1,i),nerr)
         if (nerr .ne. 0) write(0,*) '**Trouble writing ',
     &      file_name(1:nblen(file_name)),'.'
      enddo
      stop
      end

      subroutine in(lnum,line)
      character line*(*)
      do
         lnum = lnum + 1
         read(8,'(a)',iostat=ios) line
         if (ios.ne.0)
     *      write(*,'(a,i4)') '**Early end of input, line ',lnum
         if (line(1:1).ne.'#' .or. ios.ne.0) exit
      enddo
c     print*,line(1:nblen(line))
      end

      integer function nblen(string)
      character*(*) string

      do i=len(string),2,-1
         if (string(i:i) .ne. ' ') exit
      enddo
      nblen = i
      end

      integer function nclen(string)
c     Returns length of non-comment part of string:
c     | or # introduce a comment
      character*(*) string

      ib = index(string,'|')
      is = index(string,'#')
      if (ib.eq.0 .and. is.eq.0) then
         nclen = len(string)
      else if (ib.gt.0 .and. is.gt.0) then
         nclen = min(ib,is)-1
      else
         nclen = max(ib,is)-1
      endif
      end
