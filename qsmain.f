c     Usage:
c        qsmain [ -prog | -noprog ] <input>
c
c     where <input> is the name of an input file.  This must be present on the
c     command line.
c
c     Use -prog to show progress bar during calculation [default]
c
c     Use -noprog to suppress progress bar during calculation

      program qseis
      implicit none
c
      include 'qsglobal.h'
c
c     work space
c
      integer i,ji,jk,istp,nssel,runtime,nblen
      double precision pi,srate
      logical grnexist
      integer time
      character arg*80
c
c     read input file file
c
      print *,'######################################################'
      print *,'#                                                    #'
      print *,'#               Welcome to the program               #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#        QQQ     SSSS    EEEEE    III     SSSS       #'
      print *,'#       Q   Q   S        E         I     S           #'
      print *,'#       Q Q Q    SSS     EEEE      I      SSS        #'
      print *,'#       Q  QQ       S    E         I         S       #'
      print *,'#        QQQQ   SSSS     EEEEE    III    SSSS        #'
      print *,'#                                                    #'
      print *,'#                  (Version 2006)                    #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#                      by                            #'
      print *,'#                 Rongjiang Wang                     #'
      print *,'#              (wang@gfz-potsdam.de)                 #'
      print *,'#                                                    #'
      print *,'#           GeoForschungsZentrum Potsdam             #'
      print *,'#           Last modified: October, 2006             #'
      print *,'######################################################'
      print *,'                          '
      ji = 0
      inputfile = ' '
      oprog = .true.
      do i=1,iargc()
         if (i.le.ji) cycle
         call getarg(i,arg)
         if (arg(1:1) .ne. '-') then
            inputfile = arg
         elseif (arg .eq. '-noprog') then
            oprog = .false.
         elseif (arg .eq. '-prog') then
            oprog = .true.
         else
            jk = nblen(arg)
            write(*,'(3a)') '**Invalid option: ',arg(1:jk),', skipped'
         endif
      enddo
      if (inputfile.eq.' ') stop '**Missing input file'
      i = nblen(inputfile)
      write(*,'(a,a)')' The input data file is ',inputfile(1:i)
      runtime=time()
c
      pi=4.d0*datan(1.d0)
c
      open(10,file=inputfile,status='old')
      call qsgetinp(10,srate,nssel)
      close(10)
c
      if(nssel-ssel(7).gt.0)then
        grnexist=.false.
        call qswvint(srate)
        call qsmultis(grnexist)
        iexist=0
        do istp=1,7
          do i=1,4
            if(fsel(i,istp).eq.1)then
              call qsfftinv(i,istp)
            endif
          enddo
        enddo
      else
        grnexist=.true.
        call qsmultis(grnexist)
      endif
c
      runtime=time()-runtime
      write(*,'(a)')' #############################################'
      write(*,'(a)')' #                                           #'
      write(*,'(a)')' #      End of computations with qseis06     #'
      write(*,'(a)')' #                                           #'
      write(*,'(a,i10,a)')' #       Run time: ',runtime,
     +                                           ' sec            #'
      write(*,'(a)')' #############################################'
1001  format(2i7,E12.4,a)
1002  format(i4,a,E12.4,a,$)
1003  format(E12.5,$)
1004  format(2E12.4,$)
1005  format(2E12.4)
 500  stop
      end

      integer function nblen(string)
      character*(*) string

      do i=len(string),2,-1
         if (string(i:i) .ne. ' ') exit
      enddo
      nblen = i
      end
