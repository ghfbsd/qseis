      character lab*16
      lab = '5 Hz'
      do i=1,50
         call spin(i,50,lab(1:4))
         call sleep(1)
      enddo
      end

      subroutine spin(ino,imx,lab)
      character lab*(*)
      parameter (llen=80, lst=4)
      character str*(lst), line*(llen)
      data icnt/0/, str/'|/-\'/

      il = llen - len(lab)
      is = float(ino)/float(imx) * il
      id = il - is + 1
      ic = mod(icnt,4) + 1
      do i=1,il
         if (i.le.is) line(i:i) = '*'
         if (i.eq.is+1) line(i:i) = str(ic:ic)
         if (i.gt.is+1) line(i:i) = '.'
      enddo
      line(il+1:) = lab
      icnt = icnt + 1
      write(*,'(a,a,$)') line,char(13)
      if (ino.ge.imx) write(*,*)
      end
