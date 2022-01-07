FC=gfortran
FFLAGS=-g -fbounds-check -Wunused
FFLAGS=-O2

LSAC = /usr/local/lib/libsacio-64.a

OBJ= bessj.o bessj0.o bessj1.o caxcb.o cdgemp.o cmemcpy.o four1.o \
	getdata.o qsam2ve.o qsbsj.o qsfftinv.o qsgetinp.o qshksh.o qskern.o \
	qslayer.o qsmain.o qsmultis.o qspsv.o qsqmodel.o qssh.o qssource.o \
	qssublay.o qsve2am.o qswavelet.o qswaveno.o qswvint.o taper.o wavelet.o

qsam2ve.o qsbsj.o qsfftinv.o qsgetinp.o qshksh.o qskern.o qslayer.o \
	qsmain.o qsmultis.o qspsv.o qsqmodel.o qssh.o qssource.o qssublay.o \
	qsve2am.o qswavelet.o qswaveno.o qswvint.o: qsglobal.h

qsmain: $(OBJ)
	$(FC) -o qsmain $(OBJ)

qseissac: qseissac.o tokens.o gcdist.o
	$(FC) -o qseissac qseissac.o tokens.o gcdist.o ${LSAC}

clean: ; /bin/rm -f *.o

distclean: ; /bin/rm -f *.o qsmain qseissac *.t[rtzv] *.{seis,mex}{,-2}.[rtz]
