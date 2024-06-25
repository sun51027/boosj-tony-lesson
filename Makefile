#
#  Makefile 
#  000725 by Stefan Haberhauer NEC-ESS 
#

# Suffixes rule
.SUFFIXES :
.SUFFIXES : .o .f 

# F90 compiler
#F90 = sxf90
#F90 = f90
#F90=ifort
F90=gfortran
FC=${F90}
# options to invalidate environment variable for F90
#CLEARFLAGS = -clear
CLEARFLAGS = 

# F90 user compile options
#FFLOAT = -float0 -ew -sx4
#FFLOAT = -float0 -ew -sx5
#float2 FFLOAT = -float2 -sx4
#FHOPT = -Chopt
#FBASE = -P stack -f0 -eb \
#        -Wf,-pvctl noassume loopcnt=100000 vr256 vwork=stack vworksz=4M noverrchk novlchk matmul
#FFLAGS = ${FFLOAT} ${FBASE} 
#FFLAGS = -w -r8 -i4 -O3 -xT
FFLAGS = -fdefault-real-8 -O3
# linker
LD = ${F90}

# link options
#LDFLAGS = ${FFLOAT}
#LDFLAGS = ${FFLAGS}
#LDFLAGS = ${FFLAGS} -Wl,-framework -Wl,accelerate
# link library
#LIBRARY = -leispack_64 -lblas_64
##float2 LIBRARY = -leispack -lblas

EXEC = tprbal.x

MODULES = tprbap_modules.f

FSRCS = bophys_ap.f \
	bootps.f \
	cospol.f \
	eqinvm_ap.f \
	extint_ap.f \
	intext.f \
	lamcal_ap.f \
	lamnew.f \
	lgikvm_ap.f \
	materp_ap.f \
	mercap.f \
	metric_ap.f \
	mtaskb_ap.f \
	mtaskl_ap.f \
	ploteq.f \
	sgedi.f  \
	sgefa.f  \
	sgesl.f  \
	tprall_bap.f \
	trgfun.f \
	veqrec_ap.f \
	vforce_ap.f \
	vmtobo_ap.f \
	xminv.f 

OBJS = 	bophys_ap.o \
	bootps.o \
	cospol.o \
	eqinvm_ap.o \
	extint_ap.o \
	intext.o \
	lamcal_ap.o \
	lamnew.o \
	lgikvm_ap.o \
	materp_ap.o \
	mercap.o \
	metric_ap.o \
	mtaskb_ap.o \
	mtaskl_ap.o \
	ploteq.o \
	sgedi.o  \
	sgefa.o  \
	sgesl.o  \
	tprall_bap.o \
	tprbap_modules.o \
	trgfun.o \
	veqrec_ap.o \
	vforce_ap.o \
	vmtobo_ap.o \
        functions90_it.o \
	xminv.o 

${EXEC} : ${OBJS} 
	${LD} ${CLEARFLAGS} -o ${EXEC} ${LDFLAGS} ${OBJS} ${LIBRARY}
#
# General rule
.o.o :
	${F90} ${CLEARFLAGS} -c ${FINL} ${FFLAGS} $<

# High optimization
#lamnew.o: lamnew.o
# 	${F90} ${CLEARFLAGS} -c ${FHOPT} ${FFLAGS} $<
#
#sgefa.o: sgefa.o
#	${F90} ${CLEARFLAGS} -c ${FHOPT} ${FFLAGS} $<
# 
materp_ap.o: materp_ap.f tprbap_modules.o
	${F90} ${CLEARFLAGS} -c ${FHOPT} ${FFLAGS} $<

eqinvm_ap.o: eqinvm_ap.f tprbap_modules.o
	${F90} ${CLEARFLAGS} -c ${FHOPT} ${FFLAGS} $<

veqrec_ap.o: veqrec_ap.f tprbap_modules.o
	${F90} ${CLEARFLAGS} -c ${FHOPT} ${FFLAGS} $<

ploteq.o: ploteq.f tprbap_modules.o
	${F90} ${CLEARFLAGS} -c ${FHOPT} ${FFLAGS} $<

vmtobo_ap.o: vmtobo_ap.f 
	${F90} ${CLEARFLAGS} -c ${FHOPT} ${FFLAGS} $<

mtaskb_ap.o: mtaskb_ap.f tprbap_modules.o
	${F90} ${CLEARFLAGS} -c ${FFLAGS} $<

tprall_bap.o: tprall_bap.f tprbap_modules.o
	${F90} ${CLEARFLAGS} -c ${FFLAGS} $<

newcase:
	@rm eqinvm_ap.o materp_ap.o mtaskb_ap.o ploteq.o tprall_bap.o veqrec_ap.o
	make

clean :
	@rm -f ${EXEC} ${OBJS} *.L


# Dependencies of FORTRAN source files
#eqinvm_ap.f: parbal.inc tprcom.bal
#materp_ap.f: parbal.inc tprcom.bal
#mtaskb_ap.f: parbal.inc tprcom.bal
#ploteq.f: parbal.inc tprcom.bal
#veqrec_ap.f: parbal.inc tprcom.bal
