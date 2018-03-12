# off-tree build into $(DESTDIR)
DESTDIR  = build/

#Uncomment below to specify default options
MPICC       =   mpicc
MPICXX      = mpic++

#
# Optimized defaults for icc
#OPTIMIZE =  -fopenmp -O3 -ipo -g -Wall -xHost -no-prec-div -fp-model fast -wd11021
#MS try1
#OPTIMIZE =  -fopenmp -O3 -g -Wall -xHost -no-prec-div -fp-model fast -wd11021
# MS try2
OPTIMIZE =  -fopenmp -O3 -g -Wall -ffast-math -march=native

#GSL_INCL = $(shell pkg-config --cflags gsl)
#We don't want to add -lm here on icc
#GSL_LIBS = $(filter-out -lm,$(shell pkg-config --libs gsl))
GSL_INCL = -I/usr/local/include/gsl
GSL_LIBS = -L/usr/local/gsl-2.4/lib -lgsl -lgslcblas


# MS
#-DPETAPM_ORDER=1
#-DTOPNODEFACTOR=5.0


#--------------------------------------- Basic operation mode of code
OPT += -DDENSITY_INDEPENDENT_SPH
#OPT += -DLIGHTCONE                       # write a lightcone on the fly; in development
#OPT += VALGRIND  # allow debugging with valgrind, disable the GADGET memory allocator.

# flags shall that always be there they need to be cleaned up
OPT += -DOPENMP_USE_SPINLOCK
OPT += -DSPH_GRAD_RHO  # calculate grad of rho in SPH, required for Krumholtz & Gnedin H2 SFR

#--------------------------------------- SFR/feedback model
# Star formation master switch. Also enables the Wind model
OPT	+=  -DSFR

#-------------------------------------- AGN stuff
OPT	+=  -DBLACK_HOLES             # enables Black-Holes (master switch)

#-------------------------------------------- Things for special behaviour
OPT	+=  -DNO_ISEND_IRECV_IN_DOMAIN     #sparse MPI_Alltoallv do not use ISEND IRECV
