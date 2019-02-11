include ${FSLCONFDIR}/default.mk

PROJNAME = fabber_cest

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_BOOST} -I..
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L../fabber_core

FSLVERSION= $(shell cat ${FSLDIR}/etc/fslversion | head -c 1)
ifeq ($(FSLVERSION), 5) 
  NIFTILIB = -lfslio -lniftiio 
else 
  NIFTILIB = -lNewNifti
endif

LIBS = -lutils -lnewimage -lmiscmaths -lprob -lnewmat ${NIFTILIB} -lznz -lz -ldl

XFILES = fabber_cest

# Forward models
OBJS =  fwdmodel_cest.o spline_interpolator.o

# For debugging:
OPTFLAGS = -ggdb
#OPTFLAGS =

# Pass Git revision details
GIT_SHA1:=$(shell git describe --dirty)
GIT_DATE:=$(shell git log -1 --format=%ad --date=local)
CXXFLAGS += -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_DATE="\"${GIT_DATE}\""

#
# Build
#

all:	${XFILES} libfabber_models_cest.a

# models in a library
libfabber_models_cest.a : ${OBJS}
	${AR} -r $@ ${OBJS}

# fabber built from the FSL fabbercore library including the models specifieid in this project
fabber_cest : fabber_client.o ${OBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ $< ${OBJS} -lfabbercore -lfabberexec ${LIBS}

# DO NOT DELETE
