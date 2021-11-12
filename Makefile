include ${FSLCONFDIR}/default.mk

PROJNAME = fabber_cest
LIBS = -lfsl-fabberexec -lfsl-fabbercore -lfsl-newimage \
       -lfsl-miscmaths -lfsl-NewNifti -lfsl-utils \
       -lfsl-cprob -lfsl-znz -ldl
XFILES = fabber_cest
SOFILES = libfsl-fabber_models_cest.so

# Forward models
OBJS =  fwdmodel_cest.o spline_interpolator.o

# For debugging:
#OPTFLAGS = -ggdb

# Pass Git revision details
GIT_SHA1:=$(shell git describe --dirty)
GIT_DATE:=$(shell git log -1 --format=%ad --date=local)
CXXFLAGS += -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_DATE="\"${GIT_DATE}\""

#
# Build
#

all: ${XFILES} ${SOFILES}

# models in a library
libfsl-fabber_models_cest.so : ${OBJS}
	${CXX} ${CXXFLAGS} -shared -o $@ $^ ${LDFLAGS}

# fabber built from the FSL fabbercore library including the models specifieid in this project
fabber_cest : fabber_client.o | libfsl-fabber_models_cest.so
	${CXX} ${CXXFLAGS} -o $@ $< -lfsl-fabber_models_cest ${LDFLAGS}
