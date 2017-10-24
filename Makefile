include ${FSLCONFDIR}/default.mk

PROJNAME = fabber_cest

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_BOOST} -I..
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L../fabber_core

LIBS = -lutils -lnewimage -lmiscmaths -lprob -lnewmat -lfslio -lniftiio -lznz -lz -ldl

XFILES = fabber_cest

# Forward models
OBJS =  fwdmodel_cest.o spline_interpolator.o

# For debugging:
OPTFLAGS = -ggdb
#OPTFLAGS =

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
