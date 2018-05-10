############################################################
#
# $Id$
#
# Copyright (c) INRIA 2012
#
# AUTHOR:
# Gregoire Malandain (gregoire.malandain@inria.fr)
# 
#
#



ETCDIR=../make
include ${ETCDIR}/GNUmakefile.vars # variables par defaut



#
# some variables
#

CREATE_DIRS=TRUE
# DO_NOT_CHECK_DISTANT_OBJS=TRUE
ALWAYS_UPDATE_DEPEND=FALSE



#
# compilation, linking
#

CC=gcc

CPPFLAGS += -D_$(MYCPU)_

# CFLAGS += -fbounds-checking
CFLAGS += -ansi
CFLAGS += -Wall
CFLAGS += -fsigned-char 
CFLAGS += -fsigned-bitfields

# CSTYLE est passe a la compilation et a l'edition de liens
CSTYLE += -fopenmp

OPT=-O

IFLAGS += -I../libbasic/
LDPATH += -L../libbasic/lib/$(INTERNECPU)$(INTERNEOS)
LDLIBS += -lbasic

IFLAGS += -I../libio/
LDPATH += -L../libio/lib/$(INTERNECPU)$(INTERNEOS)
LDLIBS += -lio

LDLIBS += -lm -lz



#
# dependencies
#



#
# files
#

FILES = bal-behavior.c \
	bal-block-tools.c \
	bal-block.c \
	bal-blockmatching-param-tools.c \
	bal-blockmatching-param.c \
	bal-blockmatching.c \
	bal-estimator.c \
	bal-field-tools.c \
	bal-field.c \
	bal-image-tools.c \
	bal-image.c \
	bal-lineartrsf-tools.c \
	bal-lineartrsf.c \
	bal-matrix.c \
	bal-pyramid.c \
	bal-stddef.c \
	bal-transformation-tools.c \
	bal-transformation.c \
	bal-vectorfield.c \
\
	bal-tests.c



#
# targets
#

EXECS = test-copy \
	test-pyramid \
	test-pyramid-image \
	applyTrsf \
	applyTrsfToPoints \
        build-pyramid-image \
	copyTrsf \
	createGrid \
	createTrsf \
	createVectorTrsf \
	cropImage \
	invTrsf \
	printTrsf \
	composeTrsf \
	blockmatching

valgrind : 
	${CC} -g ${CSTYLE} ${IFLAGS} -o ${BINDIR}/valgrind-blockmatching \
	${FILES} blockmatching.c ${LDPATH} ${LDLIBS}

all-valgrind : 
	${CC} -g ${IFLAGS} -o ${BINDIR}/valgrind-test-copy \
	${FILES} test-copy.c ${LDPATH} ${LDLIBS}
	${CC} -g ${IFLAGS} -o ${BINDIR}/valgrind-applyTrsf \
	${FILES} applyTrsf.c ${LDPATH} ${LDLIBS}
	${CC} -g ${IFLAGS} -o ${BINDIR}/valgrind-copyTrsf \
	${FILES} copyTrsf.c ${LDPATH} ${LDLIBS}
	${CC} -g ${IFLAGS} -o ${BINDIR}/valgrind-createGrid \
	${FILES} createGrid.c ${LDPATH} ${LDLIBS}
	${CC} -g ${IFLAGS} -o ${BINDIR}/valgrind-createRandomTrsf \
	${FILES} createRandomTrsf.c ${LDPATH} ${LDLIBS}
	${CC} -g ${IFLAGS} -o ${BINDIR}/valgrind-createVectorTrsf \
	${FILES} createVectorTrsf.c ${LDPATH} ${LDLIBS}
	${CC} -g ${IFLAGS} -o ${BINDIR}/valgrind-invTrsf \
	${FILES} invTrsf.c ${LDPATH} ${LDLIBS}
	${CC} -g ${IFLAGS} -o ${BINDIR}/valgrind-printTrsf \
	${FILES} printTrsf.c ${LDPATH} ${LDLIBS}
	${CC} -g ${IFLAGS} -o ${BINDIR}/valgrind-composeTrsf \
	${FILES} composeTrsf.c ${LDPATH} ${LDLIBS}
	${CC} -g ${IFLAGS} -o ${BINDIR}/valgrind-blockmatching \
	${FILES} blockmatching.c ${LDPATH} ${LDLIBS}


LIB = libblockmatching.a

exec :
	$(MAKE) $(EXECS)

DATE=`date +"%Y-%m-%d_%H-%M-00"`
BAKCOPY=BAK_${DATE}
DEVCOPY=DEV_${DATE}

bak :
	@echo ${BAKCOPY}
	@mkdir -p ${BAKCOPY}
	@cp -p *.c *.h  ${BAKCOPY}/
	@sed -e "s#-I\.\.#-I\.\./\.\.#g" -e "s#-L\.\.#-L\.\./\.\.#g" \
		-e "s#ETCDIR=\.\.#ETCDIR=\.\./\.\.#g" GNUmakefile > ${BAKCOPY}/GNUmakefile
	@echo "done"

dev :
	@echo ${DEVCOPY}
	@mkdir -p ${DEVCOPY}
	@cp -p *.c *.h  ${DEVCOPY}/
	@sed -e "s#-I\.\.#-I\.\./\.\.#g" -e "s#-L\.\.#-L\.\./\.\.#g" \
		-e "s#ETCDIR=\.\.#ETCDIR=\.\./\.\.#g" GNUmakefile > ${DEVCOPY}/GNUmakefile
	@echo "done"



#
# rules
#

include ${ETCDIR}/GNUmakefile.rules 
