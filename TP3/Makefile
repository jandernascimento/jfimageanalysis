#######################################
# Makefile PBM                        #
#                                     #
# E.B.                                #
#######################################


#PROG = ppmtopgm pxmtopxm pgmtopgm
PROG = imagetransform

all : $(PROG)

# Variables pour la compilation des fichiers
CC        =  gcc
CFLAGS    =  -g -Wall -std=c99
CPPFLAGS  =  -DDEBUG
LDFLAGS   =  -g -lm

pxmtopxm : pxmtopxm.o Util.o
ppmtopgm : ppmtopgm.o Util.o

imagetransform : imagetransform.o Util.o

clean :
	@rm -f *.o

cleanall : clean
	@rm -f $(PROG)

