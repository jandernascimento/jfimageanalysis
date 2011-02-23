#######################################
# Makefile PBM                        #
#                                     #
# E.B.                                #
#######################################


PROG = ppmtopgm pxmtopxm pgmtopgm

all : $(PROG)

# Variables pour la compilation des fichiers
CC        =  gcc
CFLAGS    =  -g -Wall
CPPFLAGS  =  -DDEBUG
LDFLAGS   =  -g -lm

pxmtopxm : pxmtopxm.o Util.o
ppmtopgm : ppmtopgm.o Util.o

pgmtopgm : pgmtopgm.o Util.o

clean :
	@rm -f *.o

cleanall : clean
	@rm -f $(PROG)

