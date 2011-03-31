#######################################
# Makefile PBM                        #
#                                     #
# E.B.                                #
#######################################


PROG = showregion #pxmtopxm pgmtopgm

all : $(PROG)

# Variables pour la compilation des fichiers
CC        =  gcc
CFLAGS    =  -g -Wall -std=c99
CPPFLAGS  =  -DDEBUG
LDFLAGS   =  -g -lm

#pxmtopxm : pxmtopxm.o Util.o
showregion : showregion.o Util.o

#pgmtopgm : pgmtopgm.o Util.o

clean : .PHONY
	@rm -f *.o
	@rm showregion

cleanall : clean
	@rm -f $(PROG)

