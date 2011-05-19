#######################################
# Makefile PBM                        #
#                                     #
# E.B.                                #
#######################################


PROG = tp6  #pxmtopxm pgmtopgm

all : $(PROG)

# Variables pour la compilation des fichiers
CC        =  gcc
CFLAGS    =  -g -Wall -std=c99
CPPFLAGS  =  -DDEBUG
LDFLAGS   =  -g -lm

#pxmtopxm : pxmtopxm.o Util.o
tp6 : tp6.o Util.o 

clean : 
	@rm -f *.o
	@rm tp6

.PHONY: clean

cleanall : clean
	@rm -f $(PROG)

