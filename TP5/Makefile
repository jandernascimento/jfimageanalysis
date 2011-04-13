#######################################
# Makefile PBM                        #
#                                     #
# E.B.                                #
#######################################


PROG = tp5  #pxmtopxm pgmtopgm

all : $(PROG)

# Variables pour la compilation des fichiers
CC        =  gcc
CFLAGS    =  -g -Wall -std=c99
CPPFLAGS  =  -DDEBUG
LDFLAGS   =  -g -lm

tp5 : tp5.o Util.o 

clean : 
	@rm -f *.o
	@rm tp5

.PHONY: clean

cleanall : clean
	@rm -f $(PROG)

