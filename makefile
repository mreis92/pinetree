CC = gcc

DEBUG = -g -Wall -Wextra -pedantic -fdiagnostics-show-option #-Wno-long-long -Wno-return-type -Wno-parentheses
OPT = -O3
OMP = -fopenmp

RNALIB = ViennaRNA-2.1.8/lib 
RNAHEADERS = ViennaRNA-2.1.8/H

CUSTOM_LIBS = -lRNA
LIBS = -lm $(CUSTOM_LIBS)
EXE = pinetree
EXE_EVO = pinetree_evo

FILES = annotation.c dataset.c fasta.c pinetree_utils.c strmap.c util.c   
PFILES = pinetree.c accessibility.c align.c RNAup_cmdl.c RNAplfold_cmdl.c RNAup.c RNAplfold.c 
EFILES = model.c evo.c pinetree_evo.c

VER=1.0

pinetree: 
	export OMP_NUM_THREADS=4
	$(CC) $(OPT) $(OMP) $(FILES) $(PFILES) -o $(EXE) -I$(RNAHEADERS) -L$(RNALIB) $(LIBS)
	
evo: 
	export OMP_NUM_THREADS=4
	$(CC) $(OPT) $(OMP) $(FILES) $(EFILES) -o $(EXE_EVO) -L. -lm

debug:
	export OMP_NUM_THREADS=2
	$(CC) $(DEBUG) $(OMP) $(FILES) $(PFILES) -o $(EXE) -I$(RNAHEADERS) -L$(RNALIB) $(LIBS)
	
debug_evo:
	export OMP_NUM_THREADS=2
	$(CC) $(DEBUG) $(OMP) $(FILES) $(EFILES) -o $(EXE_EVO) -L. -lm

all: pinetree evo
	
clean:
	rm -f $(EXE) $(EXE_EVO) *~ *.o

pack-pinetree:
	mkdir -p pinetree-$(VER)
	cp *.c *.h *.py *.sh makefile fasta36 pinetree-$(VER)
	tar cvf pinetree-$(VER)-src.tar pinetree-$(VER)
	gzip -f --best pinetree-$(VER)-src.tar
	rm -rf pinetree-$(VER)



