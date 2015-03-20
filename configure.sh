wget http://www.tbi.univie.ac.at/RNA/download.php?id=viennarna-2.1.8 -O ViennaRNA.tar.gz &
wget http://faculty.virginia.edu/wrpearson/fasta/CURRENT/fasta-36.3.7a.tar.gz -O FASTA.tar.gz

tar -zxvf ViennaRNA.tar.gz
rm ViennaRNA.tar.gz

tar -zxvf FASTA.tar.gz
rm FASTA.tar.gz

wait 

(cd ViennaRNA-2.1.8; ./configure; make) &
(cd fasta-36.3.7a/src; make -f ../make/Makefile.linux_sse2 all; mv ../bin/fasta36 ../../ ; rm -rf ../../fasta-36.3.7a) &

wait

make all