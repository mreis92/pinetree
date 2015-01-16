#!/bin/sh

TFLAG=false
MFLAG=false
OUTPUT=output

NPROC=20

USAGE="	USAGE: bash pinetree.sh -t path/to/transcript/file -m path/to/mirna/file\n
		optional flags:\n
		-n number of processors\n
		-h help\n
		-o path/to/output/folder (one will be created if it doesn't exist)"

while getopts ":t:m:n:o:h" opt; do
  case $opt in
    t)
      TRANSCRIPT=$OPTARG
      TFLAG=true
      ;;
	m)
      MIRNA=$OPTARG
      MFLAG=true
      ;;
	n)
      NPROC=$OPTARG
      ;;
	o)
      OUTPUT=$OPTARG
      ;;
	h)
      echo -e $USAGE
      exit 1
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

if [ $TFLAG == false ] || [ $MFLAG == false ] ; then
	echo -e $USAGE
	exit 1
fi

mkdir -p $OUTPUT

FASTA_OUT=$OUTPUT/results.txt
FASTA_PROC=$OUTPUT/processed_results.csv
	
./fasta36 -T $NPROC -q -n -E 150 $MIRNA $TRANSCRIPT 1 &> $FASTA_OUT 

python process_alignments.py $TRANSCRIPT $MIRNA $FASTA_OUT > $FASTA_PROC

./pinetree -C CONFIG -t $TRANSCRIPT -m $MIRNA -a $FASTA_PROC -o $OUTPUT/pinetree 

rm $FASTA_OUT
rm $FASTA_PROC






