#!/bin/sh

TFLAG=false
MFLAG=false
OUTPUT=output

NPROC=20

USAGE="	USAGE: bash pinetree.sh -t path/to/transcript/file -m path/to/mirna/file\n
		optional flags:\n
		-n number of processors\n
		-h help\n
		-A accessibility mode\n
		-N normalize scores between 0 and 1 (cannot be used with accessibility)\n
		-H enable human readable output\n
		-o path/to/output/folder (one will be created if it doesn't exist)\n
		-a path/to/annotation/file"

while getopts ":t:m:n:o:a:ANhH" opt; do
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
    a)
      ANNOTATION="-A $OPTARG"
      ;;
    N)
      NORMALIZATION="-N"
      ;;
    A)
	  ACCESSIBILITY="-z"
	  ;;
	H)
	  HUMAN_OUTPUT="-H"
	  ;;
	h)
      echo -e $USAGE
      exit 1
      ;;
    \?)
      echo -e "Invalid option: -$OPTARG\n" >&2
      exit 1
      ;;
    :)
      echo -e "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

if [ -n "$NORMALIZATION" ] && [ -n "$ACCESSIBILITY" ] ; then
	echo -e $USAGE
	exit 1
fi

if [ $TFLAG == false ] || [ $MFLAG == false ] ; then
	echo -e $USAGE
	exit 1
fi

mkdir -p $OUTPUT

FASTA_OUT=$OUTPUT/results.txt
FASTA_PROC=$OUTPUT/processed_results.csv

START=$(date +"%s")

./fasta36 -T $NPROC -q -n -E 150 $MIRNA $TRANSCRIPT 1 &> $FASTA_OUT 
python process_alignments.py $TRANSCRIPT $MIRNA $FASTA_OUT > $FASTA_PROC

./pinetree -C CONFIG -t $TRANSCRIPT -m $MIRNA -a $FASTA_PROC -o $OUTPUT/pinetree \
$NORMALIZATION $ACCESSIBILITY $ANNOTATION -s $START $HUMAN_OUTPUT






