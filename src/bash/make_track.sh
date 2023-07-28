#!/bin/bash
# ------------------------------------------------------------------
# [Author] Eugenio Mattei
#          Make BigWig tracks from fragment file
# ------------------------------------------------------------------

VERSION=0.1.0
SUBJECT=task-make-track
USAGE="Usage: make_track -o [outfile] chrom_sizes fragment"

# --- Options processing -------------------------------------------
if [ $# == 0 ]; then
    echo $USAGE
    exit 1;
fi

while getopts ":o:vh" optname
  do
    case "$optname" in
      "v")
        echo "Version $VERSION"
        exit 0;
        ;;
      "o")
        echo "-o argument: $OPTARG"
        output_file=$OPTARG
        ;;
      "h")
        echo $USAGE
        exit 0;
        ;;
      "?")
        echo "Unknown option $OPTARG"
        exit 0;
        ;;
      ":")
        echo "No argument value for option $OPTARG"
        exit 0;
        ;;
      *)
        echo "Unknown error while processing options"
        exit 0;
        ;;
    esac
  done

shift $(($OPTIND - 1))

if [ -z "$2" ]; then
    echo $USAGE
    exit 1;
fi


# --- Locks -------------------------------------------------------
LOCK_FILE=/tmp/$SUBJECT.lock
if [ -f "$LOCK_FILE" ]; then
   echo "Script is already running"
   exit
fi

trap "rm -f $LOCK_FILE" EXIT
touch $LOCK_FILE

TMPBED=$(mktemp)
TMPGRAPH=$(mktemp)

fragment_file=$1
chrom_sizes=$2

# --- Body --------------------------------------------------------
#  SCRIPT LOGIC GOES HERE
pigz -c -d -p 8 $fragment_file | \
awk -v OFS="\t" '{print $1,$2-4,$2+4"\n"$1,$3-4,$3+4}' | \
sort --parallel=4 -k1,1 -k2,2n > $TMPFILE

insertion_number=$(wc -l < $TEMPFILE)
scale_factor=$(bc <<< "scale=6;10000000/$(echo $insertion_number)")

bedtools merge -i $TEMPFILE -c 1 -o count | \
awk -v scaling=$scale_factor -v OFS="\t" '{$4=$4*scaling; print $0}' | \
sort --parallel=8 -k1,1 -k2,2n > $TMPGRAPH

bedGraphToBigWig $TMPGRAPH $chrom_sizes $output_file
# -----------------------------------------------------------------