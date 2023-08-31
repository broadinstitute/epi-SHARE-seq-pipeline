#!/bin/bash
# ------------------------------------------------------------------
# [Author] Eugenio Mattei
#          Make BigWig tracks from fragment file
# ------------------------------------------------------------------

VERSION=0.1.0
SUBJECT=task-make-track
USAGE="Usage: make_track -d [library_size]-o [outfile] chrom_sizes fragment"

# --- Options processing -------------------------------------------
if [ $# == 0 ]; then
    echo $USAGE
    exit 1;
fi

while getopts ":o:d:vh" optname
  do
    case "$optname" in
      "v")
        echo "Version $VERSION"
        exit 0;
        ;;
      "d")
        echo "Library size for scaling: $OPTARG"
        library_size=$OPTARG
        ;;
      "o")
        echo "Final bigwig will be written to: $OPTARG"
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

chrom_sizes=$1
fragment_file=$2

if [ -z "$library_size" ];then
    echo "Computing library size"
    library_size=$(pigz -dc -p 8 $fragment_file | awk '{count+=1}END{print count*2}')
fi
echo "Library size is: $library_size"

# --- Body --------------------------------------------------------
#  SCRIPT LOGIC GOES HERE
pigz -c -d -p 8 $fragment_file | \
awk -v OFS="\t" '{print $1,$2-4,$2+4"\n"$1,$3-4,$3+4}' | \
sort --parallel=4 -k1,1 -k2,2n > $TMPBED

scale_factor=$(bc <<< "scale=6;10000000/$(echo $library_size)")

bedtools merge -i $TMPBED -c 1 -o count | \
awk -v scaling=$scale_factor -v OFS="\t" '{$4=$4*scaling; print $0}' | \
sort --parallel=8 -k1,1 -k2,2n > $TMPGRAPH

bedGraphToBigWig $TMPGRAPH $chrom_sizes $output_file

rm $TMPBED $TMPGRAPH
# -----------------------------------------------------------------