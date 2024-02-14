#!/bin/bash

usage (){
echo "  -1 <reads_1>        short read pair 1                               [          required ]"
echo "  -i <reads_2>        short read pair 2                               [          required ]"
echo "  -l <long_reads>     long reads file                                 [          required ]"
echo "  -d <draft>          draft assembly                                  [          required ]"
echo "  -B <long_read_map>  mapping of long reads to draft                  [ default:      none]"
echo "  -t <threads>        the number of threads to use.                   [ default:        1 ]"
echo "  -T <tempdir>   directory to store intermediate files                [ default:    temp/ ]"
echo "  -h             display this help and exit"
exit 1
}

reads1=""
reads2=""
longreads=""
longbam=""
draft=""
threads="1"
sortmem="10G"
kmerlen="17"
while getopts ":k:i:t:fT:o:l:pm:h" opt; do
  case $opt in
    1)
        reads1="$OPTARG"
        ;;
    2)
        reads2="$OPTARG"
        ;;
    l)
        longreads="$OPTARG"
        ;;
    d)
        draft="$OPTARG"
        ;;
    B)
        longbam="$OPTARG"
        ;;
    t)
        threads="$OPTARG"
        ;;
    T)
        tempdir="$OPTARG"
        ;;
    h)
        usage
        ;;
    \?) 
        echo "Invalid option -$OPTARG" >&2
        ;;
  esac
done

if [ "$reads1" == "" ]; then
    echo "Option -1 <reads_1> needed."
    exit 1
else
    echo "Short read file 1: $reads1"
fi

if [ "$reads2" == "" ]; then
    echo "Option -2 <reads_1> needed."
    exit 1
else
    echo "Short read file 2: $reads2"
fi

if [ "$longreads" == "" ]; then
    echo "Option -l <long reads> needed."
    exit 1
else
    echo "Long reads: $longreads"
fi

if [ "$draft" == "" ]; then
    echo "Option -d <draft assemblys> needed."
    exit 1
else
    echo "Draft: $longreads"
fi

echo "Using $threads threads"
echo "Using $tempdir as temporary directory."
mkdir -p $tempdir

cp suk $tempdir
cp hypo $tempdir
cp run_overlap.sh $tempdir
cp run_scaffold.sh $tempdir
cd $tempdir

if [ "$longbam" == "" ]; then
    echo "Mapping long reads to draft"
    minimap2 -ax map-ont -t 40 $draft $longreads | samtools view -bS | samtools sort -@ $threads -m $sortmem -o long_align.bam
else
    echo "Long mapping: $longbam"
    ln -s $longbam long_align.bam
fi

echo "[STEP 1] Getting solid kmers"
echo $reads1 > shorts.txt
echo $reads2 >> shorts.txt
./suk -k 17 -i @shorts.txt -t 40 -e
cd ..

echo "[STEP 2] Scanning misjoin"
python scan_misjoin.py $draft long_align.bam misjoin.fa

echo "[STEP 3] Finding overlaps"
./run_overlap.sh -k SUK_k17.bv -i misjoin.fa -l $longreads -t $threads -o overlap -T overlap_temp

echo "[STEP 4] Realignment for polishing"
minimap2 -I 64G -ax map-ont -t 40 overlap.fa $longreads | samtools view -bS | samtools sort -@ 10 -m 10G -o overlap_long.bam
minimap2 -I 64G -ax sr -t 40 overlap.fa $reads1 $reads2 | samtools view -bS | samtools sort -@ 10 -m 10G -o overlap_short.bam

echo "[STEP 5] Polishing"
./hypo -d $tempdir/overlap.fa -s 3g -B overlap_long.bam -C 60 -b overlap_short.bam -r @shorts.txt -c 100 -t $threads

echo "[STEP 6] Scaffolding"
./run.sh -k SUK_k17.bv -i hypo_overlap.fa -l l.fq.gz -t 40
