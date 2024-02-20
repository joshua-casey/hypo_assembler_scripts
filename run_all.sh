#!/bin/bash

usage (){
echo "  -1 <reads_1>        short read pair 1                               [          required ]"
echo "  -2 <reads_2>        short read pair 2                               [          required ]"
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
tempdir="temp/"
while getopts "1:2:l:d:B:t:T:h" opt; do
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

if [ "$longbam" == "" ]; then
    echo "Mapping long reads to draft"
    minimap2 -ax map-ont -t 40 $draft $longreads | samtools view -bS | samtools sort -@ $threads -m $sortmem -o $tempdir/long_align.bam
else
    echo "Long read mapping: $longbam"
    ln -s $longbam $tempdir/long_align.bam
fi

echo "[STEP 1] Getting solid kmers"
echo $reads1 > $tempdir/shorts.txt
echo $reads2 >> $tempdir/shorts.txt
./suk -k 17 -i @"$tempdir"/shorts.txt -t 40 -e
mv SUK_k17.bv $tempdir/SUK_k17.bv

echo "[STEP 2] Scanning misjoin"
./find_misjoin $draft $tempdir/long_align.bam $tempdir/misjoin.fa

echo "[STEP 3] Finding overlaps"
./run_overlap.sh -k $tempdir/SUK_k17.bv -i $tempdir/misjoin.fa -l $longreads -t $threads -o $tempdir/overlap -T $tempdir/overlap_temp

echo "[STEP 4] Realignment for polishing"
minimap2 -I 64G -ax map-ont -t $threads $tempdir/overlap.fa $longreads | samtools view -bS | samtools sort -@ 10 -m 10G -o $tempdir/overlap_long.bam
minimap2 -I 64G -ax sr -t $threads $tempdir/overlap.fa $reads1 $reads2 | samtools view -bS | samtools sort -@ 10 -m 10G -o $tempdir/overlap_short.bam

echo "[STEP 5] Polishing"
./hypo -d $tempdir/overlap.fa -s 3g -B $tempdir/overlap_long.bam -C 60 -b $tempdir/overlap_short.bam -r @"$tempdir"/shorts.txt -c 100 -t $threads -w $tempdir/hypo_wdir -o $tempdir/polished.fa

echo "[STEP 6] Scaffolding"
./run_scaffold.sh -k $tempdir/SUK_k17.bv -i $tempdir/polished.fa -l l.fq.gz -t 40 -o $tempdir/scaffold
