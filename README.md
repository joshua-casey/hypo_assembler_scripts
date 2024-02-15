### Input files:
1. asm.fa - initial flye assembly
2. l.fq.gz - long reads file
3. r1.fq.gz - short reads file
4. r2.fq.gz - short reads file
5. shorts.txt - text file containing path to r1.fq.gz and r2.fq.gz
6. l.bam - optional - l.fq.gz aligned to asm.fa

### Prereqs:
1. Python lib requirement: pysam and biopython
2. KMC3, minimap2, samtools in path

#### Automated Compilation
1. Run build_all.sh
2. All used executables will be at directory run_all/

#### Manual Compilation
Compile C++ scripts on overlap/ and scaffold/:
1. mkdir build
2. cd build
3. cmake ..
4. make
5. Put the resulting binary (i.e. find_overlap and find_scaffold) in base directory (i.e. overlap/ and scaffold/)
    i.e. cp find_overlap ../ and cp find_scaffold ../
Build hypo polisher
1. mkdir build
2. cmake ..
3. make
4. Executable will be in build/bin/hypo

### Running Steps

#### Automated Running Steps
1. Run run_all.sh on run_all/

2. All used executables will be at directory run_all/

1. Align l.fq.gz to asm.fa and sort
minimap2 -ax map-ont -t 40 asm.fa l.fq.gz | samtools view -bS | samtools sort -@ 10 -m 10G -o long_read_align.bam
Output: long_read_align.bam

2. Run suk on shorts.txt
./suk -k 17 -i @shorts.txt -t 40 -e
Output: SUK_k17.bv

3. Run scan_misjoin.py
python scan_misjoin.py asm.fa long_read_align.bam misjoin.fa
Output: misjoin.fa

4. Run overlap/run_overlap.sh
./run_overlap.sh -k SUK_k17.bv -i misjoin.fa -l l.fq.gz -t 40
Output: overlap.fa

5. Realign short and long reads to overlap.fa
minimap2 -ax map-ont -t 40 overlap.fa l.fq.gz | samtools view -bS | samtools sort -@ 10 -m 10G -o overlap_long.bam
minimap2 -ax sr -t 40 overlap.fa r1.fq.gz r2.fq.gz | samtools view -bS | samtools sort -@ 10 -m 10G -o overlap_short.bam
Output: overlap_long.bam and overlap_short.bam

5. Run hypo polisher
./hypo -d overlap.fa -s 3g -B overlap_long.bam -C 60 -b overlap_short.bam -r @shorts.txt -c 100 -t 40
Output: hypo_overlap.fa

6. Run scaffold/run_scaffold.sh
./run_scaffold.sh -k SUK_k17.bv -i hypo_overlap.fa -l l.fq.gz -t 40
Output: scaffold_1.fa and scaffold_2.fa
