#!/bin/bash
#SBATCH -p priority
#SBATCH -t 0-12:00
#SBATCH -c 4
#SBATCH --mem=40G
#SBATCH -o %j.out
#SBATCH -e %j.err 
#SBATCH--mail-type=ALL


##########################################
#                                        #
#            USAGE AND SETUP             #
#                                        #
##########################################

# USAGE
USAGE="
DESCRIPTION
This workflow is designed to find non-templated 3'-end repetitive tails in TAIL-seq sequencing data

USAGE
sbatch tail-seq.sh \$1 \$2 \$3

COMMAND LINE VARIABLES
\$1 = Sample (String. Prefix of fastq files, e.g. fastq files = ./../SRR1005384_1.fastq and ./../SRR1005384_2.fastq, Sample = SRR1005384)
\$2 = Species (Possible values: human, mouse, fly, worm, worm-orsay, or arabidopsis)
\$3 = Repeat (Possible values: combinations of A, T, C, G, and U)

EXAMPLE
sbatch tail-seq.sh SRR1005384 human UGUGUGUGUGUGUGUGUGUG

NOTE
Sometimes repeats are found in read names which causes the workflow to crash. Strip read names of unnecessary information by running this command: sed -i \"s/ .*//\" SRR1005384_2.fastq
"
if [[ $1 = "--help" || $1 = "-h" ]] || [[ -z "$1" ]]
then
	echo "$USAGE"
	exit 0
fi

# SET REPEAT AND SPECIES-SPECIFIC VARIABLES BASED ON COMMAND LINE INPUTS
# make repeat variable
REPEAT=$(echo "$3")
# covert $REPEAT to uppercase
REPEAT=$(printf '%s\n' "$REPEAT" | awk '{ print toupper($0) }')
# convert U to T
REPEAT=$(echo "$REPEAT" | sed 's/U/T/g')
# check that $REPEAT contains only standard IUPAC nucleotides
REPEAT_CHARACTERS=$(echo $REPEAT | sed -e "s/./\0\n/g" | sort -u)
for i in $(echo $REPEAT_CHARACTERS)
do
	if [[ "$i" != A ]] &&
	   [[ "$i" != T ]] &&	
	   [[ "$i" != C ]] &&
	   [[ "$i" != G ]] &&
	   [[ "$i" != U ]]
	then
    	echo "$USAGE"
    	echo "ERROR: The repeat can only be comprised of combinations of A, T, C, G and U"
    	exit 1
	fi
done

# set species-specific variables (genomes, annotation files, etc.)
if [[ $2 = "human" ]]
then
	STAR_GENOME=/n/groups/kennedy/pagano/pUG_RNAs/references/hg38/hg38_STAR_genomeDir
	FEATURE_COUNTS_ANNO=/n/groups/kennedy/pagano/pUG_RNAs/references/hg38/Homo_sapiens.GRCh38.100.with.rmsk_TE.gtf
	TAXID=9606
	MAX_INTRON=0
	MIN_INTRON=21
else
	if [[ $2 = "mouse" ]]
	then
		STAR_GENOME=/n/groups/kennedy/pagano/pUG_RNAs/references/m38/m38_STAR_genomeDir
		FEATURE_COUNTS_ANNO=/n/groups/kennedy/pagano/pUG_RNAs/references/m38/Mus_musculus.GRCm38.100.with.rmsk_TE.gtf
		TAXID=10090
		MAX_INTRON=0
		MIN_INTRON=21
	else
		if [[ $2 = "fly" ]]
		then
			STAR_GENOME=/n/groups/kennedy/pagano/pUG_RNAs/references/dm6/dm6_STAR_genomeDir
			FEATURE_COUNTS_ANNO=/n/groups/kennedy/pagano/pUG_RNAs/references/dm6/Drosophila_melanogaster.BDGP6.28.100.with.rmsk_TE.gtf
			TAXID=7227
			MAX_INTRON=250000
			MIN_INTRON=25
		else
			if [[ $2 = "worm" ]]
			then
				STAR_GENOME=/n/groups/kennedy/pagano/pUG_RNAs/references/ce11/ce11_STAR_genomeDir
				FEATURE_COUNTS_ANNO=/n/groups/kennedy/pagano/pUG_RNAs/references/ce11/Caenorhabditis_elegans.WBcel235.100.with.rmsk_TE.gtf
				TAXID=6239
				MAX_INTRON=105000
				MIN_INTRON=20
			else
			    if [[ $2 = "arabidopsis" ]]
				then
					STAR_GENOME=/n/groups/kennedy/pagano/pUG_RNAs/references/TAIR10/TAIR10_STAR_genomeDir
					FEATURE_COUNTS_ANNO=/n/groups/kennedy/pagano/pUG_RNAs/references/TAIR10/Arabidopsis_thaliana.TAIR10.47.with.rmsk_TE.gtf
					TAXID=3702
					MAX_INTRON=6000
					MIN_INTRON=40
				else
			    	if [[ $2 = "worm-orsay" ]]
					then
						STAR_GENOME=/n/groups/kennedy/pagano/pUG_RNAs/references/ce11_orsay_STAR_genomeDir
						FEATURE_COUNTS_ANNO=/n/groups/kennedy/pagano/pUG_RNAs/references/orsay_virus/Caenorhabditis_elegans.WBcel235.100.with.rmsk_TE.and.orsay.RNAs.gtf
						TAXID=6239,977912
						MAX_INTRON=105000
						MIN_INTRON=20
					else
			    		echo "$USAGE"
			    		echo "ERROR: Not a valid value for species. Possible values: human, mouse, fly, worm, or arabidopsis"
			    		exit 1
			    	fi
				fi
			fi
		fi
	fi
fi

# load modules
module load star/2.7.3a
module load samtools/1.9
module load datamash/1.2
module load java/jdk-1.8u112

# make analysis directory
if [[ ! -d "$1"_"$3"_TAIL-seq ]]
then
	mkdir "$1"_"$3"_TAIL-seq
fi

# change to analysis directory
cd "$1"_"$3"_TAIL-seq

# make REPEAT_REVCOMP variable
paste -d '\n' <(echo \>REPEAT) <(echo "$REPEAT") > repeat.fa
REPEAT_REVCOMP=$(/n/groups/kennedy/pagano/modules/seqtk/seqtk seq -r repeat.fa | grep -v ">")
rm repeat.fa

##########################################
#                                        #
#   CANDIDATE GATHERING AND ALIGNMENT    #
#                                        #
##########################################

# determine number of PE reads
TOTAL_PE_READS=$(echo $(cat ./../"$1"_2.fastq | wc -l)/4 | bc)

# remove duplicates?

# grab reads with repeat
grep -h -i --no-group-separator $REPEAT_REVCOMP -B 1 -A 2 ./../"$1"_2.fastq > "$1"_"$REPEAT"_candidates.fastq

# determine number of reads with repeat
READS_WITH_REPEAT=$(echo $(cat "$1"_"$REPEAT"_candidates.fastq | wc -l)/4 | bc)

# compile list of candidate read names
awk 'NR % 4 == 1' "$1"_"$REPEAT"_candidates.fastq | awk -F " " '{print $1}' | sed 's/^@//' > "$1"_"$REPEAT"_candidates_read_names

# align reads to genome
mkdir STAR
STAR --runMode alignReads \
--runThreadN 4 \
--outFilterScoreMinOverLread 0 \
--outFilterMatchNminOverLread 0.3 \
--genomeDir $STAR_GENOME \
--alignIntronMax $MAX_INTRON \
--alignIntronMin $MIN_INTRON \
--outSAMtype BAM Unsorted \
--outFileNamePrefix ./STAR/"$1"_"$REPEAT"_candidates_ \
--readFilesIn "$1"_"$REPEAT"_candidates.fastq

# sort alignments
samtools sort --threads 3 ./STAR/"$1"_"$REPEAT"_candidates_Aligned.out.bam -o "$1"_"$REPEAT"_candidates_sorted.bam

# remove secondary aligments
samtools view -h --threads 3 -F 256 "$1"_"$REPEAT"_candidates_sorted.bam -O SAM -o "$1"_"$REPEAT"_candidates_primary.sam 

# compile list of mapped candidates read names
cut -f 1 "$1"_"$REPEAT"_candidates_primary.sam | grep -v "^@" | sort > "$1"_"$REPEAT"_candidates_mapped_read_names_.2

# compile list of unmapped candidates read names
sort "$1"_"$REPEAT"_candidates_read_names | comm - "$1"_"$REPEAT"_candidates_mapped_read_names_.2 | cut -f 1 | sort | awk '$1 != ""' > "$1"_"$REPEAT"_candidates_unmapped_read_names_.2

# create read pair .1 identifiers for mapped and unmapped reads
sed 's/.2$/.1/' "$1"_"$REPEAT"_candidates_mapped_read_names_.2 > "$1"_"$REPEAT"_candidates_mapped_read_names_.1 &
sed 's/.2$/.1/' "$1"_"$REPEAT"_candidates_unmapped_read_names_.2 > "$1"_"$REPEAT"_candidates_unmapped_read_names_.1 &
wait


##########################################
#                                        #
#     PROCESSING UNMAPPED CANDIDATES     #
#                                        #
##########################################

# determine number of candidates that failed to map
UNMAPPED_CANDS=$(cat "$1"_"$REPEAT"_candidates_unmapped_read_names_.2 | wc -l)

# grab read 1s of unmapped candidate read 2s
/n/groups/kennedy/pagano/modules/bbmap/filterbyname.sh \
in=./../"$1"_1.fastq \
out="$1"_read_1s_of_unmapped_"$REPEAT"_candidates.fastq \
names="$1"_"$REPEAT"_candidates_unmapped_read_names_.1 \
include=t

# align reads to genome
STAR --runMode alignReads \
--runThreadN 4 \
--genomeDir $STAR_GENOME \
--alignIntronMax $MAX_INTRON \
--alignIntronMin $MIN_INTRON \
--outSAMtype BAM Unsorted \
--outFileNamePrefix ./STAR/"$1"_read_1s_of_unmapped_"$REPEAT"_candidates_ \
--readFilesIn "$1"_read_1s_of_unmapped_"$REPEAT"_candidates.fastq

########################
###FOR MAPPED READ 1s###
######################## 
# convert bam to sam
samtools view -h --threads 3 -O SAM ./STAR/"$1"_read_1s_of_unmapped_"$REPEAT"_candidates_Aligned.out.bam -o "$1"_read_1s_of_unmapped_"$REPEAT"_candidates_unfiltered.sam

# compile list of mapped read 1s from unmapped candidate read 2s
cut -f 1 "$1"_read_1s_of_unmapped_"$REPEAT"_candidates_unfiltered.sam | grep -v "^@" | sort | uniq > "$1"_mapped_read_1s_read_names_.1

# create read pair .2 identifiers for mapped and unmapped reads
sed 's/.1$/.2/' "$1"_mapped_read_1s_read_names_.1 > "$1"_mapped_read_1s_read_names_.2

# grab read 2s of mapped read 1s
/n/groups/kennedy/pagano/modules/bbmap/filterbyname.sh \
in=./../"$1"_2.fastq \
out="$1"_read_2s_of_mapped_read_1s_unfiltered.fastq \
names="$1"_mapped_read_1s_read_names_.2 \
include=t

# covert fastq to fasta
/n/groups/kennedy/pagano/modules/seqtk/seqtk seq -A "$1"_read_2s_of_mapped_read_1s_unfiltered.fastq > "$1"_read_2s_of_mapped_read_1s_unfiltered.fasta

# blast read 2s of mapped read 1s against nr/nt database
blastn -db /n/groups/kennedy/pagano/databases/blastdb_nt/nt -query "$1"_read_2s_of_mapped_read_1s_unfiltered.fasta -out "$1"_best_blastn_target_for_read_2s_of_mapped_read_1s_unfiltered.txt -evalue 1e-05 -max_target_seqs 1 -subject_besthit -taxids "$TAXID" -outfmt "6 std scomnames stitle qcovs qseq"

# establish false positive cutoff
REP_LEN=${#REPEAT}
READ_LENGTH=$(/n/groups/kennedy/pagano/modules/seqtk/seqtk comp "$1"_read_2s_of_mapped_read_1s_unfiltered.fastq | datamash count 2 mean 2 sstdev 2 | cut -f 2)
CUTOFF=$(echo \(\($READ_LENGTH-$REP_LEN\)*100\)/$READ_LENGTH | bc)

# create list of reads with blast hits
cut -f 1 "$1"_best_blastn_target_for_read_2s_of_mapped_read_1s_unfiltered.txt | sort | uniq > "$1"_unmapped_read_2s_with_mapped_read_1s_and_blast_hits

# identify false positives
> "$1"_read_2s_of_mapped_read_1s_false_positives_.2
for i in $(cat "$1"_unmapped_read_2s_with_mapped_read_1s_and_blast_hits)
do 
	grep -m 1 -F "$i" "$1"_best_blastn_target_for_read_2s_of_mapped_read_1s_unfiltered.txt > blastn_temp_"$i"
	paste <(cut -f 1-15 blastn_temp_"$i") <(awk '{print $NF}' blastn_temp_"$i" | sed 's/-//g') > blastn_temp_"$i"_hyphens_removed
	mv blastn_temp_"$i"_hyphens_removed blastn_temp_"$i"
	QUERY_SIZE=$(echo $(awk '{print $NF}' blastn_temp_"$i" | wc -c) - 1 | bc)
	QUERY_COV=$(echo $QUERY_SIZE\*100\/$READ_LENGTH | bc)
	paste <(cat blastn_temp_"$i") <(echo $QUERY_COV) > blastn_temp_"$i"_query_coverage_added
	mv blastn_temp_"$i"_query_coverage_added blastn_temp_"$i"
	if [[ $(awk '{print $NF}' blastn_temp_"$i") -le $CUTOFF ]]
	then
		if [[ $(grep "$REPEAT_REVCOMP" blastn_temp_"$i" | wc -l) = 0 ]]
		then 
			cat blastn_temp_"$i" >> "$1"_read_2s_of_mapped_read_1s_blastn_hits_with_templated_hits_removed.txt
		else
			echo $i >> "$1"_read_2s_of_mapped_read_1s_false_positives_.2
		fi	
	else
		echo $i >> "$1"_read_2s_of_mapped_read_1s_false_positives_.2
	fi	
	rm blastn_temp_"$i"
done

# create read pair .1 identifiers
sed 's/.2$/.1/' "$1"_read_2s_of_mapped_read_1s_false_positives_.2 | sort > "$1"_read_2s_of_mapped_read_1s_false_positives_.1

# filter false positives from read name list
sort "$1"_mapped_read_1s_read_names_.1 | comm - "$1"_read_2s_of_mapped_read_1s_false_positives_.1 | cut -f 1 | sort | awk '$1 != ""' > "$1"_mapped_read_1s_read_names_false_positives_removed_.1

# create read pair .2 identifiers
sed 's/.1$/.2/' "$1"_mapped_read_1s_read_names_false_positives_removed_.1 > "$1"_mapped_read_1s_read_names_false_positives_removed_.2
cat "$1"_mapped_read_1s_read_names_false_positives_removed_.1 "$1"_mapped_read_1s_read_names_false_positives_removed_.2 | sort > "$1"_mapped_read_1s_read_names_false_positives_removed_.interleaved

# create an interleaved fastq of single-end (read 1) mapped repeat candidates
/n/groups/kennedy/pagano/modules/bbmap/filterbyname.sh \
in=./../"$1"_1.fastq \
in2=./../"$1"_2.fastq \
out="$1"_"$REPEAT"_candidates_READ1_MAPPED_READ2_UNMAPPED_interleaved.fastq \
names="$1"_mapped_read_1s_read_names_false_positives_removed_.interleaved \
include=t

# create a fastq of single-end (read 1) mapped repeat candidates
/n/groups/kennedy/pagano/modules/bbmap/filterbyname.sh \
in=./../"$1"_1.fastq \
out="$1"_"$REPEAT"_candidates_READ1_MAPPED_READ2_UNMAPPED_1.fastq \
names="$1"_mapped_read_1s_read_names_false_positives_removed_.1 \
include=t

# align reads to genome
STAR --runMode alignReads \
--runThreadN 4 \
--genomeDir $STAR_GENOME \
--alignIntronMax $MAX_INTRON \
--alignIntronMin $MIN_INTRON \
--outSAMtype BAM Unsorted \
--outFileNamePrefix ./STAR/"$1"_read_1s_of_unmapped_"$REPEAT"_candidates_false_positives_removed_ \
--readFilesIn "$1"_"$REPEAT"_candidates_READ1_MAPPED_READ2_UNMAPPED_1.fastq

# sort bam, index, and convert to sam
samtools sort --threads 3 ./STAR/"$1"_read_1s_of_unmapped_"$REPEAT"_candidates_false_positives_removed_Aligned.out.bam -O BAM -o "$1"_"$REPEAT"_candidates_READ1_MAPPED_READ2_UNMAPPED.bam
samtools index "$1"_"$REPEAT"_candidates_READ1_MAPPED_READ2_UNMAPPED.bam
samtools view --threads 3 -O SAM "$1"_"$REPEAT"_candidates_READ1_MAPPED_READ2_UNMAPPED.bam -o "$1"_"$REPEAT"_candidates_READ1_MAPPED_READ2_UNMAPPED.sam

# determine number of READ1s_MAPPED_READ2s_UNMAPPED
READ1s_MAPPED_READ2s_UNMAPPED=$(cat "$1"_mapped_read_1s_read_names_false_positives_removed_.1 | wc -l)

##########################
###FOR UNMAPPED READ 1s###
########################## 
# compile list of unmapped read 1s from unmapped candidate read 2s
sort "$1"_"$REPEAT"_candidates_unmapped_read_names_.1 | comm - "$1"_mapped_read_1s_read_names_.1 | cut -f 1 | sort | awk '$1 != ""' > "$1"_unmapped_read_1s_read_names_.1

# create read pair .2 identifiers
sed 's/.1$/.2/' "$1"_unmapped_read_1s_read_names_.1 > "$1"_unmapped_read_1s_read_names_.2
cat "$1"_unmapped_read_1s_read_names_.1 "$1"_unmapped_read_1s_read_names_.2 | sort > "$1"_unmapped_read_1s_read_names_.interleaved

# create an interleaved fastq of single-end (read 1) unmapped repeat candidates
/n/groups/kennedy/pagano/modules/bbmap/filterbyname.sh \
in=./../"$1"_1.fastq \
in2=./../"$1"_2.fastq \
out="$1"_"$REPEAT"_candidates_READ1_UNMAPPED_READ2_UNMAPPED_interleaved.fastq \
names="$1"_unmapped_read_1s_read_names_.interleaved \
include=t

# covert fastq to fasta
/n/groups/kennedy/pagano/modules/seqtk/seqtk seq -A "$1"_"$REPEAT"_candidates_READ1_UNMAPPED_READ2_UNMAPPED_interleaved.fastq > "$1"_"$REPEAT"_candidates_READ1_UNMAPPED_READ2_UNMAPPED_interleaved.fasta

# blast read 2s of mapped read 1s against nr/nt database
blastn -db /n/groups/kennedy/pagano/databases/blastdb_nt/nt -query "$1"_"$REPEAT"_candidates_READ1_UNMAPPED_READ2_UNMAPPED_interleaved.fasta -out "$1"_best_blastn_target_for_UNMAPPED_read_1_and_read_2_candidates.txt -evalue 1e-05 -max_target_seqs 1 -subject_besthit -outfmt "6 std scomnames stitle qcovs qseq"

# determine number of READ1s_UNMAPPED_READ2s_UNMAPPED
READ1s_UNMAPPED_READ2s_UNMAPPED=$(cat "$1"_unmapped_read_1s_read_names_.1 | wc -l)


##########################################
#                                        #
#      PROCESSING MAPPED CANDIDATES      #
#                                        #
##########################################

# determine number of candidates that mapped
MAPPED_CANDS=$(cat "$1"_"$REPEAT"_candidates_mapped_read_names_.2 | wc -l)

# separate forward and reverse alignments
samtools view --threads 3 -F 16 "$1"_"$REPEAT"_candidates_primary.sam -o "$1"_primary_forward_alignments.sam & 
samtools view --threads 3 -f 16 "$1"_"$REPEAT"_candidates_primary.sam -o "$1"_primary_reverse_alignments.sam &
wait

# select soft-clipped alignments
samtools view -H "$1"_"$REPEAT"_candidates_sorted.bam > "$1"_soft_clipped_forward_alignments.sam &
samtools view -H "$1"_"$REPEAT"_candidates_sorted.bam > "$1"_soft_clipped_reverse_alignments.sam &
wait
awk '$6 ~ /S/' "$1"_primary_forward_alignments.sam >> "$1"_soft_clipped_forward_alignments.sam &
awk '$6 ~ /S/' "$1"_primary_reverse_alignments.sam >> "$1"_soft_clipped_reverse_alignments.sam &
wait

# convert sam to fastq
samtools fastq --threads 3 "$1"_soft_clipped_forward_alignments.sam > "$1"_soft_clipped_forward_alignments.fastq &
samtools fastq --threads 3 "$1"_soft_clipped_reverse_alignments.sam > "$1"_soft_clipped_reverse_alignments.fastq &
wait

# grab read names (to be used in for loops)
cut -f 1 "$1"_soft_clipped_forward_alignments.sam | grep -v ^@ > "$1"_soft_clipped_forward_alignments_read_names &
cut -f 1 "$1"_soft_clipped_reverse_alignments.sam | grep -v ^@ > "$1"_soft_clipped_reverse_alignments_read_names &
wait

# extract and expand cigar scores
paste <(grep -v ^@ "$1"_soft_clipped_forward_alignments.sam | cut -f 1) <(grep -v ^@ "$1"_soft_clipped_forward_alignments.sam | cut -f 6 | sed 's/[A-Z]/&\t/g') > "$1"_soft_clipped_forward_alignments_cigar_scores &
paste <(grep -v ^@ "$1"_soft_clipped_reverse_alignments.sam | cut -f 1) <(grep -v ^@ "$1"_soft_clipped_reverse_alignments.sam | cut -f 6 | sed 's/[A-Z]/&\t/g') > "$1"_soft_clipped_reverse_alignments_cigar_scores &
wait

# create a fastq of 5' clippings that are greater than or equal to the length of the repeat
# in TAIL-seq, read 2 is the reverse read from the 3'-end of the RNA, so we're looking for non-templated repeats at the 5' end of this read.
# cigar score is in reference to the reference genome (Watson strand).
	#for forward alignments, S in field 2 indicates 5' soft clipping.
	#for reverse alignments, S in the last field indicates 5' soft clipping.
> "$1"_forward_alignments_5p-end_read_clips.fastq
> "$1"_reverse_alignments_5p-end_read_clips.fastq

REP_LEN=${#REPEAT}

for i in $(cat "$1"_soft_clipped_forward_alignments_read_names)
do
	grep -m 1 -F -w "$i" "$1"_soft_clipped_forward_alignments_cigar_scores > forward_alignment
	if [[ $(awk '{print $2}' forward_alignment | sed 's/[0-9]//g') = S ]]
	then
		fivep_clip=$(awk '{print $2}' forward_alignment | sed 's/[S]//g')
		if [[ $fivep_clip -ge "$REP_LEN" ]]
		then
			grep -m 1 -F -w @"$i" -A 3 "$1"_soft_clipped_forward_alignments.fastq > "$i".fastq
			awk 'FNR>=1 && FNR<=1' "$i".fastq >> "$1"_forward_alignments_5p-end_read_clips.fastq
			awk 'FNR>=2 && FNR<=2' "$i".fastq | cut -c 1-"$fivep_clip" >> "$1"_forward_alignments_5p-end_read_clips.fastq
			awk 'FNR>=3 && FNR<=3' "$i".fastq >> "$1"_forward_alignments_5p-end_read_clips.fastq
			awk 'FNR>=4 && FNR<=4' "$i".fastq | cut -c 1-"$fivep_clip" >> "$1"_forward_alignments_5p-end_read_clips.fastq
		fi
	fi
	if [[ -f "$i".fastq ]]
	then
		rm "$i".fastq
	fi
	rm forward_alignment
done &
for i in $(cat "$1"_soft_clipped_reverse_alignments_read_names)
do
	grep -m 1 -F -w "$i" "$1"_soft_clipped_reverse_alignments_cigar_scores > reverse_alignment
	if [[ $(awk '{print $NF}' reverse_alignment | sed 's/[0-9]//g') = S ]]
	then
		fivep_clip=$(awk '{print $NF}' reverse_alignment | sed 's/[S]//g')
		if [[ $fivep_clip -ge "$REP_LEN" ]]
		then
			grep -m 1 -F -w @"$i" -A 3 "$1"_soft_clipped_reverse_alignments.fastq > "$i".fastq
			awk 'FNR>=1 && FNR<=1' "$i".fastq >> "$1"_reverse_alignments_5p-end_read_clips.fastq
			awk 'FNR>=2 && FNR<=2' "$i".fastq | cut -c 1-"$fivep_clip" >> "$1"_reverse_alignments_5p-end_read_clips.fastq
			awk 'FNR>=3 && FNR<=3' "$i".fastq >> "$1"_reverse_alignments_5p-end_read_clips.fastq
			awk 'FNR>=4 && FNR<=4' "$i".fastq | cut -c 1-"$fivep_clip" >> "$1"_reverse_alignments_5p-end_read_clips.fastq
		fi
	fi
	if [[ -f "$i".fastq ]]
	then
		rm "$i".fastq
	fi
	rm reverse_alignment
done &
wait

# combine forward alignment and reverse alignment clippings
cat "$1"_forward_alignments_5p-end_read_clips.fastq "$1"_reverse_alignments_5p-end_read_clips.fastq > "$1"_all_alignments_5p-end_read_clips.fastq

# grab clippings with repeat
grep -i --no-group-separator "$REPEAT_REVCOMP" -B 1 -A 2 "$1"_all_alignments_5p-end_read_clips.fastq > "$1"_all_alignments_5p-end_read_clips_containing_"$REPEAT".fastq

# compile a list of read names of clippings with repeat
awk 'NR % 4 == 1' "$1"_all_alignments_5p-end_read_clips_containing_"$REPEAT".fastq | sed 's/^@//' > "$1"_5p-end_read_clips_containing_"$REPEAT"_read_names_.2

# grab full length reads
/n/groups/kennedy/pagano/modules/bbmap/filterbyname.sh \
in=./../"$1"_2.fastq \
out="$1"_5p-end_"$REPEAT"_clips_full_length_read_2s_unfiltered.fastq \
names="$1"_5p-end_read_clips_containing_"$REPEAT"_read_names_.2 \
include=t

# covert fastq to fasta
/n/groups/kennedy/pagano/modules/seqtk/seqtk seq -A "$1"_5p-end_"$REPEAT"_clips_full_length_read_2s_unfiltered.fastq > "$1"_5p-end_"$REPEAT"_clips_full_length_read_2s_unfiltered.fasta

# blast full-length read 2s of clippings with repeat against nr/nt database
blastn -db /n/groups/kennedy/pagano/databases/blastdb_nt/nt -query "$1"_5p-end_"$REPEAT"_clips_full_length_read_2s_unfiltered.fasta -out "$1"_best_blastn_target_for_5p-end_"$REPEAT"_clips_full_length_read_2s_unfiltered.txt -evalue 1e-05 -max_target_seqs 1 -subject_besthit -taxids "$TAXID" -outfmt "6 std scomnames stitle qcovs qseq"

# establish false positive cutoff
REP_LEN=${#REPEAT}
READ_LENGTH=$(/n/groups/kennedy/pagano/modules/seqtk/seqtk comp "$1"_5p-end_"$REPEAT"_clips_full_length_read_2s_unfiltered.fastq | datamash count 2 mean 2 sstdev 2 | cut -f 2)
CUTOFF=$(echo \(\($READ_LENGTH-$REP_LEN\)*100\)/$READ_LENGTH | bc)

# create list of reads with blast hits
cut -f 1 "$1"_best_blastn_target_for_5p-end_"$REPEAT"_clips_full_length_read_2s_unfiltered.txt | sort | uniq > "$1"_mapped_read_2s_with_5p-end_"$REPEAT"_clips_and_blast_hits

# identify false positives
> "$1"_candidate_read_2s_false_positives_.2
for i in $(cat "$1"_mapped_read_2s_with_5p-end_"$REPEAT"_clips_and_blast_hits)
do 
	grep -m 1 -F "$i" "$1"_best_blastn_target_for_5p-end_"$REPEAT"_clips_full_length_read_2s_unfiltered.txt > blastn_temp_"$i"
	paste <(cut -f 1-15 blastn_temp_"$i") <(awk '{print $NF}' blastn_temp_"$i" | sed 's/-//g') > blastn_temp_"$i"_hyphens_removed
	mv blastn_temp_"$i"_hyphens_removed blastn_temp_"$i"
	QUERY_SIZE=$(echo $(awk '{print $NF}' blastn_temp_"$i" | wc -c) - 1 | bc)
	QUERY_COV=$(echo $QUERY_SIZE\*100\/$READ_LENGTH | bc)
	paste <(cat blastn_temp_"$i") <(echo $QUERY_COV) > blastn_temp_"$i"_query_coverage_added
	mv blastn_temp_"$i"_query_coverage_added blastn_temp_"$i"
	if [[ $(awk '{print $NF}' blastn_temp_"$i") -le $CUTOFF ]]
	then
		if [[ $(grep "$REPEAT_REVCOMP" blastn_temp_"$i" | wc -l) = 0 ]]
		then 
			cat blastn_temp_"$i" >> "$1"_read_2s_with_5p-end_"$REPEAT"_clips_blastn_hits_with_templated_hits_removed.txt
		else
			echo $i >> "$1"_candidate_read_2s_false_positives_.2
		fi	
	else
		echo $i >> "$1"_candidate_read_2s_false_positives_.2
	fi	
	rm blastn_temp_"$i"
done

# sort false positive read names
sort "$1"_candidate_read_2s_false_positives_.2 > temp
mv temp "$1"_candidate_read_2s_false_positives_.2

# filter false positives from read name list
sort "$1"_5p-end_read_clips_containing_"$REPEAT"_read_names_.2 | comm - "$1"_candidate_read_2s_false_positives_.2 | cut -f 1 | sort | awk '$1 != ""' > "$1"_5p-end_read_clips_containing_"$REPEAT"_read_names_false_positives_removed_.2

# create read pair .1 identifiers
sed 's/.2$/.1/' "$1"_5p-end_read_clips_containing_"$REPEAT"_read_names_false_positives_removed_.2 > "$1"_5p-end_read_clips_containing_"$REPEAT"_read_names_false_positives_removed_.1

# grab reads
/n/groups/kennedy/pagano/modules/bbmap/filterbyname.sh \
in=./../"$1"_1.fastq \
out="$1"_5p-end_read_clips_containing_"$REPEAT"_full_length_read_1.fastq \
names="$1"_5p-end_read_clips_containing_"$REPEAT"_read_names_false_positives_removed_.1 \
include=t &
/n/groups/kennedy/pagano/modules/bbmap/filterbyname.sh \
in=./../"$1"_2.fastq \
out="$1"_5p-end_read_clips_containing_"$REPEAT"_full_length_read_2.fastq \
names="$1"_5p-end_read_clips_containing_"$REPEAT"_read_names_false_positives_removed_.2 \
include=t &
wait

# align candidate paired-reads to genome
# outFilterMatchNminOverLread needs to be lower here 
STAR --runMode alignReads \
--runThreadN 4 \
--outFilterScoreMinOverLread 0 \
--outFilterMatchNminOverLread 0.2 \
--outFilterMultimapNmax 200 \
--genomeDir $STAR_GENOME \
--alignIntronMax $MAX_INTRON \
--alignIntronMin $MIN_INTRON \
--outSAMtype BAM Unsorted \
--outSAMunmapped Within KeepPairs \
--outSAMorder Paired \
--outFileNamePrefix ./STAR/"$1"_paired_reads_of_"$REPEAT"_candidates_ \
--outReadsUnmapped Fastx \
--readFilesIn "$1"_5p-end_read_clips_containing_"$REPEAT"_full_length_read_1.fastq "$1"_5p-end_read_clips_containing_"$REPEAT"_full_length_read_2.fastq

# sort alignments
samtools sort --threads 3 -O BAM ./STAR/"$1"_paired_reads_of_"$REPEAT"_candidates_Aligned.out.bam -o "$1"_"$REPEAT"_candidates_ALL_PAIRS.bam

# select properly paired alignments
samtools view -h -b -f 2 "$1"_"$REPEAT"_candidates_ALL_PAIRS.bam -o "$1"_"$REPEAT"_candidates_PROPERLY_PAIRED.bam & 
samtools view -h -b -F 2 "$1"_"$REPEAT"_candidates_ALL_PAIRS.bam -o "$1"_"$REPEAT"_candidates_NOT_PROPERLY_PAIRED.bam &
wait 

# index bam files
samtools index "$1"_"$REPEAT"_candidates_ALL_PAIRS.bam &
samtools index "$1"_"$REPEAT"_candidates_PROPERLY_PAIRED.bam & 
samtools index "$1"_"$REPEAT"_candidates_NOT_PROPERLY_PAIRED.bam &
wait 

# convert bam to sam
samtools view -O SAM "$1"_"$REPEAT"_candidates_PROPERLY_PAIRED.bam -o "$1"_"$REPEAT"_candidates_PROPERLY_PAIRED.sam &
samtools view -O SAM "$1"_"$REPEAT"_candidates_NOT_PROPERLY_PAIRED.bam -o "$1"_"$REPEAT"_candidates_NOT_PROPERLY_PAIRED.sam &
wait 

# compile list of read names and strip read 1 identifier
cut -f 1 "$1"_"$REPEAT"_candidates_PROPERLY_PAIRED.sam | sort | uniq > "$1"_"$REPEAT"_candidates_PROPERLY_PAIRED_read_names_.1 &
cut -f 1 "$1"_"$REPEAT"_candidates_NOT_PROPERLY_PAIRED.sam | sort | uniq > "$1"_"$REPEAT"_candidates_NOT_PROPERLY_PAIRED_read_names_.1 &
wait

# create read pair .2 identifiers
sed 's/.1$/.2/' "$1"_"$REPEAT"_candidates_PROPERLY_PAIRED_read_names_.1 > "$1"_"$REPEAT"_candidates_PROPERLY_PAIRED_read_names_.2 &
sed 's/.1$/.2/' "$1"_"$REPEAT"_candidates_NOT_PROPERLY_PAIRED_read_names_.1 > "$1"_"$REPEAT"_candidates_NOT_PROPERLY_PAIRED_read_names_.2 &
wait

cat "$1"_"$REPEAT"_candidates_PROPERLY_PAIRED_read_names_.1 "$1"_"$REPEAT"_candidates_PROPERLY_PAIRED_read_names_.2 | sort > "$1"_"$REPEAT"_candidates_PROPERLY_PAIRED_read_names_.interleaved &
cat "$1"_"$REPEAT"_candidates_NOT_PROPERLY_PAIRED_read_names_.1 "$1"_"$REPEAT"_candidates_NOT_PROPERLY_PAIRED_read_names_.2 | sort > "$1"_"$REPEAT"_candidates_NOT_PROPERLY_PAIRED_read_names_.interleaved &
wait

# create an interleaved fastq of properly-paired and not properly-paired mapped repeat candidates
/n/groups/kennedy/pagano/modules/bbmap/filterbyname.sh \
in=./../"$1"_1.fastq \
in2=./../"$1"_2.fastq \
out="$1"_"$REPEAT"_candidates_PROPERLY_PAIRED_interleaved.fastq \
names="$1"_"$REPEAT"_candidates_PROPERLY_PAIRED_read_names_.interleaved \
include=t &
/n/groups/kennedy/pagano/modules/bbmap/filterbyname.sh \
in=./../"$1"_1.fastq \
in2=./../"$1"_2.fastq \
out="$1"_"$REPEAT"_candidates_NOT_PROPERLY_PAIRED_interleaved.fastq \
names="$1"_"$REPEAT"_candidates_NOT_PROPERLY_PAIRED_read_names_.interleaved \
include=t &
wait

# determine number of properly-paired and not properly-paired reads that mapped
PROPERLY_PAIRED=$(cat "$1"_"$REPEAT"_candidates_PROPERLY_PAIRED_read_names_.1 | wc -l)
NOT_PROPERLY_PAIRED=$(cat "$1"_"$REPEAT"_candidates_NOT_PROPERLY_PAIRED_read_names_.1 | wc -l)
NON_TEMPLATED=$(echo $(echo "$PROPERLY_PAIRED")+$(echo "$NOT_PROPERLY_PAIRED") | bc)


##########################################
#                                        #
#    MERGE DATA AND RUN FEATURECOUNTS    #
#                                        #
##########################################

# combine outputs from unmapped and mapped analyses
# fastq files
cat "$1"_"$REPEAT"_candidates_READ1_MAPPED_READ2_UNMAPPED_interleaved.fastq "$1"_"$REPEAT"_candidates_PROPERLY_PAIRED_interleaved.fastq "$1"_"$REPEAT"_candidates_NOT_PROPERLY_PAIRED_interleaved.fastq > "$1"_"$REPEAT"_candidates_ALL_interleaved.fastq
# alignments
samtools merge -O BAM "$1"_"$REPEAT"_candidates_ALL.bam "$1"_"$REPEAT"_candidates_READ1_MAPPED_READ2_UNMAPPED.bam "$1"_"$REPEAT"_candidates_PROPERLY_PAIRED.bam "$1"_"$REPEAT"_candidates_NOT_PROPERLY_PAIRED.bam
samtools index "$1"_"$REPEAT"_candidates_ALL.bam
samtools view --threads 3 -O SAM "$1"_"$REPEAT"_candidates_ALL.bam -o "$1"_"$REPEAT"_candidates_ALL.sam

# count number of aligned reads per gene
/n/groups/kennedy/pagano/modules/subread-2.0.0/bin/featureCounts \
-p \
-t gene \
-M \
-O \
--fraction \
-a "$FEATURE_COUNTS_ANNO" \
-o "$1"_"$REPEAT"_candidates_counts \
"$1"_"$REPEAT"_candidates_PROPERLY_PAIRED.bam "$1"_"$REPEAT"_candidates_NOT_PROPERLY_PAIRED.bam "$1"_"$REPEAT"_candidates_READ1_MAPPED_READ2_UNMAPPED.bam

# remove features with 0 counts
(head -n 2 "$1"_"$REPEAT"_candidates_counts && awk '$7 > 0 || $8 > 0 || $9 > 0' "$1"_"$REPEAT"_candidates_counts | grep -v ^Geneid) > "$1"_"$REPEAT"_candidates_counts_hits
mv "$1"_"$REPEAT"_candidates_counts_hits "$1"_"$REPEAT"_candidates_counts

# remove chr, position, and strand information for TEs
awk 'BEGIN {FS="\t"; OFS="\t"} {if($2 ~ /;/) {print $1 "\t" "." "\t" "." "\t" "." "\t" "." "\t" $6 "\t" $7 "\t" $8 "\t" $9} else print $0}' "$1"_"$REPEAT"_candidates_counts > temp
mv temp "$1"_"$REPEAT"_candidates_counts
mv "$1"_"$REPEAT"_candidates_counts "$1"_"$REPEAT"_candidates_counts.txt


##########################################
#                                        #
#   PRINTING SUMMARY FILE AND CLEAN-UP   #
#                                        #
##########################################

# print summary file
> "$1"_TAIL-seq_summary.txt
printf -- "%s\n" "TAIL-SEQ LOG FILE" >> "$1"_TAIL-seq_summary.txt
printf -- "%s\n" "[`date`] " "" >> "$1"_TAIL-seq_summary.txt

printf -- "%s\n" "==== COMMAND-LINE-INPUTS ====" >> "$1"_TAIL-seq_summary.txt
printf -- "%s\n" "Sample: "$1"" >> "$1"_TAIL-seq_summary.txt
printf -- "%s\n" "Species: "$2"" >> "$1"_TAIL-seq_summary.txt
printf -- "%s\n" "Repeat: "$3" " "" >> "$1"_TAIL-seq_summary.txt

printf -- "%s\n" "======== STATISTICS =========" >> "$1"_TAIL-seq_summary.txt
printf -- "%s\n" "Total reads processed: $(printf "%'.f\n" $TOTAL_PE_READS)" >> "$1"_TAIL-seq_summary.txt
printf -- "%s\n" "3'-end reads with "$REPEAT": $(printf "%'.f\n" $READS_WITH_REPEAT)" >> "$1"_TAIL-seq_summary.txt
printf -- "%s\t" "" ""-MAPPED" "3\'-ENDS:" "$(printf "%'.f\n" $MAPPED_CANDS)"" >> "$1"_TAIL-seq_summary.txt
printf -- "%s\n" " " >> "$1"_TAIL-seq_summary.txt
printf -- "%s\t" "" "" ""-NON-TEMPLATED:" "$(printf "%'.f\n" $NON_TEMPLATED)"" >> "$1"_TAIL-seq_summary.txt
printf -- "%s\n" " " >> "$1"_TAIL-seq_summary.txt
printf -- "%s\t" "" "" "" ""-PROPERLY" "PAIRED:" "$(printf "%'.f\n" $PROPERLY_PAIRED)"" >> "$1"_TAIL-seq_summary.txt
printf -- "%s\n" " " >> "$1"_TAIL-seq_summary.txt
printf -- "%s\t" "" "" "" ""-NOT" "PROPERLY" "PAIRED:" "$(printf "%'.f\n" $NOT_PROPERLY_PAIRED)"" >> "$1"_TAIL-seq_summary.txt
printf -- "%s\n" " " >> "$1"_TAIL-seq_summary.txt
printf -- "%s\t" "" ""-UNMAPPED" "3\'-ENDS:" "$(printf "%'.f\n" $UNMAPPED_CANDS)"" >> "$1"_TAIL-seq_summary.txt
printf -- "%s\n" " " >> "$1"_TAIL-seq_summary.txt
printf -- "%s\t" "" "" ""-MATE" "MAPPED" "AND" "NON-TEMPLATED:" "$(printf "%'.f\n" $READ1s_MAPPED_READ2s_UNMAPPED)"" >> "$1"_TAIL-seq_summary.txt
printf -- "%s\n" " " >> "$1"_TAIL-seq_summary.txt
printf -- "%s\t" "" "" ""-MATE" "UNMAPPED:" "$(printf "%'.f\n" $READ1s_UNMAPPED_READ2s_UNMAPPED)"" >> "$1"_TAIL-seq_summary.txt

# make results directory and move final files into it
mkdir results
mv "$1"_"$REPEAT"_candidates_READ1_MAPPED_READ2_UNMAPPED_interleaved.fastq ./results
mv "$1"_"$REPEAT"_candidates_READ1_MAPPED_READ2_UNMAPPED.sam ./results
mv "$1"_"$REPEAT"_candidates_READ1_MAPPED_READ2_UNMAPPED.bam ./results
mv "$1"_"$REPEAT"_candidates_READ1_MAPPED_READ2_UNMAPPED.bam.bai ./results
mv "$1"_"$REPEAT"_candidates_READ1_UNMAPPED_READ2_UNMAPPED_interleaved.fastq ./results
mv "$1"_best_blastn_target_for_UNMAPPED_read_1_and_read_2_candidates.txt ./results
mv "$1"_"$REPEAT"_candidates_PROPERLY_PAIRED_interleaved.fastq ./results
mv "$1"_"$REPEAT"_candidates_PROPERLY_PAIRED.sam ./results
mv "$1"_"$REPEAT"_candidates_PROPERLY_PAIRED.bam ./results
mv "$1"_"$REPEAT"_candidates_PROPERLY_PAIRED.bam.bai ./results
mv "$1"_"$REPEAT"_candidates_NOT_PROPERLY_PAIRED_interleaved.fastq ./results
mv "$1"_"$REPEAT"_candidates_NOT_PROPERLY_PAIRED.sam  ./results
mv "$1"_"$REPEAT"_candidates_NOT_PROPERLY_PAIRED.bam ./results
mv "$1"_"$REPEAT"_candidates_NOT_PROPERLY_PAIRED.bam.bai ./results
mv "$1"_"$REPEAT"_candidates_ALL_interleaved.fastq ./results
mv "$1"_"$REPEAT"_candidates_ALL.sam ./results
mv "$1"_"$REPEAT"_candidates_ALL.bam ./results
mv "$1"_"$REPEAT"_candidates_ALL.bam.bai ./results
mv "$1"_"$REPEAT"_candidates_counts.txt ./results
mv "$1"_TAIL-seq_summary.txt ./results

#########