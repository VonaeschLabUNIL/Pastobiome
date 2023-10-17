#Pipeline for JOHI metagenomic pre-processing and processing
# Include cleaning with fastP, mOTUs2, MetaPhlan, Humann3 and RGI loops

#FatsP cleaning
mkdir CLEANED_READS
mkdir REPORTS_CLEANED_READS

# Install fast
# conda install -c bioconda fastp

# Paired-end:
for prefix in $(ls *.fastq.gz | sed -E 's/_R[12]_001[.]fastq.gz//' | uniq)
do
echo "performing cleaning on files :" "${prefix}_R1_001.fastq.gz" "${prefix}_R2_001.fastq.gz"
fastp -i "${prefix}_R1_001.fastq.gz" -I "${prefix}_R2_001.fastq.gz" -o CLEANED_READS/"${prefix}_FASTP_R1_001.fastq.gz"  -O CLEANED_READS/"${prefix}_FASTP_R2_001.fastq.gz" --report_title "${prefix}_fastp_report" --thread 38 -j REPORTS_CLEANED_READS/"${prefix}_fastp".json -h REPORTS_CLEANED_READS/"${prefix}_fastp".html
done

multiqc ./REPORTS_CLEANED_READS/ --ignore-symlinks --outdir ./REPORTS_CLEANED_READS --filename MULTIQC_ALL_SAMPLE_REPORT --fullnames --title MULTIQC_ALL_SAMPLE_REPORT

cd ./CLEANED_READS

mkdir MOTUS_ANALYSIS
mkdir MOTUS_ANALYSIS/MOTUS_MERGE
for level in $(cat LEVEL_LIST.txt)
do
  mkdir MOTUS_ANALYSIS/MOTUS_${level}_ABOND
  mkdir MOTUS_ANALYSIS/MOTUS_${level}_READS


  for prefix in $(ls *.fastq.gz | sed -E 's/_R[12]_001[.]fastq.gz//' | uniq)
  do
   echo "Analyses abondance du" "${prefix}" "au niveau" "${level}"
    motus profile -s ${prefix}_R1_001.fastq.gz -t 38 -o MOTUS_ANALYSIS/MOTUS_${level}_ABOND/${prefix}_${level}_ABOND.motus -n ${prefix}_ABOND -k ${level} -p -q
    echo "Comptage de reads du" "${prefix}" "au niveau" "${level}"
    motus profile -s ${prefix}_R1_001.fastq.gz -t 12 -o MOTUS_ANALYSIS/MOTUS_${level}_READS/${prefix}_${level}_READS.motus -n ${prefix}_READS -k ${level} -p -q -c 
    echo "Tes fichiers récapitulatif sont prêts"
  done
  motus merge -d MOTUS_ANALYSIS/MOTUS_${level}_ABOND -o MOTUS_ANALYSIS/MOTUS_MERGE/${level}_MERGE_ABONDsim.csv
  motus merge -d MOTUS_ANALYSIS/MOTUS_${level}_READS -o MOTUS_ANALYSIS/MOTUS_MERGE/${level}_MERGE_READS.motus
done

#MetaPhlAn
mkdir METAPHLAN_ANALYSIS/
mkdir METAPHLAN_ANALYSIS/PROFILES_STANDARD/
mkdir METAPHLAN_ANALYSIS/PROFILES_BIOM/
mkdir METAPHLAN_ANALYSIS/BOWTIE_OUT/
for prefix in $(ls *.fastq.gz | sed -E 's/_FASTP_R[12]_001[.]fastq.gz//' | uniq)
do
    echo "Performing metaphlan3 analysis on sample": "${prefix}_FASTP_R1_001.fastq.gz"
    metaphlan "${prefix}_FASTP_R1_001.fastq.gz" --nproc 26 --add_viruses --input_type fastq -o METAPHLAN_ANALYSIS/${file}_metaphlan3_abundance.txt --no_map
done

#To merge the profiled samples to put them all together :
merge_metaphlan_tables.py METAPHLAN_ANALYSIS/*_metaphlan3_abundance.txt > METAPHLAN_ANALYSIS/JOHI_merged_abundance_table.txt

#Humann3
for prefix in $(ls *.fastq.gz | sed -E 's/_FASTP_R[12]_001[.]fastq.gz//' | uniq)
do
	echo "Humann analysis on sample": "${prefix}_FASTP_R1_001.fastq.gz"
	humann -i "${prefix}_FASTP_R1_001.fastq.gz" --threads 26 --output-basename "${prefix}" -o HUMANN_OUTPUT_JOHI_SUBSET --output-max-decimals 2 --remove-temp-output --o-log HUMANN_OUTPUT_JOHI_SUBSET/"${prefix}".log --taxonomic-profile JOHI_merged_abundance_table.txt
done

#combining results
humann_join_tables -i hmp_subset -o hmp_subset_genefamilies.tsv --file_name genefamilies

# Normalizing to control for different sequencing depth across the samples.
humann_renorm_table -i hmp_subset_genefamilies.tsv -o hmp_subset_genefamilies-cpm.tsv --units cpm

mkdir Barplots
mkdir Barplots/log_scale
#mkdir Barplots/linear_scale

for level in $(cat Pathways_humann.txt)
do
   echo "Barplot of pathway" "${level}"
   humann_barplot --input path_abundance_country.txt --focal-metadata Country --last-metadata Country --output Barplots/log_scale/plot_log_${level}.png --focal-feature ${level} --sort sum metadata --scaling logstack
   echo "${level}" "done log scale"
   #humann_barplot --input path_abundance_country.txt --focal-metadata Country --last-metadata Country --output Barplots/linear_scale/plot_lin_${level}.png --focal-feature ${level} --sort sum metadata
   #echo  "${level}" "done linear scale"
done

# RGI loop for metagenomic short reads or genomic short reads

mkdir RGI_ANALYSIS

echo "========Loading CARD reference data======="
rgi clean --local
wget https://card.mcmaster.ca/latest/data
tar -xvf data ./card.json

echo "======Loading into local directory======"
rgi load --card_json ./card.json --local

echo "======CARD Canonical annotations======"
# Adjust version number
rgi card_annotation -i ./card.json > card_annotation.log 2>&1

ech "======Versions======"
rm card_database_v*_all.fasta
data_version=`echo card_database_v*.fasta | sed 's/.*card_database_v\(.*\).fasta/\1/'`
echo "$cmd data_version: $data_version"

echo "======Load databse======"
rgi load -i ./card.json --card_annotation card_database_v${data_version}.fasta --local

echo "====== Loop ======"

for prefix in $(ls *.fastq.gz | sed -E 's/_FASTP_R[12]_001[.]fastq.gz//' | uniq)
  do
    echo "======Analyze RGI of " "${prefix}" "======"
    rgi bwt -1 ${prefix}_FASTP_R1_001.fastq.gz -2 ${prefix}_FASTP_R2_001.fastq.gz -n 36 -o RGI_ANALYSIS/${prefix}_output --clean --local
    echo "======" "${prefix}" " is done======"
  done
