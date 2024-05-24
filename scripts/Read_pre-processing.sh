#########################################################################
#########################################################################
#########################################################################
#########################################################################

# Code for the pre-processing of sequencing reads from a Surface-Plex assay
# Including the following steps:
# 1. Generate an index with STAR from a custom .fasta file including specific antigen barcodes
# 2. Trim reads using cutadapt
# 3. Extract Unique Molecular Identifiers (UMIs) using UMI-tools
# 4. Align reads to the generated index using STAR
# 5. Remove samples that do not exist based on a list (optional)
# 6. Convert SAM files to BAM format using SAMtools
# 7. Sort and index BAM files using SAMtools
# 8. Deduplicate reads using UMIs
# 9. Convert deduplicated BAM files to text format and count unique reads

#########################################################################
#########################################################################
#########################################################################
#########################################################################

# 1. Generate an index with STAR from custom .fasta file including specific antigen barcodes

index="path/to/index/dir"

mkdir -p $index

ml STAR/2.7.9a

# Generate the genome index using STAR
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $index --genomeFastaFiles $index/antigen_bc.fasta --genomeSAindexNbases 2

# 2. Trim reads using cutadapt

fastqs="path/to/fastq/dir"

cd fastqs

ml cutadapt/4.1

trimmed_output="path/to/trimmed_output/dir"

mkdir -p $trimmed_output

# Loop through each fastq file and trim reads to just the antigen barcode + UMI
for fq in $fastqs/*/*; do
    dir_name=$(basename $(dirname $fq))
    base_name=$(basename $fq)
    mkdir -p $trimmed_output/$dir_name
    cutadapt -j 8 -u 31 -u -3 -o $trimmed_output/$dir_name/$base_name $fq
done

# 3. Extract Unique Molecular Identifiers (UMIs) using UMI-tools

ml UMI-tools/1.0.1

umi_output="path/to/umi_output/dir"

mkdir -p $umi_output

# Loop through each trimmed fastq file and extract UMIs
for fq in $trimmed_output/*/*; do
    dir_name=$(basename $(dirname $fq))
    base_name=$(basename $fq)
    mkdir -p $umi_output/$dir_name
    umi_tools extract -I $fq --bc-pattern=NNNNNNNNNN --3prime --extract-method=string --stdout $umi_output/$dir_name/$base_name
done

# 4. Align reads to the generated index using STAR

output="path/to/aligned_output/dir"

mkdir -p $output

ml STAR/2.7.9a

# Loop through each UMI-extracted fastq file and align reads
for fq in $umi_output/*/*.fastq.gz; do
    dir_name=$(basename $(dirname $fq))
    base_name=$(basename $fq)
    STAR --runThreadN 8 --readFilesIn $fq --readFilesCommand gunzip -c --genomeDir $index --outFileNamePrefix $output/$dir_name/$base_name
done

# 5. Remove samples that do not exist based on a list of well barcodes

cd $output/aligned/dir

awk -F'\t' '{print $1}' path/to/Remove_barcodes.txt | while read -r barcode; do
  rm -f ${barcode}_*.sam
done

# 6. Convert SAM files to BAM format using SAMtools

bam_output="path/to/mapped_bam_output/dir"

mkdir -p $bam_output

ml SAMtools/1.16.1

cd $output

dir_list=("1C" "1D" "2C" "2D" "3C" "3D" "4C" "4D")

# Loop through each directory and convert SAM to BAM
for dir in "${dir_list[@]}"; do
  if [ -d "$dir" ]; then
    for i in $dir/*.sam; do
      dir_name=$(basename $(dirname $i))
      base_name=$(basename $i)
      mkdir -p $bam_output/$dir_name
      out_bam=$(echo $base_name | sed 's/.fastq.gzAligned.out.sam/.bam/')
      samtools view -@ 8 -b $i > $bam_output/$dir_name/$out_bam
    done
  else
    echo "Directory $dir does not exist. Skipping."
  fi
done

# 7. Sort and index BAM files using SAMtools

sorted_bam_output="path/to/sorted_bam_output/dir"

mkdir -p $sorted_bam_output

cd $bam_output
for i in */*; do
    dir_name=$(basename $(dirname $i))
    base_name=$(basename $i)
    mkdir -p $sorted_bam_output/$dir_name
    samtools sort -@ 8 $i -o $sorted_bam_output/$dir_name/$base_name
done

cd $sorted_bam_output
for i in */*; do
    samtools index -@ 8 $i
done

# 8. Deduplicate reads using UMIs

deduplicated_output="path/to/deduplicated_bam_output/dir"

mkdir -p $deduplicated_output

# Loop through each sorted BAM file and deduplicate reads
for i in $sorted_bam_output/*/*.bam; do
    dir_name=$(basename $(dirname $i))
    base_name=$(basename $i)
    mkdir -p $deduplicated_output/$dir_name
    umi_tools dedup -I $i --method=unique -S $deduplicated_output/$dir_name/$base_name
done

# 9. Convert deduplicated BAM files to text format and count unique reads

cd $deduplicated_output
for i in */*.bam; do
    dir_name=$(basename $(dirname $i))
    base_name=$(basename $i)
    out_txt=$(echo $base_name | sed 's/.bam/.txt/')
    samtools view -h $i > $dir_name/$out_txt
done

# Loop through each text file and count unique reads
for file in */*.txt; do
    
    dir_name=$(basename $(dirname $file))
    base_name=$(basename $file)
    new_name="${dir_name}_$(echo $base_name | sed 's/_R1.txt/.txt/')"

    # Process the input file and save to the new output file
    cat $file | cut -f3 | uniq -c | column -t | tail -n +6 | awk '{gsub(/ +/, "\t")}1' > $new_name
done
