#!/bin/bash
set -e -x -o pipefail

main() {

    # Download the inputs to the worker's local storage
    mark-section "Downloading input files"
    dx-download-all-inputs

    # Extract the gzipped fastq files and adjust the variables associated with
    # them
    gunzip "$mate1_fastqgz_path"
    mate1_fastq_path="${mate1_fastqgz_path%.gz}"
    mate1_fastq_name="${mate1_fastqgz_name%.gz}"

    gunzip "$mate2_fastqgz_path"
    mate2_fastq_path="${mate2_fastqgz_path%.gz}"
    mate2_fastq_name="${mate2_fastqgz_name%.gz}"

    # Extract the tarball containing the HISAT2 reference
    mark-section "Extracting reference tarball"
    tar xf "$hisat2_index_targz_path"

    # Get the index basename to use when calling hisat2
    index_basename=${hisat2_index_targz_name%.tar.gz}/genome

    # Run hisat2
    output_sam_name="${mate1_fastq_prefix}.sam"
    output_bam_name="${mate1_fastq_prefix}.bam"
    mark-section "Running HISAT2"
    hisat2 --dta -x "$index_basename" -1 mate1.fastq -2 mate2.fastq -S "$output_sam_name"
    mkdir -p out/aligned_bam
    samtools view -b "$output_sam_name" -o out/aligned_bam/"$output_bam_name"

    # Upload the resulting file and associate it with the aligned_sam output
    mark-section "Uploading BAM results"
    dx-upload-all-outputs
    mark-success
}
