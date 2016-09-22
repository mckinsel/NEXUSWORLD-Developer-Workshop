import os
import subprocess
import traceback

import dxpy

@dxpy.entry_point("main")
def main(hisat2_index_targz, mate1_fastq, mate2_fastq):

    # First, download all the input files to local storage
    inputs_dict = dxpy.download_all_inputs()

    # Second, unzip the inputs fastqs and update paths and names
    subprocess.check_call(["gunzip", inputs_dict["mate1_fastqgz_path"]])
    inputs_dict["mate1_fastq_path"] = inputs_dict["mate1_fastqgz_path"].replace("*.gz", "")
    inputs_dict["mate1_fastq_name"] = inputs_dict["mate1_fastqgz_name"].replace("*.gz", "")

    subprocess.check_call(["gunzip", inputs_dict["mate2_fastqgz_path"]])
    inputs_dict["mate2_fastq_path"] = inputs_dict["mate2_fastqgz_path"].replace("*.gz", "")
    inputs_dict["mate2_fastq_name"] = inputs_dict["mate2_fastqgz_name"].replace("*.gz", "")

    # Third, extract the index tarball
    subprocess.check_call(["tar", "xf", inputs_dict["hisat2_index_targz_path"]])

    # Fourth, figure out what the basename of the hisat2 reference index is
    # This depends on the basename following the pattern used in the indexes
    # distributed by the authors of HISAT2, that is the index in grch37.tar.gz
    # will extract to grch37/genome*
    index_basename = os.path.join(inputs_dict["hisat2_index_targz_name"][:-len(".tar.gz")],
                                  "genome")

    # Prepare the hisat2 command and run it.
    output_sam_name = inputs_dict["mate1_fastqgz_prefix"] + ".sam"
    output_bam_name = inputs_dict["mate1_fastqgz_prefix"] + ".bam"
    hisat2_cmd_template = ("hisat2 --dta -x {index_basename} -1 {mate1_fastq} "
                           "-2 {mate2_fastq} -S {hisat2_output_sam}")
    hisat2_cmd = hisat2_cmd_template.format(
        index_basename=index_basename,
        mate1_fastq=inputs_dict["mate1_fastq_path"],
        mate2_fastq=inputs_dict["mate2_fastq_path"],
        hisat2_output_sam=output_sam_name)
    try:
        subprocess.check_call(hisat2_cmd, shell=True)
    except subprocess.CalledProcessError as exc:
        traceback.print_exc()
        raise dxpy.AppError(
            "Error while running HISAT2: {e}. Consult the log for more information".format(
                e=exc.message))
    
    # Convert the SAM to a BAM
    subprocess.check_call(["samtools", "view", "-b", output_sam_name, "-o", output_bam_name])

    # Upload the output SAM file.
    uploaded_dxfile = dxpy.upload_local_file("hisat2_output.sam")

    # Return the ID of the uploaded SAM file associated with the "aligned_sam"
    # field in the outputSpec in dxapp.json.
    return {"aligned_sam": dxpy.dxlink(uploaded_dxfile.get_id())}
