import os
import subprocess

import dxpy

@dxpy.entry_point("main")
def main(hisat2_index_targz, mate1_fastq, mate2_fastq):

    # First, download all the input files to local storage
    dxpy.download_dxfile(hisat2_index_targz, "hisat2_index.tar.gz")
    dxpy.download_dxfile(mate1_fastq, "mate1.fastq")
    dxpy.download_dxfile(mate2_fastq, "mate2.fastq")


    # Second, extract the index tarball
    proc = subprocess.Popen(["tar", "xf", "hisat2_index.tar.gz"])
    proc.wait()

    # You could also use the tarfile module for this
    # import tarfile
    # tar = tarfile.open("hisat2_index.tar.gz")
    # tar.extractall()
    # tar.close()

    # Third, figure out what the basename of the hisat2 reference index is
    # This depends on the basename following the pattern used in the indexes
    # distributed by the authors of HISAT2, that is the index in grch37.tar.gz
    # will extract to grch37/genome*
    index_basename = os.path.join(dxpy.DXFile(hisat2_index_targz).name[:-len(".tar.gz")],
                                  "genome")

    # Prepare the hisat2 command and run it.
    hisat2_cmd_template = ("hisat2 --dta -x {index_basename} -1 {mate1_fastq} "
                           "-2 {mate2_fastq} -S {hisat2_output_sam}")
    hisat2_cmd = hisat2_cmd_template.format(
        index_basename=index_basename,
        mate1_fastq="mate1.fastq",
        mate2_fastq="mate2.fastq",
        hisat2_output_sam="hisat2_output.sam")
    subprocess.check_call(hisat2_cmd, shell=True)

    # Upload the output SAM file.
    uploaded_dxfile = dxpy.upload_local_file("hisat2_output.sam")

    # Return the ID of the uploaded SAM file associated with the "aligned_sam"
    # field in the outputSpec in dxapp.json.
    return {"aligned_sam": dxpy.dxlink(uploaded_dxfile.get_id())}
