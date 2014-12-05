import os
import hashlib

def clean_up_dir(dirname):
    """Remove all the fastq files and intermediate bams/headers"""

    for filename in os.listdir(dirname):
        if(filename.endswith(".fq") 
            or filename.endswith(".paired.bam")
            or filename.endswith(".sam")):
            filepath = os.path.join(dirname, filename)
            os.remove(filepath)
            assert(not(os.path.exists(filepath)))

def make_new_dir(path):
    """Create a new directory if it doesn't exist"""

    if(not(os.path.isdir(path))):
        os.makedirs(path)

    return path

def calc_md5sum(bam_filename):
    """What it says"""

    md5worker = hashlib.md5()
    bam_file = open(bam_filename)
    block_size = 128
    while True:
        data = bam_file.read(block_size)
        if not data:
            break
        md5worker.update(data)

    return md5worker.hexdigest()
