#!/usr/bin/env python3

import click
import os
import subprocess
from pathlib import Path
from contextlib import contextmanager
import tempfile
import shutil
import re
import logging
import timeit

logger = logging.getLogger(__name__)

@contextmanager
def cwd(path):
    '''
    # Usage:
    from mycontext import cwd
    with cwd(other_dir):
        # do something in the other_dir
 
 
    # Demo:
    from mycontext import cwd
    import os
 
    logger(os.getcwd())
    with cwd('ansible'):
        logger(os.getcwd())
    logger(os.getcwd())
    '''
 
    oldpwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(oldpwd)

@contextmanager
def tmpdir():
    '''
    https://code-maven.com/python-context-tools
    # Usage:
    from  import tmpdir
    with tmpdir() as dd:
        # store files in dd
        # the whole directory will be gone when the 'with' statement ends
 
 
    # Demo
    from mycontext import tmpdir
    import os
 
    with tmpdir() as d:
        logger(d)
        os.system("touch " + d + "/abc")
        os.system("ls -l /tmp/")
        os.system("ls -l " + d)
    '''

    dd = tempfile.mkdtemp()
    try:
        yield dd
    finally:
        shutil.rmtree(dd)
# bwa index ${contigs}

# bwa mem -M -t ${task.cpus} ${contigs} ${short_reads} | samtools view -bS -| samtools sort > alignments.bam

# samtools index alignments.bam

# pilon --genome ${contigs} --frags alignments.bam --changes \
# --output ${sample_id}_pilon --fix all


def execute(cmd):
    logger.info(f"Executing command: {cmd}")
    subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)


def build_index():
    logger.info(f"Building index...")
    cmd = f"bwa index contigs.fasta"
    execute(cmd)


def map(threads):
    logger.info("Mapping....")
    aligment_bam = "mapped.sorted.bam"
    cmd = f"bwa mem -M -t {threads} contigs.fasta R1.fq.gz R2.fq.gz | samtools sort -@ {threads} -O bam -o mapped.sorted.bam"
    execute(cmd)
    execute("samtools index mapped.sorted.bam")
    return(aligment_bam)


def pilon(aligment, prefix = None):
    polished_asemblies = f"{prefix}_pilon"
    changes = f"{polished_asemblies}.changes"
    cmd = f"_JAVA_OPTIONS=-Xmx24g pilon --genome contigs.fasta --frags {aligment} --changes --output {polished_asemblies} --fix bases"
    execute(cmd)
    return(f"{polished_asemblies}.fasta", changes)


def fix_header():
    # file_name = fasta.split(".")[0]
    pat = re.compile("_pilon")
    with open("final.polished.fasta", "w") as out_fasta:
        with open("polished_pilon.fasta", "r") as in_fasta:
            for line in in_fasta:
                if line.startswith(">"):
                    fixed = pat.sub("", line)
                    out_fasta.write(fixed)
                else:
                    out_fasta.write(line)


@click.command()
@click.option("-t", "--threads", default=8, help="number of threads", show_default=True)
@click.option("-n", "--iteration", default=10, help="number of iteration", show_default=True)
@click.option("-p", "--prefix", default="polished", help="prefix", show_default=True)
# @click.option("--keep_tmp",is_flag=True, help="prefix", show_default=True)
@click.argument("R1")
@click.argument("R2")
@click.argument("contigs")
def polish(r1, r2, contigs, iteration, threads, prefix):
    """
    Polishing assemblies using Pilon
    """
    _start = timeit.default_timer()
    assert Path(r1).exists(), f"{r1} is not valid"
    assert Path(r2).exists(), f"{r2} is not valid"
    assert Path(contigs).exists(), f"{contigs} is not valid"

    # with cwd(tmpdir):
    _cwd = os.getcwd()
    tmpdir = Path(f"{_cwd}/_tmp")
    
    if tmpdir.exists():
        shutil.rmtree(str(tmpdir))
    
    tmpdir.mkdir()

    r1_path = Path(r1).absolute()
    r2_path = Path(r2).absolute()
    contigs_path = Path(contigs).absolute()
    logger.info(
        f"Paths:\nR1: {r1_path}\nR2: {r2_path}\nContigs: {contigs_path}")
    
    with cwd(str(tmpdir)):
        
        logger.info(f"Current working dir: {_cwd}")

        for _file in ["R1.fq.gz", "R2.fq.gz", "contigs.fasta"]:
            _file_path = Path(_file)
            if _file_path.exists():
                _file_path.unlink()
        
        Path("R1.fq.gz").symlink_to(r1_path)
        Path("R2.fq.gz").symlink_to(r2_path)


        j = 0
        changes = None
        for i in range(iteration):
            logger.info(f"Polishing round {i+1}:")
            if Path('contigs.fasta').exists() and i > 0:
                Path('contigs.fasta').unlink()
            Path('contigs.fasta').symlink_to(contigs_path)
            build_index()
            aligment = map(threads)
            contigs_path, changes = pilon(aligment, prefix)
            logger.info(f"Contigs: {contigs_path}")
            changes_content = Path(changes).read_text()
            if len(changes_content) == 0:
                logger.info(f"No more changes. Polishing should be stopped now.")
                fix_header()
                shutil.move("final.polished.fasta", _cwd)
                break
            else:
                changes_count = len(changes_content.split('\n'))
                logger.info(f"Changes left: {changes_count}")
                j += 1
        if (j==iteration) and changes is not None:
            shutil.move(changes,_cwd)

    _stop = timeit.default_timer()
    # shutil.rmtree(tmpdir)
    duration = _stop - _start
    logger.info(f"Duration: {duration}")


if __name__ == "__main__":
    logging.basicConfig(format='%(asctime)s - %(message)s',
                        datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)
    polish()
