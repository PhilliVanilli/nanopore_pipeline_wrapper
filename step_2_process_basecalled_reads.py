import os
import sys
import subprocess
import argparse
import pathlib
import shutil
from itertools import groupby
import pandas as pd
from Bio import SeqIO


__author__ = 'Colin Anthony'


def try_except_exit_on_fail(cmd):
    try:
        subprocess.call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print(e)
        sys.exit("exiting")


def try_except_continue_on_fail(cmd):
    try:
        subprocess.call(cmd, shell=True)
        return True
    except subprocess.CalledProcessError as e:
        print(e)
        return False


def py3_fasta_iter(fasta_name):
    """
    modified from Brent Pedersen: https://www.biostars.org/p/710/#1412
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(str(fasta_name), 'r')
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header_str = header.__next__()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())
        yield (header_str, seq)


def gather_fastqs(fastq_path, run_name, max_len, min_len):

    fastq_outpath = fastq_path.parent
    fastq_name = f"{run_name}_all.fastq"
    output_fastq = pathlib.Path(fastq_outpath, fastq_name)

    all_fastqs = list(fastq_path.glob("*.fastq"))
    summary_file_location = pathlib.Path(fastq_path, "sequencing_summary.txt")

    with open(output_fastq, 'w') as handle:
        for fastq in all_fastqs:
            for record in SeqIO.parse(open(fastq), "fastq"):
                seq_len = len(record.seq)
                if seq_len > max_len or seq_len < min_len:
                    continue
                else:
                    SeqIO.write([record], handle, "fastq")

    if output_fastq.is_file():
        return summary_file_location
    else:
        return False


def main(project_path, sample_names, reference, make_index, ref_start, ref_end, min_len, max_len):
    # set the primer_scheme directory
    script_folder = pathlib.Path(__file__).absolute().parent
    primer_scheme_dir = pathlib.Path(script_folder, "primer-schemes")

    # get folder paths
    project_path = pathlib.Path(project_path).absolute()
    run_name = project_path.parts[-1]
    fast5_dir = pathlib.Path(project_path, "fast5")
    fastq_dir = pathlib.Path(project_path, "fastq")
    sample_names = pathlib.Path(sample_names).absolute()
    demultipled_folder = pathlib.Path(project_path, "demultiplexed")
    sample_folder = pathlib.Path(project_path, "samples")

    # set dir to project dir so that output is written in correct place by external tools
    os.chdir(project_path)

    # set the reference genome
    reference_scheme = \
        {"ChikAsianECSA_V1": pathlib.Path(primer_scheme_dir, "ChikAsianECSA", "V1", "ChikAsianECSA.reference.fasta"),
         "ZikaAsian_V1":  pathlib.Path(primer_scheme_dir, "ZikaAsian", "V1", "ZikaAsian.reference.fasta"),
         "ZaireEbola_V1": pathlib.Path(primer_scheme_dir, "ZaireEbola", "V2", "ZaireEbola.reference.fasta"),
         "ZaireEbola_V2": pathlib.Path(primer_scheme_dir, "ZaireEbola", "V2", "ZaireEbola.reference.fasta"),
         "LassaL_V1": pathlib.Path(primer_scheme_dir, "LassaL", "V1", "LassaL.reference.fasta"),
         "LassaS_V1": pathlib.Path(primer_scheme_dir, "LassaS", "V1", "LassaS.reference.fasta")}

    chosen_ref_scheme = str(reference_scheme[reference])
    chosen_ref_scheme_bed_file = chosen_ref_scheme.replace(".reference.fasta", ".scheme.bed")
    ref_name, ref_seq = list(py3_fasta_iter(chosen_ref_scheme))[0]
    ref_name = ref_name.split()[0]
    if not ref_start:
        ref_start = 1
    if not ref_end:
        ref_end = len(ref_seq)

    print(f"\nReference is {chosen_ref_scheme}\n")
    print(f"\nPrimer bed file is {chosen_ref_scheme_bed_file}\n")

    # gather all fastq's into one file and all sequencing_summary.txt's into one file
    print(f"\nrunning: collecting all fastq files into one file")
    summary_file_path = gather_fastqs(fastq_dir, run_name, max_len, min_len)
    if not summary_file_path:
        print("gathering fastq files failed, concatenated fastq file was not found")
        sys.exit("exiting")

    # demultiplex with porchop
    print(f"\nrunning: porechop demultiplexing")
    master_reads_file = pathlib.Path(project_path, run_name + "_all.fastq")
    if not master_reads_file.is_file():
        print(f"could not find the concatenated fastq, {master_reads_file}")
        sys.exit("exiting")
    demultiplex_cmd = f"porechop --format fastq --verbosity 2 -i {master_reads_file} " \
                      f"--discard_middle " \
                      f"--require_two_barcodes " \
                      f"--barcode_threshold 80 " \
                      f"--threads 4 " \
                      f"--check_reads 10000 " \
                      f"--barcode_diff 5 " \
                      f"--barcode_dir {demultipled_folder} " \
                      f"> {master_reads_file}.demultiplexreport.txt"

    try_except_exit_on_fail(demultiplex_cmd)

    # add run name to each demultiplexed file
    for file in list(demultipled_folder.glob("*.fastq")):
        path = file.parents[0]
        name = file.name
        if len(name) == 10:
            new_name = f"{run_name}_{name}"
            new_file = pathlib.Path(path, new_name)
            os.rename(str(file), new_file)

    # index concatenated fastq with nanopolish
    print(f"\nrunning: nanopolish index on fast5/fastq files")
    nanopolish_index_cmd = f"nanopolish index -s {summary_file_path} -d {fast5_dir} {master_reads_file}"
    try_except_exit_on_fail(nanopolish_index_cmd)

    # concatenated demultiplexed files for each sample and setup sample names and barcode combinations
    sample_names_df = pd.read_csv(sample_names, sep=None, engine="python")
    sample_names_df['barcode_1'] = sample_names_df['barcode_1'].apply(lambda x: run_name + "_" + x + ".fastq")
    sample_names_df['barcode_2'] = sample_names_df['barcode_2'].apply(lambda x: run_name + "_" + x + ".fastq")
    sample_names_dict = sample_names_df.set_index('sample_name').T.to_dict(orient='list')
    all_sample_files = []
    for sample_name, [barcode_1, barcode_2] in sample_names_dict.items():
        sample_dir = pathlib.Path(sample_folder, sample_name)
        if not sample_dir.exists():
            pathlib.Path(sample_dir).mkdir(mode=0o777, parents=True, exist_ok=True)
        barcode_1_file = pathlib.Path(demultipled_folder, barcode_1)
        barcode_2_file = pathlib.Path(demultipled_folder, barcode_2)
        cat_outfile = pathlib.Path(sample_dir, f"{sample_name}.fastq")
        all_sample_files.append(cat_outfile)
        cat_cmd = f"cat {str(barcode_1_file)} {str(barcode_2_file)} > {cat_outfile}"
        run = try_except_continue_on_fail(cat_cmd)
        if not run:
            print("missing one or more demultiplexed files for this sample")
            continue

    if make_index:
        make_index_cmd = f"bwa index {chosen_ref_scheme}"
        try_except_exit_on_fail(make_index_cmd)

    for sample_fastq in all_sample_files:
        if not sample_fastq.is_file():
            print(f"could not find the concatenated sample fastq file: {sample_fastq}\nskipping sample")
            continue
        sample_name = pathlib.Path(sample_fastq).stem
        sample_folder = pathlib.Path(sample_fastq).parent
        sam_name = pathlib.Path(sample_folder, sample_name + "_mapped.sam")
        bam_file = pathlib.Path(sample_folder, sample_name + "_mapped.bam")
        bam_file_sorted = pathlib.Path(sample_folder, sample_name + ".sorted.bam")
        trimmed_bam_file = pathlib.Path(sample_folder, sample_name + ".sorted.primerclipped.bam")
        rename_trimmed_bam_file = pathlib.Path(sample_folder, sample_name + ".primertrimmed.sorted.bam")
        vcf_file = pathlib.Path(sample_folder, sample_name + "_polished.vcf")
        nanopolish_cons_file = pathlib.Path(sample_folder, sample_name + "_consensus_nanopolish.fasta")
        bcftools_vcf_file = pathlib.Path(sample_folder, sample_name + "_bcftools.vcf")
        bcftools_cons_file = pathlib.Path(sample_folder, sample_name + "_consensus_bcftools.fasta")
        os.chdir(sample_folder)

        # run read mapping using bwa
        print(f"\nrunning: bwa read mapping")
        bwa_cmd = f"bwa mem -t 4 -x ont2d {chosen_ref_scheme} {sample_fastq} > {sam_name}"
        run = try_except_continue_on_fail(bwa_cmd)
        if not run:
            continue

        # convert sam to bam
        print(f"\nrunning: sam to bam conversion")
        sam_bam_cmd = f"samtools view -bS {sam_name} > {bam_file}"
        run = try_except_continue_on_fail(sam_bam_cmd)
        if not run:
            continue

        # sort bam file
        print(f"\nrunning: sorting bam file")
        sort_sam_cmd = f"samtools sort -T {sample_name} {bam_file} -o {bam_file_sorted}"
        run = try_except_continue_on_fail(sort_sam_cmd)
        if not run:
            continue

        # index bam file for bamclipper
        print(f"\nrunning: indexing bam file")
        index_bam_cmd = f"samtools index {bam_file_sorted}"
        run = try_except_continue_on_fail(index_bam_cmd)
        if not run:
            continue

        # remove primer sequences with bamclipper
        print(f"\nrunning: trim primer sequences from bam file")

        trim_primer = f"bamclipper.sh -b {bam_file_sorted} " \
            f"-p {chosen_ref_scheme_bed_file} -n 4 -u 1 -d 1"
        run = try_except_continue_on_fail(trim_primer)
        if not run:
            continue

        # run nanopolish
        print(f"\nrunning: nanopolish variant calling")
        nanopolish_cmd_v11 = f"nanopolish variants " \
            f"--fix-homopolymers --snps -o {vcf_file} -w '{ref_name}:{ref_start}-{ref_end}' -t 4 --ploidy=1 -v " \
            f"-r {master_reads_file} -b {trimmed_bam_file} -g {chosen_ref_scheme} --min-candidate-frequency=0.5" \
            f"--min-candidate-depth=10 --max-haplotypes=10000 "
        print(f'{ref_name}:{ref_start}-{ref_end}')
        print(nanopolish_cmd_v11)
        run = try_except_continue_on_fail(nanopolish_cmd_v11)
        if not run:
            continue

        # make nanopolish consensus
        print(f"\nrunning: making consensuses sequence from bam and variant call (vcf) files")
        consensus_cmd = f"nanopolish vcf2fasta --skip-checks -g {chosen_ref_scheme} {vcf_file} > {nanopolish_cons_file}"
        run = try_except_continue_on_fail(consensus_cmd)
        if not run:
            continue

        # make bcftools consensus
        min_base_qual = 30 # default=13
        p_val_of_variant = 0.2 # default=0.5

        bcf_vcf_cmd = f"bcftools mpileup --min-BQ {min_base_qual} -Ou -f {chosen_ref_scheme} {trimmed_bam_file} " \
            f"| bcftools call -c -p {p_val_of_variant} --ploidy 1 -mv -Oz -o {bcftools_vcf_file}"
        bcf_index_cmd = f"bcftools index {bcftools_vcf_file}"
        bcf_cons_cmd = f"bcftools consensus -H A -f {chosen_ref_scheme} {bcftools_vcf_file} > {bcftools_cons_file}"
        run = try_except_continue_on_fail(bcf_vcf_cmd)
        if not run:
            continue
        run = try_except_continue_on_fail(bcf_index_cmd)
        if not run:
            continue
        run = try_except_continue_on_fail(bcf_cons_cmd)
        if not run:
            continue

        # make artic-ebov consensus
        shutil.copyfile(trimmed_bam_file, rename_trimmed_bam_file)
        cons_file_script = pathlib.Path(script_folder, "margin_cons.py")

        set_min_depth = 10 # default=20
        set_min_qual = 30 # default=200

        cons_cmd = f"python {cons_file_script} -r {chosen_ref_scheme} -v {vcf_file} -b {rename_trimmed_bam_file} " \
            f"-n {sample_name} -d {set_min_depth} -q {set_min_qual} "
        run = try_except_continue_on_fail(cons_cmd)
        if not run:
            continue

        # plot depth and quality for sample
        plot_file_script = pathlib.Path(script_folder, "plot_depths_qual.py")
        plot_cmd = f"python {plot_file_script} -r {chosen_ref_scheme} -v {vcf_file} -b {rename_trimmed_bam_file} " \
            f"-n {sample_name}"
        run = try_except_continue_on_fail(plot_cmd)
        if not run:
            continue

        print(f"Completed processing sample: {sample_name}")

    print("sample processing completed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process raw nanopore reads to fasta consensus sequences",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-in", "--project_path", default=argparse.SUPPRESS, type=str,
                        help="The path to the directory containing the 'fast5' and 'fastq' folders ", required=True)
    parser.add_argument("-s", "--sample_names", default=argparse.SUPPRESS, type=str,
                        help="The csv file with the sample names and corresponding barcode combinations"
                             "This file should have three columns: barcode_1,barcode_2,sample_name", required=True)
    parser.add_argument("-r", "--reference", type=str, default="ChikAsianECSA_V1",
                        help="The reference genome and primer scheme to use",
                        choices=["ChikAsianECSA_V1", "ZikaAsian_V1", "ZaireEbola_V1", "ZaireEbola_V2", "LassaL_V1",
                                 "LassaS_V1"], required=False)
    parser.add_argument("-m", "--make_index", default=False, action="store_true",
                        help="index the reference if not done previously", required=False)
    parser.add_argument("-rs", "--reference_start", default=1, type=int,
                        help="The start coordinate of the reference sequence for read mapping", required=False)
    parser.add_argument("-re", "--reference_end", default=False, type=int,
                        help="The end coordinate of the reference sequence for read mapping. Default = full length",
                        required=False)
    parser.add_argument("-mi", "--min_len", default=300, type=int, help="The minimum read length allowed",
                        required=False)
    parser.add_argument("-ma", "--max_len", default=700, type=int, help="The minimum read length allowed",
                        required=False)

    args = parser.parse_args()

    project_path = args.project_path
    sample_names = args.sample_names
    reference = args.reference
    make_index = args.make_index
    reference_start = args.reference_start
    reference_end = args.reference_end
    min_len = args.min_len
    max_len = args.max_len

    main(project_path, sample_names, reference, make_index, reference_start, reference_end, min_len, max_len)
