import os
import sys
import subprocess
import argparse
import pathlib
from itertools import groupby
import shutil
import pandas as pd


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


def main(project_path, run_name, sample_names, reference, primer_scheme_dir, make_index, ref_start, ref_end):

    project_path = pathlib.Path(project_path).absolute()
    fast5_dir = pathlib.Path(project_path, "fast5")
    fastq_dir = pathlib.Path(project_path, "fastq")
    sample_names = pathlib.Path(sample_names).absolute()
    demultipled_folder = pathlib.Path(project_path, "demultiplexed")
    sample_folder = pathlib.Path(project_path, "samples")
    # set dir to project dir
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
    print(f"\nrunning: artic gather to collect all fastq files into one file")
    gather_cmd = f"artic gather --guppy --min-length 100 --max-length 700 --prefix {run_name} {fastq_dir}"
    #try_except_exit_on_fail(gather_cmd)

    # demultiplex with porchop
    print(f"\nrunning: porechop demultiplexing\n")
    master_reads_file = pathlib.Path(project_path, run_name + "_all.fastq")
    if not master_reads_file.is_file():
        print(f"could not find the concatenated fastq, {master_reads_file}")
        sys.exit("exiting")
    demultiplex_cmd = f"porechop --verbosity 2 --untrimmed -i {master_reads_file} " \
                      f"--discard_middle " \
                      f"--require_two_barcodes " \
                      f"--barcode_threshold 80 " \
                      f"--threads 4 " \
                      f"--check_reads 10000 " \
                      f"--barcode_diff 5 " \
                      f"--barcode_dir {demultipled_folder} " \
                      f"> {str(master_reads_file)}.demultiplexreport.txt"

    #try_except_exit_on_fail(demultiplex_cmd)

    for file in list(demultipled_folder.glob("*.fastq")):
        path = file.parents[0]
        name = file.name
        if len(name) == 10:
            new_name = f"{run_name}_{name}"
            new_file = pathlib.Path(path, new_name)
            os.rename(str(file), new_file)

    # index concatenated fastq with nanopolish
    print(f"\nrunning: nanopolish index on fast5/fastq files\n")
    nanopolish_index_cmd = f"nanopolish index -d {fast5_dir} {master_reads_file}"
    #try_except_exit_on_fail(nanopolish_index_cmd)

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
        sample_folder = pathlib.Path(sample_fastq).parents[0]
        sam_name = pathlib.Path(sample_folder, sample_name + "_mapped.sam")
        bam_file = pathlib.Path(sample_folder, sample_name + "_mapped.bam")
        bam_file_sorted = pathlib.Path(sample_folder, sample_name + ".sorted.bam")
        trimmed_bam_file = pathlib.Path(sample_folder, sample_name + ".sorted.primerclipped.bam")
        rename_trimmed_bam_file = pathlib.Path(sample_folder, sample_name + ".primertrimmed.sorted.bam")
        vcf_file = pathlib.Path(sample_folder, sample_name + ".vcf")
        consensus_file = pathlib.Path(sample_folder, sample_name + "_consensus.fasta")

        os.chdir(sample_folder)

        # run read mapping using bwa
        print(f"\nrunning: bwa read mapping\n")
        bwa_cmd = f"bwa mem -t 4 -x ont2d {chosen_ref_scheme} {sample_fastq} > {sam_name}"
        run = try_except_continue_on_fail(bwa_cmd)
        if not run:
            continue

        # convert sam to bam
        print(f"\nrunning: sam to bam conversion\n")
        sam_bam_cmd = f"samtools view -bS {sam_name} > {bam_file}"
        run = try_except_continue_on_fail(sam_bam_cmd)
        if not run:
            continue

        # sort bam file
        print(f"\nrunning: sorting bam file\n")
        sort_sam_cmd = f"samtools sort -T {sample_name} {bam_file} -o {bam_file_sorted}"
        run = try_except_continue_on_fail(sort_sam_cmd)
        if not run:
            continue

        # index bam file for bamclipper
        print(f"\nrunning: indexing bam file\n")
        index_bam_cmd = f"samtools index {bam_file_sorted}"
        run = try_except_continue_on_fail(index_bam_cmd)
        if not run:
            continue

        # remove primer sequences with bamclipper
        print(f"\nrunning: trim primer sequences from bam file\n")

        trim_primer = f"bamclipper.sh -b {bam_file_sorted} " \
            f"-p {chosen_ref_scheme_bed_file} -n 4 -u 1 -d 1"
        run = try_except_continue_on_fail(trim_primer)
        if not run:
            continue

        # run nanopolish
        print(f"\nrunning: nanopolish variant calling\n")
        nanopolish_cmd = f"nanopolish variants --min-flanking-sequence 10 --fix-homopolymers --max-rounds=100 " \
            f"--progress -t 4 --min-flanking-sequence=30 --ploidy 1 " \
            f"--reads {master_reads_file} -o {vcf_file} -b {trimmed_bam_file} -g {chosen_ref_scheme} " \
            f"-w '{ref_name}:{ref_start}-{ref_end}' --consensus {consensus_file}"
        run = try_except_continue_on_fail(nanopolish_cmd)
        if not run:
            continue

        # copy renamed trimmed bam file
        shutil.copyfile(trimmed_bam_file, rename_trimmed_bam_file)

        # extract variant calls
        print(f"\nrunning: extracting variant calls to tab file\n")
        extract_vcf = f"vcfextract {sample_name} > {sample_name}.variants.tab"
        run = try_except_continue_on_fail(extract_vcf)
        if not run:
            continue

        # make consensus
        print(f"\nrunning: making consensus sequence from bam and variant call (vcf) files\n")
        consensus_cmd = f"margin_cons {chosen_ref_scheme} {vcf_file} {rename_trimmed_bam_file} > " \
            f"{sample_name}_their_consensus.fasta"
        run = try_except_continue_on_fail(consensus_cmd)
        if not run:
            continue

        sed_cmd = f"sed -i '1d' {sample_name}_their_consensus.fasta"
        run = try_except_continue_on_fail(sed_cmd)
        if not run:
            continue
        print(f"Completed processing sample: {sample_name}")
        input("enter")
    print("sample processing completed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process indexed, demultiplexed nanopore reads to consensus sequences",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-in", "--project_path", default=argparse.SUPPRESS, type=str,
                        help="The path to the directory containing the 'fast5' and 'fastq' folders ", required=True)
    parser.add_argument("-n", "--run_name", default=argparse.SUPPRESS, type=str,
                        help="The path to where the sequence_summary file will be written", required=True)
    parser.add_argument("-s", "--sample_names", default=argparse.SUPPRESS, type=str,
                        help="The csv file with the sample names and corresponding barcode combinations"
                             "This file should have three columns: barcode_1,barcode_2,sample_name", required=True)
    parser.add_argument("-r", "--reference", type=str, default="ChikAsianECSA_V1",
                        help="The reference genome and primer scheme to use",
                        choices=["ChikAsianECSA_V1", "ZikaAsian_V1", "ZaireEbola_V1", "ZaireEbola_V2", "LassaL_V1",
                                 "LassaS_V1"], required=False)
    parser.add_argument("-p", "--primer_scheme_dir", type=str, default=argparse.SUPPRESS,
                        help="The path to the 'primer_scheme' folder", required=True)
    parser.add_argument("-m", "--make_index", default=False, action="store_true",
                        help="The path to the 'primer_scheme' folder", required=False)
    parser.add_argument("-rs", "--reference_start", default=1, type=int,
                        help="The start coordinate of the reference sequence for read mapping", required=False)
    parser.add_argument("-re", "--reference_end", default=False, type=int,
                        help="The end coordinate of the reference sequence for read mapping. Default = full length",
                        required=False)

    args = parser.parse_args()

    project_path = args.project_path
    run_name = args.run_name
    sample_names = args.sample_names
    reference = args.reference
    primer_scheme_dir = args.primer_scheme_dir
    make_index = args.make_index
    reference_start = args.reference_start
    reference_end = args.reference_end

    main(project_path, run_name, sample_names, reference, primer_scheme_dir, make_index, reference_start, reference_end)
