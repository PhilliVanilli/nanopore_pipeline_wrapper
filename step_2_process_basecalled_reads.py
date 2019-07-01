import os
import sys
import subprocess
import argparse
import pathlib
import shutil
import collections
import datetime
import math
from itertools import groupby
import pandas as pd
from Bio import SeqIO


__author__ = 'Colin Anthony'


class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass


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


def fasta_to_dct(file_name):
    """
    :param file_name: The fasta formatted file to read from.
    :return: a dictionary of the contents of the file name given. Dictionary in the format:
    {sequence_id: sequence_string, id_2: sequence_2, etc.}
    """
    dct = collections.defaultdict(str)
    my_gen = py3_fasta_iter(file_name)
    for i, (k, v) in enumerate(my_gen):
        # resolve for duplicate sequence headers
        new_key = k.replace(" ", "_") + str(i).zfill(4)
        dct[new_key] = v.upper()

    return dct


def gather_fastqs(fastq_path, run_name, max_len, min_len):

    fastq_outpath = fastq_path.parent
    fastq_name = f"{run_name}_all.fastq"
    output_fastq = pathlib.Path(fastq_outpath, fastq_name)

    all_fastqs = list(fastq_path.glob("*.fastq"))

    with open(output_fastq, 'w') as handle:
        for fastq in all_fastqs:
            for record in SeqIO.parse(open(fastq), "fastq"):
                seq_len = len(record.seq)
                if seq_len > max_len or seq_len < min_len:
                    continue
                else:
                    SeqIO.write([record], handle, "fastq")

    if output_fastq.is_file():
        return True
    else:
        return False


def d_freq_lists(dna_list, n):
    """

    :param dna_list: (list) a alist of DNA sequences
    :param n: (int) length of alignment (ie number of positions)
    :return: (dict) a dictionary of the frequency for each base, for each site in the alignment
    """
    counts_dict = {'A': [0]*n, 'C': [0]*n, 'G': [0]*n, 'T': [0]*n, '-': [0]*n, "N": [0]*n}
    bases = ["A", "C", "G", "T", "-"]
    depth_dict = {'non_gap': [0] * n, 'gap': [0] * n}
    total = len(dna_list)
    for seq in dna_list:
        for index, dna in enumerate(seq):
            if dna.upper() not in bases:
                counts_dict["N"][index] += 1
            else:
                counts_dict[dna][index] += 1
            if dna == "-":
                depth_dict["gap"][index] += 1
            else:
                depth_dict["non_gap"][index] += 1

    for base, countslist in counts_dict.items():
        for i, cnt in enumerate(countslist):
            frq = round((cnt/total*100), 4)
            countslist[i] = frq
        counts_dict[base] = countslist

    return counts_dict, depth_dict


def consensus_maker(d, min_depth):
    """
    Create a consensus sequence from an alignment
    :param d: (dict) dictionary of an alignment (key = seq name (str): value = aligned sequence (str))
    :param min_depth: (int) the minimum depth required to call a base in the consensus (otherwise called as "!"
    :return: (str) the consensus sequence
    """
    seq_list = []
    for names, seq in d.items():
        seq_list.append(seq)
    if not seq_list:
        raise IndexError

    seq_length = len(seq_list[0])
    master_profile, depth_profile = d_freq_lists(seq_list, seq_length)

    consensus = ""
    degen = {('A', 'G'): 'R', ('C', 'T'): 'Y', ('A', 'C'): 'M', ('G', 'T'): 'K', ('C', 'G'): 'S', ('A', 'T'): 'W',
             ('A', 'C', 'T'): 'H', ('C', 'G', 'T'): 'B', ('A', 'C', 'G'): 'V', ('A', 'G', 'T'): 'D',
             ('A', 'C', 'G', 'T'): 'N'}

    for position in range(seq_length):
        if depth_profile["non_gap"][position] <= min_depth:
            consensus += str("N")
        else:
            dct = {base: master_profile[base][position] for base in ['A', 'C', 'G', 'T', 'N']}
            # get the highest frequency value
            max_freq = max(dct.values())
            # get the base with the highest frequency value
            base_with_max_freq = max(dct, key=dct.get)
            # if multiple bases share the max frequency make a list of them for degeneracy code lookup
            most_freq_bases = list(sorted(base for base in ['A', 'C', 'G', 'T'] if dct[base] == max_freq))
            if len(most_freq_bases) == 1:
                consensus += str(base_with_max_freq)
            else:
                most_freq_bases = tuple(most_freq_bases)
                consensus += str(degen[most_freq_bases])

    return consensus


def rename_fasta(fasta_file_name_path, sample_name, cons_type):
    fasta_d = fasta_to_dct(fasta_file_name_path)
    os.unlink(fasta_file_name_path)
    with open(fasta_file_name_path, 'w') as fh:
        for seq_name, seq in fasta_d.items():
            new_name = f"{sample_name}_{cons_type}"
            fh.write(f">{new_name}\n{seq}\n")


def cat_sample_names(barcode, run_name):
    if barcode != '':
        file_name = f"{run_name}_{barcode}.fastq"
    else:
        file_name = " "

    return file_name


def main(project_path, sample_names, reference, ref_start, ref_end, min_len, max_len, min_depth, run_step,
         rerun_step_only, msa_cons_only):

    # set the primer_scheme directory
    script_folder = pathlib.Path(__file__).absolute().parent
    primer_scheme_dir = pathlib.Path(script_folder, "primer-schemes")

    # get folder paths
    project_path = pathlib.Path(project_path).absolute()
    run_name = project_path.parts[-1]
    fast5_dir = pathlib.Path(project_path, "fast5")
    fastq_dir = pathlib.Path(project_path, "fastq")
    sequencing_summary_file = pathlib.Path(fastq_dir, "sequencing_summary.txt")
    sample_names = pathlib.Path(sample_names).absolute()
    demultipled_folder = pathlib.Path(project_path, "demultiplexed")
    sample_folder = pathlib.Path(project_path, "samples")
    master_reads_file = pathlib.Path(project_path, run_name + "_all.fastq")
    time_stamp = str('{:%Y-%m-%d_%H_%M}'.format(datetime.datetime.now()))
    log_file = pathlib.Path(project_path, f"{time_stamp}_{run_name}_log_file.txt")
    with open(log_file, "w") as handle:
        handle.write(f"# start of pipeline run for project: {run_name}\n")

    # set dir to project dir so that output is written in correct place by external tools
    os.chdir(project_path)

    # set the reference genome
    reference_scheme = \
        {"ChikECSA_V1_800": pathlib.Path(primer_scheme_dir, "ChikECSA800", "V1", "ChikECSA800.reference.fasta"),
         "ChikAsian_V1_400": pathlib.Path(primer_scheme_dir, "ChikAsian400", "V1", "ChikAsian400.reference.fasta"),
         "ZikaAsian_V1_400": pathlib.Path(primer_scheme_dir, "ZikaAsian400", "V1", "ZikaAsian400.reference.fasta"),
         }

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
    with open(log_file, "a") as handle:
        handle.write(f"\nReference is {chosen_ref_scheme}\nPrimer bed file is {chosen_ref_scheme_bed_file}\n")

    if run_step == 1:
        # gather all fastq's into one file and all sequencing_summary.txt's into one file
        print(f"\nrunning: collecting all fastq files into one file")
        with open(log_file, "a") as handle:
            handle.write(f"\nrunning: collecting all fastq files into one file")
        gather_fastq_status = gather_fastqs(fastq_dir, run_name, max_len, min_len)
        if not gather_fastq_status:
            print("gathering fastq files failed, concatenated fastq file was not found")
            with open(log_file, "a") as handle:
                handle.write("gathering fastq files failed, concatenated fastq file was not found\nexiting")
            sys.exit("exiting")
        if not rerun_step_only:
            run_step = 2
        else:
            sys.exit("Run step only completed, exiting")

    if run_step == 2:
        # demultiplex with porchop
        print(f"\nrunning: porechop demultiplexing")
        with open(log_file, "a") as handle:
            handle.write(f"\nrunning: porechop demultiplexing")

        if not master_reads_file.is_file():
            print(f"could not find the concatenated fastq, {master_reads_file}")
            with open(log_file, "a") as handle:
                handle.write(f"could not find the concatenated fastq, {master_reads_file}\nexiting")
            sys.exit("exiting")

        # get number of sequences and split into subfiles if too large
        max_fastq_size_porechop = 1000000
        grep_cmd = f"grep -c '^@' {master_reads_file}"
        try:
            fastq_entries = int(subprocess.check_output(grep_cmd, shell=True).decode(sys.stdout.encoding).strip())
        except subprocess.CalledProcessError as e:
            with open(log_file, "a") as handle:
                handle.write(f"could not grep the fastq file, error was:\n{e}\nexiting")
            sys.exit("exiting")

        num_of_chunks = math.ceil(fastq_entries / max_fastq_size_porechop)
        if num_of_chunks > 1.0:
            with open(log_file, "a") as handle:
                handle.write(f"\nexperiment run fastq file too large for porechop\nchunking file into parts\n")
            # file is possibly too large for system memory to run porechop, chunk into parts
            split_into_parts_by_x_lines = int(fastq_entries / num_of_chunks) * 4

            split_cmd = f"split --additional-suffix _temp_chunk.fastq -l {split_into_parts_by_x_lines} " \
                f"{master_reads_file} {str(master_reads_file.stem).replace('all', 'all_')}"

            try_except_exit_on_fail(split_cmd)

            # for each chunk, do porechop
            colleted_temp_folders = []
            for file in pathlib.Path(project_path).glob(f"{master_reads_file.stem}*temp_chunk.fastq"):
                with open(log_file, "a") as handle:
                    handle.write(f"\nrunning porechop on chunk: {file}\n")
                tmp_demix_folder = pathlib.Path(demultipled_folder, file.stem)
                tmp_demix_folder.mkdir(mode=0o777, parents=True, exist_ok=True)
                colleted_temp_folders.append(tmp_demix_folder)
                demultiplex_cmd = f"porechop --format fastq --verbosity 2 -i {file} " \
                                  f"--discard_middle " \
                                  f"--require_two_barcodes " \
                                  f"--barcode_threshold 80 " \
                                  f"--threads 8 " \
                                  f"--check_reads 10000 " \
                                  f"--barcode_diff 5 " \
                                  f"--barcode_dir {tmp_demix_folder} " \
                                  f"> {file}.demultiplexreport.txt"

                with open(log_file, "a") as handle:
                    handle.write(f"\n{demultiplex_cmd}\n")

                try_except_exit_on_fail(demultiplex_cmd)

            # collect file path for each barcode from each demultiplexed chunk
            collect_demultiplexed_files = collections.defaultdict(list)
            for folder in pathlib.Path(demultipled_folder).glob(f"*temp_chunk"):
                for file in folder.glob("*.fastq"):
                    barcode_id = file.stem
                    collect_demultiplexed_files[barcode_id].append(file)

            # concatenate barcode file from each chunk into one file
            with open(log_file, "a") as handle:
                handle.write(f"\nconcatenating porechop demultiplexed files from each chunk\n")
            for barcode, file_list in collect_demultiplexed_files.items():
                files_to_cat = " ".join([str(x) for x in file_list])
                barcode_outfile = pathlib.Path(demultipled_folder, f"{barcode}.fastq")
                cat_cmd = f"cat {files_to_cat} > {barcode_outfile}"
                try_except_exit_on_fail(cat_cmd)

            # remove each chunked file and the temp demultiplex folder and files from each chunk
            with open(log_file, "a") as handle:
                handle.write(f"\nremoving chunked files and temp demultiplexed files from each chunk\n")
            # for folder in colleted_temp_folders:
                # shutil.rmtree(str(folder))
            for file in pathlib.Path(project_path).glob(f"{master_reads_file.stem}*temp_chunk.fastq"):
                os.unlink(str(file))

        else:
            # file size within limit, no need to chunk
            demultiplex_cmd = f"porechop --format fastq --verbosity 2 -i {master_reads_file} " \
                              f"--discard_middle " \
                              f"--require_two_barcodes " \
                              f"--barcode_threshold 80 " \
                              f"--threads 8 " \
                              f"--check_reads 10000 " \
                              f"--barcode_diff 5 " \
                              f"--barcode_dir {demultipled_folder} " \
                              f"> {master_reads_file}.demultiplexreport.txt"

            with open(log_file, "a") as handle:
                handle.write(f"\n{demultiplex_cmd}\n")

            try_except_exit_on_fail(demultiplex_cmd)

        # add run name to each demultiplexed file
        for file in list(demultipled_folder.glob("*.fastq")):
            path = file.parents[0]
            name = file.name
            if len(name) == 10:
                new_name = f"{run_name}_{name}"
                new_file = pathlib.Path(path, new_name)
                os.rename(str(file), new_file)

        if not rerun_step_only:
            run_step = 3
        else:
            sys.exit("Run step only completed, exiting")

    if run_step == 3:
        # index concatenated fastq with nanopolish
        print(f"\nrunning: nanopolish index on fast5/fastq files")
        with open(log_file, "a") as handle:
            handle.write(f"\nrunning: nanopolish index on fast5/fastq files\n")
            if not sequencing_summary_file.is_file():
                handle.write(f"\nSequencing summary file not found")
                nanopolish_index_cmd = f"nanopolish index -d {fast5_dir} {master_reads_file} "
            else:
                nanopolish_index_cmd = f"nanopolish index -s {sequencing_summary_file} -d {fast5_dir} " \
                    f"{master_reads_file} "
        try_except_exit_on_fail(nanopolish_index_cmd)
        if not rerun_step_only:
            run_step = 4
        else:
            sys.exit("Run step only completed, exiting")

    if run_step == 4:
        # concatenated demultiplexed files for each sample and setup sample names and barcode combinations
        sample_names_df = pd.read_csv(sample_names, sep=None, keep_default_na=False, na_values=['NA'], engine="python")
        sample_names_df['barcode_1'] = sample_names_df['barcode_1'].apply(lambda x: cat_sample_names(x, run_name))
        sample_names_df['barcode_2'] = sample_names_df['barcode_2'].apply(lambda x: cat_sample_names(x, run_name))
        sample_names_dict = sample_names_df.set_index('sample_name').T.to_dict(orient='list')

        for sample_name, [barcode_1, barcode_2] in sample_names_dict.items():
            sample_dir = pathlib.Path(sample_folder, sample_name)
            if not sample_dir.exists():
                pathlib.Path(sample_dir).mkdir(mode=0o777, parents=True, exist_ok=True)

            # allow for case where only one barcode was specified per sample.
            barcode_1_file = pathlib.Path(demultipled_folder, barcode_1)
            if barcode_2 == " ":
                barcode_2_file = ""
            else:
                barcode_2_file = pathlib.Path(demultipled_folder, barcode_2)

            cat_outfile = pathlib.Path(sample_dir, f"{sample_name}.fastq")
            cat_cmd = f"cat {str(barcode_1_file)} {str(barcode_2_file)} > {cat_outfile}"
            print(cat_cmd)
            run = try_except_continue_on_fail(cat_cmd)
            if not run:
                print("missing one or more demultiplexed files for this sample")
                with open(log_file, "a") as handle:
                    handle.write("\nmissing one or more demultiplexed files for this sample\n")
                continue

        if not rerun_step_only:
            run_step = 5
        else:
            sys.exit("Run step only completed, exiting")

    if run_step == 5:
        make_index_cmd = f"bwa index {chosen_ref_scheme}"
        try_except_exit_on_fail(make_index_cmd)

        all_sample_files = pathlib.Path(sample_folder).glob("*/*.fastq")
        for sample_fastq in all_sample_files:
            if not sample_fastq.is_file():
                print(f"could not find the concatenated sample fastq file: {sample_fastq}\nskipping sample")
                with open(log_file, "a") as handle:
                    handle.write(f"could not find the concatenated sample fastq file: {sample_fastq}\nskipping sample")
                continue
            sample_name = pathlib.Path(sample_fastq).stem
            sample_folder = pathlib.Path(sample_fastq).parent
            sam_name = pathlib.Path(sample_folder, sample_name + "_mapped.sam")
            bam_file = pathlib.Path(sample_folder, sample_name + "_mapped.bam")
            bam_file_sorted = pathlib.Path(sample_folder, sample_name + ".sorted.bam")
            trimmed_sam_file = pathlib.Path(sample_folder, sample_name + ".sorted.primerclipped.sam")
            trimmed_bam_file = pathlib.Path(sample_folder, sample_name + ".sorted.primerclipped.bam")
            sorted_trimmed_bam_file = pathlib.Path(sample_folder, sample_name + ".sorted.primerclipped_sorted.bam")
            rename_trimmed_bam_file = pathlib.Path(sample_folder, sample_name + ".primertrimmed.sorted.bam")
            vcf_file = pathlib.Path(sample_folder, sample_name + "_polished.vcf")
            nanopolish_cons_file = pathlib.Path(sample_folder, sample_name + "_consensus_nanopolish.fasta")
            bcftools_vcf_file = pathlib.Path(sample_folder, sample_name + "_bcftools.vcf")
            bcftools_cons_file = pathlib.Path(sample_folder, sample_name + "_consensus_bcftools.fasta")
            msa_fasta = pathlib.Path(sample_folder, sample_name + "_msa_from_bam_file.fasta")
            msa_cons = pathlib.Path(sample_folder, sample_name + "_msa_consensus.fasta")
            artic_cons_file = pathlib.Path(sample_folder, f"{sample_name}_consensus_artic.fasta")
            all_consensus_sequences = pathlib.Path(sample_folder, sample_name + "_all_consensus.fasta")

            os.chdir(sample_folder)

            with open(log_file, "a") as handle:
                handle.write(f"\n\n________________\nStarting processing sample: {sample_name}\n\n________________\n")

            # run read mapping using bwa
            print(f"\nrunning: bwa read mapping")
            bwa_cmd = f"bwa mem -t 4 -x ont2d {chosen_ref_scheme} {sample_fastq} -o {sam_name} 2>&1 | tee -a {log_file}"
            with open(log_file, "a") as handle:
                handle.write(f"\nrunning: bwa read mapping\n")
                handle.write(f"{bwa_cmd}\n")
            run = try_except_continue_on_fail(bwa_cmd)
            if not run:
                continue

            # convert sam to bam
            print(f"\nrunning: sam to bam conversion")
            sam_bam_cmd = f"samtools view -bSh {sam_name} -o {bam_file} 2>&1 | tee -a {log_file}"
            with open(log_file, "a") as handle:
                handle.write(f"\nrunning: sam to bam conversion\n")
                handle.write(f"{sam_bam_cmd}\n")
            run = try_except_continue_on_fail(sam_bam_cmd)
            if not run:
                continue

            # sort bam file
            print(f"\nrunning: sorting bam file")
            sort_sam_cmd = f"samtools sort -T {sample_name} {bam_file} -o {bam_file_sorted} 2>&1 | tee -a {log_file}"
            with open(log_file, "a") as handle:
                handle.write(f"\nrunning: sorting bam file\n")
                handle.write(f"{sort_sam_cmd}\n")
            run = try_except_continue_on_fail(sort_sam_cmd)
            if not run:
                continue

            # index bam file
            print(f"\nrunning: indexing bam file")
            index_bam_cmd = f"samtools index {bam_file_sorted} 2>&1 | tee -a {log_file}"
            with open(log_file, "a") as handle:
                handle.write(f"\nrunning: indexing bam file\n")
                handle.write(f"{index_bam_cmd}\n")
            run = try_except_continue_on_fail(index_bam_cmd)
            if not run:
                continue

            # remove primer sequences with custom script
            print(f"\nrunning: trim primer sequences from bam file")
            trim_script = pathlib.Path(script_folder, "clip_primers_from_bed_file.py")
            trim_primer = f"python {trim_script} -in {bam_file_sorted} -o {trimmed_sam_file} " \
                f"-b {chosen_ref_scheme_bed_file} 2>&1 | tee -a {log_file}"
            with open(log_file, "a") as handle:
                handle.write(f"\nrunning: soft clipping primer sequences from bam file\n")
                handle.write(f"{trim_primer}\n")
            run = try_except_continue_on_fail(trim_primer)
            if not run:
                continue

            # convert sam to bam
            print(f"\nrunning: sam to bam conversion of trimmed file")
            sam_bam_cmd = f"samtools view -bS {trimmed_sam_file} -o {trimmed_bam_file} 2>&1 | tee -a {log_file}"
            with open(log_file, "a") as handle:
                handle.write(f"\nrunning: sam to bam conversion\n")
                handle.write(f"{sam_bam_cmd}\n")
            run = try_except_continue_on_fail(sam_bam_cmd)
            if not run:
                continue

            # sort bam file
            print(f"\nrunning: sorting bam file")
            sort_sam_cmd = f"samtools sort -T {sample_name} {trimmed_bam_file} -o {sorted_trimmed_bam_file} " \
                f"2>&1 | tee -a {log_file}"
            with open(log_file, "a") as handle:
                handle.write(f"\nrunning: sorting bam file\n")
                handle.write(f"{sort_sam_cmd}\n")
            run = try_except_continue_on_fail(sort_sam_cmd)
            if not run:
                continue

            # index trimmed bam file
            print(f"\nrunning: indexing bam file")
            index_bam_cmd = f"samtools index {sorted_trimmed_bam_file} 2>&1 | tee -a {log_file}"
            with open(log_file, "a") as handle:
                handle.write(f"\nrunning: indexing bam file\n")
                handle.write(f"{index_bam_cmd}\n")
            run = try_except_continue_on_fail(index_bam_cmd)
            if not run:
                continue
            if not msa_cons_only:
                # run nanopolish
                print(f"\nrunning: nanopolish variant calling")
                nanopolish_cmd_v11 = f"nanopolish variants  -r --snps -o {vcf_file} " \
                    f"-w '{ref_name}:{ref_start}-{ref_end}' -t 4 --ploidy=1 -v " \
                    f"-r {master_reads_file} -b {sorted_trimmed_bam_file} -g {chosen_ref_scheme} " \
                    f"--min-candidate-frequency=0.3" \
                    f"--min-candidate-depth=10 --max-haplotypes=1000000"
                with open(log_file, "a") as handle:
                    handle.write(f"\nrunning: nanopolish variant calling using:\n")
                    handle.write(f"{nanopolish_cmd_v11}\n")

                run = try_except_continue_on_fail(nanopolish_cmd_v11)
                if not run:
                    continue

                # make nanopolish consensus
                print(f"\nrunning: making consensuses sequence from nanopolish")
                consensus_cmd = f"nanopolish vcf2fasta --skip-checks -g {chosen_ref_scheme} {vcf_file} > " \
                    f"{nanopolish_cons_file}"
                with open(log_file, "a") as handle:
                    handle.write(f"\nrunning: making consensuses sequence from from nanopolish using:\n")
                    handle.write(f"{consensus_cmd}\n")
                run = try_except_continue_on_fail(consensus_cmd)
                if not run:
                    continue

                # rename the fasta header to the sample name
                rename_fasta(nanopolish_cons_file, sample_name, "nanopolish_cons")

                # make bcftools consensus
                print(f"\nrunning: making consensuses sequence from bcftools")
                min_base_qual = 30  # default=13
                p_val_of_variant = 0.2  # default=0.5
                bcf_vcf_cmd = f"bcftools mpileup --min-BQ {min_base_qual} -Ou -f {chosen_ref_scheme} " \
                    f"{sorted_trimmed_bam_file} | bcftools call -c -p {p_val_of_variant} --ploidy 1 -v -Oz " \
                    f"-o {bcftools_vcf_file} 2>&1 | tee -a {log_file}"
                bcf_index_cmd = f"bcftools index {bcftools_vcf_file} 2>&1 | tee -a {log_file}"
                bcf_cons_cmd = f"bcftools consensus -H A -f {chosen_ref_scheme} {bcftools_vcf_file} " \
                    f"-o {bcftools_cons_file} 2>&1 | tee -a {log_file}"
                with open(log_file, "a") as handle:
                    handle.write(f"\nrunning: making consensuses sequence from bcftools:\n")
                    handle.write(f"{bcf_vcf_cmd}\n\n{bcf_index_cmd}\n\n{bcf_cons_cmd}\n")
                run = try_except_continue_on_fail(bcf_vcf_cmd)
                if not run:
                    continue
                run = try_except_continue_on_fail(bcf_index_cmd)
                if not run:
                    continue
                run = try_except_continue_on_fail(bcf_cons_cmd)
                if not run:
                    continue

                # rename the fasta header to the sample name
                rename_fasta(bcftools_cons_file, sample_name, "bcftools_cons")

                # make artic-ebov consensus
                print(f"\nrunning: making consensuses sequence from artic_ebov method")
                shutil.copyfile(sorted_trimmed_bam_file, rename_trimmed_bam_file)
                cons_file_script = pathlib.Path(script_folder, "margin_cons.py")

                set_min_depth = 100  # default=200
                set_min_qual = 200  # default=200

                cons_cmd = f"python {cons_file_script} -r {chosen_ref_scheme} -v {vcf_file} -b {rename_trimmed_bam_file} " \
                    f"-n {sample_name} -d {set_min_depth} -q {set_min_qual}"
                with open(log_file, "a") as handle:
                    handle.write(f"\nrunning: making consensuses sequence from artic_ebov method:\n")
                    handle.write(f"{cons_cmd}\n")
                run = try_except_continue_on_fail(cons_cmd)
                if not run:
                    continue

                # rename the fasta header to the sample name
                rename_fasta(artic_cons_file, sample_name, "artic_ebov_nanopolish_vcf_cons")

            # convert bam file to a mutli fasta alignment
            print(f"\nrunning: making consensuses sequence from bam to MSA with jvarkit")

            sam4web = pathlib.Path(script_folder, "jvarkit", "dist", "sam4weblogo.jar")
            msa_from_bam = f"java -jar {sam4web} -r '{ref_name}:{ref_start}-{ref_end}' -o {msa_fasta} " \
                f"{sorted_trimmed_bam_file} 2>&1 | tee -a {log_file}"
            print(msa_from_bam)

            with open(log_file, "a") as handle:
                handle.write(f"\nrunning: making consensuses sequence from bam to MSA with jvarkit\n")
                handle.write(f"{msa_from_bam}\n")
            run = try_except_continue_on_fail(msa_from_bam)
            if not run:
                continue

            # convert multi fasta alignment to consensus sequence
            fasta_msa_d = fasta_to_dct(msa_fasta)

            # set minimum depth for calling a position in the consensus sequence
            try:
                cons = consensus_maker(fasta_msa_d, min_depth)
            except IndexError as e:
                with open(log_file, "a") as handle:
                    handle.write(f"\nNo MSA made from Bam file\nno reads may have been mapped\n{e}\n")
            else:
                with open(msa_cons, 'w') as handle:
                    handle.write(f">{sample_name}_bam_msa_consensus\n{cons}\n")

            if not msa_cons_only:
                # add all consensus seqs into one file
                concat_consensus_cmd = f"cat {chosen_ref_scheme} {nanopolish_cons_file} {bcftools_cons_file} {msa_cons} " \
                    f"{artic_cons_file} > {all_consensus_sequences}"
                run = try_except_continue_on_fail(concat_consensus_cmd)
                if not run:
                    pass

            # plot depth and quality for sample
            plot_file_script = pathlib.Path(script_folder, "plot_depths_qual.py")
            plot_cmd = f"python {plot_file_script} -r {chosen_ref_scheme} -v {vcf_file} -b {rename_trimmed_bam_file} " \
                f"-n {sample_name} 2>&1 | tee -a {log_file}"
            run = try_except_continue_on_fail(plot_cmd)
            if not run:
                continue

            print(f"Completed processing sample: {sample_name}")
            with open(log_file, "a") as handle:
                handle.write(f"\n\n________________\nCompleted processing sample: {sample_name}\n\n________________\n")

    print("sample processing completed")
    with open(log_file, "a") as handle:
        handle.write(f"\nsample processing completed\n\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process raw nanopore reads to fasta consensus sequences",
                                     formatter_class=Formatter)

    parser.add_argument("-in", "--project_path", default=argparse.SUPPRESS, type=str,
                        help="The path to the directory containing the 'fast5' and 'fastq' folders ", required=True)
    parser.add_argument("-s", "--sample_names", default=argparse.SUPPRESS, type=str,
                        help="The csv file with the sample names and corresponding barcode combinations\n"
                             "This file should have three columns: barcode_1,barcode_2,sample_name", required=True)
    parser.add_argument("-r", "--reference", type=str, default="ChikAsianECSA_V1",
                        help="The reference genome and primer scheme to use",
                        choices=["ChikAsian_V1_400", "ChikECSA_V1_800", "ZikaAsian_V1_400"], required=False)
    parser.add_argument("-rs", "--reference_start", default=1, type=int,
                        help="The start coordinate of the reference sequence for read mapping", required=False)
    parser.add_argument("-re", "--reference_end", default=False, type=int,
                        help="The end coordinate of the reference sequence for read mapping. Default = full length",
                        required=False)
    parser.add_argument("-mi", "--min_len", type=int, help="The minimum read length allowed:\n = 300 for 400bp amplicon"
                                                           " design\n = 700 for 800bp amplicon design", required=True)
    parser.add_argument("-ma", "--max_len", type=int, help="The maximum read length allowed:\n = 500 for 400bp amplicon"
                                                           " design\n = 900 for 800bp amplicon design", required=True)
    parser.add_argument("-d", "--min_depth", type=int, default=100, help="The minimum coverage to call a position in "
                                                                         "the MSA to consensus", required=False)
    parser.add_argument("--run_step", default=1, type=int, required=False,
                        help="Only rerun the specified step:\n"
                             "--run_step 1 = gather fastqs\n"
                             "--run_step 2 = demultiplex\n"
                             "--run_step 3 = index master fastq file\n"
                             "--run_step 4 = concatenate demultiplexed files into sample files\n"
                             "--run_step 5 = run read mapping and all the variant calling steps on each sample\n")

    parser.add_argument("--rerun_step_only", default=False, action="store_true",
                        help="Only rerun the specified step", required=False)
    parser.add_argument("-m", "--msa_cons_only", default=False, action="store_true",
                        help="Only do MSA to consensus sequence for variant calling", required=False)

    args = parser.parse_args()

    project_path = args.project_path
    sample_names = args.sample_names
    reference = args.reference
    reference_start = args.reference_start
    reference_end = args.reference_end
    min_len = args.min_len
    max_len = args.max_len
    min_depth = args.min_depth
    run_step = args.run_step
    rerun_step_only = args.rerun_step_only
    msa_cons_only = args.msa_cons_only

    main(project_path, sample_names, reference, reference_start, reference_end, min_len, max_len, min_depth, run_step,
         rerun_step_only, msa_cons_only)
