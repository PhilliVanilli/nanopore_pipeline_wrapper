import os, time
import sys
import argparse
import re
import pathlib
import datetime
import pandas as pd
import shutil
import csv
from src.misc_functions import try_except_continue_on_fail
from src.misc_functions import try_except_exit_on_fail
from src.misc_functions import py3_fasta_iter
from src.misc_functions import cat_sample_names
from src.misc_functions import filter_length
from src.misc_functions import fasta_to_dct
from basecall_guppy import main as gupppy_basecall
from demultiplex_guppy import main as guppy_demultiplex


__author__ = 'Colin Anthony & Philippe Selhorst'


class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass


def main(project_path, reference, ref_start, ref_end, min_len, max_len, min_depth, run_step,
         rerun_step_only, basecall_mode, msa_cons, artic, cpu_threads, gpu_threads, gpu_buffers, use_gaps, use_bwa,
         guppy_path, real_time):

    # set the primer_scheme directory
    script_folder = pathlib.Path(__file__).absolute().parent
    primer_scheme_dir = pathlib.Path(script_folder, "primer-schemes")

    # get folder paths
    project_path = pathlib.Path(project_path).absolute()
    plot_folder = pathlib.Path(project_path, "seq_depth_plots")
    if os.path.exists(plot_folder):
        shutil.rmtree(plot_folder)
    plot_folder.mkdir(mode=0o777, parents=True, exist_ok=True)
    run_name = project_path.parts[-1]
    fast5_dir = pathlib.Path(project_path, "fast5")
    fastq_dir = pathlib.Path(project_path, "fastq")

    seq_summary_file_name = ""
    for file in project_path.glob('sequencing_summary*.txt'):
        seq_summary_file_name = file
    seq_summary_file = pathlib.Path(seq_summary_file_name).resolve()

    sample_names = pathlib.Path(project_path, "sample_names.csv")
    if not sample_names:
        sys.exit("Could not find sample_names.csv in project folder")
    demultiplexed_folder = pathlib.Path(project_path, "demultiplexed")
    sample_folder = pathlib.Path(project_path, "samples")

    print(f"\nProject folder is {project_path}")

    # master_reads_file = pathlib.Path(project_path, run_name + "_all.fastq")
    time_stamp = str('{:%Y-%m-%d_%H_%M}'.format(datetime.datetime.now()))
    log_file = pathlib.Path(project_path, f"{time_stamp}_{run_name}_log_file.txt")
    log_file_final = pathlib.Path(project_path, f"{time_stamp}_{run_name}_log_file_final.txt")

    with open(log_file, "w") as handle:
        handle.write(f"# start of pipeline run for project: {run_name}\n")

    now = datetime.datetime.now()
    date_time = now.strftime("%d/%m/%Y, %H:%M:%S")
    print(f"\nstart time = {date_time}\n")
    with open(log_file, "a") as handle:
        handle.write(f"\nstart time = {date_time}\n")

    # set dir to project dir so that output is written in correct place by external tools
    os.chdir(project_path)

    # set the reference genome
    reference_scheme = \
        {"ChikECSA_V1_800": pathlib.Path(primer_scheme_dir, "ChikECSA800", "V1", "ChikECSA800.reference.fasta"),
         "ChikAsian_V1_400": pathlib.Path(primer_scheme_dir, "ChikAsian400", "V1", "ChikAsian400.reference.fasta"),
         "ZikaAsian_V1_400": pathlib.Path(primer_scheme_dir, "ZikaAsian400", "V1", "ZikaAsian400.reference.fasta"),
         "SARS2_V1_800": pathlib.Path(primer_scheme_dir, "SARS2_800", "V1", "SARS2_800.reference.fasta"),
         "SARS2_V1_400": pathlib.Path(primer_scheme_dir, "SARS2_400", "V1", "SARS2_400.reference.fasta"),
         "RSVA_V1_3000": pathlib.Path(primer_scheme_dir, "RSVA_3000", "V1", "RSVA_3000.reference.fasta"),
         "RSVB_V1_3000": pathlib.Path(primer_scheme_dir, "RSVB_3000", "V1", "RSVB_3000.reference.fasta"),
         "DENV1_V1_400": pathlib.Path(primer_scheme_dir, "DENV1_400", "V1", "DENV1_400.reference.fasta"),
         "DENV2_V1_400": pathlib.Path(primer_scheme_dir, "DENV2_400", "V1", "DENV2_400.reference.fasta")
         }

    chosen_ref_scheme = str(reference_scheme[reference])
    chosen_ref_scheme_bed_file = chosen_ref_scheme.replace(".reference.fasta", ".scheme.bed")
    scheme_name = reference.replace("_V1_", "_")
    ref_name, ref_seq = list(py3_fasta_iter(chosen_ref_scheme))[0]
    ref_name = ref_name.split()[0]
    if not ref_start or ref_start == 0:
        ref_start = 1
    if not ref_end or ref_end > len(ref_seq):
        ref_end = len(ref_seq)
    # reference_slice = f'{ref_name}:{ref_start}-{ref_end}'
    print(f"\nReference is {chosen_ref_scheme}")
    print(f"\nPrimer bed file is {chosen_ref_scheme_bed_file}\n")
    with open(log_file, "a") as handle:
        handle.write(f"\nReference is {chosen_ref_scheme}\nPrimer bed file is {chosen_ref_scheme_bed_file}\n")

    if run_step == 0:
        run = gupppy_basecall(fast5_dir, guppy_path, fastq_dir, gpu_threads, basecall_mode, real_time, reference, script_folder)
        faildir = pathlib.Path(fastq_dir, "fail")
        shutil.rmtree(faildir)
        if run and not rerun_step_only:
            run_step = 1
        elif run and rerun_step_only:
            sys.exit("Run step only completed, exiting")
        else:
            sys.exit("Basecalling failed")

    if run_step == 1:
        # demultiplex
        print(f"\nrunning: demultiplexing")
        with open(log_file, "a") as handle:
            handle.write(f"\nrunning: demultiplexing")
        if not list(fastq_dir.glob("*.fastq*")):
            fastq_dir = pathlib.Path(fastq_dir, "pass")
            if not list(fastq_dir.glob("*.fastq*")):
                print(f"No fastq files found in {str(fastq_dir)} or {str(fastq_dir.parent)}")
                sys.exit("fastq files not found")
        run = guppy_demultiplex(fastq_dir, guppy_path, demultiplexed_folder, cpu_threads, gpu_buffers, gpu_threads)
        if run and not rerun_step_only:
            run_step = 2
        elif run and rerun_step_only:
            sys.exit("demultiplexing completed, exiting")
        else:
            sys.exit("demultiplexing failed")

    if run_step == 2:

        pre_existing_files = list(demultiplexed_folder.glob("*.fastq"))
        if pre_existing_files:
            print("Found existing files in top level of demultiplex folder.\nThese files will be deleted")
            for file in pre_existing_files:
                os.unlink((str(file)))

        for folder in demultiplexed_folder.glob("barcode*"):
            search = list(pathlib.Path(folder).glob("*.fastq"))
            if not search:
                print(f"no files in folder\nskipping folder: {folder}\n")
                continue
            if len(search) > 1:
                barcode_number = pathlib.Path(search[0]).parent.parts[-1]
                concat_outfile = f"cat_barcode_{barcode_number}.fastq"
                cat_cmd = f"cat "
                for file in search:
                    cat_cmd += f"{str(file)} "
                cat_cmd += f" > {concat_outfile}"
                try_except_exit_on_fail(cat_cmd)
                new_name = pathlib.Path(demultiplexed_folder, f"{run_name}_{barcode_number}.fastq")
                filtered_file = filter_length(concat_outfile, new_name, max_len, min_len)

                os.unlink(str(concat_outfile))
                if not filtered_file:
                    print(f"no sequences in file after length filtering for {concat_outfile}\n")

                # sed_syntax = r"\t/\n"
                # bash_cmd = f"cat {concat_outfile} | paste - - - - | awk 'length($2)  >= {min_len} && length($2) <= {max_len}' | sed 's/{sed_syntax}/g' > {str(new_name)}"
                # print(bash_cmd)
                # seqmagick_cmd = f"seqmagick quality-filter --min-mean-quality 0 " \
                #                 f"--min-length {min_len} --max-length {max_len} " \
                #                 f"{concat_outfile} {new_name} "
                # vsearch_cmd = f"vsearch --fastq_filter {concat_outfile} -fastq_maxlen {max_len} " \
                #               f"--fastq_qmax 100 --fastq_minlen {min_len} --fastqout {new_name}"
                # try_except_exit_on_fail(bash_cmd)

            else:
                file = pathlib.Path(search[0])
                barcode_number = file.parent.parts[-1]
                new_name = pathlib.Path(demultiplexed_folder, f"{run_name}_{barcode_number}.fastq")

                filtered_file = filter_length(file, new_name, max_len, min_len)

                if not filtered_file:
                    print(f"no sequences in file after length filtering for {file}\n")

                # sed_syntax = r"\t/\n"
                # bash_cmd = f"cat {file} | paste - - - - | awk 'length($2)  >= {min_len} && length($2) <= {max_len}' | sed 's/{sed_syntax}/g' > {str(new_name)}"
                # print(bash_cmd)
                # seqmagick_cmd = f"seqmagick quality-filter --min-mean-quality 0 " \
                #                 f"--min-length {min_len} --max-length {max_len} " \
                #                 f"{file} {new_name} "
                # vsearch_cmd = f"vsearch --fastq_filter {file} -fastq_maxlen {max_len} --fastq_minlen {min_len} " \
                #               f"--fastq_qmax 100 --fastqout {new_name}"
                # try_except_exit_on_fail(bash_cmd)

        if not rerun_step_only:
            run_step = 3
        elif rerun_step_only:
            sys.exit("filer demultiplexed files and rename them completed, exiting")
        else:
            sys.exit("filtering and renaming demultiplexed files failed")

    # if run_step == 3 and not msa_cons:
    #     # index concatenated fastq with nanopolish
    #     print(f"\nrunning: nanopolish index on fast5/fastq files")
    #     with open(log_file, "a") as handle:
    #         handle.write(f"\nrunning: nanopolish index on fast5/fastq files\n")
    #         if not sequencing_summary_file.is_file():
    #             handle.write(f"\nSequencing summary file not found")
    #             nanopolish_index_cmd = f"nanopolish index -d {fast5_dir} {master_reads_file} "
    #         else:
    #             nanopolish_index_cmd = f"nanopolish index -s {sequencing_summary_file} -d {fast5_dir} " \
    #                 f"{master_reads_file} "
    #     try_except_exit_on_fail(nanopolish_index_cmd)
    #     if not rerun_step_only:
    #         run_step = 4
    #     else:
    #         sys.exit("Run step only completed, exiting")

    if run_step == 3:
        # concatenated demultiplexed files for each sample and setup sample names and barcode combinations
        print("collecting demultiplexed files into sample.fastq files based on specified sample barcode combinations\n")
        with open(log_file, "a") as handle:
            handle.write(f"\ncollecting demultiplexed files into sample.fastq files based on specified sample "
                         f"barcode combinations\n")

        sample_names_df = pd.read_csv(sample_names, sep=None, keep_default_na=False, na_values=['NA'], engine="python")
        sample_names_df['barcode_1'] = sample_names_df['barcode_1'].apply(lambda x: cat_sample_names(x, run_name))
        sample_names_df['barcode_2'] = sample_names_df['barcode_2'].apply(lambda x: cat_sample_names(x, run_name))
        sample_names_dict = sample_names_df.set_index('sample_name').T.to_dict(orient='list')

        for sample_name, [barcode_1, barcode_2] in sample_names_dict.items():
            sample_dir = pathlib.Path(sample_folder, sample_name)
            if not sample_dir.exists():
                pathlib.Path(sample_dir).mkdir(mode=0o777, parents=True, exist_ok=True)

            # allow for case where only one barcode was specified per sample.
            barcode_1_file = pathlib.Path(demultiplexed_folder, barcode_1)
            if barcode_2 == " ":
                barcode_2_file = ""
            else:
                barcode_2_file = pathlib.Path(demultiplexed_folder, barcode_2)
            if artic:
                artic_folder = pathlib.Path(sample_dir, "art")
                if os.path.exists(artic_folder):
                    shutil.rmtree(artic_folder)
                artic_folder.mkdir(mode=0o777, parents=True, exist_ok=True)
                cat_outfile = pathlib.Path(sample_dir, f"art/{sample_name}.fastq")
                cat_cmd = f"cat {str(barcode_1_file)} {str(barcode_2_file)} > {cat_outfile}"
                print(cat_cmd)
                run = try_except_continue_on_fail(cat_cmd)
                if not run:
                    print("missing one or more demultiplexed files for this sample")
                    with open(log_file, "a") as handle:
                        handle.write("\nmissing one or more demultiplexed files for this sample\n")
                    continue
            if msa_cons:
                msa_folder = pathlib.Path(sample_dir, "msa")
                if os.path.exists(msa_folder):
                    shutil.rmtree(msa_folder)
                msa_folder.mkdir(mode=0o777, parents=True, exist_ok=True)
                if artic:
                    all_art_sample_files = pathlib.Path(sample_folder).glob("*/art/*.fastq")
                    for filepath in all_art_sample_files:
                        destinationpath = pathlib.Path(str(filepath).replace('art','msa'))
                        shutil.copy(filepath, destinationpath)
                else:
                    cat_outfile = pathlib.Path(sample_dir, f"msa/{sample_name}.fastq")
                    cat_cmd = f"cat {str(barcode_1_file)} {str(barcode_2_file)} > {cat_outfile}"
                    print(cat_cmd)
                    run = try_except_continue_on_fail(cat_cmd)
                    if not run:
                        print("missing one or more demultiplexed files for this sample")
                        with open(log_file, "a") as handle:
                            handle.write("\nmissing one or more demultiplexed files for this sample\n")
                        continue

        if not rerun_step_only:
            run_step = 4
        else:
            sys.exit("Run step only completed, exiting")

    if run_step == 4:

        print("Running variant calling on samples")
        with open(log_file, "a") as handle:
            handle.write(f"\nRunning variant calling on samples\n")

        # get number of samples
        number_samples = len(next(os.walk(sample_folder))[1])
        print(f"\nNumber of samples = {str(number_samples)}\n")

        # make bwa index if needed
        if use_bwa:
            make_index_cmd = f"bwa index {chosen_ref_scheme}"
            with open(log_file, "a") as handle:
                handle.write(f"\n{make_index_cmd}\n")
            try_except_exit_on_fail(make_index_cmd)

        # make variable for project file containing all samples' consensus sequences
        project_name = project_path.parts[-1]
        all_samples_consens_seqs = pathlib.Path(project_path, project_name + "_all_samples.fasta")

        # initialize the file, and add reference to all consensus file
        with open(all_samples_consens_seqs, 'w') as fh:
            fh.write(f">{ref_name}\n{ref_seq}\n")
        p = pathlib.Path(project_path, project_name + '_mapping.csv')
        with open(p, 'w') as fh:
            fh.close()

        max_threads = cpu_threads
        used_threads = 0
        msa_threads = 1
        artic_threads = 2
        # delete pre existing files
        for file in pathlib.Path(sample_folder).glob("*/*/*.*"):
            if not str(file).endswith(".fastq"):
                os.remove(file)

        log_file_art_temp = pathlib.Path(project_path, f"{time_stamp}_{run_name}_log_file_art_temp.txt")
        log_file_msa_temp = pathlib.Path(project_path, f"{time_stamp}_{run_name}_log_file_msa_temp.txt")
        log_file_msa = pathlib.Path(project_path, f"{time_stamp}_{run_name}_log_file_msa.txt")
        log_file_art = pathlib.Path(project_path, f"{time_stamp}_{run_name}_log_file_art.txt")

        if msa_cons:

            print(f"\n________________\n\nStarting MSA processing samples\n________________\n")
            print(f"min_depth = {min_depth}")
            with open(log_file_msa_temp, "a") as handle:
                handle.write(
                    f"\n________________\n\nStarting MSA processing samples\n________________\n")
                handle.write(f"\nmin_depth = {min_depth}\n")

            all_msa_sample_files = pathlib.Path(sample_folder).glob("*/msa/*.fastq")
            sample_no = 0
            for sample_fastq in all_msa_sample_files:
                sample_no += 1

                # get fastq path and name
                sample_path = pathlib.Path(sample_fastq).parent
                sample_name = pathlib.Path(sample_fastq).stem
                log_file_msa_sample = pathlib.Path(project_path, f"{time_stamp}_{sample_name}_log_file_msa_sample.txt")
                os.chdir(sample_path)

                # check free threads
                finished_threads = msa_threads*len(list(pathlib.Path(sample_folder).glob("*/msa/finished.txt")))
                free_threads = max_threads + finished_threads - used_threads
                print("\nfree_threads = " + str(free_threads))
                print('\n' + 'waiting for free threads')
                while free_threads < msa_threads:
                    time.sleep(5)
                    finished_threads = msa_threads*len(list(pathlib.Path(sample_folder).glob("*/msa/finished.txt")))
                    free_threads = max_threads + finished_threads - used_threads

                # check if fastq is present
                file_present = list(sample_path.glob("*.fastq"))
                if not file_present:
                    print(
                        f"\nCould not find the concatenated sample fastq file in msa folder: {sample_fastq}\nskipping sample")
                    with open(log_file_msa_sample, "a") as handle:
                        handle.write(
                            f"\nCould not find the concatenated sample fastq file in msa folder: {sample_fastq}\nskipping sample")
                    continue
                print(f"\n------->Running majority consensus pipeline for {sample_no} st/nd sample {sample_name} in new window\n")
                with open(log_file_msa_sample, "a") as handle:
                    handle.write(
                        f"\n\n------->Running majority consensus pipeline for {sample_no} st/nd sample {sample_name} in new window\n")

                # start majority consensus pipeline in new window
                majority_cmd = f"python ~/targetseq/msa_consensus.py -in {sample_fastq} -pf {plot_folder} -lf {log_file_msa_sample} " \
                               f"{use_bwa} -rs {chosen_ref_scheme} -bf {chosen_ref_scheme_bed_file} " \
                               f"-t {msa_threads} -d {min_depth} {use_gaps} -ac {all_samples_consens_seqs}"
                print(majority_cmd)
                try_except_continue_on_fail(f"gnome-terminal -- /bin/sh -c 'conda run -n target {majority_cmd} && touch finished.txt'")
                used_threads += msa_threads

        if artic:
            print(f"\n________________\n\nStarting ART processing samples\n________________\n")
            with open(log_file_art_temp, "a") as handle:
                handle.write(
                    f"\n\n________________\n\nStarting ART processing samples\n________________\n\n")

            all_art_sample_files = pathlib.Path(sample_folder).glob("*/art/*.fastq")
            sample_no = 0
            for sample_fastq in all_art_sample_files:
                sample_no += 1
                # get fastq path and name
                sample_path = pathlib.Path(sample_fastq).parent
                sample_name = pathlib.Path(sample_fastq).stem
                log_file_art_sample = pathlib.Path(project_path, f"{time_stamp}_{sample_name}_log_file_art_sample.txt")
                os.chdir(sample_path)

                # check free threads
                finished_threads = msa_threads*len(list(pathlib.Path(sample_folder).glob("*/msa/finished.txt")))+ artic_threads*len(list(pathlib.Path(sample_folder).glob("*/art/finished.txt")))
                free_threads = max_threads + finished_threads - used_threads
                print("\nfree_threads = " + str(free_threads))
                print('\n' + 'waiting for free threads')
                while free_threads < artic_threads:
                    time.sleep(5)
                    finished_threads = msa_threads*len(list(pathlib.Path(sample_folder).glob("*/msa/finished.txt")))+ artic_threads*len(list(pathlib.Path(sample_folder).glob("*/art/finished.txt")))
                    free_threads = max_threads + finished_threads - used_threads

                # check if fastq is present
                file_present = list(sample_path.glob("*.fastq"))
                if not file_present:
                    print(f"\nCould not find the concatenated sample fastq file in art folder: {sample_fastq}\nskipping sample")
                    with open(log_file_art_sample, "a") as handle:
                        handle.write(
                            f"\nCould not find the concatenated sample fastq file in art folder: {sample_fastq}\nskipping sample")
                    continue
                print(f"\n\n------->Running artic pipeline for {sample_no} st/nd sample {sample_name} in new window\n")
                with open(log_file_art_sample, "a") as handle:
                    handle.write(
                        f"\n------->Running artic pipeline for {sample_no} st/nd sample {sample_name} in new window\n\n")

                # start artic pipeline in new window
                artic_cmd = f"artic minion --normalise 500 --threads {artic_threads} --scheme-directory ~/artic-ncov2019/primer_schemes " \
                            f"--read-file {sample_fastq} --fast5-directory {fast5_dir} " \
                            f"--sequencing-summary {seq_summary_file} {scheme_name} {sample_name} " \
                            f"2>&1 | tee -a {log_file_art_sample}"
                print(artic_cmd)
                try_except_continue_on_fail(
                    f"gnome-terminal -- /bin/sh -c 'conda run -n artic-ncov2019 {artic_cmd} && touch finished.txt'")
                used_threads += artic_threads

            # write consensus to master consensus file when all artic files made
            print('\n' + 'waiting for art samples to be completed' + '\n')
            all_last_files = list(pathlib.Path(sample_folder).glob("*/art/finished.txt"))
            while len(all_last_files) < number_samples:
                time.sleep(5)
                all_last_files = list(pathlib.Path(sample_folder).glob("*/art/finished.txt"))
            else:
                all_cons_files = list(pathlib.Path(sample_folder).glob("*/art/*.consensus.fasta"))
                for cons_file in all_cons_files:
                    artic_d = fasta_to_dct(cons_file)
                    with open(all_samples_consens_seqs, 'a') as fh:
                        for name, seq in artic_d.items():
                            newname = re.sub("/ARTIC.*", "_art", name)
                            fh.write(f">{newname}\n{seq.replace('-', '')}\n")

        # construct log file and align all samples consens seqs
        if msa_cons:
            all_last_files = list(pathlib.Path(sample_folder).glob("*/msa/finished.txt"))
            print('\n' + 'waiting for msa samples to be completed' + '\n')
            while len(all_last_files) < number_samples:
                time.sleep(5)
                all_last_files = list(pathlib.Path(sample_folder).glob("*/msa/finished.txt"))

        # concat all log files
        os.chdir(project_path)

        if msa_cons:
            loglist_msa = []
            for path in pathlib.Path(project_path).glob("*_log_file_msa_sample.txt"):
                loglist_msa.append(str(path))
            sep = " "
            string_msa = sep.join(loglist_msa)
            cat_cmd = f"cat {str(log_file_msa_temp)} {string_msa} > {log_file_msa}"
            try_except_continue_on_fail(cat_cmd)
            for path in list(pathlib.Path(project_path).glob("*_log_file_msa_sample.txt")):
                os.remove(path)
            os.remove(log_file_msa_temp)


        if artic:
            loglist_art = []
            for path in pathlib.Path(project_path).glob("*_log_file_art_sample.txt"):
                loglist_art.append(str(path))
            sep = " "
            string_art = sep.join(loglist_art)
            cat_cmd = f"cat {str(log_file_art_temp)} {string_art} > {log_file_art}"
            try_except_continue_on_fail(cat_cmd)
            for path in list(pathlib.Path(project_path).glob("*_log_file_art_sample.txt")):
                os.remove(path)
            os.remove(log_file_art_temp)

        cat_cmd = f"cat {str(log_file)} {str(log_file_msa)} {str(log_file_art)} > {log_file_final}"
        try_except_continue_on_fail(cat_cmd)
        os.remove(log_file)
        if msa_cons:
            os.remove(log_file_msa)
        if artic:
            os.remove(log_file_art)

        print("aligning consensus sequence from all samples\n")
        tmp_file = pathlib.Path(project_path, "temp_aligned_file.fasta")
        mafft_cmd = f"mafft --thread -1 --auto {str(all_samples_consens_seqs)} > {str(tmp_file)}"
        print(mafft_cmd)
        run = try_except_continue_on_fail(mafft_cmd)
        if not run:
            print(f"could not align {all_samples_consens_seqs}")
            sys.exit("exiting")
        else:
            all_samples_consens_seqs.unlink()
            os.rename(str(tmp_file), str(all_samples_consens_seqs))

        # calculate & collect seq stats
        ref_name, ref_seq = list(py3_fasta_iter(chosen_ref_scheme))[0]
        ref_length = len(ref_seq)
        seqstats_outfile = pathlib.Path(project_path, f"{run_name}_seqstats.csv")
        all_consensus_d = fasta_to_dct(all_samples_consens_seqs)
        ref_lookup_name = list(all_consensus_d.keys())[0]
        aligned_ref = all_consensus_d[ref_lookup_name]
        sample_folder = pathlib.Path(project_path, "samples")
        del all_consensus_d[ref_lookup_name]
        with open(seqstats_outfile, 'w') as fh:
            fh.write("sample_name,genome_coverage,mean_depth\n")
            for v_name, v_seq in all_consensus_d.items():
                seqname = v_name[0:-11]
                mean_depth = ""
                if '_msa_' in v_name:
                    depth_file = csv.reader(open(f"{sample_folder}/{seqname}/msa/{seqname}_depth.csv", "r"))
                    for k, v in depth_file:
                        mean_depth = v
                seq_coverage = 0
                for i, base in enumerate(v_seq.upper()):
                    if base != "-" and base != "N" and aligned_ref[i] != "-":
                        seq_coverage += 1
                percent_coverage = round((seq_coverage / ref_length) * 100, 2)
                fh.write(f"{v_name},{percent_coverage},{mean_depth}\n")

    # run nanoplot
    os.chdir(sample_folder)
    nanoplotcmd= "for file in ./*/art/*.fastq;do folder=${file/.fastq/_nanoplot};nfolder=${folder/\/*\///}; NanoPlot --fastq $file -o ../nanoplot/$nfolder --tsv_stats;done"
    try_except_exit_on_fail(nanoplotcmd)

    # print end time
    now = datetime.datetime.now()
    date_time = now.strftime("%d/%m/%Y, %H:%M:%S")
    print(f"\nend time = {date_time}\n\n")
    with open(log_file_final, "a") as handle:
        handle.write(f"\nend time = {date_time}\n\n")

    print("sample processing completed\n")
    with open(log_file_final, "a") as handle:
        handle.write(f"\nsample processing completed\n\n")

    # compress fast5 files
    os.chdir(project_path)
    targzpath = pathlib.Path(project_path.parent, run_name + ".tar")
    fast5_dir_name = fast5_dir.parts[-1]
    seq_summary_file_name = pathlib.Path(seq_summary_file).name
    tarcmd = f"tar -cf {targzpath} {fast5_dir_name} {seq_summary_file_name}"
    try_except_exit_on_fail(tarcmd)
    zipcmd = f"pigz -7 -p 16 {targzpath}"
    try_except_exit_on_fail(zipcmd)
    print(tarcmd)
    with open(log_file_final, "a") as handle:
        handle.write(f"\n{tarcmd}\n\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process raw nanopore reads to fasta consensus sequences",
                                     formatter_class=Formatter)

    parser.add_argument("-in", "--project_path", default=argparse.SUPPRESS, type=str,
                        help="The path to the directory containing the 'fast5' and 'fastq' folders ", required=True)
    parser.add_argument("-r", "--reference", type=str, help="The reference genome and primer scheme to use",
                        choices=["ChikAsian_V1_400", "ChikECSA_V1_800", "ZikaAsian_V1_400", "SARS2_V1_800", "SARS2_V1_400", "RSVA_V1_3000", "RSVB_V1_3000", "DENV1_V1_400", "DENV2_V1_400"], required=True)
    parser.add_argument("-rs", "--reference_start", default=1, type=int,
                        help="The start coordinate of the reference sequence for read mapping", required=False)
    parser.add_argument("-re", "--reference_end", default=False, type=int,
                        help="The end coordinate of the reference sequence for read mapping. Default = full length",
                        required=False)
    parser.add_argument("-mi", "--min_len", type=int, default=700,
                        help="The minimum read length allowed:\n = 300 for 400bp amplicon design"
                                                             "\n = 700 for 800bp amplicon design", required=False)
    parser.add_argument("-ma", "--max_len", type=int, default=900,
                        help="The maximum read length allowed:\n = 500 for 400bp amplicon design"
                             "                                \n = 900 for 800bp amplicon design", required=False)
    parser.add_argument("-d", "--min_depth", type=int, default=100, help="The minimum coverage to call a position in "
                                                                         "the MSA to consensus", required=False)
    parser.add_argument("--run_step", default=0, type=int, required=False,
                        help="Run the pipeline starting at this step:\n"
                             "--run_step 0 = basecall reads with Guppy\n"
                             "--run_step 1 = demultiplex with Guppy\n"
                             "--run_step 2 = size filer and rename demultiplexed fastq file\n"
                             "--run_step 3 = concatenate demultiplexed files into sample files\n"
                             "--run_step 4 = run read mapping and all the variant calling steps on each sample\n")

    parser.add_argument("--run_step_only", default=False, action="store_true",
                        help="Only run the step specified in 'run_step'", required=False)
    parser.add_argument("-b", "--basecall_mode", default=1, choices=[0, 1], type=int,
                        help="0 = basecall in r10.4.1 kit14 mode\n"
                             "1 = basecall in r9.4.1 mode\n", required=False)
    parser.add_argument("-m", "--msa", default=False, action="store_true",
                        help="Generate consensus from MSA", required=False)
    parser.add_argument("-a", "--art", default=False, action="store_true",
                        help="Generate consensus with Artic pipeline", required=False)
    parser.add_argument("-c", "--cpu_threads", type=int, default=16, choices=range(0, 16),
                        help="The number of cpu threads to use for bwa, nanopolish etc...", required=False)
    parser.add_argument("-g", "--gpu_threads", type=int, default=8,
                        help="The number of gpu threads to use ...", required=False)
    parser.add_argument("-gb", "--gpu_buffers", type=int, default=15,
                        help="The number of gpu buffers to use for demultiplexing", required=False)
    parser.add_argument("--use_gaps", default='', action="store_const", const='-ug',
                        help="use gap characters when making the consensus sequences", required=False)
    parser.add_argument("--use_bwa", default='', action="store_const", const='-bwa',
                        help="use bwa instead of minimap2 to map reads to reference", required=False)
    parser.add_argument("-p", "--guppy_path", default=argparse.SUPPRESS, type=str,
                        help="The path to the guppy executables eg: '.../ont-guppy/bin/'", required=True)
    parser.add_argument("-rt", "--real_time", default=False, action="store_true",
                        help="start basecalling fast5 files in batches during sequencing", required=False)

    args = parser.parse_args()

    project_path = args.project_path
    reference = args.reference
    reference_start = args.reference_start
    reference_end = args.reference_end
    min_len = args.min_len
    max_len = args.max_len
    min_depth = args.min_depth
    run_step = args.run_step
    run_step_only = args.run_step_only
    basecall_mode = args.basecall_mode
    msa_cons = args.msa
    artic = args.art
    cpu_threads = args.cpu_threads
    gpu_threads = args.gpu_threads
    gpu_buffers = args.gpu_buffers
    use_gaps = args.use_gaps
    use_bwa = args.use_bwa
    guppy_path = args.guppy_path
    real_time = args.real_time

    main(project_path, reference, reference_start, reference_end, min_len, max_len, min_depth, run_step,
         run_step_only, basecall_mode, msa_cons, artic, cpu_threads, gpu_threads, gpu_buffers, use_gaps, use_bwa,
         guppy_path, real_time)
