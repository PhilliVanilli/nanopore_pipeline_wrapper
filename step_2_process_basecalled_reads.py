import os, time
import sys
import argparse
import re
import pathlib
import datetime
import pandas as pd
import shutil
from src.misc_functions import try_except_continue_on_fail
from src.misc_functions import try_except_exit_on_fail
from src.misc_functions import py3_fasta_iter
from src.misc_functions import cat_sample_names
from src.misc_functions import filter_length
from src.misc_functions import fasta_to_dct
from basecall_guppy import main as gupppy_basecall
from demultiplex_guppy import main as guppy_demultiplex
from all_samples_summary import main as sample_summary


__author__ = 'Colin Anthony'


class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass


def main(project_path, reference, ref_start, ref_end, min_len, max_len, min_depth, run_step,
         rerun_step_only, basecall_mode, msa_cons, cpu_cores, gpu_cores, gpu_buffers, use_gaps, use_bwa,
         guppy_path, real_time):

    threads = cpu_cores
    # set threads
    # threads = int()
    # if msa_cons:
    #     threads = cpu_cores - 8
    # else:
    #     threads = cpu_cores

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
    # sequencing_summary_file = pathlib.Path(fastq_dir, "sequencing_summary.txt")
    sample_names = pathlib.Path(project_path, "sample_names.csv")
    if not sample_names:
        sys.exit("Could not find sample_names.csv in porject folder")
    demultiplexed_folder = pathlib.Path(project_path, "demultiplexed")
    sample_folder = pathlib.Path(project_path, "samples")
    print(sample_folder)
    # master_reads_file = pathlib.Path(project_path, run_name + "_all.fastq")
    time_stamp = str('{:%Y-%m-%d_%H_%M}'.format(datetime.datetime.now()))
    log_file = pathlib.Path(project_path, f"{time_stamp}_{run_name}_log_file.txt")

    with open(log_file, "w") as handle:
        handle.write(f"# start of pipeline run for project: {run_name}\n")

    now = datetime.datetime.now()
    date_time = now.strftime("%d/%m/%Y, %H:%M:%S")
    print(f"\nstart time = {date_time}\n\n")
    with open(log_file, "a") as handle:
        handle.write(f"\nstart time = {date_time}\n\n")

    # set dir to project dir so that output is written in correct place by external tools
    os.chdir(project_path)

    # set the reference genome
    reference_scheme = \
        {"ChikECSA_V1_800": pathlib.Path(primer_scheme_dir, "ChikECSA800", "V1", "ChikECSA800.reference.fasta"),
         "ChikAsian_V1_400": pathlib.Path(primer_scheme_dir, "ChikAsian400", "V1", "ChikAsian400.reference.fasta"),
         "ZikaAsian_V1_400": pathlib.Path(primer_scheme_dir, "ZikaAsian400", "V1", "ZikaAsian400.reference.fasta"),
         "SARS2_V1_800": pathlib.Path(primer_scheme_dir, "SARS2_800", "V1", "SARS2_800.reference.fasta"),
         "SARS2_V1_400": pathlib.Path(primer_scheme_dir, "SARS2_400", "V1", "SARS2_400.reference.fasta"),
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
    print(f"\nReference is {chosen_ref_scheme}\n")
    print(f"\nPrimer bed file is {chosen_ref_scheme_bed_file}\n")
    with open(log_file, "a") as handle:
        handle.write(f"\nReference is {chosen_ref_scheme}\nPrimer bed file is {chosen_ref_scheme_bed_file}\n")

    if run_step == 0:
        run = gupppy_basecall(fast5_dir, guppy_path, fastq_dir, gpu_cores, basecall_mode, real_time, reference, script_folder)
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
        run = guppy_demultiplex(fastq_dir, guppy_path, demultiplexed_folder, threads, gpu_buffers, gpu_cores)
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

            cat_outfile = pathlib.Path(sample_dir, f"{sample_name}.fastq")
            cat_cmd = f"cat {str(barcode_1_file)} {str(barcode_2_file)} > {cat_outfile}"
            print(cat_cmd)
            run = try_except_continue_on_fail(cat_cmd)
            if not run:
                print("missing one or more demultiplexed files for this sample")
                with open(log_file, "a") as handle:
                    handle.write("\nmissing one or more demultiplexed files for this sample\n")
                continue
        for fastq in demultiplexed_folder.glob('*.fastq'):
            os.remove(str(fastq))
        if not rerun_step_only:
            run_step = 4
        else:
            sys.exit("Run step only completed, exiting")

    if run_step == 4:
        print("Running variant calling on samples")

        with open(log_file, "a") as handle:
            handle.write(f"\nRunning variant calling on samples\n")
        if use_bwa:
            make_index_cmd = f"bwa index {chosen_ref_scheme}"
            with open(log_file, "a") as handle:
                handle.write(f"\n{make_index_cmd}\n")

            try_except_exit_on_fail(make_index_cmd)

        all_sample_files = pathlib.Path(sample_folder).glob("*/*.fastq")
        number_samples = len(list(sample_folder.glob('*/*.fastq')))

        # make variable for project file containing all samples' consensus sequences
        project_name = project_path.parts[-1]
        all_samples_consens_seqs = pathlib.Path(project_path, project_name + "_all_samples.fasta")

        # initialize the file, and add reference to all consensus file
        with open(all_samples_consens_seqs, 'w') as fh:
            fh.write(f">{ref_name}\n{ref_seq}\n")
        p = pathlib.Path(project_path, project_name + '_mapping.csv')
        with open(p, 'w') as fh:
            fh.close()

        samples_run = 1
        old_number_png_files = 0
        for sample_fastq in all_sample_files:
            # get folder paths
            sample_folder = pathlib.Path(sample_fastq).parent
            sample_name = pathlib.Path(sample_fastq).stem
            os.chdir(sample_folder)
            seq_summary_file_name = ""
            for file in project_path.glob('sequencing_summary*.txt'):
                seq_summary_file_name = file
            seq_summary_file = pathlib.Path(seq_summary_file_name).resolve()
            artic_folder = pathlib.Path(sample_folder, "artic")
            if os.path.exists(artic_folder):
                shutil.rmtree(artic_folder)
            artic_folder.mkdir(mode=0o777, parents=True, exist_ok=True)

            # check if fastq is present
            if not sample_fastq.is_file():
                print(f"\nCould not find the concatenated sample fastq file: {sample_fastq}\nskipping sample")
                with open(log_file, "a") as handle:
                    handle.write(
                        f"\nCould not find the concatenated sample fastq file: {sample_fastq}\nskipping sample")
                continue
            print(f"\n________________\n\nStarting processing sample: {sample_name}\n________________\n")
            with open(log_file, "a") as handle:
                handle.write(
                    f"\n________________\n\nStarting processing sample: {sample_name}\n________________\n")

            # start artic pipeline in new window
            print(f"\n------->Running Artic's pipeline in new window\n")
            with open(log_file, "a") as handle:
                handle.write(
                    f"\n------->Running Artic's pipeline in new window\n\n")

            artic_cmd = f"artic minion --normalise 400 --threads {threads} --scheme-directory ~/artic-ncov2019/primer_schemes " \
                        f"--read-file {sample_fastq} --fast5-directory {fast5_dir} " \
                        f"--sequencing-summary {seq_summary_file} {scheme_name} {sample_name} " \
                        f"2>&1 | tee -a {log_file}"
            print(artic_cmd)
            try_except_continue_on_fail(
                f"gnome-terminal -- /bin/sh -c 'conda run -n artic-ncov2019 {artic_cmd}'")

            last_file_made = pathlib.Path(sample_folder, sample_name + ".muscle.out.fasta")
            while pathlib.Path.exists(last_file_made) == False:
                time.sleep(5)
            else:
                time.sleep(2)
                all_files = os.listdir(sample_folder)

                # write consensus to master consensus file
                artic_cons_file = pathlib.Path(sample_folder, sample_name + ".consensus.fasta")
                artic_d = fasta_to_dct(artic_cons_file)
                with open(all_samples_consens_seqs, 'a') as fh:
                    for name, seq in artic_d.items():
                        newname = re.sub("/ARTIC.*", "_art", name)
                        fh.write(f">{newname}\n{seq.replace('-', '')}\n")

                for filename in all_files:
                    if os.path.isfile(filename) and not filename.endswith('.fastq'):
                        file = os.path.join(sample_folder, filename)
                        shutil.move(file, artic_folder)

            # start majority consensus pipeline in new window
            if msa_cons:

                print(f"\n\n------->Running majority consensus pipeline in new window\n")
                with open(log_file, "a") as handle:
                    handle.write(
                        f"\n\n------->Running majority consensus pipeline in new window\n")

                majority_cmd = f"python ~/nanopore_pipeline_wrapper/msa_consensus.py -in {sample_fastq} -pf {plot_folder} -lf {log_file} " \
                               f"{use_bwa} -rs {chosen_ref_scheme} -bf {chosen_ref_scheme_bed_file} " \
                               f"-t {threads} -d {min_depth} {use_gaps} -ac {all_samples_consens_seqs}"

                print(majority_cmd)
                try_except_continue_on_fail(f"gnome-terminal -- /bin/sh -c 'conda run -n nanop {majority_cmd}'")

                # open(f"{sample_name}_msa_from_bam_file.fasta", "w+")
                last_file_made_2 = pathlib.Path(sample_folder, sample_name + "_msa_from_bam_file.fasta")
                while pathlib.Path.exists(last_file_made_2) == False:
                    time.sleep(5)
                else:
                    if samples_run + 1 <= number_samples:
                        print(f"\ncontinuing with sample {samples_run + 1}\n")

                # keep threads balanced
                number_png_files = len(list(plot_folder.glob('*_sequencing_depth.png')))
                print(f'{number_png_files} png files created')
                difference = number_png_files - old_number_png_files
                old_number_png_files = number_png_files
                threads = threads - 1 + difference
                samples_run += 1

        # run sample summary as soon as all sequencing_depth.png files made
        number_pngs = len(list(plot_folder.glob('*_sequencing_depth.png')))
        if number_pngs < number_samples:
            print('waiting for all msa to be completed')
        while number_pngs < number_samples:
            time.sleep(10)
            number_pngs = len(list(plot_folder.glob('*_sequencing_depth.png')))
        else:
            sample_summary(project_path, all_samples_consens_seqs, chosen_ref_scheme, run_name)


    now = datetime.datetime.now()
    date_time = now.strftime("%d/%m/%Y, %H:%M:%S")
    print(f"\nend time = {date_time}\n\n")
    with open(log_file, "a") as handle:
        handle.write(f"\nend time = {date_time}\n\n")

    print("sample processing completed\n")
    with open(log_file, "a") as handle:
        handle.write(f"\nsample processing completed\n\n")

    # compress fast5 files
    targzpath = pathlib.Path(project_path.parent, run_name + ".tar.gz")
    tarcmd = f"tar cf - {fast5_dir} | pigz -7 -p 16  > {targzpath}"
    try_except_exit_on_fail(tarcmd)
    print(tarcmd)
    with open(log_file, "a") as handle:
        handle.write(f"\n{tarcmd}\n\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process raw nanopore reads to fasta consensus sequences",
                                     formatter_class=Formatter)

    parser.add_argument("-in", "--project_path", default=argparse.SUPPRESS, type=str,
                        help="The path to the directory containing the 'fast5' and 'fastq' folders ", required=True)
    parser.add_argument("-r", "--reference", type=str, help="The reference genome and primer scheme to use",
                        choices=["ChikAsian_V1_400", "ChikECSA_V1_800", "ZikaAsian_V1_400", "SARS2_V1_800", "SARS2_V1_400", "DENV1_400", "DENV2_400"], required=True)
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
                        help="0 = basecall in fast mode\n"
                             "1 = basecall in high accuracy mode\n", required=False)
    parser.add_argument("-m", "--msa", default=False, action="store_true",
                        help="Generate consensus from MSA", required=False)
    parser.add_argument("-c", "--cpu_cores", type=int, default=15, choices=range(0, 16),
                        help="The number of threads to use for bwa, nanopolish etc...", required=False)
    parser.add_argument("-g", "--gpu_cores", type=int, default=8,
                        help="The number of gpu threads to use ...", required=False)
    parser.add_argument("-gb", "--gpu_buffers", type=int, default=16,
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
    cpu_cores = args.cpu_cores
    gpu_cores = args.gpu_cores
    gpu_buffers = args.gpu_buffers
    use_gaps = args.use_gaps
    use_bwa = args.use_bwa
    guppy_path = args.guppy_path
    real_time = args.real_time

    main(project_path, reference, reference_start, reference_end, min_len, max_len, min_depth, run_step,
         run_step_only, basecall_mode, msa_cons, cpu_cores, gpu_cores, gpu_buffers, use_gaps, use_bwa,
         guppy_path, real_time)
