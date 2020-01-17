import os
import sys
import argparse
import pathlib
import datetime
import pandas as pd
from src.misc_functions import try_except_continue_on_fail
from src.misc_functions import try_except_exit_on_fail
from src.misc_functions import py3_fasta_iter
from src.misc_functions import cat_sample_names
from src.misc_functions import filter_length
from basecall_guppy import main as gupppy_basecall
from demultiplex_guppy import main as guppy_demultiplex
from analyse_sample import main as sample_analysis
from all_samples_summary import main as sample_summary

__author__ = 'Colin Anthony'


class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass


def main(project_path, sample_names, reference, ref_start, ref_end, min_len, max_len, min_depth, run_step,
         rerun_step_only, basecall_mode, msa_cons_only, threads, gpu_cores, gpu_buffers, use_gaps, use_minmap2,
         guppy_path):

    # set the primer_scheme directory
    script_folder = pathlib.Path(__file__).absolute().parent
    primer_scheme_dir = pathlib.Path(script_folder, "primer-schemes")

    # get folder paths
    project_path = pathlib.Path(project_path).absolute()
    plot_folder = pathlib.Path(project_path, "seq_depth_plots")
    plot_folder.mkdir(mode=0o777, parents=True, exist_ok=True)
    run_name = project_path.parts[-1]
    fast5_dir = pathlib.Path(project_path, "fast5")
    fastq_dir = pathlib.Path(project_path, "fastq")
    sequencing_summary_file = pathlib.Path(fastq_dir, "sequencing_summary.txt")
    sample_names = pathlib.Path(sample_names).absolute()
    demultiplexed_folder = pathlib.Path(project_path, "demultiplexed")
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
    if not ref_start or ref_start == 0:
        ref_start = 1
    if not ref_end or ref_end > len(ref_seq):
        ref_end = len(ref_seq)
    reference_slice = f'{ref_name}:{ref_start}-{ref_end}'
    print(f"\nReference is {chosen_ref_scheme}\n")
    print(f"\nPrimer bed file is {chosen_ref_scheme_bed_file}\n")
    with open(log_file, "a") as handle:
        handle.write(f"\nReference is {chosen_ref_scheme}\nPrimer bed file is {chosen_ref_scheme_bed_file}\n")

    if run_step == 0:
        run = gupppy_basecall(fast5_dir, guppy_path, fastq_dir, gpu_cores, basecall_mode)
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
        if not list(fastq_dir.glob("*.fastq")):
            fastq_dir = pathlib.Path(fastq_dir, "pass")
        
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

                if not filtered_file:
                    print(f"no sequences in file after length filtering for {concat_outfile}\n")
                # bash_cmd = f"cat {concat_outfile} | paste - - - - | awk 'length($2)  >= {min_len} && length($2) <= {max_len}' | sed 's/\t/\n/g' > {str(new_name)}"
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

                # bash_cmd = f"cat {file} | paste - - - - | awk 'length($2)  >= {min_len} && length($2) <= {max_len}' | sed 's/\t/\n/g' > {str(new_name)}"
                # print(bash_cmd)
                # seqmagick_cmd = f"seqmagick quality-filter --min-mean-quality 0 " \
                #                 f"--min-length {min_len} --max-length {max_len} " \
                #                 f"{file} {new_name} "
                # vsearch_cmd = f"vsearch --fastq_filter {file} -fastq_maxlen {max_len} --fastq_minlen {min_len} " \
                #               f"--fastq_qmax 100 --fastqout {new_name}"
                # try_except_exit_on_fail(bash_cmd)

            # os.rmdir(str(folder))

        if not rerun_step_only and not msa_cons_only:
            run_step = 3
        elif not rerun_step_only and msa_cons_only:
            run_step = 4
        elif rerun_step_only:
            sys.exit("filer demiltiplexed files and rename them completed, exiting")
        else:
            sys.exit("filtering and renaming demultiplexed files failed")

    if run_step == 3 and not msa_cons_only:
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
    else:
        run_step = 4

    if run_step == 4:
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

        if not rerun_step_only:
            run_step = 5
        else:
            sys.exit("Run step only completed, exiting")

    if run_step == 5:
        print("Running variant calling on samples")
        with open(log_file, "a") as handle:
            handle.write(f"\nRunning variant calling on samples\n")
        if not use_minmap2:
            make_index_cmd = f"bwa index {chosen_ref_scheme}"
            with open(log_file, "a") as handle:
                handle.write(f"\n{make_index_cmd}\n")

            try_except_exit_on_fail(make_index_cmd)

        all_sample_files = pathlib.Path(sample_folder).glob("*/*.fastq")

        # make variable for project file containing all samples' consensus sequences
        project_name = project_path.parts[-1]
        all_samples_consens_seqs = pathlib.Path(project_path, project_name + "_all_samples.fasta")

        # initialize the file, and add reference to all consensus file
        with open(all_samples_consens_seqs, 'w') as fh:
            fh.write(f">{ref_name}\n{ref_seq}\n")

        for sample_fastq in all_sample_files:
            if not sample_fastq.is_file():
                print(f"could not find the concatenated sample fastq file: {sample_fastq}\nskipping sample")
                with open(log_file, "a") as handle:
                    handle.write(f"could not find the concatenated sample fastq file: {sample_fastq}\nskipping sample")
                continue
            run = sample_analysis(sample_fastq, plot_folder, log_file, use_minmap2, chosen_ref_scheme,
                                  chosen_ref_scheme_bed_file, threads, msa_cons_only, min_depth, use_gaps,
                                  all_samples_consens_seqs, reference_slice, master_reads_file)
            if not run:
                continue

        # align the master consensus file
        sample_summary(project_path, all_samples_consens_seqs, chosen_ref_scheme, run_name)

    print("sample processing completed\n")
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
    parser.add_argument("-mi", "--min_len", type=int, default=700,
                        help="The minimum read length allowed:\n = 300 for 400bp amplicon design"
                                                             "\n = 700 for 800bp amplicon design", required=False)
    parser.add_argument("-ma", "--max_len", type=int, default=900,
                        help="The maximum read length allowed:\n = 500 for 400bp amplicon design"
                             "                                \n = 900 for 800bp amplicon design", required=False)
    parser.add_argument("-d", "--min_depth", type=int, default=100, help="The minimum coverage to call a position in "
                                                                         "the MSA to consensus", required=False)
    parser.add_argument("--run_step", default=1, type=int, required=False,
                        help="Run the pipeline starting at this step:\n"
                             "--run_step 0 = basecall reads with Guppy\n"
                             "--run_step 1 = demultiplex with Guppy\n"
                             "--run_step 2 = size filer and rename demultiplexed fastq file\n"
                             "--run_step 3 = index master fastq file\n"
                             "--run_step 4 = concatenate demultiplexed files into sample files\n"
                             "--run_step 5 = run read mapping and all the variant calling steps on each sample\n")

    parser.add_argument("--run_step_only", default=False, action="store_true",
                        help="Only run the step specified in 'run_step'", required=False)
    parser.add_argument("-b", "--basecall_mode", default=0, choices=[0, 1], type=int,
                        help="0 = basecall in fast mode\n"
                             "1 = basecall in high accuracy mode\n", required=False)
    parser.add_argument("-m", "--msa_cons_only", default=False, action="store_true",
                        help="Only do MSA to consensus sequence for variant calling", required=False)
    parser.add_argument("-t", "--threads", type=int, default=8,
                        help="The number of threads to use for bwa, nanopolish etc...", required=False)
    parser.add_argument("-g", "--gpu_cores", type=int, default=4,
                        help="The number of gpu threads to use ...", required=False)
    parser.add_argument("-gb", "--gpu_buffers", type=int, default=16,
                        help="The number of gpu buffers to use for demultiplexing", required=False)
    parser.add_argument("--use_gaps", default=False, action="store_true",
                        help="use gap characters when making the consensus sequences", required=False)
    parser.add_argument("--use_minmap2", default=False, action="store_true",
                        help="use bwa instead of minimap2 to map reads to reference", required=False)
    parser.add_argument("-p", "--guppy_path", default=argparse.SUPPRESS, type=str,
                        help="The path to the guppy executables eg: '.../ont-guppy/bin/'", required=True)

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
    run_step_only = args.run_step_only
    basecall_mode = args.basecall_mode
    msa_cons_only = args.msa_cons_only
    threads = args.threads
    gpu_cores = args.gpu_cores
    gpu_buffers = args.gpu_buffers
    use_gaps = args.use_gaps
    use_minmap2 = args.use_minmap2
    guppy_path = args.guppy_path

    main(project_path, sample_names, reference, reference_start, reference_end, min_len, max_len, min_depth, run_step,
         run_step_only, basecall_mode, msa_cons_only, threads, gpu_cores, gpu_buffers, use_gaps, use_minmap2,
         guppy_path)
