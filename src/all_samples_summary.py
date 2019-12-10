import argparse
import pathlib
import os
import sys
from src.misc_functions import try_except_continue_on_fail
from src.misc_functions import fasta_to_dct
from src.misc_functions import py3_fasta_iter
__author__ = 'Colin Anthony'


class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass


def main(project_path, all_samples_consens_seqs, chosen_ref_scheme, run_name):

    print("aligning consensus sequence from all samples\n")
    tmp_file = pathlib.Path(project_path, "temp_aligned_file.fasta")
    mafft_cmd = f"mafft {str(all_samples_consens_seqs)} > {str(tmp_file)}"

    ref_name, ref_seq = list(py3_fasta_iter(chosen_ref_scheme))[0]
    print(mafft_cmd)
    run = try_except_continue_on_fail(mafft_cmd)
    if not run:
        print(f"could not align {all_samples_consens_seqs}")
        sys.exit("exiting")
    else:
        all_samples_consens_seqs.unlink()
        os.rename(tmp_file, str(all_samples_consens_seqs))

        # calculate coverage
        ref_length = len(ref_seq)
        coverage_outfile = pathlib.Path(project_path, f"{run_name}_genome_coverage.csv")
        all_consensus_d = fasta_to_dct(all_samples_consens_seqs)
        ref_lookup_name = list(all_consensus_d.keys())[0]
        aligned_ref = all_consensus_d[ref_lookup_name]
        del all_consensus_d[ref_lookup_name]
        with open(coverage_outfile, 'w') as fh:
            fh.write("sample_name,genome_coverage\n")
            for v_name, v_seq in all_consensus_d.items():
                seq_coverage = 0
                for i, base in enumerate(v_seq.upper()):
                    if base != "-" and base != "N" and aligned_ref[i] != "-":
                        seq_coverage += 1
                percent_coverage = round((seq_coverage / ref_length) * 100, 2)
                fh.write(f"{v_name},{percent_coverage}\n")

    print("done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script is meant to post process the consensus sequences from the'
                                                 'nanopore_pipeline, collecting consensus sequences from all samples'
                                                 'aligning them together with the mapping reference, and outputting '
                                                 '% coverage stats for each sample',
                                     formatter_class=Formatter)

    parser.add_argument('-in', '--project_path', type=str, default=None, required=True,
                        help='The path to the project folder')
    parser.add_argument('-a', '--all_samples_consens_seqs', type=str, default=None, required=True,
                        help='The path name of the file to contain consensus sequences from all samples')
    parser.add_argument('-r', '--chosen_ref_scheme', type=str, default=None, required=True,
                        help='The path and name of the reference fasta file')
    parser.add_argument('-n', '--run_name', type=str, default=None, required=True,
                        help='The name of the sequencing run')

    args = parser.parse_args()
    project_path = args.project_path
    all_samples_consens_seqs = args.all_samples_consens_seqs
    chosen_ref_scheme = args.chosen_ref_scheme
    run_name = args.run_name

    main(project_path, all_samples_consens_seqs, chosen_ref_scheme, run_name)
