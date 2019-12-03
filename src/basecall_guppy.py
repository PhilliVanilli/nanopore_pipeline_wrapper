import argparse
import pathlib
from src.misc_functions import try_except_exit_on_fail


__author__ = 'Colin Anthony'


class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass


def main(inpath, guppy_path, outpath, threads, bascall_mode):
    # force absolute file paths
    inpath = pathlib.Path(inpath).absolute()
    outpath = pathlib.Path(outpath).absolute()
    guppy_path = pathlib.Path(guppy_path).absolute()
    guppy_basecaller = pathlib.Path(guppy_path, "guppy_basecaller")
    cuda_device = "CUDA:0"
    gpu_threads = 4
    config_option = ["dna_r9.4.1_450bps_fast.cfg ", "dna_r9.4.1_450bps_hac.cfg"]
    config = config_option[bascall_mode]
    # add arg to keep fastq size to 4000 seqs in case something fails, easier to find where it went wrong and fix
    guppy_basecall_cmd = f"{str(guppy_basecaller)} -i {inpath} -r -s {outpath} -t {threads} -c {config} " \
                         f"--device {cuda_device} --gpu_runners_per_device {gpu_threads}  --num_callers 4 -x 'auto'" \
                         f"--compress_fastq --records_per_fastq 4000 --qscore_filtering 7 "

    try_except_exit_on_fail(guppy_basecall_cmd)

    print("done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='',
                                     formatter_class=Formatter)

    parser.add_argument('-in', '--inpath', type=str, default=None, required=True,
                        help='The path to the fast5 folder')
    parser.add_argument('-p', '--guppy_path', type=str, default=None, required=True,
                        help='The path to guppy exexutable')
    parser.add_argument('-o', '--outpath', type=str, default=None, required=True,
                        help='The path for the outfile')
    parser.add_argument("-t", "--threads", type=int, default=8,
                        help="The number of threads to use for basecalling", required=False)
    parser.add_argument('-b', '--bascall_mode', type=int, choices=[0, 1], default=0, required=False,
                        help='0 = Fast mode\n'
                             '1 = high accuracy mode')
    args = parser.parse_args()
    inpath = args.inpath
    guppy_path = args.guppy_path
    outpath = args.outpath
    threads = args.threads
    bascall_mode = args.bascall_mode

    main(inpath, guppy_path, outpath, threads, bascall_mode)
