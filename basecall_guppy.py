import argparse
import pathlib
from src.misc_functions import try_except_continue_on_fail


__author__ = 'Colin Anthony'


class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass


def main(inpath, guppy_path, outpath, gpu_threads, bascall_mode):
    # force absolute file paths
    inpath = pathlib.Path(inpath).absolute()
    outpath = pathlib.Path(outpath).absolute()
    guppy_path = pathlib.Path(guppy_path).absolute()
    guppy_basecaller = pathlib.Path(guppy_path, "guppy_basecaller")
    cuda_device = "CUDA:0"
    config_option = ["dna_r9.4.1_450bps_fast.cfg ", "dna_r9.4.1_450bps_hac.cfg"]
    config = config_option[bascall_mode]
    if gpu_threads == 0:
        gpu_settings = ""
    else:
        gpu_settings = f"--gpu_runners_per_device {gpu_threads}  --num_callers 4 -x 'auto' " #--device {cuda_device}

    # add arg to keep fastq size to 4000 seqs in case something fails, easier to find where it went wrong and fix
    guppy_basecall_cmd = f"{str(guppy_basecaller)} -i {inpath} -r -s {outpath} -c {config} " \
                         f"--compress_fastq --records_per_fastq 4000 --qscore_filtering 7 " \
                         f"{gpu_settings}"

    run = try_except_continue_on_fail(guppy_basecall_cmd)

    if run:
        print("basecalling completed\n")
    else:
        print("basecalling failed")

    return run


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='',
                                     formatter_class=Formatter)

    parser.add_argument('-in', '--inpath', type=str, default=None, required=True,
                        help='The path to the fast5 folder')
    parser.add_argument('-p', '--guppy_path', type=str, default=None, required=True,
                        help='The path to guppy exexutable')
    parser.add_argument('-o', '--outpath', type=str, default=None, required=True,
                        help='The path for the outfile')
    parser.add_argument("-g", "--gpu_threads", type=int, default=4,
                        help="The number of gpu threads to use, use '0' if no gpu available", required=False)
    parser.add_argument('-b', '--bascall_mode', type=int, choices=[0, 1], default=0, required=False,
                        help='0 = Fast mode\n'
                             '1 = high accuracy mode')
    args = parser.parse_args()
    inpath = args.inpath
    guppy_path = args.guppy_path
    outpath = args.outpath
    gpu_threads = args.gpu_threads
    bascall_mode = args.bascall_mode

    main(inpath, guppy_path, outpath, gpu_threads, bascall_mode)
