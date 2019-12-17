import argparse
import pathlib
from src.misc_functions import try_except_continue_on_fail


__author__ = 'Colin Anthony'


class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass


def main(inpath, guppy_path, outpath, threads, gpu_threads):
    # force absolute file paths
    inpath = pathlib.Path(inpath).absolute()
    outpath = pathlib.Path(outpath).absolute()
    guppy_path = pathlib.Path(guppy_path).absolute()
    guppy_demultiplexer = pathlib.Path(guppy_path, "guppy_barcoder")
    if gpu_threads == 0:
        gpu_settings = ''
    else:
        gpu_settings = f"-x 'auto'"

    guppy_demux_cmd = f"{str(guppy_demultiplexer)} -i {inpath} -s {outpath} -t {threads} --compress_fastq" \
                      f"--require_barcodes_both_ends --trim_barcodes --num_barcoding_buffers 2 --records_per_fastq 0" \
                      f"{gpu_settings}"

    run = try_except_continue_on_fail(guppy_demux_cmd)

    if run:
        print("demultiplexing completed\n")
    else:
        print("demultiplexing failed")

    return run


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='',
                                     formatter_class=Formatter)

    parser.add_argument('-in', '--inpath', type=str, default=None, required=True,
                        help='The path to the fastq folder')
    parser.add_argument('-p', '--guppy_path', type=str, default=None, required=True,
                        help='The path to guppy exexutable')
    parser.add_argument('-o', '--outpath', type=str, default=None, required=True,
                        help='The path for the outfile')
    parser.add_argument("-t", "--threads", type=int, default=4,
                        help="The number of threads to use for demultiplexing, bwa, nanopolish etc...", required=False)
    parser.add_argument("-g", "--gpu_threads", type=int, default=4,
                        help="The number of gpu threads to use, use '0' if no gpu available", required=False)
    args = parser.parse_args()
    inpath = args.inpath
    guppy_path = args.guppy_path
    outpath = args.outpath
    threads = args.threads
    gpu_threads = args.gpu_threads

    main(inpath, guppy_path, outpath, threads, gpu_threads)
