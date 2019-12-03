import argparse
import pathlib
from src.misc_functions import try_except_exit_on_fail
# from tensorflow.python.client import device_lib


__author__ = 'Colin Anthony'


class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass


def main(infile, guppy_path, outpath, threads):
    # force absolute file paths
    infile = pathlib.Path(infile).absolute()
    outpath = pathlib.Path(outpath).absolute()
    guppy_path = pathlib.Path(guppy_path).absolute()
    # gpu_devices = device_lib.list_local_devices()
    guppy_demultiplexer = pathlib.Path(guppy_path, "guppy_barcoder")
    gpu_device = ""
    gpu_threads = 4
    guppy_demux_cmd = f"{str(guppy_demultiplexer)} -i {infile} -s {outpath}  -t {threads} " \
                      f"--require_barcodes_both_ends --trim_barcodes " \
                      f"--num_barcoding_buffers 2 --gpu_runners_per_device {gpu_threads} -x 'auto'"

    try_except_exit_on_fail(guppy_demux_cmd)

    print("done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='',
                                     formatter_class=Formatter)

    parser.add_argument('-in', '--infile', type=str, default=None, required=True,
                        help='The path and name of the infile')
    parser.add_argument('-p', '--guppy_path', type=str, default=None, required=True,
                        help='The path to guppy exexutable')
    parser.add_argument('-o', '--outpath', type=str, default=None, required=True,
                        help='The path for the outfile')
    parser.add_argument("-t", "--threads", type=int, default=8,
                        help="The number of threads to use for demultiplexing, bwa, nanopolish etc...", required=False)

    args = parser.parse_args()
    infile = args.infile
    guppy_path = args.guppy_path
    outpath = args.outpath
    threads = args.threads

    main(infile, guppy_path, outpath, threads)
