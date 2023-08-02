import argparse
import pathlib
from src.misc_functions import try_except_continue_on_fail


__author__ = 'Colin Anthony'


class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass


def main(inpath, guppy_path, outpath, threads, gpu_buffers, use_gpu):
    # force absolute file paths
    inpath = pathlib.Path(inpath).absolute()
    outpath = pathlib.Path(outpath).absolute()
    guppy_path = pathlib.Path(guppy_path).absolute()
    guppy_demultiplexer = pathlib.Path(guppy_path, "guppy_barcoder")
    if not use_gpu:
        gpu_settings = ''
    else:
        gpu_settings = f"-x 'auto'"

    guppy_demux_cmd = f"{str(guppy_demultiplexer)} -i {inpath} -s {outpath} -t {threads} " \
                      f"--require_barcodes_both_ends --enable_trim_barcodes --num_barcoding_buffers {gpu_buffers} " \
                      f"--records_per_fastq 0 {gpu_settings}"

    run = try_except_continue_on_fail(guppy_demux_cmd)

    if run:
        print("demultiplexing completed\n")
    else:
        print("demultiplexing failed")

    return run


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Script to wrapp guppy demultiplex commands',
                                     formatter_class=Formatter)

    parser.add_argument('-in', '--inpath', type=str, default=None, required=True,
                        help='The path to the fastq folder')
    parser.add_argument('-p', '--guppy_path', type=str, default=None, required=True,
                        help='The path to guppy exexutable')
    parser.add_argument('-o', '--outpath', type=str, default=None, required=True,
                        help='The path for the outfile')
    parser.add_argument("-t", "--threads", type=int, default=4,
                        help="The number of threads to use for demultiplexing", required=False)
    parser.add_argument("-b", "--gpu_buffers", type=int, default=16,
                        help="The number of gpu buffers to use for demultiplexing", required=False)
    parser.add_argument("-g", "--use_gpu", action='store_flase', default=True,
                        help="The use this flaf to turn off gpu usage", required=False)

    args = parser.parse_args()
    inpath = args.inpath
    guppy_path = args.guppy_path
    outpath = args.outpath
    threads = args.threads
    gpu_buffers = args.gpu_buffers
    use_gpu = args.use_gpu

    main(inpath, guppy_path, outpath, threads, gpu_buffers, use_gpu)
