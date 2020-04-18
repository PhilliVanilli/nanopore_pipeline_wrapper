import argparse
import pathlib
import matplotlib.pyplot as plt
import pandas as pd
import re

__author__ = 'Colin Anthony'


class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass


def plot_depth(x_vals, depth_vals, sample_name, outfile):

    fig, ax = plt.subplots()
    ax.set_ylabel('Sequencing depth')
    ax.set_xlabel('Sequence position')
    ax.set_title(sample_name + " sequencing depth")

    # ax.set_xticks(x_labs)
    ax.tick_params(axis='x', rotation=90)

    plt.plot(x_vals, depth_vals)

    w = 8
    h = 4
    f = plt.gcf()
    f.set_size_inches(w, h)

    plt.savefig(outfile, ext="png", dpi=300, facecolor="white", bbox_inches="tight")


def plot_qual(x_vals, qual_vals, sample_name, outfile):

    fig, ax = plt.subplots()
    ax.set_ylabel('Sequencing quality')
    ax.set_xlabel('Sequence position')
    ax.set_title(sample_name + " Quality scores")

    last_used_index = 0
    ax.tick_params(axis='x', rotation=90)

    plt.plot(x_vals, qual_vals)
    w = 8
    h = 4
    f = plt.gcf()
    f.set_size_inches(w, h)

    plt.savefig(outfile, ext="png", dpi=300, facecolor="white", bbox_inches="tight")


def main(infile, outpath):
    # force absolute file paths
    infile = pathlib.Path(infile).absolute()
    outpath = pathlib.Path(outpath).absolute()
    sample_name = infile.stem
    sample_name = re.sub("_bcftools_.*", "", sample_name)
    depth_outfile = pathlib.Path(outpath, sample_name + "_vcf_depth_plot.png")
    qual_outfile = pathlib.Path(outpath, sample_name + "_vcf_qual_plot.png")
    data = pd.read_csv(infile, sep=',', header=0, parse_dates=True)

    # plot_depth(data['position'], data['seq_depth'], sample_name, depth_outfile)
    plot_qual(data['position'], data['quality_score'], sample_name, qual_outfile)

    print("done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='',
                                     formatter_class=Formatter)

    parser.add_argument('-in', '--infile', type=str, default=None, required=True,
                        help='The path and name of the infile')
    parser.add_argument('-o', '--outpath', type=str, default=None, required=True,
                        help='The path for the outfile')

    args = parser.parse_args()
    infile = args.infile
    outpath = args.outpath

    main(infile, outpath)
