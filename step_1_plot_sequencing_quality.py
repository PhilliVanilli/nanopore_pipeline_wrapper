import argparse
import pathlib
from pycoQC.pycoQC import pycoQC
from plotly.offline import plot


__author__ = 'Colin Anthony'


def main(project_path, seq_summary):
    """
    script to generate QC plots for nanopore run
    :param project_path: (str) path to the project folder (output files will be written here
    :param seq_summary: (str) path and name of the sequencing_summart.txt file
    :return:
    """
    project_path = pathlib.Path(project_path).absolute()
    seq_summary = str(pathlib.Path(seq_summary).absolute())

    plot_path = pathlib.Path(project_path, "QC_plots")
    plot_path.mkdir(mode=0o777, parents=True, exist_ok=True)

    p = pycoQC(seq_summary)

    fig_summary = p.summary()
    summary_fig = str(pathlib.Path(plot_path, "fig_summary.html"))
    plot(fig_summary, filename=summary_fig, show_link=False)

    fig_read_length = p.reads_len_1D()
    read_length_fig = str(pathlib.Path(plot_path, "fig_read_length.html"))
    plot(fig_read_length, filename=read_length_fig, show_link=False)

    fig_read_qual = p.reads_qual_1D()
    read_qual_fig = str(pathlib.Path(plot_path, "fig_read_qual.html"))
    plot(fig_read_qual, filename=read_qual_fig, show_link=False)

    fig_len_qual = p.reads_len_qual_2D()
    len_qual_fig = str(pathlib.Path(plot_path, "fig_len_qual.html"))
    plot(fig_len_qual, filename=len_qual_fig, show_link=False)

    fig_output_over_time = p.output_over_time()
    output_over_time_fig = str(pathlib.Path(plot_path, "fig_output_over_time.html"))
    plot(fig_output_over_time, filename=output_over_time_fig, show_link=False)

    fig_len_over_time = p.len_over_time()
    len_over_time_fig = str(pathlib.Path(plot_path, "fig_len_over_time.html"))
    plot(fig_len_over_time, filename=len_over_time_fig, show_link=False)

    fig_qual_over_time = p.qual_over_time()
    qual_over_time_fig = str(pathlib.Path(plot_path, "fig_qual_over_time.html"))
    plot(fig_qual_over_time, filename=qual_over_time_fig, show_link=False)

    fig_channels_activity = p.channels_activity()
    channels_activity_fig = str(pathlib.Path(plot_path, "fig_channels_activity.html"))
    plot(fig_channels_activity, filename=channels_activity_fig, show_link=False)

    print("Plotting QC stats complete")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process raw nanopore reads to consensus sequences",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-in", "--project_path", default=argparse.SUPPRESS, type=str,
                        help="The path to the directory containing the 'fast5' and 'fastq' folders. "
                             "Output interactive html plots will be written to a 'QC_plots' folder here", required=True)
    parser.add_argument("-s", "--seq_summary", default=argparse.SUPPRESS, type=str,
                        help="The path and name of the sequencing_summary.txt file", required=True)

    args = parser.parse_args()
    project_path = args.project_path
    seq_summary = args.seq_summary

    main(project_path, seq_summary)