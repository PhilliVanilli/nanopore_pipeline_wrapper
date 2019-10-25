import argparse
from collections import defaultdict
from Bio import SeqIO


__author__ = 'Colin Anthony'


def fastq_to_dict(fastq):
    fastq_d = defaultdict(list)
    for record in SeqIO.parse(open(fastq), "fastq"):
        fastq_d[record.name] = [str(record.seq), record.letter_annotations["phred_quality"]]

    return fastq_d


def main(fastq):
    fastq_to_dict(fastq)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Trim alignments from an amplicon scheme.')
    parser.add_argument('-in', '--infile', help='Input fastq filename and path')


    args = parser.parse_args()
    infile = args.infile

    main(infile)
