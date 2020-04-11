import argparse
import pathlib
import vcfpy


__author__ = 'Colin Anthony'


class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass


def main(infile, outpath):
    # force absolute file paths
    infile = pathlib.Path(infile).absolute()
    outpath = pathlib.Path(outpath).absolute()
    cons_outfile = pathlib.Path(outpath, infile.stem + "_vcf_consensus.fasta")
    depth_qual_outfile = pathlib.Path(outpath, infile.stem + "_depth_qual.csv")
    headers = "position,seq_depth,quality_score\n"
    sample_name = infile.stem.replace("_bcftools", "")

    ref_seq = ""
    cons_seq = ""
    first = True
    csv_string = ""
    ref_name = ""
    for record in vcfpy.Reader.from_path(infile):

        if not record.is_snv():
            continue
        pos = record.POS

        if first:
            ref_name = record.CHROM
            cons_seq += "N" * (pos - 1)
            ref_seq += "N" * (pos - 1)
            first = False

        ref_base = record.REF
        ref_seq += ref_base
        alt_subs = record.ALT
        qual =record.QUAL
        inf = record.INFO
        depth = inf["DP"]

        csv_string += f"{pos},{depth},{qual}\n"

        if len(alt_subs) >= 1:
            alt_base = record.ALT[0].value
        else:
            alt_base = ref_base

        cons_seq += alt_base

    with cons_outfile.open("w") as fh:
        fh.write(f">{sample_name}\n{cons_seq}\n")

    with depth_qual_outfile.open('w') as fh:
        fh.write(headers)
        fh.write(csv_string)

    with open("ref.fasta", "w") as fh:
        fh.write(f">{ref_name}\n{ref_seq}\n")

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
