#!/usr/bin/env python
from Bio import SeqIO
import sys
import vcf
import subprocess
from collections import defaultdict
import os.path
import argparse
import pathlib
# This script was modified from the Artic-EBOV repo http://artic.network/ebov/ebov-it-setup.html


def collect_depths(bamfile):
    if not os.path.exists(bamfile):
        raise SystemExit("bamfile %s doesn't exist" % (bamfile,))

    print((bamfile, sys.stderr))

    p = subprocess.Popen(['samtools', 'depth', bamfile], stdout=subprocess.PIPE)
    out, err = p.communicate()
    depths = defaultdict(dict)
    for ln in out.decode('utf-8').split("\n"):
        if ln:
            contig, pos, depth = ln.split("\t")
            depths[contig][int(pos)] = int(depth)

    return depths


def report(r, status, allele, vcffile):
    idfile = os.path.basename(vcffile).split(".")[0]
    print("%s\t%s\tstatus\t%s" % (idfile, r.POS, status), sys.stderr)
    print("%s\t%s\tdepth\t%s" % (idfile, r.POS, r.INFO.get('TotalReads', ['n/a'])), sys.stderr)
    print("%s\t%s\tbasecalledfrac\t%s" % (idfile, r.POS, r.INFO.get('BaseCalledFraction', ['n/a'])), sys.stderr)
    print("%s\t%s\tsupportfrac\t%s" % (idfile, r.POS, r.INFO.get('SupportFraction', ['n/a'])), sys.stderr)
    print("%s\t%s\tallele\t%s" % (idfile, r.POS, allele), sys.stderr)
    print("%s\t%s\tref\t%s" % (idfile, r.POS, r.REF), sys.stderr)


def main(reference, vcffile, bamfile, depth_threshold, sample_name):

    masked_positions = []
    bamfile = pathlib.Path(bamfile).absolute()
    reference = pathlib.Path(reference).absolute()
    vcffile = pathlib.Path(vcffile).absolute()
    outfile = pathlib.Path(bamfile.parent, sample_name + "_consensus_artic.fasta")
    depths = collect_depths(bamfile)
    seq = list(SeqIO.parse(open(str(reference)), "fasta"))[0]
    cons = list(seq.seq)

    for n, c in enumerate(cons):
        try:
            depth = depths[seq.id][n+1]
        except KeyError:
            depth = 0

        if depth < depth_threshold:
            cons[n] = 'N'

    for mask in masked_positions:
        cons[mask-1] = 'N'

    sett = set()
    vcf_reader = vcf.Reader(open(str(vcffile), 'r'))
    for record in vcf_reader:
        if record.ALT[0] != '.':
            # variant call
            if record.POS in masked_positions:
                report(record, "masked_manual", "n", vcffile)
                continue

            # commented out: primers removed with bamclipper
            # if 'PRIMER' in record.INFO:
            #     report(record, "primer_binding_site", "n")
            #     cons[record.POS-1] = 'N'
            #     continue

            # support = float(record.INFO['SupportFraction'])
            total_reads = int(record.INFO['TotalReads'])
            qual = record.QUAL

            ref = record.REF
            alt = str(record.ALT[0])

            if len(alt) > len(ref):
                print(f"Skipping insertion at position: {record.POS}", sys.stderr)
                continue

            if qual >= 200 and total_reads >= depth_threshold:
                if len(ref) > len(alt):
                    print(f"N-masking confident deletion at {record.POS}", sys.stderr)
                    for n in range(len(ref)):
                        cons[record.POS-1+n] = '-'
                    continue

                report(record, "variant", alt, vcffile)
                sett.add(record.POS)
                if len(ref) > len(alt):
                    print(("deletion", sys.stderr))
                    continue

                if len(alt) > len(ref):
                    print(("insertion", sys.stderr))
                    continue

                cons[record.POS-1] = str(alt)
            elif len(ref) > len(alt):
                continue
            else:
                report(record, "low_qual_variant", "n", vcffile)
                cons[record.POS-1] = 'N'
                continue

    with open(outfile, 'w') as handle:
        handle.write(f">{sample_name}\n{''.join(cons)}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="make consensus from vcf and bam file",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-r", "--reference", type=str, default=argparse.SUPPRESS,
                        help="The reference genome and primer scheme to use", required=True)
    parser.add_argument("-v", "--vcf_file", type=str, default=argparse.SUPPRESS,
                        help="The path and name of the vcf file", required=True)
    parser.add_argument("-b", "--bam_file", default=argparse.SUPPRESS, type=str,
                        help="The path and name of the sorted, trimmed bam file", required=True)
    parser.add_argument("-n", "--sample_name", type=str, default=argparse.SUPPRESS,
                        help="The sample name", required=True)
    parser.add_argument("-d", "--depth_threshold", type=int, default=10,
                        help="The minimum coverage allowed to call variants as real", required=True)

    args = parser.parse_args()

    bam_file = args.bam_file
    reference = args.reference
    vcf_file = args.vcf_file
    sample_name = args.sample_name
    depth_threshold = args.depth_threshold

    main(reference, vcf_file, bam_file, sample_name, depth_threshold)
