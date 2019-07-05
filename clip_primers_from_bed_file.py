import argparse
import sys
from copy import copy
import csv
from operator import itemgetter
import pysam


"""
Written by Nick Loman as part of the ZiBRA pipeline (zibraproject.org)
edited by Colin Anthony
"""


def check_still_matching_bases(s):
    for flag, length in s.cigartuples:
        if flag == 0:
            return True
    return False


def read_bed_file(primerset_bed):

    out_bedfile = []
    with open(primerset_bed) as csvfile:
        reader = csv.DictReader(csvfile, dialect='excel-tab')
        for row in reader:
            primer_name = row['Primer_ID']
            orientation = primer_name.split("_")
            if "LEFT" in orientation:
                row['direction'] = "+"
            elif "RIGHT" in orientation:
                row['direction'] = "-"
            else:
                print("could not parse primer name\nexpected either LEFT or RIGHT in name, using '_' as a separator")
                raise ValueError

            row['end'] = int(row['end'])
            row['start'] = int(row['start'])
            out_bedfile.append(row)

    return out_bedfile


def trim(cigar, s, start_pos, end):
    if not end:
        pos = s.pos
    else:
        pos = s.reference_end

    eaten = 0
    while 1:
        # chomp stuff off until we reach pos
        if end:
            flag, length = cigar.pop()
        else:
            flag, length = cigar.pop(0)

        if flag == 0:
            # match
            eaten += length
            if not end:
                pos += length
            else:
                pos -= length
        if flag == 1:
            # insertion to the ref
            eaten += length
        if flag == 2:
            # deletion to the ref
            if not end:
                pos += length
            else:
                pos -= length
            pass
        if flag == 4:
            eaten += length
        if not end and pos >= start_pos and flag == 0:
            break
        if end and pos <= start_pos and flag == 0:
            break

    extra = abs(pos - start_pos)

    if extra:
        if flag == 0:
            if end:
                cigar.append((0, extra))
            else:
                cigar.insert(0, (0, extra))
            eaten -= extra

    if not end:
        s.pos = pos - extra

    if end:
        cigar.append((4, eaten))
    else:
        cigar.insert(0, (4, eaten))
    # oldcigarstring = s.cigarstring
    s.cigartuples = cigar
    # print(s.query_name, oldcigarstring[0:50], s.cigarstring[0:50])


def find_primer(bed, ref_pos, direction):
    closest = min([(abs(p['start'] - ref_pos), p['start'] - ref_pos, p) for p in bed if p['direction'] == direction],
                  key=itemgetter(0))
    return closest


def is_correctly_paired(p1, p2):
    name1 = p1[2]['Primer_ID']
    name2 = p2[2]['Primer_ID']

    name1 = name1.replace('_LEFT', '')
    name2 = name2.replace('_RIGHT', '')

    return name1 == name2


def main(infile, outfile, bedfile):

    bed = read_bed_file(bedfile)
    infile = pysam.AlignmentFile(infile, "rb")
    outfile_trimmed = pysam.AlignmentFile(outfile, "wh", template=infile)
    suppl_out = outfile + "_excluded_sequences_as_supplamentary.sam"
    marked_supplamentary = pysam.AlignmentFile(suppl_out, "wh", template=infile)
    primer_mismatch_file = outfile + "_excluded_as_primer_mismatched.sam"
    marked_primer_missmatch = pysam.AlignmentFile(primer_mismatch_file, "wh", template=infile)
    total = 0
    unmapped = 0
    missmatched = 0
    suppl = 0

    for s in infile:
        total += 1
        cigar = copy(s.cigartuples)

        # logic - if alignment start site is _before_ but within X bases of  a primer site, trim it off
        if s.is_unmapped:
            # print(f"{s.query_name} skipped as unmapped")
            unmapped += 1
            continue

        if s.is_supplementary:
            # print(f"{s.query_name} skipped as supplementary")
            marked_supplamentary.write(s)
            suppl += 1
            continue

        p1 = find_primer(bed, s.reference_start, '+')
        p2 = find_primer(bed, s.reference_end, '-')

        if not is_correctly_paired(p1, p2):
            # print("mismatched primer pair. primers matched:", p1[2]['Primer_ID'], p2[2]['Primer_ID'],
            #       "this is probably two amplicons ligated together")
            marked_primer_missmatch.write(s)
            missmatched += 1
            continue

        # if the alignment starts before the end of the primer, trim to that position
        try:
            primer_position = p1[2]['end']
            # print(f"end {primer_position}\nref start {s.reference_start}")
            if s.reference_start < primer_position:
                trim(cigar, s, primer_position, 0)

            primer_position = p2[2]['start']
            # print(f"end {primer_position}\nref end {s.reference_end}")
            if s.reference_end > primer_position:
                trim(cigar, s, primer_position, 1)

        except Exception as e:
            print("problem %s" % (e,), sys.stderr)
            pass

        if not check_still_matching_bases(s):
            continue

        outfile_trimmed.write(s)

    print(f"Total: {total}\nUnmapped: {unmapped} ({round(unmapped/total*100, 2)}%)\nSupplamentary: {suppl} ({round(suppl/total*100, 2)}%)\n"
          f"Mismatched primers: {missmatched} ({round(missmatched/total*100, 2)}%)")
    print("Finished soft clipping bam file")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Trim alignments from an amplicon scheme.')
    parser.add_argument('-in', '--infile', help='Input BAM filename and path')
    parser.add_argument('-o', '--outfile', help='output SAM filename and path')
    parser.add_argument('-b', '--bedfile', help='BED file containing the amplicon scheme')

    args = parser.parse_args()
    infile = args.infile
    outfile = args.outfile
    bedfile = args.bedfile

    main(infile, outfile, bedfile)
