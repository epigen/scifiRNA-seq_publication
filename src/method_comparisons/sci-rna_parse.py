#! /usr/bin/env python

"""
Convert sci-RNA-seq files from BAM files as downloaded from SRA
(interleaved BAM files) to the format used at CeMM: a single end BAM file with
barcodes as tags.
"""

import time
import pysam
from argparse import ArgumentParser

print(f"# {time.asctime()} - Program start.")

parser = ArgumentParser()
parser.add_argument(dest="input_bam", type=str)
parser.add_argument(dest="output_bam", type=str)

args = parser.parse_args()

print(f"# {time.asctime()} - Parsed arguments.")

input_bam = pysam.AlignmentFile(args.input_bam, check_sq=False)
output_bam = pysam.AlignmentFile(args.output_bam, "wb", template=input_bam)

end = False
rec = dict()
for i, read in enumerate(input_bam):
    if i % 10000 == 0:
        print(f"# {time.asctime()} - Line {i}.")
    # make sure records are correctly interleaved
    if read.is_read2 and end:
        raise

    if read.is_read1:
        # the round2 barcode has already been demultiplexed, gets the BC tag
        rec['r1'] = read.get_tag("RG")
        # handle two round or three round experiment
        if "+" in rec['r1']:
            rec['r3'], rec['BC'] = rec['BC'].split("+")
        # read 1 first 8bp contains the UMI, this gets the RX tag
        rec['RX'] = read.seq[:8]
        # read 1 remaining sequence has the round1, this gets the r2 tag
        rec['r2'] = read.seq[8:]
        end = False

    if read.is_read2:
        output_read = pysam.AlignedSegment()
        # for the three round experiments, remove r1 barcode from read name
        if "," in read.qname:
            output_read.qname = ",".join(read.qname.split(",")[1:])
        else:
            output_read.qname = read.qname
        # use sequence and quality of read2 as the final read
        output_read.seq = read.seq
        output_read.qual = read.qual
        output_read.set_tags(sorted(rec.items()))
        output_read.is_paired = False
        output_bam.write(output_read)
        rec = dict()
output_bam.close()

print(f"# {time.asctime()} - Program end.")

"""
# Read scheme as downloaded
# # two round
1 77    CCTTGTCTCATCAGAATG      AAAAAEEEEEEEEEEEEE      RG:Z:TTGATTGGCG
1 141   AGCCTGAGCAACAGAGTGACTCTGTCGAAAGAAAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAAGGAAAGAAGGAAAAGAAAGAAAAGAAAGGAAGAAAAGGAAGAA

# # three round
ATTCAGAA,M00511:3:000000000-B2Y5T:1:1101:11056:2211 77  CTGTCTATTATGCCGAGT RG:Z:ATTCAGAAGG+AGGCGAGAGC
ATTCAGAA,M00511:3:000000000-B2Y5T:1:1101:11056:2211 141 ATACGGAGAAAATTAGCATGGTGAAGCATTTCATATTTTGCAGTCACGGGAAGGTCATTTTACTGTTTGCTGGCTAGCTCCGTCATCGTTTGAATCAAAACAAAATGGCTGTCCCCTTCTATTATACTTTRG:Z:ATTCAGAAGG+AGGCGAGAGC
"""

"""
# # two round mid-throughput
1   77  *   0   0   *   *   0   0   GCAGTTATTTAACGTACC  A6AAAEAEEEAEEEEEEE  RG:Z:AGTTGAATCC
1   141 *   0   0   *   *   0   0   GCCTATGTCTTGTTGAGTGAAAAGAAAATTTCCAGTGTTCAGTCCATTGTCC    /AAAAEAEEEA<EEAEEEAEE//AEAEEE/AE<AE/E<EAA//EEAAE<AEA    RG:Z:AGTTGAATCC
"""

"""
# # sci-plex species mixing experiment
R1: 8M10r1
R2: 52T
<T>RG: 10r2
1  77      *  0  0  *  *  0  0  ACCCGTCTCGAACGCCGG      AAAAAEEEEE66EEEEEE      RG:Z:ATTATGCAAG
1  141     *  0  0  *  *  0  0  ATTTTAAACTTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN    AAAAAEEEEEE#########################################      RG:Z:ATTATGCAAG

"""

"""
# For sci-RNA-seq v3 (Cao-Spielman 2018) Supplementary Table S11
sc_ligation_RT_1: /5Phos/CAGAGCNNNNNNNNTCCTACCAGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT (9)
sc_ligation_1: GCTCTGA CAATCAAGT/ideoxyU/ACGACGCTCTTCCGATCT ACTTGATTG T (9)
P5-1: AATGATACGGCGACCACCGAGATCTACAC CTCCATCGAG ACACTCTTTCCCTACACGACGCTCTTCCGATCT (10)
P7-1: CAAGCAGAAGACGGCATACGAGATccgaatccgaGTCTCGTGGGCTCGG (10)
"""
