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
# args = parser.parse_args([
#     '/home/arendeiro/sci-rna/data/external/GSE110823/SRR6750057.bam',
#     '/home/arendeiro/sci-rna/data/external/GSE110823/splitseq_300.annotated.bam'])
args = parser.parse_args()

print(f"# {time.asctime()} - Parsed arguments.")

input_bam = pysam.AlignmentFile(args.input_bam, check_sq=False)
output_bam = pysam.AlignmentFile(args.output_bam, "wb", template=input_bam)

# read[1 :10] -> UMI
# read[10:18] -> r3
# read[48:56] -> r2
# read[86:94] -> r1

end = False
rec = dict()
for i, read in enumerate(input_bam):
    if i % 100000 == 0:
        print(f"# {time.asctime()} - Line {i}.")
    # make sure records are correctly interleaved
    if read.is_read1 and end:
        raise

    if read.is_read1:
        output_read = pysam.AlignedSegment()
        output_read.qname = read.qname
        # use sequence and quality of read1 as the final read
        output_read.seq = read.seq
        output_read.qual = read.qual
        output_read.is_paired = False

    if read.is_read2:
        end = False

        # The following quality filtering is the same as in the yjzhang pipeline:
        # https://github.com/yjzhang/split-seq-pipeline/blob/a5fab4e220cdf1911ac40519539d1e44462c6cd7/split_seq_pipeline/tools.py#L37
        # Basically, drop reads with any base with quality under 20
        quality = True
        if any([x < 20 for x in read.query_qualities[:9]]):
            quality = False
        elif (
            any([x < 20 for x in read.query_qualities[10:18]]) or  # r3
            any([x < 20 for x in read.query_qualities[48:56]]) or  # r2
            any([x < 20 for x in read.query_qualities[86:94]])     # r1
        ):
            quality = False
        if not quality:
            rec = dict()
            continue

        # the round2 barcode has already been demultiplexed, gets the BC tag
        rec['r1'] = read.get_tag("RG")
        # handle two round or three round experiment
        if "+" in rec['r1']:
            rec['r3'], rec['BC'] = rec['BC'].split("+")
        # read 2 first 10bp contains the UMI, this gets the RX tag
        # since last base is often N, I skip it
        rec['RX'] = read.seq[:9]
        # read 2, base-pair 10 to 18 has the round3, this gets the r3 tag
        rec['r3'] = read.seq[10:18]
        # read 2, base-pair 48 to 56 has the round2, this gets the r2 tag
        # since last two bases are N's, I skip them
        rec['r2'] = read.seq[48:56]
        # read 2, base-pair 86 to 92 has the round1, this gets the r1 tag
        # since last two bases are N's, I skip them
        rec['r1'] = read.seq[86:94]
        output_read.set_tags(sorted(rec.items()))
        output_bam.write(output_read)
        rec = dict()
output_bam.close()

print(f"# {time.asctime()} - Program end.")

"""
# read[1:10]->UMI, read[10:18]->r3, read[48:56]->r2, read[86:94]->r1
TTGCTTTCCN CATACCAA GTGGCCGATGTTTCGANTCGGCGTACNACT TGGCTTNN ANAANCGNTCTTANAAAGCCAGAGCATNAT NAACNNNN
CTCATACATN TAGGATGA GTGGCCGATGTTTCGCNTCGGCGTACNACT ATTGAGNN ANCCNCGNGCTTGNGAGGCCAGAGCATNCG NAGGNNNN
CGTCTAGACN CCAGTTCA GTGGCCGATGTTTCGCNTCGGCGTACNACT GTCGTANN ANCCNCGNGCTTGNGAGGCCAGAGCATNCG CCTANNCN
CGCCTTGTTN GACTAGTA GTGGCCGATGTTTCGCNTCGGCGTACNACT ACACGANN ANCCNCGNGCTTGNGAGGACAGAGCATNCG NGTANNNN
TAGCATATAN TCCGTCTA GTGGCCGATGTTTCGCNTCGGCGTACNACT CGGATTNN ANCCNCGNGCTTGNGAGGCCAGAGCATNCG TCTTNNCN
TTGGGATGAN GATAGACA GTGGCCGATGTTTCGCNTCGGCGTACNACT CGACACNN ANCCNCGNGCTTGNGAGGCCAGAGCATNCG TCTTNNCN
AGGTCAGATN CATAGCGA CGTGGCAGATGTTTCGNATCGGCGTANGAC TAGATCNN ANTCNACNAGCTTNAGAGGCCAGAGCANTC GAGTNNAN
GTATCGTCCN GCCTCATA GTGGCCGATGTTTCGCNTCGGCGTACNACT AAACATNN ANCANCGNGCTTGNGAGGCCAGAGCATNCG CATCNNGN
ACTCGTCCGN GTGTTCTA GTGGCCGATGTTTCGCNTCGGCGTACNACT TGAAGANN ANCCNCGNGCTTGNGAGGCCAGAGCATNCG ATTGNNGN
CCTTTACTCN GCTCGGTA GTGGCCGATGTTTCGCNTCGGTACGANTCC ACCACANN CNACNTGNTTGAGNGGCCAGAGCATTCNAA TCCGNNTN
TCATTCTACN AGTACAAG GTGGCCGATGTTTCGCNTCGGCGTACNACT GTGTTCNN ANCCNCGNGCTTGNGAGGCCAGAGCATNCG GCTGNNTN
TATATATTAN ACAGCAGA GTGGCCGATGTTTCGCNTCGGCGTACNACT AGTGGTNN ANCCNCGNGCTTGNGAGGCTTTAAAATNCA GCAGNNTN
TAATTTATAN AGTGGTCA GTGGCCGATGTGTCGCNTCGGCGTACNACT AGTCACNN ANCCNCGNGCTTGNGATGCCAGAGCATNCG ATAGNNAN
CTGCGTGTAN CGAATACA GTGGCCGATGTTTCGCNTCGGCGTACNACT TTCACGNN ANCCNCGNGCTTGNGAGGCCAGAGCATNCG ACACNNCN
CTTAGCTTTN CAAGGAGC GTGGCCGATGTTTCGCNTCGGCGTACNACT AGCCATNN ANACNCGNGCTTGNGAGGCCAGAGCATNAG TGAANNGN
TTTCTTTTGN AACTCACC GTGGCCGATGTTTCGCNTCGGCGTACNACT CATACCNN ANCCNCGNGCTTGNGAGGCCAGAGCATNCG CAAGNNGN
CTTTTAAACN AGATGTAC GTGGCCGATGTTTCGCNTCGGCGTACGACT GAGTTANN ANACNCGNTCTTGNGAGGCCAGAGCATNCG CATANNAN
TGCTCCCCTN CTGGCATA GTGGCAGATGTTTCGCNTCGGCGTACNACT GTCGTANN ANCCNCGNGCTTGNGAGGCCAGAGCATNCG TCTTNNCN
TTTATGCTAN ACACAGAA GTGGCCGATGTTTCGCNTCGGCGTACGACT TGAAGANN ANAANCCNTCTTGNACGGCCAGAGCATNAG ACCANNTN
GACGGGTAGN ACAGCAGA GTGGCCGATGTTTCGCNTCGGCGTACGACT GCGAGTNN ANCCNCGNGCTTGNGAGGCCAGAGCATNCG ACAGNNGN
TTCTAACGGN TCCGTCTA GGGGGCGAGGGGTCGCNTCGGCGTACGAAT CTCAATNN ANCCNCGNGCTTGNGAGGCCAGAGCATNCG AACTNNCN
ATGATTTGAN TCAGACTA GTGGCCGATGTTTCGCNTCGGCGTACGACT AACAACNN ANCCNCGNGCTTGNGAGGCCAGAGCATNCG AATGNNTN
ACGAACTTAN CGGATTGC GTGGCCGATGTGTCGCNTCGGCGTACGACT CAAGGANN ANCCNCGNGCTTGNGAGGCCAGAGCATNCG CTAANNTN
TCTCGCTTTN TTCACGCA GTGGCCGCTGTTTCGCNTCGGCGTACGACT ACGTATNN ANCCNCGNGCTTGNGAGGCCAGAGCATNCG AACGNNAN
TGTACCTTTN TGGTGGTA GTGGCCGATGTTTCGCNTCGGCGTACGACT CTCAATNN ANACNCGNGATTGNGAGGCCAGAGCATNAG CGAANNTN
AATTCGGGTN GCTAACGA GTGGCCGATGTTTCGCNTCGGCGTACGACT ACACAGNN ANCCNCGNGCTTGNGAGGCCAGAGCATNCG AACTNNCN
AACAGTAGTN CTGAGCCA GTGGCCGATGTTTCGCNTCGGCGTACGACT TCCGTCNN ANACNCGNGCTTGNGAGGCCAGAGCATNAG AGATNNCN
TATGTTCAGN CGGATTGC GTGGCCGATGTTTCGCNTCGGCGTGACTCC GACAACNN CNACNTGNTTGAGNGGCCAGAGCATTCNAT CATGNNTN
GCAATTTTCN CACCGAGA GTGGCCGATGGGGCGCNTCGGCGTACGACT TGGCACNN ANCCNCGNGCTTGNGAGGCCAGAGCATNCG CGCCNNTN
TTTACTCATN CGAACTTA GTGGCCGATGTTTCGCNTCGGCGTACGACT CCTCCANN ANACNCGNGCTTGNGAGGCCAGAGCATNCG CATCNNGN
GTCGGTTAAN CATACCAA GTGGCCGATGTTTCGCNTCGGCGTACGACT GATAGANN ANACNCGNTCTTGNGGCCAGAGCATTCNTC TCCTNNTN
TTCTTGTTTN CACCTTAC GTGGCCGATGTTTCGCNTCGGCGTACGACT ATCCTGNN ANCCNCGNGCTTGNGAGGCCAGAGCATNCG TGAANNGN
AGGTGAAGTN AACCGAGA GTGGCCGATGTTTCGCNTCGGCGTACGACT TGGAACNN ANCCNCGNGCTTGNGAGGCCAGAGCATNCG ATTGNNGN
ATACATTCGN TGGAACAA GTGGCCGATGTTTCGCNTCGGCGTACGACT CTCAATNN ANCCNCGNGCTTGNGAGGCCAGAGCATNCG GATGNNTN
TGAAATAGAN ATCCTGTA GTGGCCGATGTTTCGCNTCGGCGTACGACT AAACATNN ANCCNCGNGCTGTNTCTTATACACATCNGA CGCTNNCN
GCATGTGGGG ACAGATTC GTGGCCGATGTTTCGCNTCGGCGTACGACT GAATCTNN ANACNCGTGCTTGNGAGGCCAGAGCATNAG CAAGNNTN
GTTACCGTTT CCTCTATC GTGGCCGATGGGTCGCNTCGGCGTACGACT CCTAATNN TNTCNCTTATACANATCTGACGCTGCCNAC GAGCNNTN
AATACCTATN ATTGAGGA GTGGCCGATGTTTCGCNTCGGCGTACGACT AAGGACNN ANCCNCGNGCTTGNGAGGCCAGAGCATNCG ACAGNNTN
ATCTATGAGG TTCACGCA GTGGCCGATGTTTCGCNTCGGCGTACGACT TGGTGGNN ANCANGTGCTTGANAGGCCAGAGCATTNGA ACGCNNAN
CGTAGTAGTT GCTAACGA GTGGCCGATGTTTCGCNTCGGCGTACGACT CAGCGTNN ANCCNCGTGCTTGNGAGGCCAGAGCATNCG CCGTNNGN
"""
