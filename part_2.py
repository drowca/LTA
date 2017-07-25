## TODO:
## this part should be passed by options
chrom = 'TNFb_full'
reference_file = "data/TNFb_ref.fas"
vcffile = 'data/TNF.second.vcf'
phase_pfx = 'data/TNFb_res'
phase_in_fn = phase_pfx + '.out'
frame_positions_fn = 'data/start_stop.txt'
pos_start = 0
pos_stop = 0

refseq=''
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
for seq_record in SeqIO.parse(reference_file, "fasta",  IUPAC.ambiguous_dna):
    if seq_record.id != chrom:
        continue
    refseq = seq_record

varlist = dict()
import vcf
vcf_reader=vcf.Reader(open(vcffile, 'r'))
for record in vcf_reader:
    if record.CHROM != chrom:
        continue
    varlist[record.POS] = {'alleles': record.alleles,\
                            'start': record.affected_start,\
                            'end': record.affected_end}

def read_section(filename, section):
    """Read section of file where divisions are made with lines
    BEGIN Section_Name; END Section_Name. Specifically for PHASE
    output file"""
    with open(filename) as file:
        recording = False
        for line in file:
            if "BEGIN " + section in line:
                recording = True
            elif "END " + section in line:
                recording = False
            elif recording:
                yield line

poslist = []
for line in open(phase_in_fn,  "r"):
    if line.startswith('Positions of loci'):
        tmp = line.partition(':')
        poslist = [int(x) for x in tmp[2].split()]

# TODO: this part should be made into try-except clause
if pos_start == 0 and pos_stop == 0:
    pos_stop = None

alleles = dict()
for line in read_section(phase_in_fn, 'LIST_SUMMARY'):
    a = line.split()
    allele_seq = refseq
    # WARNING: this won't work with multidigit allele symbols in PHASE output!!!
    varmap = dict(zip(poslist, ''.join(str(x) for x in a[1:-1])))
    for pos in reversed(sorted(poslist)):
        substitute = Seq(str(varlist[pos]['alleles'][int(varmap[pos])-1]), IUPAC.ambiguous_dna)
        startN = pos - 1 if varlist[pos]['start'] == pos else varlist[pos]['start']
        endN = varlist[pos]['end']
        allele_seq = allele_seq[:startN] + substitute + allele_seq[endN:]
    # TODO check pos_start < pos_stop in case of RC
    alleles[a[0]] = allele_seq[pos_start:pos_stop]

from Bio.SeqRecord import SeqRecord
genotypes_writeout = []
for line in read_section(phase_in_fn, 'BESTPAIRS_SUMMARY'):
    aid,  tmp,  haps = line.partition(':')
    hap1,  tmp,  hap2 = haps.lstrip(' (').rstrip(')\n').partition(',')
    genotypes_writeout.append(SeqRecord(alleles[hap1].seq,  id=aid + '-1',  description=chrom))
    genotypes_writeout.append(SeqRecord(alleles[hap2].seq,  id=aid + '-2',  description=chrom))

SeqIO.write(genotypes_writeout, chrom + ".fas", "fasta")
