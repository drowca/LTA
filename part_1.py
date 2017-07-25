## TODO:
## this part should be passed by options
chrom = 'TNFb_full'
reference_file = "data/TNFb_ref.fna"
vcffile = 'data/TNF.second.vcf'
phase_pfx = 'data/' + chrom
phase_in_fn = phase_pfx + '.inp'
phase_known_fn = phase_pfx + '.known'

locilist = dict()
genotypes = dict()
ids = set()

import vcf
vcf_reader=vcf.Reader(open(vcffile, 'r'))
for record in vcf_reader:
    if record.CHROM != chrom:
        continue
    locilist[record.POS] = 'S' if record.is_snp else 'M'
    for sample in record.samples:
        ids.add(sample.sample)
        x = sample.data.GT
        if not sample.sample in genotypes.keys():
            genotypes[sample.sample] = dict()
        if x != '.':
            genotypes[sample.sample][record.POS] = [int(a) + 1 for a in x.split('/')]
        else:
            genotypes[sample.sample][record.POS] = ['?', '?'] if record.is_snp else ['-1', '-1']

if not len(locilist):
    exit()

phase_in = open(phase_in_fn,  'w')
#phase_known = open(phase_known_fn,  'w')
print(len(ids),  file=phase_in)
print(len(locilist),  file=phase_in)
print("P %s" % " ".join([str(x) for x in sorted(locilist.keys())] ),  file=phase_in)
print(" ".join(locilist[y] for y in sorted(locilist.keys()) ),  file=phase_in)
for id in ids:
    print(id,  file=phase_in)
    for j in [0, 1]:
        print(" ".join(str(genotypes[id][i][j]) for i in sorted(locilist.keys())),  file=phase_in)
phase_in.close()
