# LTA
Two-part python script associated with an article "Intronic variants of the lymphotoxin alpha (LTA) gene, within the MHC complex, affect parasite susceptibility in a free-living mammal"

This is a prototype of wrapper scripts allowing usage of [PHASE](http://stephenslab.uchicago.edu/phase/download.html) to retrive haplotypes from partially phased vcf file. 

The script uses [Biopython](http://biopython.org) to handle sequences and [PyVCF](https://pypi.python.org/pypi/PyVCF) for VCF format.

## Usage
In the beginning the input file for PHASE is prepared from provided unphased VCF. Then one should run PHASE manually using generated file. Eventually variant sequences are produced by superimposing phased sets of variants onto reference sequence. At current stage only FASTA format is supported. All parameters are hardcoded - please edit headers of source files in order to use script s with your data.
