from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

genInstall.install('GRCh37', rsync=False, bash=True)

matrices = matGen.SigProfilerMatrixGeneratorFunc("test", "GRCh37", "datasets/unzipped/21BRCA/21BRCA/21BRCA_vcf",plot=True, exome=False, bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)

