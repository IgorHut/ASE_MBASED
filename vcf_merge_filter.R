library(VariantAnnotation)
library(SummarizedExperiment)
library(GenomicRanges)

# Load in two neede .vcf files vep.vcf and vep.phased.vcf

vcf_vep = readVcf("~/Documents/GitHub/ASE_MBASED/G20479.vep.vcf")
vcf_phased =  readVcf("~/Documents/GitHub/ASE_MBASED/G20479.vep.phased.vcf")

vcf_vep<- vcf_vep[isSNV(vcf_vep)]
dim(vcf_vep)
vcf_phased <- vcf_phased[isSNV(vcf_phased)]
dim(vcf_phased)

vcf_vep
vcf_phased
info(vcf_vep)
rowRanges(vcf_vep)
rowRanges(vcf_phased)


head(geno(vcf_vep)$GT)

head(geno(vcf_phased)$GT)

sum(rownames(geno(vcf_vep)$GT) != rownames(geno(vcf_phased)$GT))
