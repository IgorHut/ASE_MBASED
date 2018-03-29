library(VariantAnnotation)
library(SummarizedExperiment)
library(GenomicRanges)

# Load in two neede .vcf files vep.vcf and vep.phased.vcf

vcf_vep = readVcf("~/Documents/GitHub/ASE_MBASED/G20479.vep.vcf")
vcf_phased =  readVcf("~/Documents/GitHub/ASE_MBASED/G20479.vep.phased.gz")

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

vcf_phased <- vcf_phased[as.vector(geno(vcf_phased)$GT!='1|1')]
vcf_vep<-vcf_vep[rownames(vcf_phased)]

dim(vcf_phased)
dim(vcf_vep)


geno(vcf_vep)@listData$PH_info = geno(vcf_phased)$GT
vcf_vep
geno(vcf_vep)
geno(vcf_vep)$PH_info

## Current 'info' fields.
rownames(info(header(vcf_vep)))

newInfo <- DataFrame(Number=1, Type="String",
                     Description="Phasing information",
                     row.names="PH_info")
info(header(vcf_vep)) <- rbind(info(header(vcf_vep)), newInfo)
geno(header(vcf_vep)) <- rbind(geno(header(vcf_vep)), newInfo)

rownames(info(header(vcf_vep)))


writeVcf(vcf_vep, "vcf_with_ph_info.vcf")
vcf_prob = readVcf("vcf_with_ph_info.vcf")
vcf_prob
geno(vcf_prob)$PH_info

header(vcf_vep)@header@listData$FORMAT@rownames <- 
  c(header(vcf_vep)@header@listData$FORMAT@rownames, "PH_info")

header(vcf_vep)@header@listData$FORMAT@nrows <- integer(6)

header(vcf_vep)@header@listData$FORMAT@listData$Number <- 
  c(header(vcf_vep)@header@listData$FORMAT@listData$Number, 1)

header(vcf_vep)@header@listData$FORMAT@listData$Type <- 
  c(header(vcf_vep)@header@listData$FORMAT@listData$Type, "String")

header(vcf_vep)@header@listData$FORMAT@listData$Description <- 
  c(header(vcf_vep)@header@listData$FORMAT@listData$Description, "Phasing information")

geno(header(vcf_vep))