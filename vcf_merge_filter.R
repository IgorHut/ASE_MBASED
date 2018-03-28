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

#Packaging of working df
head(as.data.frame(rowRanges(vcf_vep)$REF)$x)
head(as.data.frame(rowRanges(vcf_vep)$ALT)$value)
head(unlist(info(vcf_vep)$AC))
head(info(vcf_vep)$AN)
head(as.vector(unname(geno(vcf_phased)$GT)))

mamatata_df = data.frame(Ph1 = as.data.frame(rowRanges(vcf_vep)$REF)$x, 
                         Ph2 =  as.data.frame(rowRanges(vcf_vep)$ALT)$value,
                         Ph1_count = unlist(info(vcf_vep)$AC),
                         Ph2_count = info(vcf_vep)$AN - unlist(info(vcf_vep)$AC),
                         phasing_info = as.vector(unname(geno(vcf_phased)$GT)))
rownames(mamatata_df) = names(rowRanges(vcf_phased))
head(mamatata_df)

#Swap Ph1, Ph2, Ph1_count, Ph2_count for phasing_info = 1|0
mamatata_df[mamatata_df$phasing_info == "1|0", c("Ph1", "Ph2", "Ph1_count", "Ph2_count")] <- 
  mamatata_df[mamatata_df$phasing_info == "1|0", c("Ph2", "Ph1", "Ph2_count", "Ph1_count")]

head(mamatata_df)
