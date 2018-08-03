suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(SummarizedExperiment)
  library(GenomicRanges) })
args <- commandArgs(trailingOnly = TRUE)

fileIn=args[1]
fileInPh=args[2]
fileOutPh = substr(fileInPh,1, nchar(fileInPh)-3)

# Load in two neede .vcf files vep.vcf and vep.phased.vcf
vcf_vep = readVcf(fileIn)
vcf_phased = readVcf(fileInPh)

geno(vcf_vep)@listData$PH_info = geno(vcf_phased)$GT

# Make the new info fields.
newInfo <- DataFrame(Number=1, Type="String",
                     Description="Phasing information",
                     row.names="PH_info")

# Add the new info fields
geno(header(vcf_vep)) <- rbind(geno(header(vcf_vep)), newInfo)

#Write VCF 
writeVcf(vcf_vep, fileOutPh)
