

geno(vcf_vep)@listData$PH_info = geno(vcf_phased)$GT
vcf_vep
geno(vcf_vep)
geno(vcf_vep)$PH_info

writeVcf(vcf_vep, "vcf_with_ph_info.vcf")


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
