# MBASED Allele-Specific Expression analysis report
### The script reads an annotated vcf file and performs MBASED ASE analysis
suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(MBASED)
  library(BiocParallel)
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)

set.seed(988482)


fileInput=args[1]    #'G20479.vep.vcf'
fileOutput=args[2]

noSim=as.integer(args[4])

####-----Create empty outputs if VCF empty, corrupted, or no SNPs in genes 
dfG = data.frame(gene_name = as.character(), 
                 majorAlleleFrequency = as.character(), 
                 pValueASE = as.character(), 
                 pValueHeterogeneity = as.character())

dfM = data.frame(gene_name = as.character(), 
                 seqnames = as.character(), 
                 start = as.character(), 
                 mutation_id = as.character(), 
                 RefAlt = as.character(), 
                 count = as.character(), 
                 VAF = as.character(), 
                 phasing_info = as.character(), 
                 geneConsequence = as.character()) 


vcf = tryCatch(readVcf(fileInput), error =  function(cond) { message("Empty or corrupted VCF file"); return(0)})

if(is.integer(vcf)){
  write.table(dfG, file=paste0(fileOutput,".mbased_genes.tsv"), sep="\t", row.names = F, quote = FALSE)
  write.table(dfM, file=paste0(fileOutput,".mbased_mut.tsv"),   sep="\t", row.names = F, quote = FALSE)
  message("Corrupted VCF file")
  quit()
}

if(sum(dim(vcf))==0){
  write.table(dfG, file=paste0(fileOutput,".mbased_genes.tsv"), sep="\t", row.names = F, quote = FALSE)
  write.table(dfM, file=paste0(fileOutput,".mbased_mut.tsv"),   sep="\t", row.names = F, quote = FALSE)
  message("Empty VCF file")
  quit()
} 

#------Filter SNPs 

vcf<- vcf[isSNV(vcf)]
vcf <- vcf[lapply(info(vcf)$CSQ,length)>0] # info from VEP 

if(length(vcf)==0){
  message ("There are no SNPs within gene regions. Therefore, no further analysis could be employed.")
  write.table(dfG, file=paste0(fileOutput,".mbased_genes.tsv"), sep="\t", row.names = F, quote = FALSE)
  write.table(dfM, file=paste0(fileOutput,".mbased_mut.tsv"),   sep="\t", row.names = F, quote = FALSE)
  quit()
}

#----Check if VCF is phased and contains 0|1 and 1|0 

phased = 'PH_info' %in% names(geno(vcf))
phases <- rep(0, dim(vcf)[1])
if(phased){
  phases <- as.vector(unname(geno(vcf)$PH_info))
  if(sum(c("0|1", "1|0") %in% levels(as.factor(phases)))>0){ # if any "0|1", "1|0"
    message ("VCF Phased")
    vcf <- vcf[as.vector(geno(vcf)$PH_info!='1|1')]
    phases <- as.vector(unname(geno(vcf)$PH_info))
  }else{
    phased = FALSE
  }
}else{
  message ("VCF not phased, MBASED phasing performed")
}

# Find position of GENE names in vep annotated vcf

vepD <- info(header(vcf))[c("CSQ"),]$Description
vepD <- strsplit(vepD,'Format:', fixed=TRUE)[[1]][2]
vepD <- trimws(vepD, ("b"))
vepD <- strsplit(vepD,'|', fixed=TRUE)[[1]]
pos  <- match("SYMBOL", vepD)
posID  <- match("Gene", vepD)
posCon  <- match("Consequence", vepD)

mutInGene <- function(x, p){ strsplit(x,'|', fixed=TRUE)[[1]][p] }

getGeneName <- function(x, p){
  l = as.character(lapply(x, mutInGene, p))
  paste(unique(l), sep='',collapse = ';')
}

geneList <- sapply(info(vcf)$CSQ, getGeneName, pos)
geneConsequence <- sapply(info(vcf)$CSQ, getGeneName, posCon)

##------------Phase SNVs-------------------

mamatata = data.frame(Ph1 = as.data.frame(rowRanges(vcf)$REF)$x, 
                      Ph2 = as.data.frame(rowRanges(vcf)$ALT)$value,
                      Ph1_count = sapply(geno(vcf)$AD, "[[", 1),
                      Ph2_count = sapply(geno(vcf)$AD, "[[", 2),
                      phasing_info = phases, 
                      ranges = ranges(rowRanges(vcf)), 
                      seqnames = as.character(seqnames(vcf)),
                      gene_name = geneList,
                      geneConsequence = geneConsequence,
                      stringsAsFactors = FALSE)


dz = dim(mamatata[(mamatata$Ph1_count + mamatata$Ph2_count)==0,])[1]
if(dz>0){
  mamatata = mamatata[mamatata$Ph1_count + mamatata$Ph2_count!=0,] # Total number of variants has to be bigger than 0 
  if(dim(mamatata)[1]==0){
    message ("Allelic depths for the ref and alt alleles missing")
    quit()
  }else{
    message ("Allelic depths for ", dz, " ref and alt alleles are missing")
    
  }
}

#Swap Ph1, Ph2, Ph1_count, Ph2_count for phasing_info = 1|0
if(phased){
  mamatata[mamatata$phasing_info == "1|0", c("Ph1", "Ph2", "Ph1_count", "Ph2_count")] <- mamatata[mamatata$phasing_info == "1|0", c("Ph2", "Ph1", "Ph2_count", "Ph1_count")]
}
### Same mutation in two or more genes (Conjoined genes): 
mamatata = mamatata %>% separate_rows(gene_name, sep=';')
mamatata = mamatata[order(mamatata$seqnames, mamatata$gene_name),]

##-----------MBASED--Prepare samples-------------------

mySNVs <- GRanges(seqnames = mamatata$seqnames,
                  ranges = IRanges(start = mamatata$ranges.start, 
                                   end = mamatata$ranges.end, 
                                   width = mamatata$ranges.width, 
                                   names = mamatata$ranges.names),
                  aseID   = as.character(mamatata$gene_name),
                  allele1 = as.character(mamatata$Ph1),
                  allele2 = as.character(mamatata$Ph2))

mySample <- SummarizedExperiment(assays = list(
  lociAllele1Counts = matrix(mamatata$Ph1_count, ncol=1, dimnames=list(names(mySNVs), 'mySample')),
  lociAllele2Counts = matrix(mamatata$Ph2_count, ncol=1, dimnames=list(names(mySNVs), 'mySample'))),
  #lociRhos=matrix(rep(0.9,length(mamatata$Ph1)), dimnames=list(names(mySNVs), 'mySample'))),
  rowRanges = mySNVs)

#----Run MBASED----PHASED--------

ASEres <- runMBASED(ASESummarizedExperiment=mySample,
                    isPhased=phased,
                    numSim=noSim,
                    BPPARAM = SerialParam()) #MulticoreParam or SerialParam

#------Write results------

dfG = data.frame(gene_name = rownames(assays(ASEres)$majorAlleleFrequency),
                 majorAlleleFrequency = round(as.numeric(assays(ASEres)$majorAlleleFrequency),2),
                 pValueASE = round(as.numeric(assays(ASEres)$pValueASE),3),
                 pValueHeterogeneity = round(as.numeric(assays(ASEres)$pValueHeterogeneity),2))

dfM = data.frame(mutation_id = rownames(assays(metadata(ASEres)$locusSpecificResults)$MAF),
                 maf = as.numeric(assays(metadata(ASEres)$locusSpecificResults)$MAF),
                 allele_is_major = as.logical(assays(metadata(ASEres)$locusSpecificResults)$allele1IsMajor))

# Swap back Ph1, Ph2
if(phased){
  mamatata[mamatata$phasing_info == "1|0", c("Ph1", "Ph2", "Ph1_count", "Ph2_count")] <- mamatata[mamatata$phasing_info == "1|0", c("Ph2", "Ph1", "Ph2_count", "Ph1_count")]
}

mamatata$count = paste(mamatata$Ph1_count, mamatata$Ph2_count, sep=';')
mamatata$RefAlt = paste(mamatata$Ph1, mamatata$Ph2, sep='')

mamatata = mamatata %>% select(-one_of(c('ranges.names','ranges.end','ranges.width','Ph1_count','Ph2_count','Ph1','Ph2'))) # No need for these columns
mamatata = mamatata %>% rename(start = ranges.start)

dfM = cbind(mamatata, dfM )

dfM[dfM$phasing_info == "1|0","allele_is_major"]  = !dfM[dfM$phasing_info == "1|0","allele_is_major"] # Swap allele_is_major
dfM[dfM$allele_is_major == TRUE,"maf"]  = 1 - dfM[dfM$allele_is_major == TRUE,"maf"] 
dfM = dfM %>% rename(VAF = maf)
dfM$VAF = round(dfM$VAF,2)

cols = c("seqnames", "start", "mutation_id", "RefAlt", "count", "VAF", "phasing_info", "geneConsequence")
dfM = dfM[,c("gene_name" ,cols)]

dfM = dfM %>% group_by_at(vars(cols))  %>% summarise(gene_name = paste(gene_name, collapse=";"))
dfM = dfM[order(dfM$seqnames,dfM$start),c('gene_name',cols)]

cat(paste("One-sample analysis with counts over", nrow(dfG), "transcripts\n"))

write.table(dfG, file=paste0(fileOutput,".mbased_genes.tsv"), sep="\t", row.names = F, quote = FALSE)
write.table(dfM, file=paste0(fileOutput,".mbased_mut.tsv"),   sep="\t", row.names = F, quote = FALSE)


