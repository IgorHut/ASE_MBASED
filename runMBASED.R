# Script reads annotated vcf and performs MBASED ASE analysis

suppressMessages(library(data.table))
suppressMessages(source("http://bioconductor.org/biocLite.R"))
suppressMessages(library("VariantAnnotation")) #load the package
suppressMessages(library(MBASED))

set.seed(988482)
args <- commandArgs(trailingOnly = TRUE)

noSim=1000
if(length(args)>2){
  noSim=as.integer(args[3])
}


fileInput=args[1]
fileOutput=args[2]

vcf = readVcf(fileInput)

# Filter SNP in exons

vcfSNP <- vcf[isSNV(vcf)]

exonGeneName <- function(x){
  s = x[grepl('exon',x)]
  strsplit(s[1], '|', fixed=TRUE)[[1]][5]
}


vcfSNP <- vcfSNP[sapply(info(vcfSNP)$CSQ, function(x) sum(grepl('exon',x)))>0]
geneList <- sapply(info(vcfSNP)$CSQ, exonGeneName)


##------------Phase SNVs-------------------


d = data.frame(rowData(v)['REF'], rowData(v)['ALT'], geno(v)$GT)
colnames(d) <- c('ph1', 'ph2', 'GT')

## sort columns ph1 and ph2 on GT=='0|1' or '1|0'
## pass d$ph1 as allele1 and d$ph2 as allele2 

##-----------MBASED--Prepare samples-------------------

mySNVs <- GRanges(seqnames=as.character(seqnames(vcfSNP)),
                  ranges=ranges(rowRanges(vcfSNP)),
                  aseID=geneList,
                  allele1=as.character(unlist(rowData(vcfSNP)[,'ALT'])),
                  allele2=as.character(rowData(vcfSNP)[,'REF']))



mySample <- SummarizedExperiment(assays=list(
  lociAllele1Counts=matrix(as.integer(info(vcfSNP)$AC), ncol=1, dimnames=list(names(mySNVs), 'mySample')),
  lociAllele2Counts=matrix(as.integer(info(vcfSNP)$AN-info(vcfSNP)$AC), ncol=1, dimnames=list(names(mySNVs), 'mySample'))), 
  rowRanges=mySNVs)

#----Run MBASED

ASEresults_1s <- runMBASED(
  ASESummarizedExperiment=mySample,
  isPhased=FALSE,
  numSim=noSim,
  BPPARAM = MulticoreParam())#SerialParam())

#------Write results------

df = data.frame(assays(ASEresults_1s)$majorAlleleFrequency, 
                assays(ASEresults_1s)$pValueASE, 
                assays(ASEresults_1s)$pValueHeterogeneity)


df <- cbind(newColName = rownames(df), df)

colnames(df)=c('Gene', 'majorAlleleFrequency','pValueASE','pValueHeterogeneity')
write.table(df,fileOutput, quote=FALSE, row.names = FALSE, sep=" ")
