# Script reads annotated vcf and performs MBASED ASE analysis
#' title: "MBASED Allele-Specific Expression analysis report"
#' output:
#'   html_document:
#'       toc: true
#' ---

####-------------Parse arguments

#+ options, echo = F, results = 'hide', warning = F, error = F, message = F
options(stringsAsFactors = FALSE, width = 120)

#+ collect arguments, echo = F, results = 'hide', warning = F, error = F, message = F
args <- commandArgs(TRUE)

# Parse arguments (the expected form is --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
ASEarg <- as.list(as.character(argsDF$V2))
names(ASEarg) <- argsDF$V1

#+ title, echo = FALSE
if (!is.null(ASEarg$title)) cat(ASEarg$title)

# Number of simulations
if(is.null(ASEarg$numSim)) {
  ASEarg$numSim <- 1000
}

## Number of CPUs default value is 2
if(is.null(ASEarg$ncpus2use)) {
  ASEarg$ncpus2use <- 2
}

# sigTreshold default (significance treshold)
if(is.null(ASEarg$CutOff)) {
  ASEarg$CutOff <- 0.1
}

####----------------MBASED-------------------

suppressPackageStartupMessages({
  library("knitr")
  library("gplots")
  library("VariantAnnotation")
  library("MBASED")
  library("BiocParallel")
})

set.seed(988482)



######Test


#ASEarg <- new.env()
#ASEarg$ncpus2use <- 2
#ASEarg$numSim <- 1000
#ASEarg$CutOff <- 0.1
#ASEarg$title <- "Allele_Specific_Expression"
#as.list(ASEarg)


###########



vcf = readVcf('G20479.vep.vcf')#(ASEarg$addVCF)
vcfP = readVcf('phased.vcf.gz')

# Filter SNP in exons

vcfSNP <- vcf[isSNV(vcf)]
vcfP <- vcfP[rownames(vcfSNP)]
vcfP <- vcfP[as.vector(geno(vcfP)$GT!='1|1')]
vcfSNP<-vcfSNP[rownames(vcfP)]

exonGeneName <- function(x){ s = x[grepl('exon',x)]; strsplit(s[1], '|', fixed=TRUE)[[1]][5]}


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

ASEres <- runMBASED(
  ASESummarizedExperiment=mySample,
  isPhased=FALSE,
  numSim=as.numeric(ASEarg$numSim),
  BPPARAM = MulticoreParam(workers=ASEarg$ncpus2use))#SerialParam()) workers=ASEarg$ncpus2use)

#------Write results------


df = data.frame(rownames(assays(ASEres)$majorAlleleFrequency),
                assays(ASEres)$majorAlleleFrequency,
                assays(ASEres)$pValueASE,
                p.adjust(assays(ASEres)$pValueASE, "fdr"),
                assays(ASEres)$pValueHeterogeneity,
                p.adjust(assays(ASEres)$pValueHeterogeneity, "fdr"))

df1 = data.frame(rownames(assays(metadata(ASEres)$locusSpecificResults)$MAF),
                 assays(metadata(ASEres)$locusSpecificResults)$MAF,
                 assays(metadata(ASEres)$locusSpecificResults)$allele1IsMajor)


colnames(df)=c('Gene', 'majorAlleleFrequency','pValueASE', 'adjPvalueASE','pValueHeterogeneity', 'adjPvalueHeterogeneity')
colnames(df1)=c('Allele', 'MAF', 'allele1IsMajor')


#' counts matrix size:
#+ matrix size, echo = F, warning = F, error = F, message = F
cat(paste("One-sample analysis with counts over", nrow(df), "transcripts\n"))

#' ### Table of adjusted p-values transcripts:
#+ Result table, echo = F, warning = F, error = F, message = F
head(df)

#' ### Table of locus specific results:
#+ Result table1, echo = F, warning = F, error = F, message = F
head(df1)

# Write csv file on the output of the tool
write.table(as.data.frame(df), file=paste0(ASEarg$title,".out.csv"), sep="\t", row.names = F, quote = FALSE)
#write.table(df,fileOutput, quote=FALSE, row.names = FALSE, sep=" ")

#' ### Histogram of adjusted p-values
#+ pvalhist, fig.width = 12, fig.height = 12, echo = F, warning = F, error = F, message = F
hist(as.numeric(df[,2]), main = "Histogram of ASE p-values", xlab = "Adjusted p-values", breaks = 5, col = "powderblue")

# Output workspace image.
save.image(paste0(ASEarg$title,".env.RData"))

# R data file output file
p <- rowRanges(ASEres)
save(p, file=paste0(ASEarg$title,".out.RData"))

#' # Session info
#+ session, echo = F, warning = F, error = F, message = F
sessionInfo()

