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

# Find position of GENE names in vep annotated vcf

vepD <- info(header(vcf))[c("CSQ"),]$Description
vepD <- strsplit(vepD,'Format:', fixed=TRUE)[[1]][2]
vepD <- trimws(vepD, ("b"))
vepD <- strsplit(vepD,'|', fixed=TRUE)[[1]]
pos  <- match("Gene", vepD)

vcfSNP <- vcf[sapply(info(vcfSNP)$CSQ, function(x) sum(grepl('exon',x)))>0]
mutInGene <- function(x, p){ s = x[grepl('exon',x)]; strsplit(s[1], '|', fixed=TRUE)[[1]][p]}
geneName  <- function(x, p){
    l = as.character(lapply(x, mutInGene, p))
    paste(unique(l), sep='',collapse = ';')
}

geneList <- sapply(info(vcfSNP)$CSQ, geneName, posID)

##------------Phase SNVs-------------------

mamatata = data.frame(Ph1 = as.data.frame(rowRanges(vcf)$REF)$x,
                      Ph2 =  as.data.frame(rowRanges(vcf)$ALT)$value,
                      Ph1_count = sapply(geno(vcf)$AD, "[[", 1),
                      Ph2_count = sapply(geno(vcf)$AD, "[[", 2),
                      phasing_info = as.vector(unname(geno(vcf)$PH_info)),
                      ranges=ranges(rowRanges(vcf)),
                      seqnames=as.character(seqnames(vcf)))

rownames(mamatata) = names(rowRanges(vcfP))

#Swap Ph1, Ph2, Ph1_count, Ph2_count for phasing_info = 1|0
mamatata[mamatata$phasing_info == "1|0", c("Ph1", "Ph2", "Ph1_count", "Ph2_count")] <-
mamatata[mamatata$phasing_info == "1|0", c("Ph2", "Ph1", "Ph2_count", "Ph1_count")]


##-----------MBASED--Prepare samples-------------------


mySNVs <- GRanges(seqnames=mamatata$seqnames,
                  ranges=IRanges(start=mamatata$ranges.start,
                                 end = mamatata$ranges.end,
                                 width = mamatata$ranges.width,
                                 names = mamatata$ranges.names),
                  aseID=geneList,
                  allele1=as.character(mamatata$Ph1),
                  allele2=as.character(mamatata$Ph2))


mySample <- SummarizedExperiment(assays=list(
    lociAllele1Counts=matrix(mamatata$Ph1_count, ncol=1, dimnames=list(names(mySNVsPhased), 'mySample')),
    lociAllele2Counts=matrix(mamatata$Ph2_count, ncol=1, dimnames=list(names(mySNVsPhased), 'mySample'))),
    rowRanges=mySNVs)


#----Run MBASED


ASEres <- runMBASED(
ASESummarizedExperiment=mySample,
isPhased=TRUE,
numSim=noSim,
BPPARAM = MulticoreParam())


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

