# Script reads annotated vcf and performs MBASED ASE analysis
#' title: "MBASED Allele-Specific Expression analysis report"
#' output:
#'   html_document:
#'       toc: true
#' ---
args <- commandArgs(trailingOnly = TRUE)

####----------------MBASED-------------------

suppressPackageStartupMessages({
  library("knitr")
  library("VariantAnnotation")
  library("MBASED")
  library("BiocParallel")
  library(ggplot2)
  library(plotly)
  library(rtracklayer)
  
  ###  devtools::install_github('hadley/ggplot2')
})

set.seed(988482)

fileInput=args[1]
fileOutput=args[2]
noSim=as.integer(args[3])

vcf = readVcf(fileInput)
vcf = readVcf('GSM1206240.decontaminated.genome.split.realigned.recalibrated.hc.vep_phased.vcf')
#------Filter SNPs 

vcf <- vcf[isSNV(vcf)]
vcf <- vcf[as.vector(geno(vcf)$PH_info!='1|1')]
vcf <- vcf[lapply(info(vcf)$CSQ,length)>0] 

#-------------Filter SNPs in exons

# Find position of GENE names in vep annotated vcf

vepD <- info(header(vcf))[c("CSQ"),]$Description
vepD <- strsplit(vepD,'Format:', fixed=TRUE)[[1]][2]
vepD <- trimws(vepD, ("b"))
vepD <- strsplit(vepD,'|', fixed=TRUE)[[1]]
pos  <- match("SYMBOL", vepD)
posID  <- match("Gene", vepD)

#vcf <- vcf[sapply(info(vcf)$CSQ, function(x) sum(grepl('exon',x)))>0]

mutInGene <- function(x, p){ strsplit(x, '|', fixed=TRUE)[[1]][p]}
geneName <- function(x, p){
  l = as.character(lapply(x, mutInGene, p))
  paste(unique(l), sep='',collapse = ';')
}

geneList <- sapply(info(vcf)$CSQ, geneName, posID)
geneNameList <- sapply(info(vcf)$CSQ, geneName, pos)

if(length(vcf)==0){
  message ("No alleles in Genes")
  quit()
  }

##------------Phase SNVs-------------------

mamatata = data.frame(Ph1 = as.data.frame(rowRanges(vcf)$REF)$x, 
                      Ph2 =  as.data.frame(rowRanges(vcf)$ALT)$value,
                      Ph1_count = sapply(geno(vcf)$AD, "[[", 1),
                      Ph2_count = sapply(geno(vcf)$AD, "[[", 2),
                      phasing_info = as.vector(unname(geno(vcf)$PH_info)), 
                      ranges=ranges(rowRanges(vcf)), 
                      seqnames=as.character(seqnames(vcf)),
                      geneName=geneList)

rownames(mamatata) = names(rowRanges(vcf))

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
  lociAllele1Counts=matrix(mamatata$Ph1_count, ncol=1, dimnames=list(names(mySNVs), 'mySample')),
  lociAllele2Counts=matrix(mamatata$Ph2_count, ncol=1, dimnames=list(names(mySNVs), 'mySample'))),
  rowRanges=mySNVs)


#----Run MBASED----PHASED--------

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
write.table(as.data.frame(df), file=paste0(fileOutput,".tsv"), sep="\t", row.names = F, quote = FALSE)
#write.table(df,fileOutput, quote=FALSE, row.names = FALSE, sep=" ")

#' ### Histogram of adjusted p-values
#+ pvalhist, fig.width = 12, fig.height = 12, echo = F, warning = F, error = F, message = F
hist(as.numeric(df[,2]), main = "Histogram of ASE p-values", xlab = "Adjusted p-values", breaks = 5, col = "powderblue")

# Output workspace image.
save.image(paste0(fileOutput,".env.RData"))

# R data file output file
p <- rowRanges(ASEres)
save(p, file=paste0(fileOutput,".out.RData"))


############ Read GTF and find gene positions#########

gtf1 <- readGFF("Homo_sapiens.GRCh37.87.gtf", version=2L)
gtf <- gtf1[gtf1$type=='gene',]
gtf <- gtf[ , colSums(is.na(gtf)) == 0]
g <- gtf[gtf$gene_id %in% df[,'Gene'],]

##  Ne slazu se dimenzije!!!!
d1 =  df[df[,'Gene'] %in% g$gene_id,]


geneChr <- g$seqid
geneID <- g$gene_id
geneStart <- g$start
geneEnd <- g$end
geneStrand <- g$strand
geneP <- sapply(geneID, function(x, df) df[df$Gene==x,]$adjPvalueASE, df)

data <- data.frame(geneChr, geneID, geneStart, geneEnd, geneP, geneStrand)


fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}



data$geneChr <- as.numeric(as.character(data$geneChr))
#data$geneChr <- as.numeric(levels(data$geneChr))[data$geneChr]


p <- ggplot() + scale_y_continuous(name="Chromosome", breaks = c(1:23)) +
                scale_x_continuous(name="Gene Position",labels=fancy_scientific) 
 
p <- p + geom_rect(data = data[data$geneStrand=='+',], aes(xmax = geneStart, xmin = geneEnd, ymax = geneChr + 0.4, ymin = geneChr + 0.1, fill=geneP, color=geneP)) 
#p <- p + scale_fill_gradient(low="green", high="green4") 
           
p <- p + geom_rect(data = data[data$geneStrand=='-',], aes(xmax = geneStart, xmin = geneEnd, ymax = geneChr - 0.1, ymin = geneChr - 0.4, fill=geneP, color=geneP))
#p <- p + scale_fill_gradient2(low="gray90", high="blue")
#p <- p + coord_cartesian(xlim = c(5000, 1000000)) 
p
 
 
p <- ggplotly(p)
p


#' # Session info
#+ session, echo = F, warning = F, error = F, message = F
sessionInfo()

 