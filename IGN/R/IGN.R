#' IGN: Invariable gene set based normalization for chromatin accessibility
#'
#' IGN is performed by normalizing the promoter chromatin accessibility signals 
#' for a given gene set that is unchanged in expression, usually obtained from 
#' accompanying RNA-seq data, and extrapolating to scale the genome-wide chromatin 
#' accessibility profile. This function allows users to normalize ATAC/DNase-seq signal
#' matrix based on promoter signal of invariable genes. 
#'
#' IGN is designed to normalize a set (two or more) ATAC-seq or DNase-seq datasets. 
#' It requires the corresponding gene expression profiling data (e.g., RNA-seq) from 
#' the same cell samples for invariable gene selection. If the accompanying gene expression 
#' data is not available, a given gene set from other sources can be used too, as long as 
#' the assumption of unchanged expression holds. The IGN method comprises three main steps:
#' 1) invariable gene set selection, 
#' 2) modeling promoter signal, 
#' 3) normalizing genome-wide signal 
#' 
#' At the invariable gene selection step, we first filter out genes with low expression levels in any conditions, 
#' i.e., only kept expressed genes with RPKM ≥ 1. Next, we select a set of genes with the least expression 
#' changes between conditions. Typically, this can be the top 200 genes with the smallest absolute value of 
#' log2-scaled fold changes. The size of the gene set, 200, is chosen considering a balance between a 
#' large-enough sample size for modeling and a minimal sampling from the whole transcriptome to avoid 
#' real differential expression. The selected gene set is referred to as invariable genes and is used in the subsequent steps.
#' 
#' At the final step, we apply these same coefficients to the genome-wide ATAC-seq signals 
#' for Condition B to normalize the genomic ATAC-seq profile.
#'
#' In case there are more than two sample/conditions to normalize, a base condition/sample 
#' is compared with each of the other conditions/samples using the same approach described 
#' above. Thus, the ATAC-seq profiles from all the other conditions/samples are normalized to the 
#' base condition and ATAC-seq data across all conditions/samples are normalized and comparable for downstream analysis. 
#'
#' 
#' @return \code{IGN} returns a list with two elements. 1) The normzlied targetSignal 
#'     matrix and 2) the normalize coefficient. 
#'     \describe{
#'         \item{\code{normzlied targetSignal}}{A signal matrix with the same dimension
#'         to the input targetSignal matirx. The first column (first sample's signal) 
#'         is the same as the input targetSiganl because all the other samples were normalized
#'         to be comparable with the first sample. This element is set to NA when the \code{targetSignal} 
#'         matrix is not inputted (only for estimating normalize coefficient).}
#'         \item{\code{normalize coefficient}}{A data frame recording the linear
#'         transformation coefficients of each sample. Note that the coefficient of the first sample is 1}
#'     }

#' @param promoterSignal Required. A matrix of promoter signal of genome-wide genes by samples. 
#'     Each row should represent a gene's promoter region (e.g., +/-1kb from TSS) and
#'     each column for a sample. The columns should be same to the \code{targetSignal} matrix
#'     (introduced below) and the rows should be the same as the \code{RNArpkm} matirx or
#'     contain the \code{InvariableGeneList} if inputted (introduced below). We recommend to 
#'     perform log transformation for the count data matrix to have better model building. 
#' @param targetSignal a matrix of signal of target regions by samples. Each row should
#'     represent a target genomic interval (e.g., merged ATAC-seq peaks from different
#'     conditions) and each column for a sample. The columns should be same to the
#'     \code{promoterSignal} matrix. We recommend to perform log transformation for the count
#'     data matrix to have better model building. The transformaion should be applied to 
#'     both \code{targetSignal} and \code{promoterSignal} matrix.
#' @param InvariableGeneList a list of pre-identified Invariable gene. All the genes in the list
#'     should be included in the promtoerSignal matrix. The promoter signal of these Invariable
#'     genes are assumed to follow the same distribution and are used in the model building.
#'     by default IGN take 200 or more Invariable genes. The Invariable gene number should be ≥ 50
#' @param RNArpkm Required if no \code{InvariableGeneList} inputted. A matrix of expression index (RPKM) of genome-wide genes by samples. 
#'     this parameter only takes effect when \code{InvariableGeneList} is not inputted. 
#'     Genes with sufficient expression index (RPKM higher than \code{RPKMcutoff}) in all samples/columns
#'     are selected as the candidate of Invariable genes. Users can also input other types of expression index 
#'     (e.g., signal from microarry) and customize the cutoff specified by \code{RPKMcutoff}.
#' @param RNAlogfc Required if no \code{InvariableGeneList} inputted. A matrix of fold change (log scale) in all comparison. Users can calculate the log fold change 
#'     using the \code{RNArpkm} matrix or using some external software (e.g., apeglm for log fold change shrinkage).
#'     IGN sort invariable gene candidates (i.e., genes with RPKM ≥ 1) by the absolute log fold change and the top 200 
#'     (default) genes with lowest absolute log fold change are selected as Invariable genes. The top Invariable gene number is
#'     specified by the parameter \code{topInvariableGene}. For more than one column provided, IGN will take average of the 
#'     absolute log fold change for sorting. 
#' @param RPKMcutoff cutoff for selecting sufficient expression gene. Default is 1 (for RPKM matrix)
#' @param topInvariableGene number of top invariable genes with lowest (absolute) expression fold changes. Default is 200 (genes).
#' @examples 
#' ## run IGN on testing dataset
#' library(IGN)
#' # run IGN with pre-identified gene list
#' data(ATACpromoterSignal,package="IGN")
#' data(ATACtargetSignal,package="IGN")
#' data(SDgene,package="IGN")
#' testoutList <- IGN(ATACpromoterSignal,ATACtargetSignal,InvariableGeneList=SDgene)
#' # normalized target signal matrix
#' head(testoutList[[1]])
#' # normalize coefficient
#' testoutList[[2]]
#' # run IGN with expression data
#' data(ATACpromoterSignal,package="IGN")
#' data(ATACtargetSignal,package="IGN")
#' data(RNArpkm,package="IGN")
#' data(RNAlogfc,package="IGN")
#' testoutRPKM <- IGN(ATACpromoterSignal,ATACtargetSignal,RNArpkm,RNAlogfc)
#' # normalized target signal matrix
#' head(testoutRPKM[[1]])
#' # normalize coefficient
#' testoutRPKM[[2]]
#' @keywords IGN
#' @export

IGN <- function( promoterSignal, targetSignal=NULL, RNArpkm=NULL, RNAlogfc=NULL, 
                    topInvariableGene=200, InvariableGeneList=NULL , RPKMcutoff = 1  ){
  if(is.null(InvariableGeneList)){
    #message("Identify Invariable genes based on expression data")
    if(is.null(RNArpkm)){
      stop("No RNArpkm table inputted. Fail to identify invariable gene",call.=FALSE)
    }
    if(is.null(RNAlogfc)){
      stop("No RNAlogfc table inputted. Fail to identify invariable gene",call.=FALSE)
    }
    geneIDX <- intersect(intersect(rownames(RNArpkm),rownames(RNAlogfc)),rownames(promoterSignal))
    RPKM <- as.matrix(RNArpkm[geneIDX,])
    LFC <- as.matrix(RNAlogfc[geneIDX,])
    rownames(LFC) <- geneIDX#rownames(RNAlogfc)
    colnames(LFC) <- colnames(RNAlogfc)
    highexpLFC <- as.matrix(LFC[which(apply(RPKM,1,min) >= RPKMcutoff), ])
    InvariableGene <- rownames(highexpLFC)[order(apply(abs(highexpLFC),1,mean),na.last=T)[1:topInvariableGene]]
    
  }else{
    #message("InvariableGene list inputted. Use Invariable genes defined in the InvariableGeneList")
    InvariableGene <- intersect(InvariableGeneList, rownames(promoterSignal))
    if(length(InvariableGene) < 50){
      stop(paste0("Too few invariableGene (",length(InvariableGene) ," genes) detected for model building"),call.=FALSE)
    }
  }
  
  promoter_Sig <- promoterSignal[InvariableGene,]
  allSamples <- colnames(promoter_Sig)
  #message("Model building with promoter signal of Invariable genes")
  if(is.null(targetSignal)){
    #message("No target region matrix inputted. Only estimate the coefficients from the model")
    normalizeCoeff <- c(1,0)
    for(SAMPLE in allSamples[2:length(allSamples)]){
      X <- as.numeric(promoter_Sig[,SAMPLE])
      Y <- as.numeric(promoter_Sig[,allSamples[1]])
      thismodel <- lm(Y~X)
      a <- thismodel$coefficients[2]
      b <- thismodel$coefficients[1]
      normalizeCoeff <- rbind(normalizeCoeff,c(a,b))
      normalizedData <- NULL
    }
  }else{
    if( !setequal(colnames(targetSignal), colnames(promoter_Sig))){
      stop("The samples(columns) are different between promoterSignal matrix and targetSig matrix",call.=FALSE)
    }
    normalizedData <- as.numeric(targetSignal[,allSamples[1]])
    normalizeCoeff <- c(1,0)
    for(SAMPLE in allSamples[2:length(allSamples)]){
      X <- as.numeric(promoter_Sig[,SAMPLE])
      Y <- as.numeric(promoter_Sig[,allSamples[1]])
      thismodel <- lm(Y~X)
      a <- thismodel$coefficients[2]
      b <- thismodel$coefficients[1]
      thisnormdata <- targetSignal[,SAMPLE]*a+b
      normalizedData <- cbind(normalizedData, thisnormdata)
      normalizeCoeff <- rbind(normalizeCoeff,c(a,b))
    }
    rownames(normalizedData) <- rownames(targetSignal)
    colnames(normalizedData) <- colnames(targetSignal)
  }
  rownames(normalizeCoeff) <- colnames(promoter_Sig)
  colnames(normalizeCoeff) <- c("a_slope","b_intercept")
  
  return(list(normalizedData,normalizeCoeff))
}
