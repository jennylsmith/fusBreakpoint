
#' Define Fusion Category for Harmonization Steps
#'
#' @param geneA column name which contains the gene symbol for fusion partner A
#' @param geneB column name which contains the gene symbol for fusion partner A
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' fusion.df <- data.frame(geneNameA="NUP98", geneNameB="KDM5A", Breakpoint="chr11:0000000|chr12:0000000") %>%
#'  rowwise() %>%
#'  mutate(Fusion.Category=fusion_category(geneNameA, geneNameB)) %>%
#'  ungroup()
#' }
fusion_category <- function(geneA,geneB){
  fusion <- c(geneA,geneB)[order(c(geneA,geneB))]
  fusion <- paste(fusion, collapse="-")
  return(fusion)
}


#' Order breakpoints Numerically for harmonization steps
#'
#' @param brkpt the column name containing the breakpoint position
#' @param sep the seperator between the two breakpoints
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' fusion.df <- data.frame(geneNameA="NUP98", geneNameB="KDM5A", Breakpoint="chr12:0000000|chr11:0000000") %>%
#'  rowwise() %>%
#'  mutate(Breakpoint.Ordered=order_breakpoints(Breakpoint)) %>%
#'  ungroup()
#'
#' }
order_breakpoints <- function(brkpt, sep="|"){
  split <- str_split(brkpt, pattern = paste0("\\", sep))[[1]]
  ord <- split[order(split)]
  ord <- paste(ord,collapse = sep)
  return(ord)
}


#' Create JSON file for input into St. Jude Genome Paint
#'
#' @param df dataframe for a single fusion  with the breakpoint of the fusion-genes divided into 4 columns, chromosome and position for each gene in the fusion.
#' @param fusion_genes_GR fusion_genes_GR is granges list object that contains genes that involved in the fusions in the input dataframe
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' genesGR <- transcriptsBy(Grch37.txdb, by="gene")
#' fusion.df <- data.frame(Fusion.Category="CBFB-MYH11",
#'  FusionID="My Fusion Name",
#'  chromosomes_A="chr16",
#'  pos_A=00000100,
#'  chromosomes_B="chr16",
#'  pos_B=00000200)
#'
#' formatted <- genome_paint_format(df=fusion.df, fusion_genes_GR = genesGR)
#'
#' }
genome_paint_format <- function(df,fusion_genes_GR){
  #https://docs.google.com/document/d/1owXUQuqw5hBHFERm0Ria7anKtpyoPBaZY_MCiXXf5wE/edit#
  #df is a dataframe for a single fusion  with the breakpoint of the fusion-genes divided into 4 columns, chromosome and position for each gene in the fusion.
  #fusion_genes_GR is granges list object that contains each individual gene of interest - eg NUP98 and all its associated partners.

  grA <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(df[["chromosomes_A"]]),
                                ranges = IRanges(start=df[["pos_A"]]))

  grB <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(df[["chromosomes_B"]]),
                                ranges = IRanges(start=df[["pos_B"]]))

  #Find which NUP98-Fusion gene belongs with the genomic coordinate
  breakpoint_infoA <- IRanges::subsetByOverlaps(fusion_genes_GR, grA)
  breakpoint_infoB <- IRanges::subsetByOverlaps(fusion_genes_GR, grB)

  #Define first gene
  geneA <- names(breakpoint_infoA)
  breakpoint_infoA <- unlist(breakpoint_infoA)
  chrA <- seqnames(breakpoint_infoA)
  chrA <- paste0("chr",unique(as.vector(chrA)))
  posA <- as.numeric(as.character(ranges(grA)))
  strA <- strand(breakpoint_infoA)
  strA <- unique(as.vector(strA))


  #define the second gene
  geneB <- names(breakpoint_infoB)
  breakpoint_infoB <- unlist(breakpoint_infoB)
  chrB <- seqnames(breakpoint_infoB)
  chrB <- paste0("chr",unique(as.vector(chrB)))
  posB <- as.numeric(as.character(ranges(grB)))
  strB <- strand(breakpoint_infoB)
  strB <- unique(as.vector(strB))

  #Example SV
  # chr12	12029972	12029972	{"chrB": "chr21", "dt": 5, "posB": 36417895, "sample": "SJETV039_D", "strandA": "+", "strandB": "-"}
  # chr21	36417895	36417895	{"chrA": "chr12", "dt": 5, "posA": 12029972, "sample": "SJETV039_D", "strandA": "+", "strandB": "-"}

  #Example Fusion
  #chr1  19157328  19157328  {"chrB": "chr13", "dt": 2, "geneA": "chr1", "geneB": "VWA8", "posB": 42278669, "sample": "SJRHB009_D", "strandA": "-", "strandB": "+"}
  #chr1  19157486  19157486  {"chrA": "chr13", "dt": 2, "geneA": "VWA8", "geneB": "chr1", "posA": 42278658, "sample": "SJRHB009_D", "strandA": "+", "strandB": "-"}

  #will use Rjson
  ID <- df[["FusionID"]]
  jsonA <- rjson::toJSON(list("chrB"=chrB,"dt"= "2", "geneA"=geneA, "geneB"=geneB,
                              "posB"=posB, "sample"=ID,"strandA"=strA,"strandB"=strB))
  jsonB <- rjson::toJSON(list("chrA"=chrA, "dt"="2", "geneA"=geneB,"geneB"=geneA,
                              "posA"=posA, "sample"=ID, "strandA"=strA, "strandB"=strB))

  #final format to upload to genome paint
  out <- data.frame(chr=c(chrA,chrB),
                    start=c(posA,posB),
                    end=c(posA,posB),
                    json=c(jsonA, jsonB))

  return(out)
}


#' Add Intron Rank Column to GRanges containing gene introns
#'
#' @param intron_GR GRanges list object that contains intron ranges of genes involved in the fusions in the input dataframe
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' intronGR <- intronsByTranscript(Grch37.txdb, use.names=TRUE)
#' update.intronGR <- mcol_intron_rank(intronGR)
#' }
mcol_intron_rank <- function(intron_GR){
  for (i in seq_along(intron_GR)){
    strand <- BiocGenerics::strand(intron_GR[[i]])
    strand <- unique(as.character(strand))

    n <-  nrow(S4Vectors::mcols(intron_GR[[i]]))
    if(strand == "-"){
      S4Vectors::mcols(intron_GR[[i]])$intron_rank <- paste0("intron",n:1)
    }else{
      S4Vectors::mcols(intron_GR[[i]])$intron_rank <- paste0("intron", 1:n)
    }
  }

  return(intron_GR)
}



#' Provides Exon Number from Genomic Coordinates.
#'
#' @param df dataframe for a single fusion  with the breakpoint of the fusion-genes divided into 4 columns, chromosome and position for each gene in the fusion.
#' @param exon_ranges_GR GRanges list object that contains exon ranges of transcripts involved in the fusions in the input dataframe. No duplicated gene IDs - one transcript ID per one gene ID
#' @param intron_ranges_GR GRanges list object that contains intron ranges of transcripts involved in the fusions in the input dataframe. No duplicated gene IDs - one transcript ID per one gene ID
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' my_transcript_IDs <- c("ENSTXXX","ENSTXXXXX")
#' intronGR <- intronsByTranscript(Grch37.txdb, use.names=TRUE)[my_transcript_IDs]
#' update.intronGR <- mcol_intron_rank(intronGR)
#' exonGR <- exonsBy(Grch37.txdb, by="tx", use.names=TRUE)[my_transcript_IDs]
#' fusions.df <- data.frame(chromosome_A="chr16",pos_A=00000100, chromosome_B="chr16", pos_B=00000200)
#' junctions.df <- define_exon_junctions(fusions.df,exonGR,update.intronGR)
#' }
#'
#' @import dplyr
define_exon_junctions <- function(df, exon_ranges_GR, intron_ranges_GR=NULL){
  grA <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(df[["chromosomes_A"]]),
                                ranges = IRanges(start=df[["pos_A"]]))

  grB <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(df[["chromosomes_B"]]),
                                ranges = IRanges(start=df[["pos_B"]]))

  #Find which Fusion gene exon belongs with the genomic coordinate
  exon_infoA <- IRanges::subsetByOverlaps(exon_ranges_GR, grA)
  exon_infoB <- IRanges::subsetByOverlaps(exon_ranges_GR, grB)

  #Check if any of the breakpoints actually fall within an intron range
  non_exonic <- sapply(list(exon_infoA, exon_infoB), function(x) length(x)==0)
  if(any(non_exonic) & !is.null(intron_ranges_GR)){
    to_update <- c("A","B")[which(non_exonic)]
    for(i in 1:length(to_update)){
      in_var_name <- paste0("gr", to_update)[i]
      out_var_name <- paste0("exon_info", to_update)[i]
      introns <- IRanges::subsetByOverlaps(intron_ranges_GR, get(in_var_name))
      assign(out_var_name, introns)
    }
  }

  #Reformat the metadata cols (mocols)
  exon_infoA <- IRanges::subsetByOverlaps(unlist(exon_infoA),grA)
  exon_infoA <- as.data.frame(exon_infoA) %>%
    tibble::rownames_to_column("gene") %>%
    rename_all(~paste0(.,"A"))

  exon_infoB <- IRanges::subsetByOverlaps(unlist(exon_infoB), grB)
  exon_infoB <- as.data.frame(exon_infoB) %>%
    tibble::rownames_to_column("gene") %>%
    rename_all(~paste0(.,"B"))

  #if the genomic ranges does not contain the breakpoint location, the dataframe will be empty
  no_data <- sapply(list(exon_infoA, exon_infoB), function(x) nrow(x))  == 0
  if(any(no_data)){
    to_update <- c("A","B")[which(no_data)]
    for(i in 1:length(to_update)){
      var_name <- paste0("exon_info", to_update)[i]
      obj <- get(var_name)
      dummy_data <- rep(NA,length(colnames(obj)))
      names(dummy_data) <- colnames(obj)

      obj <- obj %>%
        bind_rows(dummy_data)
      assign(var_name, obj)
    }
  }

  #rowbind the metadata cols for the breakpoint pair.
  res <- dplyr::bind_cols(df, exon_infoA,exon_infoB) %>%
    #Sample=df[["Sample"]],
    mutate(Fusion=paste(geneA,geneB, sep="-")) %>%
    unite("Exons", matches("_rank[AB]"),sep = "_",remove = F)


  #For consistency, need to add empty columns for the intron ranks
  if (!"intron_rankA" %in% colnames(res)){
    res[["intron_rankA"]] <- NA
  }

  if (!"intron_rankB" %in% colnames(res)){
    res[["intron_rankB"]] <- NA
  }


  return(res)
}





