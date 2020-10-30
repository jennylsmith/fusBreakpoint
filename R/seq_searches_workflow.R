#Jenny Smith
#Sept. 17, 2020
#purpose: examine junction sequences in RNA-seq reads.


###### Subset BAMS
#' Subset a bam for specific chromome(s) or ranges based on input breakpoint or fusion data file.
#'
#' @param Sample_ID character string of patient ID
#' @param breakpoint  character vector with breakpoint string (can be named)
#' @param fusion_data dataframe with fusion results from CICERO, TransAbyss, and STAR-Fusion concatenated
#' @param file_manifest dataframe with Sample_ID and the full file path of the BAM file.
#' @param by_chr boolean - get reads mapped to the entire chromosome containing the breakpoint or just +/- 1e6 bp from breakpoint
#' @param scan boolean - read into memory or not
#' @param Seqlens named numeric vector with chromosome lengths
#' @param outdir destination for subsetted bam files if scan = FALSE
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' my_files_df <- data.frame(Sample="PAZZUM.03A.01R", filepath=dir("/path/to/sample.bam"))
#' subset_bam("PAZZUM.03A.01R", breakpoint=list("16:1250000|16:1260000"), my_files_df)
#'}
#' @import dplyr
subset_bam <- function(Sample_ID, breakpoint=NULL,
                       fusion_data=NULL,
                       file_manifest,
                       by_chr=FALSE, scan=FALSE,Seqlens=NULL,
                       outdir=file.path("/fh/scratch/delete90/jlsmith3/subset_bams")){


  if(by_chr){
    offset <- FALSE
    if(is.null(Seqlens)){
      message("For whole chromosome subset, must provide seqlens vector.")
      return(Seqlens)
    }
  }else{
    offset <- TRUE
  }

  if(all(is.null(breakpoint), is.null(fusion_data))){
    message("Breakpoint and fusion data cannot both be NULL. Please provide either a breakpoint as character vector
    with the format of `chr:basepairs|chr:basepairs`, for example c(CBFB.MYH11=\"16:15814908|16:67116211\", KMT2A.MLLT10=\"10:21959378|11:118353210\").
    Or fusion data  is a dataframe with the fusion breakpoints.")
    return(NULL)
  }


  if(!is.null(fusion_data)){
    #Define samples and fusions
    fusion_data <- filter(fusion_data, Sample==Sample_ID)
    sample <- unique(pull(fusion_data,"Sample"))
    fusion <- paste(pull(fusion_data, "Fusion.Category"), collapse = "_")
    print(c(sample,fusion))

    #Define breakpoints
    breakpoint <- pull(fusion_data, "breakpoint.TA")
    breakpoint <- ifelse(is.na(breakpoint), pull(fusion_data, "Breakpoints.STAR"), breakpoint)

    #Define output file name
    fname <- paste0(Sample_ID,"_",fusion,".bam")
  }else{
    #Define output file name
    fname <- paste0(Sample_ID,"_subset.bam")
  }



  create_ranges <- function(x,offset=TRUE){
    if(offset){
      #plus/minus 1 million bp from the breakpoint
      pos <- as.numeric(x[2])
      pos <- IRanges::IRanges(start=pos-1e6, end=pos+1e6)
    }else{
      chr <- x[1]
      #this is sepecific to GRCh37-lite and reference from Starfusion CTAT.
      #Will need to make more generalized ASAP.
      pattern <- paste(c(paste0("chromosome ",chr,","),
                       paste0(chr,"\\s[0-9XYMT]{1,2}$")),
                       collapse = "|")
      pos <- IRanges::IRanges(start=0, end=Seqlens[grep(pattern,names(Seqlens))])
    }
    return(pos)
  }

  #Define breakpoints and basepair ranges
  breakpoint <- stringr::str_split(breakpoint,"\\|")
  if(length(breakpoint) > 1){
    seqnames <- sapply(breakpoint,
                       function(chr_pos) sapply(stringr::str_split(chr_pos,":"),function(chr) chr[1]))

    chr_ranges <- unlist(lapply(breakpoint,
                                function(chr_pos) lapply(stringr::str_split(chr_pos,":"), create_ranges, offset=offset)))

  }else{
    breakpoint <- breakpoint[[1]]
    seqnames <- sapply(stringr::str_split(breakpoint,":"),function(x) x[1])
    chr_ranges <- lapply(stringr::str_split(breakpoint,":"),create_ranges, offset=offset)
  }


  #Define location of bam file
  bam <- filter(file_manifest,Sample==Sample_ID) %>%
    pull(filepath)

  if(length(bam) > 1){
    message("File manifest has duplicates! Must have one bam file per patient ID.")
    return(filter(file_manifest,Sample==Sample_ID))
  }


  bam_file <- Rsamtools::BamFile(bam)
  if(length(bam_file$index) == 0){
    #does NOT listen to the destination parameter... could cause issues.
    Rsamtools::indexBam(bam, destination=dest)
    bam_file <- Rsamtools::BamFile(bam,asMates=TRUE)
  }

  #Define the ranges to subset
  n <- length(seqnames)
  ranges <- unlist(IRanges::IRangesList(chr_ranges))
  gr <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(seqnames,rep(1,n)),
                               ranges = ranges)
  gr <- unique(gr) #incase of intrachromosomal fusions, dont want to duplicate the same ranges

  # define parameters. Check other options to filter by with scanBamWhat()
  params <- Rsamtools::ScanBamParam(which = gr, what = Rsamtools::scanBamWhat())
  dest <- outdir

  if(scan){
    #reads in the subsetted bam file into memory
    bam_reads <- Rsamtools::scanBam(bam_file, param = params)
    return(bam_reads)
  }else{
    #saves the subsetted bam file
    Rsamtools::filterBam(bam_file, param = params,
                         destination = paste0(dest,fname),
                         indexDestination=FALSE)
  }

}

#### Create pdict for sequence searches
#' Create Biostrings pdict (pattern dictionary) from a dataframe for effecient sequence searches of RNAseq reads.
#'
#' @param theoretical_seqs_df a dataframe with computationally derived exon junction sequences
#'
#' @return
#' @export
#'
#' @examples
#'my_seqs <- data.frame(Name=as.character(1:10),Sequence=rep(c("ATCGCCCGTTA"), 10))
#'my_pdict_obj <- create_custom_pdict(theoretical_seqs_df=my_seqs)
#'
#' @import dplyr
create_custom_pdict <- function(theoretical_seqs_df){

  #computationsal Breakpoint sequences per transcript isoform to search for
  #this will need to be run 2-4x bc the nearly 10k sequences requires way too much memory
  nr <- nrow(theoretical_seqs_df)
  print(paste("There are ",nr," fusion sequences to examine."))
  if(nr > 500){
    slices <- split(1:nr, cut(1:nr, 4, labels = FALSE))
  }else{
    slices <- list(1:nr)
  }

  #Create PDict object with the computationally derived exon-exon junction sequences to search for.
  seqExonsList <- list()
  for(i in 1:length(slices)){
    seqExons <- theoretical_seqs_df %>%
      dplyr::slice(slices[[i]]) %>%
      dplyr::pull(Sequence, name=Name) %>%
      Biostrings::DNAStringSet()

    #Add reverse compliment to search
    revComp <- Biostrings::reverseComplement(seqExons)
    names(revComp) <- paste0("revComp",names(seqExons))

    #combine the sequences and the reverse complement  into one DNAstringset
    seqExons <- c(seqExons, revComp)
    seqExons <- seqExons[order(names(seqExons))]

    #convert into pdict
    seqExons <- Biostrings::PDict(seqExons)  #(3) later matchPdict can only be used with max.mismatch=0.???
    seqExonsList[[i]] <- seqExons

    rm(seqExons)
  }

  message(paste("Finished PDict processing with",sum(sapply(seqExonsList, length)),"theoretical fusion sequences."))


  return(seqExonsList)
}




######### Sequence Pattern Matching

#'  Sequence pattern matching of RNAseq reads and counts the number of breakpoint junction reads.
#'
#' @param Sample_ID character string of patient ID
#' @param pdict_obj pdict class object form biostrings.
#' @param RNAseqReads DNAstring set object with RNAseq reads to query
#' @param contigSeqs DNAstring set object with fusion algorithm (eg STAR-fusion) derived fusion contig sequences
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' my_seqs <- data.frame(Name=paste0("seq",1:10),Sequence=rep(c("ATCGCCCGTTA"), 10))
#' my_files_df <- data.frame(Sample="PAZZUM.03A.01R", filepath=dir("/path/to/sample.bam"))
#'
#' my_contigs <- data.frame(contig="AATGGCTATC", Name=paste0("seq",1)) %>% pull(contig,name=Name)
#' my_contigs <- Biostrings::DNAStringSet(my_contigs)
#'
#' my_pdict_obj <- create_custom_pdict(theoretical_seqs_df=my_seqs)
#' my_bam_reads <- subset_bam("PAZZUM.03A.01R", breakpoint=list("16:1250000|16:1260000"), my_files_df)
#'
#' evidence <- match_reads_contigs("PAZZUM.03A.01R", my_pdict_obj, my_bam_reads, contigSeqs=my_contigs)
#' }
#' @import dplyr
match_reads_contigs <- function(Sample_ID, pdict_obj, RNAseqReads, contigSeqs=NULL){

  #match RNA-seq reads (Error: cannot allocate vector of size 154.7 Gb)
  read.matches <- Biostrings::vcountPDict(pdict = pdict_obj,
                              subject = RNAseqReads,
                              max.mismatch = 0, min.mismatch = 0,
                              with.indels = FALSE)
  rownames(read.matches) <- names(pdict_obj)

  #Create a tally of the RNAseq read hits.
  #Add column for total number of RNA-seq reads searched!
  read.match.df <- data.frame(Num_Junction_Breakpoint_Reads=rowSums(read.matches)) %>%
    rownames_to_column("Theoretical_Exon_Junctions") %>%
    mutate_at(vars(Theoretical_Exon_Junctions), ~gsub("revComp","", .)) %>%
    group_by(Theoretical_Exon_Junctions) %>%
    mutate(Num_Junction_Breakpoint_Reads=sum(Num_Junction_Breakpoint_Reads)) %>%
    ungroup() %>%
    filter(!duplicated(Theoretical_Exon_Junctions)) %>%
    mutate(Num_RNAseq_Reads_Searched=length(RNAseqReads))


  #Final Results table
  rm(read.matches)
  gc()
  results <- tibble(Sample=rep(Sample_ID,
                               nrow(read.match.df))) %>%
    bind_cols(.,read.match.df)

  #Match contig sequences
  if(!is.null(contigSeqs)){
    contig.matches <- Biostrings::vcountPDict(pdict = pdict_obj,
                                  subject = contigSeqs,
                                  max.mismatch = 0, min.mismatch = 0,
                                  with.indels = FALSE)
    rownames(contig.matches) <- names(pdict_obj)
    colnames(contig.matches) <- names(contigSeqs)

    #Create a tally of the contig hits
    contig.match.df <- data.frame(Num_Contig_Matches=rowSums(contig.matches)) %>%
      rownames_to_column("Theoretical_Exon_Junctions") %>%
      mutate_at(vars(Theoretical_Exon_Junctions), ~gsub("revComp","", .)) %>%
      group_by(Theoretical_Exon_Junctions) %>%
      mutate(Num_Contig_Matches=sum(Num_Contig_Matches)) %>%
      ungroup() %>%
      filter(!duplicated(Theoretical_Exon_Junctions))

    #Also return the contig sequence itself
    contig.seq.df <- data.frame(contig.matches) %>%
      rownames_to_column("Theoretical_Exon_Junctions") %>%
      mutate_at(vars(Theoretical_Exon_Junctions), ~gsub("revComp","", .)) %>%
      #make seq and revcomp seq identical results (eg. a hit for seq or revcomp seq considered a hit)
      group_by(Theoretical_Exon_Junctions) %>%
      mutate_at(vars(matches("contig")), ~sum(.)) %>%
      ungroup() %>%
      #make into long format and substitue a hit ==1 to hit == contig.sequence
      tidyr::gather(Contig.Name, Contig.Match, matches("^contig")) %>%
      group_by(Contig.Name) %>%
      mutate_at(vars(Contig.Match), ~case_when(
        . == 0 ~ as.character(.),
        . == 1 ~ as.character(contigSeqs[[unique(Contig.Name)]]))) %>%
      group_by(Theoretical_Exon_Junctions, .add=TRUE) %>%
      filter(!duplicated(Theoretical_Exon_Junctions)) %>%
      ungroup() %>%
      tidyr::spread(Contig.Name,Contig.Match) %>%
      rename_all(~gsub("contig.", "Matching_Contig_for_Exon_Junctions_", .))

    #Final Results table
    rm(contig.matches)
    gc()
    results <- results %>%
      left_join(., contig.match.df, by="Theoretical_Exon_Junctions") %>%
      left_join(., contig.seq.df, by="Theoretical_Exon_Junctions")
  }


  return(results)

}



######## workflow
#' Workflow wrapper to extract RNAseq reads and search for custom breakpoint exon junction sequences.
#'
#' @param Sample_ID character string of patient ID
#' @param fusion_data dataframe with fusion results from CICERO, TransAbyss, and STAR-Fusion concatenated
#' @param breakpoint list class with breakpoint string
#' @param file_manifest dataframe with Sample_ID and the full file path of the BAM file.
#' @param Seqlens named numeric vector with chromosome lengths
#' @param fusion_theoretical_pdict a pdict obj with computationally derived exon junction sequences
#' @param cicero.contig a dataframe with contig sequences from CICERO output
#' @param TA.contigs a dataframe with contig sequences from TransAbyss output
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' my_seqs <- data.frame(Name=paste0("seq",1:10),Sequence=rep(c("ATCGCCCGTTA"), 10))
#' my_files_df <- data.frame(Sample="PAZZUM.03A.01R", filepath=dir("/path/to/sample.bam"))
#' chr_lengths <- c("16"=24000000)
#'
#' my_pdict_obj <- create_custom_pdict(theoretical_seqs_df=my_seqs)
#'
#' examine_breakpoint_evidence("PAZZUM.03A.01R", breakpoint=list("16:1250000|16:1260000"),
#' my_files_df,chr_lengths, my_pdict_obj, cicero_contig_df, transabyss_contig_df)
#'}
#' @import dplyr
examine_breakpoint_evidence <- function(Sample_ID,
                                        fusion_data=NULL,
                                        breakpoint=NULL,
                                        file_manifest,
                                        Seqlens,
                                        fusion_theoretical_pdict,
                                        cicero.contig,
                                        TA.contigs){

  print(paste("Processing Sample:",Sample_ID))

  #Begin subsetting the RNAseq BAM file reads
  BAM <- subset_bam(Sample_ID = Sample_ID,
                    breakpoint=breakpoint,
                    fusion_data = fusion_data,
                    file_manifest = file_manifest,
                    scan = TRUE,
                    by_chr=TRUE,
                    Seqlens=Seqlens)

  #Select only the RNAseq reads to keep in memory
  seqs <- lapply(BAM, `[[`, "seq")
  read_names <- lapply(BAM, `[[`, "qname")
  print(paste("Subsetting BAM for the Ranges:", paste(names(seqs),collapse="; " )))


  #for loop to concatentate the ranges into a single DNAstringset object.
  RNAseq <- Biostrings::DNAStringSet() #unlist() doesnt work with the seqs object for some reason, and Ive tried all lapply/sapply possibities.
  for(i in 1:length(seqs)){
    reads <- seqs[[i]]
    names(reads) <- read_names[[i]]
    RNAseq <- Biostrings::DNAStringSet(c(RNAseq,reads))
  }
  rm(BAM)

  message(paste("Finished subsetting  BAM file for RNAseq read sequences.
                There are",length(RNAseq),"RNAseq reads to search."))

  #Contigs to Search
  try(
    cicero <- cicero.contig %>%
      filter(grepl(Sample_ID,Sample)) %>%
      mutate(Name=paste0("contig.cicero.",1:nrow(.))) %>%
      pull(contig,name=Name) %>%
      Biostrings::DNAStringSet(),
    silent = T
  )

  #Contigs to Search
  try(
    TA <- TA.contigs %>%
      filter(grepl(Sample_ID,Sample)) %>%
      mutate(Name=paste0("contig.TA.",1:nrow(.))) %>%
      pull(contig, name=Name) %>%
      Biostrings::DNAStringSet(),
    silent = T
  )


  #Combine all contigs
  if(all(exists("TA") & exists("cicero"))){
    contigs <- c(cicero,TA)
    rm(TA,cicero)
  }else if(all(exists("TA") & !exists("cicero"))){
    contigs <- TA
    rm(TA)
  }else if(all(!exists("TA") & exists("cicero"))){
    contigs <- c(cicero)
    rm(cicero)
  }

  if(exists("contigs")){
    message("Finished subsetting fusion algorithm contig sequences.")
  }else{
    message("No contig sequences found.")
    contigs <- NULL
  }

  #Search for the computationally derived seqs
  final_results <- NULL
  for(i in 1:length(fusion_theoretical_pdict)){
    theoretic_txs <- match_reads_contigs(Sample_ID=Sample_ID,
                                         pdict_obj = fusion_theoretical_pdict[[i]],
                                         RNAseqReads = RNAseq,
                                         contigSeqs = contigs)
    if(is.null(final_results)){
      final_results <- theoretic_txs
    }else{
      final_results <- bind_rows(final_results,
                                 theoretic_txs)
    }
    rm(theoretic_txs) #remove large objects from memory
    message("Finished search for ", paste(i," out of", length(fusion_theoretical_pdict)),
            " parts of theoretical sequences.")
  }

  message("Finished search all theoretical sequences.")
  return(final_results)
}




#
#' Function to create a fusion category by alphabetical order
#'
#' @param geneA the column name that contains fusion partner #1
#' @param geneB the column name that contains fusion partner #2
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' my_fusions <- fusion_df %>%
#' rowwise() %>%
#' mutate(Fusion.Category=fusion_category(geneA=geneNameA, geneB=geneNameB)) %>%
#' ungroup()
#'}
fusion_category <- function(geneA,geneB){
  fusion <- c(geneA,geneB)[order(c(geneA,geneB))]
  fusion <- paste(fusion, collapse="-")
  return(fusion)
}



#
#' Function to sort breakpoints for easier matching of identical breakpoints across different datasets
#'
#' @param brkpt a column name for the breakpoint information
#' @param sep is a string for the seperator used between breakpoint junctions eg, "16:15000000|16:16000000" the sep is "|" pipe.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' my_fusions <- fusion_df %>%
#' rowwise() %>%
#' mutate(Breakpoint_Ordered=order_breakpoints(breakpoint_column)) %>%
#' ungroup()
#'}
order_breakpoints <- function(brkpt, sep="|"){
  split <- str_split(brkpt, pattern = paste0("\\", sep))[[1]]
  ord <- split[order(split)]
  ord <- paste(ord,collapse = sep)
  return(ord)
}



#
#' A function to reformat a breakpoint information for genome paint
#'
#' @param df data frame with the column names "chomosomes_A", "pos_A", "chromosomes_B", "pos_B"
#' @param fusion_genes_GR the GR ranges object with the fusion genes positions
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' library(HS.txDB)
#' library(GenomicFeatures)
#' my_fusions <- data.frame(geneA="NUP98", chomosomes_A"="11", "pos_A"="12000000,
#'  geneB="KDM5A", chromosomes_B"="12", "pos_B"="14000000")
#'
#' my_genes_GR <- transcriptsBy(txDB, by="gene")[c("NUP98", "KDM5A")]
#'
#' forGenPaint <- genome_paint_format(df=my_fusions, fusion_genes_GR=my_genes_GR)
#' }
#' @import dplyr
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


#' Function to define the exon number for each gene partner in fusion.
#'
#' @param df data frame with the column names "Sample","chomosomes_A", "pos_A", "chromosomes_B", "pos_B"
#' @param exon_ranges_GR the GRanges object with the fusion genes exon positions. One transcript ID per gene symbol.
#' @param intron_ranges_GR the GRanges object with the fusion genes intron positions.One transcript ID per gene symbol. Must have mcol of "intron_rank".
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' library(HS.txDB)
#' library(GenomicFeatures)
#'
#' my_fusions <- data.frame("Sample"="Patient1", geneA="NUP98", chomosomes_A"="11", "pos_A"="12000000,
#'  geneB="KDM5A", chromosomes_B"="12", "pos_B"="14000000")
#'
#' my_exons_GR <- exonsBy(txDB, by="tx")[c("ENST00000399788", "ENST00000439151")]
#'
#' fusion_exons <- define_exon_junctions(df=my_fusions, fusion_genes_GR=my_exons_GR)
#' }
#' @import dplyr
define_exon_junctions <- function(df, exon_ranges_GR, intron_ranges_GR=NULL){
  grA <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(df[["chromosomes_A"]]),
                                ranges = IRanges(start=df[["pos_A"]]))

  grB <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(df[["chromosomes_B"]]),
                                ranges = IRanges(start=df[["pos_B"]]))

  #Find which NUP98-Fusion gene exon belongs with the genomic coordinate
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
  exon_infoA <- as.data.frame(S4Vectors::mcols(exon_infoA)) %>%
    tibble::rownames_to_column("gene") %>%
    rename_all(~paste0(.,"A"))

  exon_infoB <- IRanges::subsetByOverlaps(unlist(exon_infoB), grB)
  exon_infoB <- as.data.frame(S4Vectors::mcols(exon_infoB)) %>%
    tibble::rownames_to_column("gene") %>%
    rename_all(~paste0(.,"B"))

  #rowbind the metadata cols for the breakpoint pair.
  res <- dplyr::bind_cols(exon_infoA,exon_infoB) %>%
    mutate(Sample=df[["Sample"]],
           Fusion=paste(geneA,geneB, sep="-")) %>%
    unite("Exons", matches("_rank[AB]"),sep = "_",remove = F)

  return(res)
}





