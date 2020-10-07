#Jenny Smith
#Sept. 17, 2020
#purpose: examine junction sequences in RNA-seq reads.


###### Subset BAMS
#' Title
#'
#' @param Sample_ID character string of patient ID
#' @param breakpoint list class with breakpoint string
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
    message("Breakpoint and fusion data cannot both be NULL. Please provide either a breakpoint in list
    with the format of chr:basepairs|chr:basepairs, for example list(CBFB.MYH11=16:15814908|16:67116211). List can be any length.Or a dataframe with the desired fusion data.")
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
    breakpoint <- stringr::str_split(breakpoint,"\\|")

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
      #this is sepecific to GRCh37-lite. Will need to make more generalized ASAP.
      pattern <- paste0("chromosome ",chr,", GRCh37 primary reference assembly") #for full annotated chrs
      pos <- IRanges::IRanges(start=0, end=Seqlens[grep(pattern,names(Seqlens))])
    }
    return(pos)
  }

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
#' Title
#'
#' @param theoretical_seqs_df a dataframe with computationally derived exon junction sequences
#'
#' @return
#' @export
#'
#' @examples
#'my_seqs <- data.frame(Name=paste0("seq",1:10),Sequence=rep(c("ATCGCCCGTTA"), 10))
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
      pull(Sequence, name=Name) %>%
      Biostrings::DNAStringSet()
    #Add reverse compliment to search
    seqExons <- c(seqExons,  magrittr::set_names(Biostrings::reverseComplement(seqExons),
                                                 paste0("revComp",names(seqExons))))
    seqExons <- seqExons[order(names(seqExons))]
    seqExons <- Biostrings::PDict(seqExons)  #(3) later matchPdict can only be used with max.mismatch=0.???
    seqExonsList[[i]] <- seqExons

    rm(seqExons)
  }

  message(paste("Finished PDict processing with",sum(sapply(seqExonsList, length)),"theoretical fusion sequences."))


  return(seqExonsList)
}




######### Sequence Pattern Matching

#' Title
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
      spread(Contig.Name,Contig.Match) %>%
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
#' Title
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



  #Begin subsetting the RNAseq BAM file reads
  BAM <- subset_bam(Sample_ID = Sample_ID,
                    breakpoint=breakpoint,
                    fusion_data = fusion_data,
                    file_manifest = file_manifest,
                    scan = TRUE,
                    by_chr=TRUE,
                    Seqlens=Seqlens)

  #Select only the RNAseq reads to keep in memory
  #This is based on a single fusion -- need to check this with interchromosomal fusions
  RNAseq <- BAM[[1]]$seq
  names(RNAseq) <- BAM[[1]]$qname
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
