

getSequence <- function(alingedSequences = allGenomes, annotationInfo = annotationInfo, sample = "crassa_BA"){
  
  contig = annotationInfo %>% 
    dplyr::filter(Sample == sample) %>%
    dplyr::select(contig)
  
  contig = contig[[1]]  
  
  start = annotationInfo %>% 
    dplyr::filter(Sample == sample) %>%
    dplyr::select(start)
  
  start = start[[1]]  
  
  
  singleSequence = alingedSequences[,c(contig,start) ]
  
  
  singleSequence = singleSequence  %>%
    dplyr::rename(Chr = contig , Position = start)
  
  return (singleSequence)
  
}

annotateSites <- function(methylationSites, annotation){
  
  genomeSites = methylationSites %>% dplyr::select(Chr,Position) %>% distinct()
  
  
  genomeSites = makeGRangesFromDataFrame(genomeSites,
                                         seqnames.field = "Chr", 
                                         start.field = "Position", 
                                         end.field = "Position")
  
  
  MethylatedIntron = findOverlaps(genomeSites,unlist(annotation$introns))
  MethylatedExons = findOverlaps(genomeSites,unlist(annotation$exons))
  MethylatedfiveUTRs = findOverlaps(genomeSites,unlist(annotation$fiveUTRs))
  MethylatedthreeUTRs = findOverlaps(genomeSites,unlist(annotation$threeUTRs))
  
  #intergenic = distanceToNearest(genomeSites,annotation$transcripts )
  
  
  #intergenicHits = data.frame(intergenic) %>% filter(distance > 0)
  #transcriptHits = data.frame(intergenic) %>% filter(distance == 0)
  
  TEProximity  =  data.frame(distanceToNearest(genomeSites,annotation$TEs ))
  TEHist =  data.frame(TEProximity) %>% filter(distance == 0)
  
  genomeSites$TEdistance =1
  genomeSites[TEProximity$queryHits]$TEdistance = TEProximity$distance
  genomeSites$transcript = 0
  genomeSites$intergenic = 0
  genomeSites$exon =  0
  genomeSites$exon[data.frame(MethylatedExons)$queryHits] = 1
  genomeSites$transcript[data.frame(MethylatedExons)$queryHits] = 1
  genomeSites$intron =  0
  genomeSites$intron[data.frame(MethylatedIntron)$queryHits] = 1
  genomeSites$transcript[data.frame(MethylatedIntron)$queryHits] = 1
  genomeSites$futr =  0
  genomeSites$futr[data.frame(MethylatedfiveUTRs)$queryHits] = 1
  genomeSites$transcript[data.frame(MethylatedfiveUTRs)$queryHits] = 1
  genomeSites$tutr =  0
  genomeSites$tutr[data.frame(MethylatedthreeUTRs)$queryHits] = 1
  genomeSites$transcript[data.frame(MethylatedthreeUTRs)$queryHits] = 1
  
  genomeSites$TE =  0
  genomeSites$TE[genomeSites$TEdistance == 0] = 1
  
  
  
  
  genomeSites.df = data.frame(genomeSites)
  
  genomeSites.df$annotation = "Intergenic"
  genomeSites.df$annotation[genomeSites.df$transcript == 1] = "CDS"
  genomeSites.df$annotation[genomeSites.df$exon == 1] = "CDS"
  genomeSites.df$annotation[genomeSites.df$intron == 1] = "Intron"
  genomeSites.df$annotation[genomeSites.df$futr == 1] = "5UTR"
  genomeSites.df$annotation[genomeSites.df$tutr == 1] = "3UTR"
  genomeSites.df$annotation[genomeSites.df$TE == 1] = "TE"
  genomeSites.df$intergenic[genomeSites.df$annotation == "Intergenic"] = 1

  
transcripts.overlap  =   findOverlaps(genomeSites,annotation$transcripts)
transcript.df = annotation$transcripts %>%
  as.data.frame()%>%
  rownames_to_column(var = "subjectHits") %>%
  mutate(subjectHits = as.integer(subjectHits)) %>%
  dplyr::select(subjectHits, width, tx_name) %>%
  inner_join(as.data.frame(transcripts.overlap)) %>%
  dplyr::select(- subjectHits) 


  genomeSites.df %>%
    dplyr::select(-width) %>%
    rownames_to_column(var = "queryHits") %>%
    mutate(queryHits = as.integer(queryHits)) %>%
    left_join(as.data.frame(transcript.df)) ->
    genomeSites.df.annotated
  
  
    
  

  genomeSites.df.annotated %>%
    dplyr::rename(Chr = "seqnames") %>%
    dplyr::rename(Position = "start") %>%
    dplyr::select(-end, -strand) ->
    genomeSites.df.annotated
  
  
  return(genomeSites.df.annotated)
  
}


getDistribution <- function(files = files ,  pattern = "crassa-BA", methylatedFileDir){
  
  
  subFiles =  files[grep(pattern = pattern, x =files )]
  samples = gsub(pattern = ".4cov.keep0meth.meth0IFnotPval.methylated", replacement = "",x = subFiles)
  
  
  samplesInfo = list()
  for(i in 1:length(subFiles)){
    file = paste(methylatedFileDir,subFiles[i], sep = "/")
    dat = read_delim(file = file,delim = " ", col_names = FALSE )
    dat$sample = samples[i]
    samplesInfo[[samples[i]]] = dat
  }
  
  
  methylationPattern = bind_rows(samplesInfo)
  methylationPattern$methylated = 0
  methylationPattern$methylated[methylationPattern$X8 < 0.05] = 1
  
  genomeSites.df.annotated
  
  
  methyaled.sites = methylationPattern %>%
    mutate(fraction = X7/X6) %>% 
    dplyr::select(X1, X2, X6, fraction, X8,methylated, sample ) %>%
    dplyr::rename(Chr = "X1" , Position = "X2", pValue = "X8",counts = "X6") %>%
    dplyr::select(Chr,Position,pValue,counts,fraction, sample,methylated)
  
  return(methyaled.sites)
}




plotMethylationDistribution <- function(methylationPattern,figureDirectory ,pattern,sampleSize = 100000){
  
  distribution2 = methylationPattern %>% 
    dplyr::select(Chr,Position,sample, fraction) %>% 
    spread(key = sample, value = fraction, fill =  0) %>% 
    dplyr::select(-Chr,-Position)
  
  
  distribution_sample = sample_n(distribution2, size = sampleSize)
  
  pairsFile = paste(figureDirectory, 
                    paste(pattern,sampleSize,"pairs.png", sep = "."),
                    sep = "/")
  pairsPlot = ggpairs(distribution_sample, lower = list(continuous = wrap("points", alpha = 0.3, size=0.1)))
  ggsave(filename = pairsFile, plot = pairsPlot)
  
}


plotPairwiseCorrelation <- function(fileDirectory, species = "crassa_BA"){
  
  genomeMethyaledSitesFile = paste(fileDirectory,
                                   paste(species,"annotated.tsv", sep = "."),
                                   sep = "/")
  methyaled.sites.annotated = read_tsv( file = genomeMethyaledSitesFile)
  
  methyaled.sites.annotated.sex = methyaled.sites.annotated %>%  
    filter(str_detect(string = sample, pattern = "sex")) %>%
    mutate(state = "sex")
  plotMethylationDistribution (figureDirectory =fileDirectory ,methylationPattern = methyaled.sites.annotated.sex, 
                               pattern = paste(species,"sex",sep = "." ))
  
  
  methyaled.sites.annotated.veg = methyaled.sites.annotated %>%  
    filter(str_detect(string = sample, pattern = "veg"))%>%
    mutate(state = "veg")
  
  plotMethylationDistribution (figureDirectory =fileDirectory ,methylationPattern = methyaled.sites.annotated.veg, 
                               pattern = paste(species,"veg",sep = "." ))
  
  
}



getAnnotation <- function(SampleName = "crassa-BA-veg", annotationInfo, annotationDir){
  
  geneFile = annotationInfo %>% 
    dplyr::filter(sample == SampleName) %>%
    dplyr::select(gene)
  
  geneFile = geneFile$gene[[1]]  
  TEFile = annotationInfo %>% 
    dplyr::filter(sample == SampleName) %>%
    dplyr::select(TE)
  
  TEFile = TEFile$TE[[1]]
  
  gtfFile <- paste(annotationDir,geneFile, sep = "/")
  
  txdb <- makeTxDbFromGFF(file=gtfFile, dataSource="ensemblgenomes",format = "gff3")
  
  TEFile <- paste(annotationDir,TEFile, sep = "/")
  
  TEdf =  read.table(file=TEFile,  sep = "\t" ,header =  FALSE) 
  colnames(TEdf) = c("seqid", "TE_start","TE_stop","TE_score","TE_ID", "TE_type","TE_ID2")
  
  TEs = makeGRangesFromDataFrame(TEdf,
                                 seqnames.field = "seqid", 
                                 start.field = "TE_start", 
                                 end.field = "TE_stop")
  
  TEs = reduce(TEs)
  
  
  
  exons = exonsBy(x =txdb,  by = c("tx", "gene"))
  introns = intronsByTranscript(x = txdb)
  fiveUTRs = fiveUTRsByTranscript(x = txdb)
  threeUTRs = threeUTRsByTranscript(x = txdb)
  transcripts = transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"))
  
  genomeSize = nrow(data.frame(unlist(coverage(unlist(exons)))))
  exonsSize = sum(data.frame(unlist(coverage(unlist(exons))))[1] >0 )
  intronsSize = sum(data.frame(unlist(coverage(unlist(introns))))[1] >0 )
  threeUTRsSize = sum(data.frame(unlist(coverage(unlist(threeUTRs))))[1] >0 )
  fiveUTRsSize = sum(data.frame(unlist(coverage(unlist(fiveUTRs))))[1] >0 )
  TEsize = sum(data.frame(unlist(coverage(TEs)))[1] >0 )
  
  intergenicSize = genomeSize - TEsize - exonsSize - intronsSize - threeUTRsSize - fiveUTRsSize
  CDSSize = exonsSize  - threeUTRsSize - fiveUTRsSize
  
  GenomeAnnotation = data.frame(annotation = c("Intergenic","CDS","Intron","5UTR","3UTR","TE"), 
                                nrOfNucleoties = c(intergenicSize,CDSSize,intronsSize,fiveUTRsSize,threeUTRsSize,TEsize))
  
  list(threeUTRs)
  
  annotations = list()
  annotations[["introns"]] = introns
  annotations[["exons"]] = exons
  annotations[["threeUTRs"]] = threeUTRs
  annotations[["fiveUTRs"]] = fiveUTRs
  annotations[["TEs"]] = TEs
  annotations[["transcripts"]] = transcripts
  annotations[["distribution"]] = GenomeAnnotation
  annotations[["species"]] = sample
  
  return(annotations)
  
}


  