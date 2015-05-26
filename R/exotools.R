#FUNCTIONS FOR EXOTOOLS WORKFLOW DEFINED HERE

#BAM IMPORT FOR GVIZ - modified from https://gist.github.com/sidderb/e485c634b386c115b2ef

BamImportMouse = function (selection) {
  require(Rsamtools)
  Input1=mouse@filename[1]
  ChIP1=mouse@filename[2]
  Input2=mouse@filename[3]
  ChIP2=mouse@filename[4]
  
  #scan ChIP 
  param = ScanBamParam(what = c("pos", "qwidth", "strand"), which = selection,
                       flag = scanBamFlag(isUnmappedQuery = FALSE))
  i1 = scanBam(Input1, param = param)[[1]]
  c1 = scanBam(ChIP1, param = param)[[1]]
  i2 = scanBam(Input2, param = param)[[1]]
  c2 = scanBam(ChIP2, param = param)[[1]]
  
  
  #Replicate 1
  i1GR = GRanges(strand=i1[["strand"]], ranges=IRanges(i1[["pos"]],
                                                       width = i1[["qwidth"]]), seqnames=seqnames(selection)[1])
  c1GR = GRanges(strand=c1[["strand"]], ranges=IRanges(c1[["pos"]],
                                                       width = c1[["qwidth"]]), seqnames=seqnames(selection)[1])
  #Replicate 2
  i2GR = GRanges(strand=i2[["strand"]], ranges=IRanges(i2[["pos"]],
                                                       width = i2[["qwidth"]]), seqnames=seqnames(selection)[1])
  c2GR = GRanges(strand=c2[["strand"]], ranges=IRanges(c2[["pos"]],
                                                       width = c2[["qwidth"]]), seqnames=seqnames(selection)[1])
  
  # calc coverage on both strands:
  cov_list = list("Input1" = coverage(ranges(i1GR), width=end(selection)),
                  "ChIP1" = coverage(ranges(c1GR), width=end(selection)),
                  "Input2" = coverage(ranges(i2GR), width=end(selection)),
                  "ChIP2" = coverage(ranges(c2GR), width=end(selection))) 
  
  pos = sort(unique(unlist(lapply(cov_list, function(y) c(start(y), end(y))))))
  
  # build final GR
  cov_gr = GRanges(seqnames = seqnames(selection)[1], 
                   ranges=IRanges(start=head(pos, -1),end=tail(pos, -1)),
                   Input1=as.numeric(cov_list[["Input1"]][head(pos, -1)])*10e6/mouse@count[1],
                   ChIP1=as.numeric(cov_list[["ChIP1"]][head(pos, -1)])*10e6/mouse@count[2],
                   Input2=as.numeric(cov_list[["Input2"]][head(pos, -1)])*10e6/mouse@count[3],
                   ChIP2=as.numeric(cov_list[["ChIP2"]][head(pos, -1)])*10e6/mouse@count[4]
  )
  return(cov_gr)
}

#BAMIMPORT HUMAN
BamImportHuman = function (selection) {
  require(Rsamtools)
  Input1=human@filename[1]
  ChIP1=human@filename[2]
  
  #scan ChIP 
  param = ScanBamParam(what = c("pos", "qwidth", "strand"), which = selection,
                       flag = scanBamFlag(isUnmappedQuery = FALSE))
  i1 = scanBam(Input1, param = param)[[1]]
  c1 = scanBam(ChIP1, param = param)[[1]]
  
  
  #Replicate 1
  i1GR = GRanges(strand=i1[["strand"]], ranges=IRanges(i1[["pos"]],
                                                       width = i1[["qwidth"]]), seqnames=seqnames(selection)[1])
  c1GR = GRanges(strand=c1[["strand"]], ranges=IRanges(c1[["pos"]],
                                                       width = c1[["qwidth"]]), seqnames=seqnames(selection)[1])
  
  # calc coverage on both strands:
  cov_list = list("Input1" = coverage(ranges(i1GR), width=end(selection)),
                  "ChIP1" = coverage(ranges(c1GR), width=end(selection)))
  
  pos = sort(unique(unlist(lapply(cov_list, function(y) c(start(y), end(y))))))
  
  # build final GR
  cov_gr = GRanges(seqnames = seqnames(selection)[1], 
                   ranges=IRanges(start=head(pos, -1),end=tail(pos, -1)),
                   Input1=as.numeric(cov_list[["Input1"]][head(pos, -1)])*10e6/human@count[1],
                   ChIP1=as.numeric(cov_list[["ChIP1"]][head(pos, -1)])*10e6/human@count[2]
  )
  return(cov_gr)
}



#plot GVIZ
plotGviz <- function(mouse_gr) {
  if (class(mouse_gr)== "character") {
    mouse_gr<-gsub(",","",mouse_gr)
    mseq<-sapply(strsplit(mouse_gr,":"),function(x) x[1])
    temp<-sapply(strsplit(mouse_gr,":"),function(x) x[2])
    mstart<-as.numeric(sapply(strsplit(temp,"-"),function(x) x[1]))
    mend<-as.numeric(sapply(strsplit(temp,"-"),function(x) x[2]))
    mgr<-GRanges(seq=mseq,IRanges(start=mstart,end=mend))
  } else {
    mgr<-mouse_gr
    mseq<-as.character(seqnames(mgr))
    mstart<-start(mgr)
    mend<-end(mgr)
  }  
  print(mgr)
  mChIP<-BamImportMouse(mgr)
  mTxTrack<-GeneRegionTrack(mouse@tx, chromosome = mseq, start = mstart, end = mend,name="UCSC")
  #mIdeoTrack <- IdeogramTrack(genome = "mm9", chromosome = mseq) 
  mChIP_Track<-DataTrack(mChIP,groups=c("Input1","ChIP1","Input2","ChIP2"),
                         legend=TRUE,genome="mm9",name="Mouse ChIP-Seq (cpm)",
                         chromosome=mseq, col=c("blue","purple","green","red"))
  
  #human
  suppressWarnings(temp<-unlist(liftOver(mgr,m2h_chain)))
  #subset temp to a single chromosome if more than one
  #temp<-keepSeqlevels(temp,names(which.max(table(seqnames(temp)))))
  #reduce temp to gene-size chunks
  temp<-reduce(temp,min.gapwidth=10000)
  #if more than one pick one with similar size to mgr
  temp<-temp[which.min(abs(IRanges::width(mgr)-IRanges::width(temp)))]
  hseq<-as.character(seqnames(temp))
  hstart<-min(start(temp))
  hend<-max(end(temp))
  hgr<-GRanges(seq=hseq,IRanges(start=hstart,end=hend))
  print(hgr)
  hChIP<-BamImportHuman(hgr)
  hTxTrack<-GeneRegionTrack(human@tx, chromosome = hseq, start = hstart, end = hend,name="UCSC")
  #hIdeoTrack <- IdeogramTrack(genome = "hg19", chromosome = hseq) 
  hChIP_Track<-DataTrack(hChIP,groups=c("Input","ChIP"),
                         legend=TRUE,genome="hg19",name="Human ChIP-Seq (cpm)",
                         chromosome=hseq, col=c("blue","green"))
  
  grid.newpage()
  ncols<-2
  nrows<-1
  pushViewport(viewport(layout = grid.layout(nrows,ncols)))
  pushViewport(viewport(layout.pos.col = 1 , layout.pos.row=1))
  plotTracks(list(mChIP_Track,mTxTrack,AxisTrack), from = mstart, to = mend, type="a", col.histogram=NA, 
             scale = 0.25, labelPos = "below",
             cex.title=1, cex.axis=1, title.width=0.4,add=TRUE)
  popViewport(1)
  
  pushViewport(viewport(layout.pos.col = 2 , layout.pos.row=1))
  plotTracks(list(hChIP_Track,hTxTrack,AxisTrack), from = hstart, to = hend, type="a", col.histogram=NA, 
             scale = 0.25, labelPos = "below",
             cex.title=1, cex.axis=1, title.width=0.4,add=TRUE)
  popViewport(1)
  
}


extractDNAfromGR <- function(gr,build=Mmusculus) {
  #initialize the column
  gr$dna<-DNAStringSet("NNNN")
  
  # look through Chromosomes
  for (i in as.character(levels(seqnames(gr)))) {
    idx<-as.logical(seqnames(gr)==i)
    gr[idx]$dna <- DNAStringSet(build[[i]],start=start(gr[idx]),end=end(gr[idx]))
  }
  return(gr)
}




pileupGR<- function(gr,fl) {
  
  size=width(gr[1])
  xlim<-c((size %/% 2)*-1,size %/% 2)
 
  
  #Widen gr by ncycle to get data for the 5' end of the bottom strand
  gr<-gr+ncycle
  sbp<-Rsamtools::ScanBamParam(which=gr)
  res<-Rsamtools::pileup(fl, scanBamParam=sbp, pileupParam=p_param)
  res[res$strand=="-","pos"]<-res[res$strand=="-","pos"]+ncycle-1
  
  #create labels from gr to match results table
  gr_labels<-paste0(seqnames(gr),":",start(gr),"-",end(gr))
  idx<-match(sapply(strsplit(as.character(res$which_label),"\\."),function(x) x[1]),gr_labels)
  #idx<-match(blah,gr_labels)
  
  #add start and strand info to results table
  res$site_strand<-as.character(strand(gr[idx]))
  res$start<-as.numeric(start(gr[idx]))
  
  res$rpos<-999
  plus_idx<-rownames(res[res$site_strand!="-",])
  res[plus_idx,]$rpos<-(res[plus_idx,]$pos-res[plus_idx,]$start-width(gr[1])%/%2)*-1
  
  minus_idx<-rownames(res[res$site_strand=="-",])
  res[minus_idx,]$rpos<-res[minus_idx,]$pos-res[minus_idx,]$start-width(gr[1])%/%2
  
  
  #assign each count to top or bottom strand
  res$rstrand<-"unknown_strand"
  if ("*" %in% res$site_strand | "+" %in% res$site_strand) {
    res[res$strand=="+" & res$site_strand!="-",]$rstrand<-"top"
    res[res$strand=="-" & res$site_strand!="-",]$rstrand<-"bottom"
  }
  
  if ("-" %in% res$site_strand) {
    res[res$strand=="-" & res$site_strand=="-",]$rstrand<-"top"
    res[res$strand=="+" & res$site_strand=="-",]$rstrand<-"bottom"
    
  }
  
  #use tapply to summarize
  data<-as.data.frame(tapply(res$count,list(res$rpos,res$rstrand),sum))
  data$pos<-rownames(data)
  
  #tidyr
  temp <- data %>% 
    gather(strand,count,top,bottom,na.rm=TRUE,convert=TRUE) %>%
    mutate(pos=as.integer(pos)*-1) %>%
    dplyr::filter(pos >= xlim[1] & pos <= xlim[2]) 
  
  temp[temp$strand=="bottom",]$count<-temp[temp$strand=="bottom",]$count *-1
  y<-max(abs(temp$count))*1.2
  p<-ggplot(temp, aes(x = pos, y = count, group=strand, fill=strand)) +
    #geom_point() + geom_line()+
    geom_area(position="dodge") + theme_classic() + xlim(xlim) +ylim(c(-y,y))
  
  return(p)
  
}


findsites <- function (pattern,seq,build=Mmusculus) {
  if (class(pattern)=="character") {pattern<-DNAString(pattern)}
  if (pattern == reverseComplement(DNAString(pattern))) {
    sites<-GRanges()
    temp<-matchPattern(DNAString(pattern),build[[seq]],fixed="subject")
    sites<-GRanges(seqnames=seq,ranges=IRanges(start=start(temp),
                                               end=end(temp)),strand="*") 
    return(sites)
  } else {
    top_sites<-GRanges()
    bot_sites<-GRanges()
    temp<-matchPattern(DNAString(pattern),build[[seq]],fixed="subject")
    temp2<-matchPattern(reverseComplement(DNAString(pattern)),build[[seq]],fixed="subject")
    if (length(temp) > 0) { 
      top_sites<-GRanges(seqnames=seq,ranges=IRanges(start=start(temp),
                                                     end=end(temp)),strand="+") } 
    if (length(temp2) > 0) { 
      bot_sites<-GRanges(seqnames=seq,ranges=IRanges(start=start(temp2),
                                                     end=end(temp2)),strand="-") } 
    return(suppressWarnings(c(top_sites,bot_sites)))
  }
}



