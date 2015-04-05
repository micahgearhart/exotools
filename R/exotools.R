p_param <- PileupParam(cycle_bins=c(0,1),
                       distinguish_nucleotides=FALSE,
                       distinguish_strands=TRUE,
                       min_nucleotide_depth=1)

pileupGR<- function(gr,fl) {
  
  size=width(gr[1])
  xlim<-c((size %/% 2)*-1,size %/% 2)
  
  #Widen gr by 50 to get data for the 5' end of hte bottom strand
  gr<-gr+50
  sbp<-ScanBamParam(which=gr)
  res<-pileup(fl, scanBamParam=sbp, pileupParam=p_param)
  res[res$strand=="-","pos"]<-res[res$strand=="-","pos"]+49
  
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
  
  # implement dplyer version of tapply
  #res %>% group_by(rpos,rstrand) %>% 
  #    summarize(count=sum(count)) %>%
  #    gather(rstrand,count,top,bottom,na.rm=TRUE,convert=TRUE) %>%
  #    mutate(rpos=as.integer(rpos)*-1) %>%
  #    mutate(count[strand=="bottom"]=count*-1) %>%
  #    dplyr::filter(rpos >= xlim[1] & rpos <= xlim[2]) %>%
  
  
  #ggplot - stacked
  #p<-ggplot(temp, aes(x = pos, y = count, colour = strand, group = strand)) +
  #geom_line() + geom_point() + theme_bw() + xlim(xlim)
  
  #ggplot - flipped  & filled 
  
  temp[temp$strand=="bottom",]$count<-temp[temp$strand=="bottom",]$count *-1
  y<-max(abs(temp$count))*1.2
  p<-ggplot(temp, aes(x = pos, y = count, group=strand, fill=strand)) +
    geom_point() + geom_line()+
    geom_area(position="dodge") + theme_bw() + xlim(xlim) +ylim(c(-y,y))
  
  return(p)
  
}


findsites <- function (pattern,seq) {
  top_sites<-GRanges()
  bot_sites<-GRanges()
  temp<-matchPattern(DNAString(pattern),Mmusculus[[seq]],fixed="subject")
  temp2<-matchPattern(reverseComplement(DNAString(pattern)),Mmusculus[[seq]],fixed="subject")
#  print(paste0(seq," Plus strand:  ",length(temp)))
#  print(paste0(seq," Minus strand:  ",length(temp2)))
  if (length(temp) > 0) { top_sites<-GRanges(seqnames=seq,ranges=IRanges(start=start(temp),end=end(temp)),strand="+") } 
  if (length(temp2) > 0) { bot_sites<-GRanges(seqnames=seq,ranges=IRanges(start=start(temp2),end=end(temp2)),strand="-") } 
  return(suppressWarnings(c(top_sites,bot_sites)))
}



