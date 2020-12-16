#set working directory
setwd("/Users/adrien/Box/Adrien_simulations/Paper1_Adrien/lastBam/Dog/")

#list all sam files
ll <- list.files(pattern=".*sam")

#loop over sam files
for(l in ll){
   start.time <- Sys.time()
   print(l)
   print(paste("Calculating read mapping metrics for", l))
   #count number of lines to skip
   system(paste("grep '^@'", l, "| wc > skip.line.txt"))
   zz <- strsplit(readLines("skip.line.txt"), " ")[[1]]
   skip.lines <- as.numeric(zz[zz!=""][1])
   ff <- fread(l, skip=skip.lines, head=F,fill=TRUE, na.strings="")

   #extract info
   dt <- data.table(UnqID=1:nrow(ff)) # unique ID to allow for multiply mapped reads
   dt[, Read:=ff$V1]
   dt[, True.Chr:=gsub(":0|:1", "", ff$V3)]
   full.string <- stri_split(ff$V4, fixed="\t")
   dt[, Mapped.Chr:=sapply(full.string, "[", 3)]
   #ss <- stri_split(sapply(full.string, "[", 1), fixed=":")
   #dt[, True.Pos:=as.numeric(sapply(ss, "[", 2))]
   #dt[, Mapped.Pos:=as.numeric(sapply(full.string, "[", 4))]
   dt[, MAPQ:=ff$V5]
   #dt[, MAPQ:=as.numeric(sapply(full.string, "[", 5))]
   #dt[, Strand:=sapply(ss, "[", 1)]
   #dt[which(Strand==0), Mapped.Pos:=as.numeric(Mapped.Pos-5)]
   #dt[which(Strand==1), Mapped.Pos:=as.numeric(Mapped.Pos+5)]
   #dt[, Cigar:=sapply(full.string, "[", 6)]
   #dt[, Seq:=sapply(full.string, "[", 10)]
   #dt[, Read.Length:=nchar(sapply(full.string, "[", 10))]
   
   #add unique identifier (for multi-mapped reads)
   
   
   #clean up environment
   #rm(list=ls()[!ls() %in% c("dt", "ll", "l", "start.time")])
   
   #estimate distance between true and mapped read position
   #dt[, Clip.Count:=stri_count(Cigar, regex="S|H")] # count number of clips
   
   #read on forward strand with no clipping
   #dt[Strand==0 & Mapped.Chr==True.Chr & Clip.Count==0, 
   #   Mapped.Dist:=Mapped.Pos-True.Pos]
   
   #read on reverse strand with no clipping
   #dt[Strand==1 & Mapped.Chr==True.Chr & Clip.Count==0, 
   #   Mapped.Dist:=Mapped.Pos-(250-Read.Length)-True.Pos]
   
   #if(dt[, any(Clip.Count>0)]){
      #reads on foward strand with clipping
   #   dt.for <- dt[Strand==0 & Mapped.Chr==True.Chr & Clip.Count>0]
   #   rgx <- sapply(stri_split(dt.for[, Cigar], 
    #                           regex="S|H"), "[", 1)
    #  spp <- rgx %in% 1:1000
     # pos.adj <- as.numeric(rgx[spp]) # adjustment due to clipping
      
   #   dt[UnqID %in% dt.for[spp, UnqID], 
    #     Mapped.Dist:=Mapped.Pos-pos.adj-True.Pos] # clipped from start of read
     # dt[UnqID %in% dt.for[!spp, UnqID], 
     #    Mapped.Dist:=Mapped.Pos-True.Pos] # clipped from end of read
      
      #reads on reverse strand with clipping
      #dt.rev <- dt[Strand==1 & Mapped.Chr==True.Chr & Clip.Count>0]
      #rgx <- sapply(stri_split(dt.rev[, Cigar], 
         #                      regex="S|H"), "[", 1)
      
      #spp <- rgx %in% 1:1000
      #pos.adj <- as.numeric(rgx[spp]) # adjustment due to clipping
      
      #dt[UnqID %in% dt.rev[spp, UnqID], 
      #   Mapped.Dist:=Mapped.Pos-(250-Read.Length)-pos.adj-True.Pos] # clipped from start of read
      #dt[UnqID %in% dt.rev[!spp, UnqID], 
      #   Mapped.Dist:=Mapped.Pos-(250-Read.Length)-True.Pos] # clipped from end of read
   #}
   
   #output file
   #out <- dt[order(True.Chr, True.Pos)]
   out <- dt[order(True.Chr)]
   out[, UnqID:=1:.N]
   out.name <- paste0("MapDist_", gsub(".sam", ".txt", gsub("FQ_Qualities_", "", l)))
   fwrite(out, out.name, row.names=F, col.names=T, sep="\t", quote=F, na="NA")
   
   end.time <- Sys.time()
   print(paste("Run time:", round(end.time-start.time, 2)))
}
  
out
