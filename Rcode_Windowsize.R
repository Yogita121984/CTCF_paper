start = seq(1,seqlengths(Hsapiens)["chr1"],window_size)
end = start -1
end = end[-1]
start = start[-length(start)]
length(start)
498501
length(end)
[1] 498501
all_signals <- data.frame( start, end)
nrow(all_signals)
all_signals["chr"] <- chr1
all_signals <- all_signals[c("chr", "start", "end")]
write.table( all_signals, file="chr1_windowsize2", col.names=T, row.names=F, quote=F, sep="\t")


#########

 for( i in 1:length( chrs ))
 {
 start = seq(1,seqlengths(Hsapiens)[chrs[i]],window_size)
 end = start -1
 end = end[-1]
 start = start[-length(start)]
 all_signals <- data.frame( start, end)
 all_signals["chr"] <- chrs[i]
 all_signals <- all_signals[c("chr", "start", "end")]
 }





#############
library (BSgenome)
library (BSgenome.Hsapiens.UCSC.hg19)
window_size=9000
chrs = paste("chr", c("1","10","11","12","13","14","15","16","17","18","19","2","20","21","22","3","4","5","6","7","8","9","X"), sep = "")
try<-data.frame()
for( i in 1:length( chrs ))
 {
 start = seq(1,seqlengths(Hsapiens)[chrs[i]],window_size)
 end = start -1
 end = end[-1]
 start = start[-length(start)]
 all_signals <- data.frame( start, end)
 all_signals["chr"] <- chrs[i]
 all_signals <- all_signals[c("chr", "start", "end")]
 try <-rbind(try, all_signals)
 }
write.table( try, file="hg19_windowsize", col.names=T, row.names=F, quote=F, sep="\t")
