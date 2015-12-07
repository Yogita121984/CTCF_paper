slideNo <- "H3k4me1" 
signalname <- "H3K4Me1_H1hesc""

distances <- as.integer( c( 10000,1000000,2000000 ))
graph_name <- paste( paste("Fig6", "H3k4me1", "H1hesc", sep="_"),"pdf", sep=".")
pdf( graph_name )
par(mfrow=c(3,1))

for( i in 1:length( distances ) )
{
   
   inputFileName <- paste( slideNo, "TotalGenes1000Samples", distances[i], sep="_" )
   chip_data <- read.table(file=inputFileName, header=T)
   infile2 <- paste(slideNo, "3Signals", sep="_")
   data_genes <- read.table( file=infile2, header=T)
      tssDistances <- c( "10 kb", "1 Mb", "2 Mb")
      mailTitile=paste(signalname, "counts - TSS ±",tssDistances[i],sep=" "  )
      lowerCR  <- quantile(log2(chip_data[ , 1]), c(0.0275, 0.975))[[1]]
      upperCR  <- quantile(log2(chip_data[ , 1]), c(0.0275, 0.975))[[2]]
      hist( log2(chip_data[, 1]), xlim=c( log2(c( min(chip_data[, 1]), (data_genes[i,1] )))), 
      main=mailTitile, xlab="")
      arrows(lowerCR, 0, lowerCR, 10 - 0.015, length = 0, col="black", lwd=5)
      arrows(upperCR, 0, upperCR, 10 - 0.015, length = 0, col="black", lwd=5)
      #for target genes
      
      text(log2(data_genes[i,1]),150, "P < 0.001",pos=2)
}
legend("topright", c("Critical Region", "CTCF binding region"), lty=c(1,1), 
col=c("black", "red"))
dev.off()
