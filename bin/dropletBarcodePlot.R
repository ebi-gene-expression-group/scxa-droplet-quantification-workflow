#!/usr/bin/env Rscript

library(ggplot2)
library(DropletUtils)
library(gridExtra)

cl <- commandArgs(trailingOnly = TRUE)

barcode_file <- cl[1]
result_matrix <- cl[2]
label <- cl[3]
output_plot_file <- cl[4]

if ( ! file.exists(barcode_file) ){
    write("Input barcodes file name does not exist", stderr())
    q(status = 1) 
}

if ( ! dir.exists(result_matrix) ){
    write("Input matrix directory does not exist", stderr())
    q(status = 1) 
}


hist_bins <- 50

# Pick a cutoff on count as per https://github.com/COMBINE-lab/salmon/issues/362#issuecomment-490160480

pick_roryk_cutoff = function(bcs){
  bcs_hist = hist(log10(bcs), plot=FALSE, n=hist_bins)
  mids = bcs_hist$mids
  vals = bcs_hist$count
  wdensity = vals * (10^mids) / sum(vals * (10^mids))
  baseline <- median(wdensity)
  
  # Find highest density in upper half of barcode distribution
  
  peak <- which(wdensity == max(wdensity[((length(wdensity)+1)/2):length(wdensity)]))
  
  # Cutoff is the point before the peak at which density falls below 2X baseline
  
  10^mids[max(which(wdensity[1:peak] < (1.5*baseline)))]
}

# Plot densities 

barcode_density_plot = function(bcs, roryk_cutoff, knee, inflection, name = 'no name') {
  bcs_hist = hist(log10(bcs), plot=FALSE, n=hist_bins)
  counts = bcs_hist$count
  mids = bcs_hist$mids
  y = counts * (10^mids) / sum(counts * (10^mids))
  qplot(y, 10^mids) + geom_point() + theme_bw() + ggtitle(name) + ylab('Count') + xlab ('Density') +
    geom_hline(aes(yintercept = roryk_cutoff, color = paste('roryk_cutoff =', length(which(bcs > roryk_cutoff)), 'cells'))) + 
    geom_hline(aes(yintercept = inflection, color = paste('dropletutils_inflection =', length(which(bcs > inflection)), 'cells'))) +
    geom_hline(aes(yintercept = knee, color = paste('dropletutils_knee =', length(which(bcs > knee)), 'cells'))) +
    scale_y_continuous(trans='log10') + theme(axis.title.y=element_blank()) + labs(color='Thresholds')
}  

# Plot a more standard barcode rank plot

barcode_rank_plot <- function(br.out, roryk_total_cutoff, knee, inflection, name='no name'){
  ggplot(data.frame(br.out), aes(x=rank, y=total)) + geom_line() + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10') + theme_bw() + 
    geom_hline(aes(yintercept = knee, color = 'dropletutils_knee')) + 
    geom_hline(aes(yintercept = inflection, color = 'dropletutils_inflection')) +
    geom_hline(aes(yintercept = roryk_total_cutoff, color = 'roryk_cutoff')) +
    ggtitle(name) + ylab('Count') + xlab('Rank') + theme(legend.position = "none")
}

# Plot the different plots and threshold statistics alongside one another

raw_barcodes <- read.delim(barcode_file, header = FALSE)
result_matrix <- read10xCounts(result_matrix)  
processed_barcode_counts <- data.frame(V1 = colData(result_matrix)$Barcode, V2=colSums(assays(result_matrix)[[1]]))
processed_barcode_counts <- processed_barcode_counts[order(processed_barcode_counts$V2, decreasing = TRUE), ]

barcode_results <- list(
    raw = raw_barcodes,
    processed = processed_barcode_counts
)

plots <- lapply(names(barcode_results), function(name){
  
  barcodes <- barcode_results[[name]]  
  
  # Get the roryk cutoff
  roryk_count_cutoff <- pick_roryk_cutoff(barcodes$V2)
  
  # Run dropletUtils' barcodeRanks to get knee etc
  br.out <- barcodeRanks(t(barcodes[,2,drop=FALSE]))
  
  dropletutils_knee <- metadata(br.out)$knee
  dropletutils_inflection <- metadata(br.out)$inflection
  
  list(
    dropletutils = barcode_rank_plot(br.out, roryk_count_cutoff, dropletutils_knee, dropletutils_inflection, name = paste(label, name)),
    roryk = barcode_density_plot(barcodes$V2, roryk_count_cutoff, dropletutils_knee, dropletutils_inflection, name = paste(label, name))
  )
})
names(plots) <- names(barcode_results)

png(width = 1000, height = 800, file=output_plot_file)
grid.arrange(plots$raw$dropletutils, plots$raw$roryk, plots$processed$dropletutils, plots$processed$roryk, nrow=2)
dev.off()
