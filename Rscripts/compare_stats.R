library(ggplot2)
library(tools)

plotcontig <- function(x) {
  basename = file_path_sans_ext(basename(x))
  t <- read.table(x, header=TRUE) # load file
  # t$RIBOSOMAL != 1 & 
  tlong <- subset(t,t$LENGTH >= 500 & t$MT != 1 & COVERAGE < 1000)
  p <- ggplot(data=tlong,aes(x=GC,y=COVERAGE,
                             color=factor(CAPSID),
                             shape=factor(floor(ORF_MEAN)),
                             size=LENGTH)) + 
    geom_point() + scale_fill_brewer() + ggtitle(basename) 
  ggsave(filename=sprintf("genome_plots/%s.pdf",basename),p,width=15)
}

files <- list.files(path="results/genome_stats_nofilter", pattern="*.sum_stat.tsv", 
                    full.names=TRUE, recursive=FALSE)

lapply(files, plotcontig)

files <- list.files(path="results/genome_stats_ref", pattern="*.sum_stat.tsv", 
                    full.names=TRUE, recursive=FALSE)

lapply(files, plotcontig)
