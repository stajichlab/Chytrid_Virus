library(ggplot2)
library(tools)
library(parallel)
library(gridExtra)
library(ggpubr)
#t <- read.table("results/genome_stats/Allomyces_arbuscula_Burma_1F.sum_stat.tsv", header=TRUE) # load file
#tlong <- subset(t,t$LENGTH >= 500 & t$MT == 0 & COVERAGE < 1000 & t$ORF_COUNT > 10)
#tlong$CAPSID <- factor(tlong$CAPSID)
#tlong$NCVOG_PERCENT <- factor( floor(tlong$NCVOG_PERCENT / 10) )
#ggplot(data=tlong,aes(x=GC,y=INTERGENIC_MEDIAN,
#                           shape=CAPSID,
#                           size = LENGTH,
#                           color=NCVOG_PERCENT)) +
#  geom_point() + scale_colour_brewer(palette="Set1") + scale_fill_brewer(palette="Set1")

plotcontig_intergenic <- function(datfile) {
  basename = file_path_sans_ext(basename(datfile))
  basename = gsub("\\.sum_stat","",basename)
  t <- read.table(datfile, header=TRUE) # load file
  # t$RIBOSOMAL != 1 &
  tlong <- subset(t,t$LENGTH >= 500 & t$MT == 0 & COVERAGE < 1000 & t$ORF_COUNT > 10)
  tlong$CAPSID <- factor(tlong$CAPSID)
  tlong$NCVOG_PERCENT <- factor( floor(tlong$NCVOG_PERCENT / 10) )
  p <- ggplot(data=tlong,aes(x=GC,y=INTERGENIC_MEDIAN,
                             shape=CAPSID,
                             size = LENGTH,
                             color=NCVOG_PERCENT)) +
    geom_point() + scale_colour_brewer(palette="Set1") + scale_fill_brewer(palette="Set1") +
    ggtitle(basename) + theme_bw()
  
    ggsave(filename=sprintf("genome_plots/per-strain/%s.intergenic.pdf",basename),p,width=8)
 p
}

plotcontig_pfam <- function(datfile) {
  basename = file_path_sans_ext(basename(datfile))
  basename = gsub("\\.sum_stat","",basename)
  t <- read.table(datfile, header=TRUE) # load file
  # t$RIBOSOMAL != 1 &
  tlong <- subset(t,t$LENGTH >= 500 & t$MT == 0 & COVERAGE < 1000 & t$ORF_COUNT > 10)
  tlong$CAPSID <- factor(tlong$CAPSID)
  tlong$NCVOG_PERCENT <- factor( floor(tlong$NCVOG_PERCENT / 10) )
  ggplot(data=tlong,aes(x=GC,y=PFAM_PERCENT,
                             shape=CAPSID,
                             color=NCVOG_PERCENT,
                             size=ORF_MEDIAN)) +
    geom_point() + scale_fill_brewer(palette="Set1") + ggtitle(basename) + theme_bw()
}

plotcontig_orfsize <- function(datfile) {
  basename = file_path_sans_ext(basename(datfile))
  basename = gsub("\\.sum_stat","",basename)
  t <- read.table(datfile, header=TRUE) # load file
  # t$RIBOSOMAL != 1 &
  tlong <- subset(t,t$LENGTH >= 500 & t$MT == 0 & COVERAGE < 1000 & t$ORF_COUNT > 10)
  tlong$CAPSID <- factor(tlong$CAPSID)
  tlong$NCVOG_PERCENT <- factor( floor(tlong$NCVOG_PERCENT / 10) )
  p <- ggplot(data=tlong,aes(x=GC,y=ORF_MEDIAN,
                             shape=CAPSID,
                             size = LENGTH,
                             color=NCVOG_PERCENT)) +
    geom_point() + scale_colour_brewer(palette="Set1") + scale_fill_brewer(palette="Set1") +
    ggtitle(basename) + theme_bw()
  p
}

files <- list.files(path="results/genome_stats", pattern="*.sum_stat.tsv",
                    full.names=TRUE, recursive=FALSE)

plots <- mclapply(files, plotcontig_intergenic)
outplotfile <- "genome_plots/cleanasm_intergenic.pdf"
multi.page <-  ggpubr::ggarrange(plotlist = plots, nrow=3,ncol=3)
ggpubr::ggexport(multi.page, filename = outplotfile,width=20,height=15)

plots <- mclapply(files, plotcontig_pfam)
outplotfile <- "genome_plots/cleanasm_pfamstat.pdf"
multi.page <-  ggpubr::ggarrange(plotlist = plots, nrow=3,ncol=3)
ggpubr::ggexport(multi.page, filename = outplotfile,width=20,height=15)

plots <- mclapply(files,plotcontig_orfsize)
outplotfile <- "genome_plots/cleanasm_orfsize.pdf"
multi.page <-  ggpubr::ggarrange(plotlist = plots, nrow=3,ncol=3)
ggpubr::ggexport(multi.page, filename = outplotfile,width=20,height=15)

files <- list.files(path="results/genome_stats_nofilter", pattern="*.sum_stat.tsv",
                    full.names=TRUE, recursive=FALSE)

plots <- mclapply(files, plotcontig_intergenic)
outplotfile <- sprintf("genome_plots/nofilter_intergenic.pdf")
multi.page <-  ggpubr::ggarrange(plotlist = plots, nrow=3,ncol=3)
ggpubr::ggexport(multi.page, filename = outplotfile,width=20,height=15)

plots <- mclapply(files, plotcontig_pfam)
outplotfile <- "genome_plots/nofilter_pfamstat.pdf"
multi.page <-  ggpubr::ggarrange(plotlist = plots, nrow=3,ncol=3)
ggpubr::ggexport(multi.page, filename = outplotfile,width=20,height=15)

plots <- mclapply(files,plotcontig_orfsize)
outplotfile <- "genome_plots/nofilter_orfsize.pdf"
multi.page <-  ggpubr::ggarrange(plotlist = plots, nrow=3,ncol=3)
ggpubr::ggexport(multi.page, filename = outplotfile,width=20,height=15)

files <- list.files(path="results/genome_stats_ref", pattern="*.sum_stat.tsv",
                    full.names=TRUE, recursive=FALSE)

plots <- mclapply(files, plotcontig_intergenic)
outplotfile <- sprintf("genome_plots/refgenomes_intergenic.pdf")
multi.page <-  ggpubr::ggarrange(plotlist = plots, nrow=3,ncol=3)
ggpubr::ggexport(multi.page, filename = outplotfile,width=20,height=15)

plots <- mclapply(files, plotcontig_pfam)
outplotfile <- "genome_plots/refgenomes_pfamstat.pdf"
multi.page <-  ggpubr::ggarrange(plotlist = plots, nrow=3,ncol=3)
ggpubr::ggexport(multi.page, filename = outplotfile,width=20,height=15)

plots <- mclapply(files,plotcontig_orfsize)
outplotfile <- "genome_plots/refgenomes_orfsize.pdf"
multi.page <-  ggpubr::ggarrange(plotlist = plots, nrow=3,ncol=3)
ggpubr::ggexport(multi.page, filename = outplotfile,width=20,height=15)

