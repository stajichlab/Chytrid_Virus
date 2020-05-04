library(ggplot2)

files <- list.files(path="results/genome_stats_nofilter", pattern="Burma_1F.sum_stat.tsv", 
                    full.names=TRUE, recursive=FALSE)

lapply(files, function(x) {
  t <- read.table(x, header=TRUE) # load file
  tlong <- subset(t,t$LENGTH >= 500 & t$RIBOSOMAL != 1 & t$MT != 1 & COVERAGE < 1000)
  p <- ggplot(data=tlong,aes(x=GC,y=COVERAGE,
                               color=factor(CAPSID),
                               shape=factor(floor(ORF_MEAN)),
                               size=LENGTH)) + 
    geom_point() + scale_fill_brewer() + ggtitle(x) 
  ggsave(filename="genome_plots_%03d.pdf",p,width=15)
})

all=read.table("results/genome_stats_nofilter/Burma_1F.sum_stat.tsv",header=T)
sumstat = subset(all,all$LENGTH >= 500 & all$RIBOSOMAL != 1 & all$MT != 1)


hist(sumstat$ORF_MEAN,100)
hist(log(sumstat$LENGTH)/log(10),100)
hist(sumstat$GC,100)
subset(all,all$GC < 40)

all=read.table("results/genome_stats_nofilter/California_12.sum_stat.tsv",header=T)
sumstat <- subset(all,all$LENGTH >= 500)

p <- ggplot(data=sumstat,aes(x=GC,y=COVERAGE,color=ORF_PER_KB,
                             size=ORF_MEAN)) + 
  geom_point() + scale_fill_brewer()
p

hist(sumstat$ORF_MEAN,100)
hist(log(sumstat$LENGTH)/log(10),100)
hist(sumstat$GC,100)
subset(all,all$GC < 40)

all=read.table("results/genome_stats_nofilter/JEL0729.sum_stat.tsv",header=T)
sumstat <- subset(all,all$LENGTH >= 500)

p <- ggplot(data=sumstat,aes(x=GC,y=COVERAGE,
                             color=factor(CAPSID),
                             shape=factor(floor(ORF_MEAN)),
                             size=LENGTH)) + 
  geom_point() + scale_fill_brewer()
p

hist(sumstat$ORF_MEAN,100)
hist(log(sumstat$LENGTH)/log(10),100)
hist(sumstat$GC,100)
subset(all,all$GC < 40)