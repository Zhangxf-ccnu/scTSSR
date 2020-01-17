library("SingleCellExperiment")
library("TSCAN")
library("scater")
library("ggplot2")
library("igraph")
source("TSCAN.R")



# run tscan
## Deng dataset
## load the imputation results
deng_pre <- readRDS('scTSSR-data/deng_pre.rds')
deng_tssr <- readRDS('scTSSR-data/deng_tssr.rds')
deng_alra <- readRDS('scTSSR-data/deng_alra.rds')
deng_saver <- readRDS('scTSSR-data/deng_saver.rds')
deng_scimpute <- readRDS('scTSSR-data/deng_scimpute.rds')
deng_scrmd <- readRDS('scTSSR-data/deng_scrmd.rds')
deng_viper <- readRDS('scTSSR-data/deng_viper.rds')
deng_drimpute <- readRDS('scTSSR-data/deng_drimpute.rds')
deng_magic <- readRDS('scTSSR-data/deng_magic.rds')
deng_scrabble <- readRDS('scTSSR-data/deng_scrabble.rds')
deng_scvi <- readRDS('scTSSR-data/deng_scvi.rds')
deng_autoimpute <- readRDS('scTSSR-data/deng_autoimpute.rds')
deng_zinbwave <- readRDS('scTSSR-data/deng_zinbwave.rds')


deng_cellLabels <- factor(colnames(deng_count_data),
                         levels=c('zygote', 'early 2-cell', 'mid 2-cell', 'late 2-cell',
                                  '4-cell', '8-cell', '16-cell', 'early blastocyst',
                                  'mid blastocyst', 'late blastocyst'))

deng.raw.tscan <- my.TSCAN(deng_pre, deng_cellLabels)
deng.tssr.tscan <- my.TSCAN(deng_tssr, deng_cellLabels)
deng.rmd.tscan <- my.TSCAN(t(exp(deng_scrmd)-1), deng_cellLabels)
deng.alra.tscan <- my.TSCAN(t(exp(deng_alra$A_norm_rank_k_cor_sc)-1), deng_cellLabels)
deng.saver.tscan <- my.TSCAN(deng_saver, deng_cellLabels)
deng.viper.tscan <- my.TSCAN(deng_viper, deng_cellLabels)
deng.scimpute.tscan <- my.TSCAN(deng_scimpute, deng_cellLabels)
deng.zinbwave.tscan <- my.TSCAN(deng_zinbwave, deng_cellLabels)
deng.autoimpute.tscan <- my.TSCAN(deng_autoimpute, deng_cellLabels)
deng.scvi.tscan <- my.TSCAN(deng_scvi, deng_cellLabels)
deng.scrabble.tscan <- my.TSCAN(deng_scrabble, deng_cellLabels)
deng.magic.tscan <- my.TSCAN(deng_magic, deng_cellLabels)
deng.drimpute.tscan <- my.TSCAN(deng_drimpute, deng_cellLabels)

save.image("scTSSR-data/deng_tscan_results.RData")


## plot trajectory
load("scTSSR-data/deng_tscan_results.RData")
count <- deng_tssr
colnames(count) <- c(1:ncol(count))
procdata <- TSCAN::preprocess(count)
lpsmclust <- TSCAN::exprmclust(procdata)
pdf(paste0('fig_trajectory_deng_tssr.pdf'), 6, 6)
plotmclust2(lpsmclust)
dev.off()



## plot POS
POS <- data.frame(method = c('Observed','scTSSR',  'ALRA','AutoImpute','DrImpute', 
                             'MAGIC','SAVER', 'scImpute','SCRABBLE','scRMD','scVI', 'VIPER', 'ZINB-WaVE'),
                 POS=c(deng.raw.tscan$POS, deng.tssr.tscan$POS, deng.alra.tscan$POS,
                       deng.autoimpute.tscan$POS, deng.drimpute.tscan$POS, deng.magic.tscan$POS,
                       deng.saver.tscan$POS,deng.scimpute.tscan$POS, deng.scrabble.tscan$POS,
                       deng.rmd.tscan$POS, deng.scvi.tscan$POS, deng.viper.tscan$POS, deng.zinbwave.tscan$POS))
POS$method <- factor(POS$method, levels = c('Observed','scTSSR',  'ALRA','AutoImpute','DrImpute', 
                                            'MAGIC','SAVER', 'scImpute','SCRABBLE','scRMD','scVI', 'VIPER', 'ZINB-WaVE'))

theme_set(theme_bw())
p <- ggplot(data = POS,
            aes(y = abs(POS),
                x = method,
                fill = method))  + scale_fill_manual(values=c("#E69F00", "#009E73","#56B4E9",  "#CC79A7", "#0072B2", '#F0E442','#D55E00',
                                                              '#3399FF', '#FF99CC', '#996600', "#CC6666", "#9999CC", "#66CC99"))  #+ scale_x_continuous(breaks=c(0.2,0.4,0.6,0.8))

p <- p + geom_line(size=2) + geom_bar(stat="identity", position=position_dodge())

p <- p + theme(axis.title.x = element_text(size = 20, face = "bold"),
               axis.title.y = element_text(size = 20, face = "bold", angle = 90),
               axis.text.x = element_text(size = 15),
               axis.text.y = element_text(size = 15),
               strip.text.x = element_text(size = 20, face = "bold"),
               strip.text.y = element_text(size = 20, face = "bold"))
p <- p + theme(axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               legend.title = element_text(size = 15, face="bold"),
               legend.text = element_text(size = 15, face="bold"))
p <- p + xlab("") + ylab("POS") + coord_cartesian(ylim = c(0.75, 0.95)) 

ggsave(p, file="POS_deng.pdf", width = 7.5, height = 5)

## plot kendall
cor <- data.frame(method = c('Observed','scTSSR',  'ALRA','AutoImpute','DrImpute', 
                             'MAGIC','SAVER', 'scImpute','SCRABBLE','scRMD','scVI', 'VIPER', 'ZINB-WaVE'),
                 cor=c(deng.raw.tscan$cor.kendall, deng.tssr.tscan$cor.kendall, deng.alra.tscan$cor.kendall,
                       deng.autoimpute.tscan$cor.kendall, deng.drimpute.tscan$cor.kendall, deng.magic.tscan$cor.kendall,
                       deng.saver.tscan$cor.kendall, deng.scimpute.tscan$cor.kendall, deng.scrabble.tscan$cor.kendall,
                       deng.rmd.tscan$cor.kendall, deng.scvi.tscan$cor.kendall, deng.viper.tscan$cor.kendall, deng.zinbwave.tscan$cor.kendall))
cor$method <- factor(cor$method, levels = c('Observed','scTSSR',  'ALRA','AutoImpute','DrImpute', 
                                            'MAGIC','SAVER', 'scImpute','SCRABBLE','scRMD','scVI', 'VIPER', 'ZINB-WaVE'))

theme_set(theme_bw())
p <- ggplot(data = cor,
            aes(y = abs(cor),
                x = method,
                fill = method))  + scale_fill_manual(values=c("#E69F00", "#009E73","#56B4E9",  "#CC79A7", "#0072B2", '#F0E442','#D55E00',
                                                              '#3399FF', '#FF99CC', '#996600', "#CC6666", "#9999CC", "#66CC99"))  

p <- p + geom_line(size=2) + geom_bar(stat="identity", position=position_dodge())

p <- p + theme(axis.title.x = element_text(size = 20, face = "bold"),
               axis.title.y = element_text(size = 20, face = "bold", angle = 90),
               axis.text.x = element_text(size = 15),
               axis.text.y = element_text(size = 15),
               strip.text.x = element_text(size = 20, face = "bold"),
               strip.text.y = element_text(size = 20, face = "bold"))
p <- p + theme(axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               legend.title = element_text(size = 15, face="bold"),
               legend.text = element_text(size = 15, face="bold"))
p <- p + xlab("") + ylab("Cor (Kendall)") + coord_cartesian(ylim = c(0.65, 0.82)) 

ggsave(p, file="kendall_cor_deng.pdf", width = 7.5, height = 5)




