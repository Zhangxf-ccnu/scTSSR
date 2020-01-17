library("scTSSR")
library("SAVER")
library("scRMD")
library("scImpute")
library("VIPER")
library("DrImpute")
library("Rmagic")
library("zinbwave")
library("SCRABBLE")
source("alra.R")


# perform different imputation methods on Deng dataset

deng_count_data <- readRDS('scTSSR-data/deng_count_data.rds')

deng_pre <- log_normalization(deng_count_data, percent=0.1, preprocess.only=TRUE)

deng_lognor <- log_normalization(deng_count_data, percent=0.1)

## scTSSR
deng_tssr <- scTSSR(deng_count_data, percent=0.1, learning_rate=0.0001, epochs=100, all_samples_as_batch=TRUE)$estimate

## ALRA
deng_alra <- alra(t(deng_lognor))

## SAVER
deng_saver <- saver(deng_pre, ncores=4)$estimate

## scimpute
scimpute(count_path = 'scTSSR-data/deng_pre.rds',
         infile = "rds",           
         outfile = "rds",          
         out_dir = "./",          
         drop_thre = 0.5,          
         ncores = 2)
deng_scimpute <- as.matrix(read.rds("scimpute_count.rds"))      

## scRMD
deng_scrmd <- rmd(t(deng_lognor))$exprs

## VIPER
deng_viper <- VIPER(as.matrix(deng_pre), report = FALSE)$imputed

## DrImpute
count_DrImpute = DrImpute::DrImpute(deng_lognor)

count_DrImpute = exp(count_DrImpute)-1    

rownames(count_DrImpute) = rownames(deng_lognor)
colnames(count_DrImpute) = colnames(deng_lognor)
count_DrImpute = as.matrix(count_DrImpute)

## MAGIC
count <- readRDS('scTSSR-data/deng_pre.rds')
count.t <- t(count)
count.normalized <- Rmagic::library.size.normalize(count.t)
count.log <- log(count.normalized + 1)

count_MAGIC <- Rmagic::magic(count.log)

count_MAGIC <- t(as.matrix(count_MAGIC))
count_MAGIC[count_MAGIC<0]<-0
count_MAGIC <- exp(count_MAGIC)-1

row.names(count_MAGIC) <- colnames(count.log)
colnames(count_MAGIC) <- rownames(count.log)
count_MAGIC <- as.matrix(count_MAGIC)

## ZINB-WaVE
m <- zinbFit(round(deng_pre))
result = imputeZeros(m, t(deng_pre))
result = t(result)

rownames(result) = rownames(deng_pre)
colnames(result) = colnames(deng_pre)
result = as.matrix(result)

## SCRABBLE
count = readRDS('scTSSR-data/deng_pre.rds')
input_count = vector('list', 2)
input_count[[1]] = count
parameter <- c(1,1e-6,1e-4)
result <- scrabble(input_count, parameter = parameter)

rownames(result) = rownames(count)
colnames(result) = colnames(count)
result = as.matrix(result)







