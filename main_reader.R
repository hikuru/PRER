rm(list=ls())

source("reader.R")

wd <- 'data/downloaded_data/'
processes_dataDir <- 'data/processed_data/'

# Parse and save all 10 cancer types
# Make 100 bootstrap samples for each cancer type

gbm<-reader(131,wd,"TCGA_GBM_RPPA_RBN-2015-02-24","gbm")
save(file = paste(processes_dataDir,'gbm', sep = ''),gbm)

blca<-reader(131,wd,"TCGA_BLCA_RPPA_RBN-2015-02-24","blca")
save(file = paste(processes_dataDir,'blca', sep = ''),blca)

hnsc<-reader(131,wd,"TCGA_HNSC_RPPA_RBN-2015-02-24","hnsc")
ssave(file = paste(processes_dataDir,'hnsc', sep = ''),hnsc)

kirc<-reader(131,wd,"TCGA_KIRC_RPPA_RBN-2015-02-24","kirc")
save(file = paste(processes_dataDir,'kirc', sep = ''),kirc)

luad<-reader(131,wd,"TCGA_LUAD_RPPA_RBN-2015-02-24","luad")
save(file = paste(processes_dataDir,'luad', sep = ''),luad)

lusc<-reader(131,wd,"TCGA_LUSC_RPPA_RBN-2015-02-24","lusc")
save(file = paste(processes_dataDir,'lusc', sep = ''),lusc)

coad<-read<-reader(131,wd,"TCGA_COAD_RPPA_RBN-2015-02-24","coad")
save(file = paste(processes_dataDir,'coad', sep = ''),coad)

ucec<-read<-reader(131,wd,"TCGA_UCEC_RPPA_RBN-2015-02-24","ucec")
save(file = paste(processes_dataDir,'ucec', sep = ''),ucec)

read<-read<-reader(131,wd,"TCGA_READ_RPPA_RBN-2015-02-24","read")
save(file = paste(processes_dataDir,'read', sep = ''),read)

brca<-reader(131,wd,"TCGA_BRCA_RPPA_RBN-2015-02-24","brca")
save(file = paste(processes_dataDir,'brca', sep = ''),brca)

ov<-reader(131,wd,"TCGA_OV_RPPA_RBN-2015-02-24","ov")
save(file = paste(processes_dataDir,'ov', sep = ''),ov)

