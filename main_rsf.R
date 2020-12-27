library(randomForestSRC)
library(survival)
library(glmnet)
library(foreach)
library(doParallel)


source(paste0(getwd(),"models.R"))

options(rf.cores=detectCores(), mc.cores=detectCores())
registerDoParallel(8)

###############################################################################################
###############################################################################################
# APPLY PRER GENERATION
###############################################################################################
###############################################################################################
# First read random walks
randomWalks <- read.table(file = "randomWalk/randomWalks.txt", header = F, sep = " ")

# Generate PRERs for each cancer type
ov_PRER <- comparator3(rppa = rppa, dd = ov$data,randomWalks = randomWalks,nnetwork)

gbm_PRER <- comparator3(rppa = rppa, dd = gbm$data,randomWalks = randomWalks,nnetwork)

brca_PRER <- comparator3(rppa = rppa, dd = brca$data,randomWalks = randomWalks)

blca_PRER <- comparator3(rppa = rppa, dd = blca$data,randomWalks = randomWalks)

kirc_PRER <- comparator3(rppa = rppa, dd = kirc$data,randomWalks = randomWalks)

lusc_PRER <- comparator3(rppa = rppa, dd = lusc$data,randomWalks = randomWalks)

luad_PRER <- comparator3(rppa = rppa, dd = luad$data,randomWalks = randomWalks)

coad_PRER <- comparator3(rppa = rppa, dd = coad$data,randomWalks = randomWalks)

ucec_PRER <- comparator3(rppa = rppa, dd = ucec$data,randomWalks = randomWalks)

hnsc_PRER <- comparator3(rppa = rppa, dd = hnsc$data,randomWalks = randomWalks)


###############################################################################################
###############################################################################################
# APPLY COX SCREEN
###############################################################################################
###############################################################################################
# Start Cox Screen to eliminate features which are not significant

ov_PRER_cox <- filt_rw(screening.cox(ov,ov_PRER,1)

gbm_PRER_cox <- filt_rw(screening.cox(gbm,gbm_PRER,1))

kirc_PRER_cox <- filt_rw(screening.cox(kirc,kirc_PRER,1))

blca_PRER_cox <- filt_rw(screening.cox(blca,blca_PRER,1))

brca_PRER_cox <- filt_rw(screening.cox(brca,brca_PRER,1))

hnsc_PRER_cox <- filt_rw(screening.cox(hnsc,hnsc_PRER,1))

luad_PRER_cox <- filt_rw(screening.cox(luad,luad_PRER,1))

coad_PRER_cox <- filt_rw(screening.cox(coad,coad_PRER,1))

ucec_PRER_cox <- filt_rw(screening.cox(ucec,ucec_PRER,1))

lusc_PRER_cox <- filt_rw(screening.cox(lusc,lusc_PRER,1))

###############################################################################################
###############################################################################################
# APPLY RANDOM SURVIVAL FOREST MODEL
###############################################################################################
###############################################################################################
# Start RSF model training for each cancer type
OV_PRER_RESULTS <- applyRSF(ov,ov_PRER,ov_PRER_cox,1)

KIRC_PRER_RESULTS <- applyRSF(kirc,kirc_PRER,kirc_PRER_cox,1)

GBM_PRER_RESULTS <- applyRSF(gbm,gbm_PRER,gbm_PRER_cox,1)

BLCA_PRER_RESULTS <- applyRSF(blca,blca_PRER,blca_PRER_cox,1)

BRCA_PRER_RESULTS <- applyRSF(brca,brca_PRER,brca_PRER_cox,1)

HNSC_PRER_RESULTS <- applyRSF(hnsc,hnsc_PRER,hnsc_PRER_cox,1)

LUAD_PRER_RESULTS <- applyRSF(luad,luad_PRER,luad_PRER_cox,1)

LUSC_PRER_RESULTS <- applyRSF(lusc,lusc_PRER,lusc_PRER_cox,1)

COAD_PRER_RESULTS <- applyRSF(coad,coad_PRER,coad_PRER_cox,1)

UCEC_PRER_RESULTS <- applyRSF(ucec,ucec_PRER,ucec_PRER_cox,1)

# Display C-index quantiles
quantile(BRCA_PRER_RESULTS$cinds)

quantile(OV_PRER_RESULTS$cinds)

quantile(KIRC_PRER_RESULTS$cinds)

quantile(HNSC_PRER_RESULTS$cinds)

quantile(LUSC_PRER_RESULTS$cinds)

quantile(LUAD_PRER_RESULTS$cinds)

quantile(GBM_PRER_RESULTS$cinds)

quantile(BLCA_PRER_RESULTS$cinds)

quantile(COAD_PRER_RESULTS$cinds)

quantile(UCEC_PRER_RESULTS$cinds)
