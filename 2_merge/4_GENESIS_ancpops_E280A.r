################################################################################
################################################################################
## Install packages
## requires R version 4
# conda create -n g_env -c bioconda bioconductor-genesis
# conda activate g_env

# I used the pruining settings from the vignette
# https://bioconductor.org/packages/devel/bioc/vignettes/GENESIS/inst/doc/pcair.html

## Start PC-Relate and PC-AiR
library(GENESIS)
library(SNPRelate)
library(GWASTools)

iDir <- "/scratch/eger/projects/E280A_SNP_arrays/E280A_340/merge_hgdp_1kg/"
tDir <- "/scratch/eger/projects/E280A_SNP_arrays/E280A_340/GENESIS_v2/"

# make tDir if it doesn't exist
if (!dir.exists(tDir)) {
  dir.create(tDir)
}

# prefix of dataset
annot1 <- "hgdp_1kg.AFR_AMR_EUR_unrelated_E280A_340.auto.hg38"

################################################################################
## convert plink to .gds file

fold <- tDir
file <- annot1

snpgdsBED2GDS(bed.fn = paste(iDir,file,".bed",sep=""), 
              bim.fn = paste(iDir,file,".bim",sep=""), 
              fam.fn = paste(iDir,file,".fam",sep=""), 
              out.gdsfn = paste(fold,file,".gds",sep=""))

################################################################################
## LD prune

# read in GDS data
gdsfile <- paste(fold,file,".gds",sep="")
gds <- snpgdsOpen(gdsfile)

# prune SNPs
snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=10e6, 
                          ld.threshold=sqrt(0.1))
pruned <- unlist(snpset, use.names=FALSE)
length(pruned) # 161,855

# save pruned
save(list = c("pruned"), file = paste(fold,file,"_pruned.Rdata",sep=""))

################################################################################
## Run KING to obtain robust and homogenous inferrence of sample relatedness

load(paste(fold,file,"_pruned.Rdata",sep=""))

# Calculate IBD coefficients by KING method of moment
king <- snpgdsIBDKING(gds) # no need to prune SNPs for KING

# get the IBD table and save it
kinship <- snpgdsIBDSelection(king)
relatedpairs <- kinship
relatedpairs <- relatedpairs[which(relatedpairs$kinship>2^(-11/2)),]
write.table(relatedpairs, paste(fold,file,"_king_relatedpairs.txt",sep=""), row.names = F, quote = F)

# convert KING to R matrix (when running KING in R)
KINGmat <- kingToMatrix(king)

# Using 1138 samples provided
# Identifying clusters of relatives...
#     334 relatives in 19 clusters; largest cluster = 213
# Creating block matrices for clusters...
# 804 samples with no relatives included
# Putting all samples together into one block diagonal matrix

# save KINGmat
save(list = c("KINGmat"), file = paste(fold,file,"_KINGmat.Rdata",sep=""))

# remove king objects that are no longer needed
rm(king)
rm(kinship)
rm(relatedpairs)

################################################################################
## Run PC-AiR

# load(paste(fold,file,"_KINGmat.Rdata",sep=""))

# PC-AiR without reference populations
mypcair <- pcair(gds,
                 kinobj = KINGmat,
                 divobj = KINGmat,
                 snp.include = pruned)

# summary of output
summary(mypcair)

# Sample Size: 1138 
# Unrelated Set: 870 Samples 
# Related Set: 268 Samples 

# save mypcair
save(list = c("mypcair"), file = paste(fold,file,"_pcair.Rdata",sep=""))

# save matrix of PCs
PCs <- mypcair$vectors
write.table(PCs, paste(fold,file,"_pcair_PCs.txt",sep=""), 
            quote = F, col.names = F)

# save The variance proportion for each principal component
write.table(mypcair$varprop[!is.na(mypcair$varprop)], 
            paste(fold,file,"_pcair_varprop.txt",sep=""), 
            quote = F, col.names = c("varprop"), row.names = F)

# save lists of related and unrelated sets
write.table(mypcair$rels, paste(fold,file,"_pcair_rels.txt",sep=""), 
            quote = F, col.names = F, row.names = F)
write.table(mypcair$unrels, paste(fold,file,"_pcair_unrels.txt",sep=""), 
            quote = F, col.names = F, row.names = F)

# close GDS data
snpgdsClose(gds)

################################################################################
## Run PC-Relate

# read in gds data
gdsfile <- paste(fold,file,".gds",sep="")
gds <- GdsGenotypeReader(filename = gdsfile)

# create a GenotypeData class object
genoData <- GenotypeData(gds)

# create genotype block data object
genoData_block <- GenotypeBlockIterator(genoData, snpInclude=pruned)

# run my pcrelate
mypcrelate <- pcrelate(genoData_block, pcs = mypcair$vectors[,1:2], 
                       training.set = mypcair$unrels)

# save mypcrelate
save(list = c("mypcrelate"), file = paste(fold,file,"_pcrelate.Rdata",sep=""))

################################################################################
## Get relatedness info / GRM from PC-Relate

relatedpairs <- mypcrelate$kinBtwn
relatedpairs <- relatedpairs[which(relatedpairs$kin>2^(-11/2)),]
inbreedcoef <- mypcrelate$kinSelf
GRM <- pcrelateToMatrix(mypcrelate)

unrelatedpairs <- mypcrelate$kinBtwn
unrelatedpairs <- unrelatedpairs[which(unrelatedpairs$kin<2^(-11/2)),]

save(list = c("GRM"), file = paste(fold,file,"_GRM.Rdata",sep=""))
write.table(unrelatedpairs, paste(fold,file,"_unrelatedpairs.txt",sep=""), row.names = F, quote = F)
write.table(relatedpairs, paste(fold,file,"_relatedpairs.txt",sep=""), row.names = F, quote = F)
write.table(inbreedcoef, paste(fold,file,"_inbreedcoef.txt",sep=""), row.names = F, quote = F)
write.table(mypcrelate$kinBtwn, paste(fold,file,"_kinship.txt",sep=""), row.names = F, quote = F)