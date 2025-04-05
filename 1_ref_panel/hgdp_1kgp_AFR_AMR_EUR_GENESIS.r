################################################################################
################################################################################
## Install packages
## requires R version 4
# conda create -n g_env -c bioconda bioconductor-genesis
# conda activate g_env

# vignette for sequence data
# https://mirrors.dotsrc.org/bioconductor-releases/3.14/bioc/vignettes/GENESIS/inst/doc/assoc_test_seq.html

# ">0.354" duplicate/MZ twin
# [0.177, 0.354] 1st-degree
# [0.0884, 0.177] 2nd-degree
# [0.0442, 0.0884] 3rd-degree

## Start PC-Relate and PC-AiR
library(GENESIS)
library(SNPRelate)
library(GWASTools)
library(data.table)

# v2: AFR, AMR, EUR n = 2252
iDir <- "/scratch/eger/projects/hgdp_1kg/"
tDir <- "/scratch/eger/projects/hgdp_1kg/GENESIS_v2/"

# prefix of dataset
annot1 <- "hgdp_1kg.hg38"

################################################################################
## subset the data
annot2 <- "hgdp_1kg.AFR_AMR_EUR.auto.geno001.maf005.hwe10-5"
command <- paste("plink", 
                 "--bfile", paste(iDir, annot1, sep=""),
                 "--keep", paste0(tDir, "hgdp_1kg.AFR_AMR_EUR_samples.txt"),
                 "--geno 0.01",
                 "--maf 0.05",
                 "--hwe 0.00001",
                 "--allow-no-sex",
                 "--make-bed",
                 "--out", paste(tDir, annot2, sep=""),
                 sep=" ")
system(command)   
# 5,750,680 variants and 2252 people
annot1 <- annot2

# add SNP names so that pruning will work
# Read the .bim file
bim_file <- paste(tDir, annot1, ".bim", sep="")
bim <- fread(bim_file, header = FALSE)

# Modify variant names to chr:pos:A2:A1 format
bim[, V2 := paste0(V1, ":", V4, ":", V6, ":", V5)]

# Save the updated .bim file
fwrite(bim, bim_file, sep = "\t", quote = FALSE, col.names = FALSE)

## prune the data
# this is a light filtering to make snpgdsLDpruning() run faster
prune <- "hgdp_1kg.AFR_AMR_EUR.auto.geno001.maf005.hwe10-5.indep_SNPs"
command <- paste("plink", 
                 "--bfile", paste(tDir, annot1, sep=""),
                 "--indep-pairwise 10000 100 0.35",
                 "--out", paste(tDir, prune, sep=""),
                 sep=" ")
system(command)

annot2 <- "hgdp_1kg.AFR_AMR_EUR.auto.geno001.maf005.hwe10-5.pruned"
command <- paste("plink", 
                 "--bfile", paste(tDir, annot1, sep=""),
                 "--extract", paste0(tDir, prune, ".prune.in"),
                 "--allow-no-sex",
                 "--make-bed",
                 "--out", paste(tDir, annot2, sep=""),
                 sep=" ")
system(command)   
# 756,053 variants and 2252 people 

################################################################################
## convert plink to .gds file

fold <- tDir
file <- annot2

snpgdsBED2GDS(bed.fn = paste(fold,file,".bed",sep=""), 
              bim.fn = paste(fold,file,".bim",sep=""), 
              fam.fn = paste(fold,file,".fam",sep=""), 
              out.gdsfn = paste(fold,file,".gds",sep=""))

################################################################################
## LD prune

# read in GDS data
gdsfile <- paste(fold,file,".gds",sep="")
gds <- snpgdsOpen(gdsfile)

# prune SNPs
snpset <- snpgdsLDpruning(gds, method="corr", 
                          slide.max.bp=10e6, 
                          ld.threshold=sqrt(0.1))

pruned <- unlist(snpset, use.names=FALSE)
length(pruned) # 339,619

# save pruned
save(list = c("pruned"), file = paste(fold,file,"_pruned.Rdata",sep=""))

################################################################################
## Run KING to obtain robust and homogenous inferrence of sample relatedness

# Calculate IBD coefficients by KING method of moment
# usually do not used pruned but this is a big dataset
# We select a subset of variants for this calculation with LD pruning.
king <- snpgdsIBDKING(gds, snp.id=pruned)

# get the IBD table and save it
kinship <- snpgdsIBDSelection(king)
relatedpairs <- kinship
relatedpairs <- relatedpairs[which(relatedpairs$kinship>2^(-11/2)),]

# check what the highest kinship coeff is
tail(relatedpairs[order(relatedpairs$kinship),])

# save the related pairs
write.table(relatedpairs, paste(fold,file,"_king_relatedpairs.txt",sep=""), 
            row.names = F, quote = F)

# convert KING to R matrix (when running KING in R)
KINGmat <- kingToMatrix(king)
# 798 samples with no relatives included

# save KINGmat
save(list = c("KINGmat"), file = paste(fold,file,"_KINGmat.Rdata",sep=""))

# remove king objects that are no longer needed
rm(king)
rm(kinship)
rm(relatedpairs)

################################################################################
## Run PC-AiR

#load(paste(fold,file,"_pruned.Rdata",sep=""))
#load(paste(fold,file,"_KINGmat.Rdata",sep=""))

# PC-AiR
# default kin.threshold = 2^(-11/2) ~ 0.022
mypcair <- pcair(gds,
                 kinobj = KINGmat,
                 divobj = KINGmat,
                 snp.include = pruned)

# summary of output
summary(mypcair)
# Unrelated Set: 1704 Samples 
# Related Set: 548 Samples 

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

# # save lists of related and unrelated sets
# write.table(mypcair$rels, paste(fold,file,"_pcair_rels.txt",sep=""), 
#             quote = F, col.names = F, row.names = F)
write.table(mypcair$unrels, paste(fold,file,"_pcair_unrels.txt",sep=""), 
            quote = F, col.names = F, row.names = F)

# close GDS data
snpgdsClose(gds)
