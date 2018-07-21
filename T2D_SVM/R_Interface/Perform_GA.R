#!/usr/bin/R
#2017/11/27
# R CMD BATCH '--args [popN] [Pco] [Pmu] [K] [modN] [Pnew] [mR_mR_corr] [miR_miR_corr] [miR_mR_corr] [mR_mR_corr_pvalue] [miR_miR_corr_pvalue]
# [miR_mR_pvalue] [Output_Filename]' --no-save --no-restore Perform_GA.R
#
#-----------------------------------------------------------------------------
# initialize parameter
#-----------------------------------------------------------------------------

args <- commandArgs(TRUE);
# population size
popN <- as.numeric(args[1])
# crossover probability
Pco <- as.numeric(args[2])
# putation probability
Pmu <- as.numeric(args[3])
# iteration times
K <- as.numeric(args[4])
# number of modules
modN <- as.numeric(args[5])
#probability of introducing fresh chromesomes
Pnew <- as.numeric(args[6])
# iteration times for Local search
LK <- as.numeric(args[7])
###read in correlation matrix
#miR_miR_corr
miR_miR_corr <- as.character(args[8])
#mR_mR_corr
mR_mR_corr <- as.character(args[9])
#miR_mR_corr
miR_mR_corr <-as.character(args[10])

###Read in p value matrix
#miR_miR_corr_pvalue
miR_miR_corr_pvalue <- as.character(args[11])
#mR_mR_corr_pvalue
mR_mR_corr_pvalue <- as.character(args[12])
#miR_mR_corr_pvalue
miR_mR_corr_pvalue <- as.character(args[13])
# File  for output
output <- as.character(args[14])

# Script for generating score matrix
#-----------------------------------------------------------------
# Functions
#-----------------------------------------------------------------

# Split strings by delimited to matrix
strsplit2 <- function (x, split, ...)
{
    x <- as.character(x)
    n <- length(x)
    s <- strsplit(x, split = split, ...)
    nc <- unlist(lapply(s, length))
    out <- matrix("", n, max(nc))
    for (i in 1:n) {
        if (nc[i])
            out[i, 1:nc[i]] <- s[[i]]
    }
    out
}

# Transform matrix with three columns to vector, the third column contains the value
# The first two column are paste by : to be names of the vector
unmatrix <- function (x)
{
    rnames <- rownames(x)
    cnames <- colnames(x)
    if (is.null(rnames))
        rnames <- paste("r", 1:nrow(x), sep = "")
    if (is.null(cnames))
        cnames <- paste("c", 1:ncol(x), sep = "")
    nmat <- outer(rnames, cnames, paste, sep = ":")
    vlist <- c(t(x))
    names(vlist) <- c(t(nmat))
    return(vlist)
}

# Transform list to a matrix with three columns
list2matrix <- function(list,names1,names2)
{
	list <- list[(list[,1] %in% names1) & (list[,2] %in% names2),]
	out <- matrix(0,nrow=length(names1),ncol=length(names2),dimnames=list(names1,names2))
	for(i in 1:dim(list)[1]){
		out[as.character(list[i,1]),as.character(list[i,2])] <- list[i,3]
	}
	out
}


#-----------------------------------------------------------------------------
# data process
#-----------------------------------------------------------------------------

# Read in corelation matrices
miR_miR_corr_raw <- as.matrix(read.table(miR_miR_corr,sep="\t",check.names=F))
mR_mR_corr_raw <- as.matrix(read.table(mR_mR_corr,sep="\t"))
miR_mR_corr_raw <- as.matrix(read.table(miR_mR_corr,sep="\t"))

# correct bugs which may existed in colnames
colnames(miR_miR_corr_raw) <- rownames(miR_miR_corr_raw)
colnames(mR_mR_corr_raw) <- rownames(mR_mR_corr_raw)
colnames(miR_mR_corr_raw) <- rownames(mR_mR_corr_raw)

# Read in p value matrices
miR_miR_corr_p <- as.matrix(read.table(miR_miR_corr_pvalue,sep="\t",check.names=F))
miR_mR_corr_p <- as.matrix(read.table(miR_mR_corr_pvalue,sep="\t"))
mR_mR_corr_p <- as.matrix(read.table(mR_mR_corr_pvalue,sep="\t"))

# Correct the scenario that different p values for the same pairs of PCC
miR_miR_corr_p[lower.tri(miR_miR_corr_p)] <- t(miR_miR_corr_p)[lower.tri(miR_miR_corr_p)]
mR_mR_corr_p[lower.tri(mR_mR_corr_p)] <- t(mR_mR_corr_p)[lower.tri(mR_mR_corr_p)]

# Set the diagonal of p value matrices to 1
diag(miR_miR_corr_p) <- 1
diag(mR_mR_corr_p) <- 1

# Read in target prediction of TFs & miRNAs
TF_gene_raw <- read.table("TF_gene.txt",sep=" ")
miR_gene_raw <- read.table("miR_gene.txt",sep=" ")

# Find the overlap genes among target prediction of miRNAs, target prediction of TFs and differential expressed genes
gene_list <- sort(intersect(intersect(unique(miR_gene_raw$V3),unique(rownames(mR_mR_corr_raw))),unique(TF_gene_raw$V3)))
TF_list <- sort(intersect(unique(TF_gene_raw$V2),gene_list))
miR_list <- sort(intersect(unique(miR_gene_raw$V2),unique(rownames(miR_mR_corr_raw))))
nTF_list <- sort(gene_list[-match(TF_list,gene_list)])
gene_list <- c(TF_list,nTF_list)

# Transform target prediction of TFs to symmetric matrix
TF_gene_inter <- list2matrix(TF_gene_raw[,c(2,3,1)],TF_list,gene_list)
# Apply similarity score cutoff for target prediction of TFs
TF_gene_inter[TF_gene_inter<0.99] <- 0

# Transform target prediction of miRNAs to symmetric matrix
miR_gene_raw[,1] <- as.numeric(miR_gene_raw[,1])/max(as.numeric(miR_gene_raw[,1]))
miR_gene_inter <- list2matrix(miR_gene_raw[,c(2,3,1)],miR_list,gene_list)

# Calculate the length of miRNAs, TFs and nTFs
miR_length <- length(miR_list)
TF_length <- length(TF_list)
nTF_length <- length(nTF_list)

N <- miR_length + TF_length + nTF_length
TF_s_index <- miR_length+1
nTF_s_index <- miR_length+TF_length+1
mR_length <- TF_length+nTF_length
rg_length <- miR_length+TF_length

# save the orignal correlation matrices
miR_miR_corr_b <- miR_miR_corr_raw
miR_mR_corr_b <- miR_mR_corr_raw
mR_mR_corr_b <- mR_mR_corr_raw

# Apply a cutoff to correlation coefficient
miR_miR_corr_raw[miR_miR_corr_p > 1e-03] <- 0
miR_mR_corr_raw[miR_mR_corr_p > 1e-03] <- 0
mR_mR_corr_raw[mR_mR_corr_p > 1e-03] <- 0

# Keep only the miRNAs/genes with predicted information
miR_miR_corr_new <- miR_miR_corr_raw[miR_list,miR_list]
miR_mR_corr_new <- miR_mR_corr_raw[miR_list,gene_list]
mR_mR_corr_new <- mR_mR_corr_raw[gene_list,gene_list]

#diag(miR_miR_corr_new) <- 0
#diag(mR_mR_corr_new) <- 0

# Any PCC between regulator and its target genes are kept regardless of the p values
miR_mR_corr_new[which(miR_gene_inter>0,arr.ind=T)] <- miR_mR_corr_b[which(miR_gene_inter>0,arr.ind=T)]
mR_mR_corr_new[which(TF_gene_inter>0,arr.ind=T)] <- mR_mR_corr_b[which(TF_gene_inter>0,arr.ind=T)]
mR_mR_corr_new[which(t(TF_gene_inter>0),arr.ind=T)] <- mR_mR_corr_b[which(t(TF_gene_inter>0),arr.ind=T)]

#miR_miR_corr <- as.matrix(miR_miR_corr_new[miR_list,miR_list])
miR_miR_corr <- as.matrix(miR_miR_corr_new[miR_list,miR_list])
miR_mR_corr <- as.matrix(miR_mR_corr_new[miR_list,gene_list])
miR_mR_corr[miR_list,TF_list] <- miR_mR_corr[miR_list,TF_list]

#######mR_mR_corr
mR_mR_corr <- as.matrix(mR_mR_corr_new[gene_list,gene_list])


miR_corrs <- cbind(miR_miR_corr,miR_mR_corr[miR_list,gene_list]/2)
# Use negative correlation only
#
# miR_corrs <- cbind(miR_miR_corr,miR_mR_corr[miR_list,gene_list]/2)
# miR_corrs[miR_corrs>0] <- 0
#

# Parameters for control the contibuition of sequence-based prediction(k1) and correaation(k2)
k1=1
k2=1

# Construct Socre matrix as below
#  -----------------------
# | miR-miR miR-TF miR-nTF|
# | TF-miR  TF-TF  TF-nTF |
# |	nTF-miR nTF-TF nTF-nTF|
#  -----------------------
Score_matrix <- as.matrix(rbind(k1*cbind(matrix(0,nrow=miR_length,ncol=miR_length), miR_gene_inter/2) + k2*abs(miR_corrs), cbind(k2*t(abs(miR_mR_corr[miR_list,gene_list]))/2 + k1*t(miR_gene_inter)/2, rbind(k1*cbind(TF_gene_inter[TF_list,TF_list],TF_gene_inter[TF_list,nTF_list]/2),k1*cbind(t(TF_gene_inter[TF_list,nTF_list])/2,matrix(0,nrow=nTF_length,ncol=nTF_length))) + k2*abs(mR_mR_corr))))

# Scale each part of score matrix to [0.5,1]
tmp <- Score_matrix[1:miR_length,1:miR_length]
Score_matrix[1:miR_length,1:miR_length] <- (1-0.5)*(tmp-min(tmp[tmp>0]))/(max(tmp)-min(tmp[tmp>0]))+0.5
tmp <- Score_matrix[TF_s_index:rg_length,1:miR_length]
Score_matrix[TF_s_index:TF_s_index:rg_length,1:miR_length] <- (1-0.5)*(tmp-min(tmp[tmp>0]))/(max(tmp)-min(tmp[tmp>0]))+0.5
tmp <- Score_matrix[1:miR_length,TF_s_index:rg_length]
Score_matrix[1:miR_length,TF_s_index:rg_length] <- (1-0.5)*(tmp-min(tmp[tmp>0]))/(max(tmp)-min(tmp[tmp>0]))+0.5
tmp <- Score_matrix[1:miR_length,nTF_s_index:N]
Score_matrix[1:miR_length,nTF_s_index:N] <- (1-0.5)*(tmp-min(tmp[tmp>0]))/(max(tmp)-min(tmp[tmp>0]))+0.5
tmp <- Score_matrix[nTF_s_index:N,1:miR_length]
Score_matrix[nTF_s_index:N,1:miR_length] <- (1-0.5)*(tmp-min(tmp[tmp>0]))/(max(tmp)-min(tmp[tmp>0]))+0.5

tmp <- Score_matrix[TF_s_index:rg_length,TF_s_index:rg_length]
Score_matrix[TF_s_index:rg_length,TF_s_index:rg_length] <- (1-0.5)*(tmp-min(tmp[tmp>0]))/(max(tmp)-min(tmp[tmp>0]))+0.5
tmp <- Score_matrix[TF_s_index:rg_length,nTF_s_index:N]
Score_matrix[TF_s_index:rg_length,nTF_s_index:N] <- (1-0.5)*(tmp-min(tmp[tmp>0]))/(max(tmp)-min(tmp[tmp>0]))+0.5
tmp <- Score_matrix[nTF_s_index:N,TF_s_index:rg_length]
Score_matrix[nTF_s_index:N,TF_s_index:rg_length] <- (1-0.5)*(tmp-min(tmp[tmp>0]))/(max(tmp)-min(tmp[tmp>0]))+0.5

tmp <- Score_matrix[nTF_s_index:N,nTF_s_index:N]
Score_matrix[nTF_s_index:N,nTF_s_index:N] <- (1-0.5)*(tmp-min(tmp[tmp>0]))/(max(tmp)-min(tmp[tmp>0]))+0.5

#####scorematrix
Score_matrix[Score_matrix < 0.5] <- 0
colnames(Score_matrix) <- rownames(Score_matrix)

# Generate information of category of names
# Example:
# Name 			Category
# hsa-miR-15a	miR
# MSIS1			TF
# HAND1			nTF
#


#####gene_category
geneCa <- rbind(cbind(miR_list,rep("miR",length(miR_list))),cbind(TF_list,rep("TF",length(TF_list))),cbind(nTF_list,rep("nTF",length(nTF_list))))
colnames(geneCa) <- c("Name","Category")


# Write score matrix, mRNA-mRNA PCC and category of names to files
#write.table(Score_matrix,"Score_matrix_even_inter.txt",sep="\t",row.names=T,col.names=T)
#write.table(geneCa,"gene_category.txt",sep="\t",row.names=T,col.names=T)
#write.table(mR_mR_corr,"mR_mR_corr.txt",sep="\t",row.names=T,col.names=T)

##########
#gene_category <- geneCa
#Score_matrix_even_inter <- Score_matrix
##########


# Perform module identification, result is outputed to R data object: **.RData
# Three file required: 	Score matrix
#						PCC matrix of mRNA-mRNA pairs
#						Category of gene names
#
# 
# Suggested parameter setup: popN (population size): 100;
#							 Pco (crossover probability): 0.7;
#							 Pmu (mutation probability): 0.001;
#							 modN (# of modules generated): 10;
#							 K (GA iteration times): 5000;
#							 LK (Local search iteration times): 1000;
#							 Pnew (probability of random immigrant): 0.01
#							 miR_per (Option of prefered miRNA percentage): 0.02
#							 TF_per (Option of prefered TF-gene percentage): 0.04
#							 nTF_per (Option of prefered nTF-gene percentage): 0.02
#-----------------------------------------------------------------------------

####
#score_matrix
#Score_matrix_even_inter <- as.character(args[8])
#mR_mR_corr
#mR_mR_corr <- as.character(args[9])
#gene_category
#gene_category <- as.character(args[10])
####
########
# Options for control the size of regulators
#miR_per <- as.numeric(args[7])*100
#TF_per <- as.numeric(args[8])*100
# Initial population setup for GA
#nTF_per <- as.numeric(args[9])*100
##########

# loading GA functions
source("GA_1.2.R")

# set parameter
#popN <- 100
#Pco <- 0.7
#Pmu <- 0.001
#K <- 500
#modN <- 10
#Pnew <- 0.01
#LK <-1000

# Options for control the size of regulators
TF_per <- 0.04*100
miR_per <- 0.02*100
# Initial population setup for GA
nTF_per <- 0.02*100

#-----------------------------Data Import-----------------------------------------------
# Loading Score Matrix
# Loading PCC matrix
# Calculate # of miRNAs/TFs/nTFs
##########
#Score_matrix <- Score_matrix_even_inter
colnames(Score_matrix) <- rownames(Score_matrix)
mR_corr <- mR_mR_corr
colnames(mR_corr) <- rownames(mR_corr)
#geneCa <- gene_category
##########
miR_list <- as.character(geneCa[which(geneCa[,"Category"]=="miR"),"Name"])
TF_list <- as.character(geneCa[which(geneCa[,"Category"]=="TF"),"Name"])
nTF_list <- as.character(geneCa[which(geneCa[,"Category"]=="nTF"),"Name"])
miR_length <- length(miR_list)
TF_length <- length(TF_list)
nTF_length <- length(nTF_list)

N <- miR_length + TF_length + nTF_length
TF_s_index <- miR_length+1
nTF_s_index <- miR_length+TF_length+1
mR_length <- TF_length+nTF_length
rg_length <- miR_length+TF_length

#-----------------------------------------------------------------------------
# Formal Run of Algorithm
#-----------------------------------------------------------------------------

# Multiple core setup
registerDoMC(cores=10)
#final <- Module_Identification(Pmu,Pco,popN,modN,K,LK,Pnew,verb=TRUE)

#os = file(output,"w")
#lapply(final,write,os, append =F)

##############################################################################################################################
##############################################################################################################################
###############################################################
###############################################################
###############################################################
#-----------------------------------------------------------------------------
# For evaluation of GA and LS, the two component are actually run separately
#-----------------------------------------------------------------------------

verb <- TRUE
print("Algorithm Start")
Modules <- matrix(logical(),nrow=modN,ncol=N)
# Matrix for store GA solutions
ga_result_out <- matrix(logical(),nrow=modN,ncol=mR_length)
ga_scores <- list()
ls_scores <- list()
m <- 1
score_sub <- as.matrix(Score_matrix)
corr <- abs(as.matrix(mR_corr))
while(m <= modN)
{
  print(paste("Module",m))
  # Perform GA to identify co-expressed gene set
  ga_result_tmp <- GA(Pmu,Pco,popN,modN,K,Pnew,corr,verb,TRUE)
  ga_result <- ga_result_tmp$pop
  ga_scores[[m]] <- ga_result_tmp$record
  # Repair GA solution
  best_chr <- module_cleanup(ga_result[1,],corr)
  # In case there no nTF genes exists, regenerate new solutions
  if(length(which(best_chr[(TF_length+1):mR_length])) < 3)
  {
    print("No nTF genes exists. Regenerate.")
  }
  else
  {
    ga_result_out[m,] <- best_chr
    # Perform LS to detect best regulators
    solution_tmp <- Local_Search(best_chr,LK,score_sub,verb,TRUE)
    solution <- solution_tmp$best
    ls_scores[[m]] <- solution_tmp$record
    Modules[m,] <- solution
    # Update correlation matrix in case duplicated gene sets are identified by GA
    corr[best_chr,best_chr] <- matrix(0,nrow=length(which(best_chr)),ncol=length(which(best_chr)))
    m <- m+1
  }
} # module end

###############################################################
###############################################################
###############################################################
###############################################################
###############################################################

#-------------------------------------------------------------------------
# Functions used in evaluation
#-------------------------------------------------------------------------

unmatrix <- function (x)
{
  rnames <- rownames(x)
  cnames <- colnames(x)
  if (is.null(rnames))
    rnames <- paste("r", 1:nrow(x), sep = "")
  if (is.null(cnames))
    cnames <- paste("c", 1:ncol(x), sep = "")
  nmat <- outer(rnames, cnames, paste, sep = ":")
  vlist <- c(t(x))
  names(vlist) <- c(t(nmat))
  return(vlist)
}

strsplit2 <- function (x, split, ...)
{
  x <- as.character(x)
  n <- length(x)
  s <- strsplit(x, split = split, ...)
  nc <- unlist(lapply(s, length))
  out <- matrix("", n, max(nc))
  for (i in 1:n) {
    if (nc[i])
      out[i, 1:nc[i]] <- s[[i]]
  }
  out
}

# Score function for module evaluation
cal_score <- function(mat,score_sub)
{

  score_m <- score_sub[mat,mat]
  score <- sum(score_m)/length(which(score_m>0))
}

special_pie <- function(x,main="",anno=TRUE)
{
  slices <- as.numeric(x)
  if(anno)
  {
    lbls <- c("PCC&Binding","PCC only","Binding only")
    pct <- format(slices/sum(slices)*100,digit=2)
    lbls <- paste(lbls,pct,sep=":")
    lbls <- paste(lbls,"%",sep="")
    pie(slices,labels=lbls,main=main,col=c("lightblue", "mistyrose", "lightcyan"))
  }
  else
  {
    pie(slices,labels="",col=c("lightblue", "mistyrose", "lightcyan"))
    legend("topright",legend=c("PCC&Binding","PCC only","Binding only"),fill=c("lightblue", "mistyrose", "lightcyan"),cex=0.7,bty="n")
  }
}

#---------------------------------------------------------------------------------------#



mod_matrix <- Modules

# Write individual interactions to file
# Scores, PCCs and nodes of each interaction are output

tmp <- matrix(0,nrow=nTF_length,ncol=mR_length)
TF_gene_inter_sub <- rbind(TF_gene_inter,tmp)
rownames(TF_gene_inter_sub) <- colnames(TF_gene_inter_sub)

# Matrix to store summary of interactions
Module_stats <- matrix("",nrow=modN,ncol=6)
colnames(Module_stats) <- c("ID","Nodes","Inter","PCC&Binding","PCC only","Binding only")

module_all <- NULL
for(i in 1:10)
{
  geneCa_sub <- geneCa[mod_matrix[i,],]
  miR_sub <- as.character(geneCa_sub[geneCa_sub[,2]=="miR",1])
  TF_sub <- as.character(geneCa_sub[geneCa_sub[,2]=="TF",1])
  nTF_sub <- as.character(geneCa_sub[geneCa_sub[,2]=="nTF",1])
  mR_sub <- c(TF_sub,nTF_sub)


  mod_list <- c(unmatrix(miR_miR_corr_new[miR_sub,miR_sub]),unmatrix(miR_mR_corr_new[miR_sub,mR_sub]),unmatrix(mR_mR_corr_new[mR_sub,mR_sub]))
  score_sub <- c(unmatrix(Score_matrix[miR_sub,miR_sub]),unmatrix(Score_matrix[miR_sub,mR_sub]),unmatrix(Score_matrix[mR_sub,mR_sub]))
  seq_pre <- c(rep(0,length(miR_sub)^2),unmatrix(miR_gene_inter[miR_sub,mR_sub]),unmatrix(TF_gene_inter_sub[mR_sub,mR_sub]))

  mod_sub_names <- strsplit2(names(mod_list),split=":")

  mod_sub <- data.frame(mod_sub_names[,1],geneCa[match(mod_sub_names[,1],geneCa[,1]),2],mod_sub_names[,2],geneCa[match(mod_sub_names[,2],geneCa[,1]),2],as.numeric(mod_list),as.numeric(seq_pre),as.numeric(score_sub))
  mod_sub_out <- mod_sub[mod_sub[,7]!=0,]
  colnames(mod_sub_out) <- c("Node1","Type1","Node2","Type2","PCC","Seq_pre","Score")
  mod_sub_out <- mod_sub_out[!duplicated(mod_sub_out[,7]),]
  module_all <- rbind(module_all,cbind(rep(i,dim(mod_sub_out)[1]),mod_sub_out))
  pcc_inter <- ifelse(mod_sub_out[,5]== 0,0,ifelse(mod_sub_out[,5]<0,-1,1))
  seq_inter <- ifelse(mod_sub_out[,6]==0,0,1)
  mod4plot <- cbind(mod_sub_out[,1:4],pcc_inter,seq_inter)
  mod4plot <- mod4plot[mod4plot[,2] !="nTF" | mod4plot[,4]!="nTF",]
  # Each module are write to individual files for Cytoscape
  write.table(mod4plot,paste(output,"_","module_",i,".txt",sep=""),row.names=F,sep="\t",quote=F)
  Module_stats[i,] <- c(i,paste(count_category(mod_matrix[i,]),collapse="/"),dim(mod_sub_out)[1],length(which(mod_sub_out[,6]!=0 & mod_sub_out[,5]!=0)),length(which(mod_sub_out[,5]!=0 & mod_sub_out[,6]==0)),length(which(mod_sub_out[,6]!=0 & mod_sub_out[,5]==0)))
}
write.table(module_all,paste(output,"_","Modules.txt"),sep="\t",row.names=F,col.names=T,quote=F)
write.table(Module_stats,paste(output,"_","Module_stats.txt"),sep="\t",row.names=F,col.names=T,quote=F)

#---------------------------------------------------------------------------------------#


#---------------------------------------------------------------------------------------#
# generate random background
#---------------------------------------------------------------------------------------#
registerDoMC(cores=10)

n <- 1000
rand_node_score <- NULL
for(i in 1:dim(mod_matrix)[1])
{
  print(i)
  cat_count <- count_category(mod_matrix[i,])
  tmp <- foreach(i = 1:n, .combine="rbind", .multicombine=TRUE, .inorder=F) %do% c(sample(1:miR_length,cat_count[1],replace=F),sample(TF_s_index:(nTF_s_index-1),cat_count[2],replace=F),sample(nTF_s_index:N,cat_count[3],replace=F))
  tmp_score <- foreach(i = 1:n, .combine="c", .multicombine=TRUE, .inorder=T) %dopar% cal_score(tmp[i,],Score_matrix)
  rand_node_score <- c(rand_node_score,tmp_score)
}
base_scores <- foreach(i = 1:modN, .combine="c", .multicombine=TRUE, .inorder=T) %dopar% cal_score(mod_matrix[i,],Score_matrix)
base_scores <- apply(mod_matrix,1,cal_score,Score_matrix)



#---------------------------------------------------------------------------------------#
# permutation test for identified modules
#---------------------------------------------------------------------------------------#

n=10000
substitue_score <- matrix(0,nrow=n,ncol=modN)
for(j in 1:modN)
{
  print(j)
  mod <- mod_matrix[j,]
  for(i in 1:n)
  {
    ind <- which(mod)
    randN <- sample(1:length(ind),1)
    rand_loc <- sample(1:length(ind),randN)
    rep <- sample((1:N)[-1*ind],randN)
    ind[rand_loc] <- rep
    substitue_score[i,j] <- cal_score(sort(ind),Score_matrix)
  }
}

p_value <- vector("numeric",modN)
for(i in 1:modN)
{
  p_value[i] <- length(which(substitue_score[,i]>base_scores[i]))/(n+1)
}


#---------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------#
# Generate table which includes miRNA/gene names for using DAVID to perform GO term and
# KEGG pathway analysis
#---------------------------------------------------------------------------------------#
out <- matrix("",nrow=10,ncol=6)
out_4GO <- NULL
for (i in 1:modN)
{
  geneCa_sub <- geneCa[mod_matrix[i,],]
  miR_sub <- as.character(geneCa_sub[geneCa_sub[,2]=="miR",1])
  TF_sub <- as.character(geneCa_sub[geneCa_sub[,2]=="TF",1])
  nTF_sub <- as.character(geneCa_sub[geneCa_sub[,2]=="nTF",1])

  out[i,1] <- i
  out[i,2] <- paste(miR_sub,collapse=",")
  out[i,3] <- paste(TF_sub,collapse=",")
  out[i,4] <- paste(nTF_sub,collapse=",")
  out[i,5] <- format(base_scores[i],digits=4)
  out[i,6] <- format(p_value[i],scientific=T)
  tmp <- c(TF_sub,nTF_sub)
  out_4GO <- rbind(out_4GO,cbind(rep(i,length(tmp)),tmp))
}
colnames(out) <- c("ID","miRNA","TF","nTF","Score","P value")
write.table(out,paste(output,"_","Table_module_list.txt"),quote=F,sep="\t",row.names=F)

write.table(out_4GO,paste(output,"_","GO_test.txt"),quote=F,sep="\t",row.names=F,col.names=F)
#---------------------------------------------------------------------------------------#

#os = file(output,"w")
#lapply(final,write,os, append =F)


