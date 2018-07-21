#!/usr/bin/env Rscript

#Created by Lingge Feng
#July 18th, 2018 
#Hua Medicine Confidential 

args = commandArgs(trailingOnly=TRUE)
### input data format: HbA1c, Age, BMI, TG, GLU(0-120min),INS(0-120min), Cp(0-120min)
### .csv format
####################
library(kernlab)
library(sampling)
library(caret)
library(Bolstad2)

load("./trained_svm.R.RData")

x <- read.csv(args[1],sep = ",", head = TRUE)
# test if there is at least one argument: if not, return an error
if (ncol(x) < 16) {
  stop("Not enough input argument ", call.=FALSE)
}
# test if there is missing value
if (nrow(na.omit(x)) != nrow(x)){
  stop("Containing missing values ", call.=FALSE)
}

tmp.org <- AUCset(x) 
tmp.sub <- tmp.org[,c("HbA1c","Age","BMI","TG","GLU","CP.0h","INS.0h","HI.BI","auc.cp","auc.glu","auc.ins","auc.g.i")]

###93% accuracy rate
TDClass(tmp.sub)

# A: LoIR_BFD 24
# B: LoBFD    71
# C: HiBFD    26
# D: LoIR     44
# E: IR       18
# F: HIR_BFD  38

