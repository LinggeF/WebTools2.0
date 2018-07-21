###CMD inputing: Rscript Data_Loading.R 'HbA1c	Age	BMI	TG	GLU	GLU.0.5	GLU.1.0h	GLU.2.0h	CP.0h	CP.0.5h	CP.1.0h	CP.2.0h	INS.0h	INS.0.5h	INS.1.0h	INS.2.0h'

### This Rscript will return a dataframe that could be use for the performPred.R

args <- commandArgs(TRUE);

HbA1c <- as.numeric(args[1])
Age <- as.numeric(args[2])
BMI <- as.numeric(args[3])
TG <- as.numeric(args[4])
GLU <- as.numeric(args[5])
GLU.0.5 <- as.numeric(args[6])
GLU.1.0h <- as.numeric(args[7])
GLU.2.0h <- as.numeric(args[8])
CP.0h <- as.numeric(args[9])
CP.0.5h <- as.numeric(args[10])
CP.1.0h <- as.numeric(args[11])
CP.2.0h <- as.numeric(args[12])
INS.0h <- as.numeric(args[13])
INS.0.5h <- as.numeric(args[14])
INS.1.0h <- as.numeric(args[15])
INS.2.0h <- as.numeric(args[16])

x <- as.data.frame(cbind(HbA1c,	Age,	BMI,	TG,	GLU,	GLU.0.5,	GLU.1.0h,	GLU.2.0h,	CP.0h,	CP.0.5h,	CP.1.0h,	CP.2.0h,	INS.0h,	INS.0.5h,	INS.1.0h,	INS.2.0h))

print(x)




