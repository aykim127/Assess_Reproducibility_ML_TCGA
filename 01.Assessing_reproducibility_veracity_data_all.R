# load required packages
require(TCGA2STAT)

##########################################################################################
# 0. Load data RNASeq RPKM data 
##########################################################################################

blca0 <- getTCGA(disease="BLCA", data.type="RNASeq", type="RPKM", clinical=TRUE)
brca0 <- getTCGA(disease="BRCA", data.type="RNASeq", type="RPKM", clinical=TRUE)
coad0 <- getTCGA(disease="COAD", data.type="RNASeq", type="RPKM", clinical=TRUE)
coadread0 <- getTCGA(disease="COADREAD", data.type="RNASeq", type="RPKM", clinical=TRUE)
esca0 <- getTCGA(disease="ESCA", data.type="RNASeq", type="RPKM", clinical=TRUE)
hnsc0 <- getTCGA(disease="HNSC", data.type="RNASeq", type="RPKM", clinical=TRUE)
kipan0 <- getTCGA(disease="KIPAN", data.type="RNASeq", type="RPKM", clinical=TRUE)
kirc0 <- getTCGA(disease="KIRC", data.type="RNASeq", type="RPKM", clinical=TRUE)
kirp0 <- getTCGA(disease="KIRP", data.type="RNASeq", type="RPKM", clinical=TRUE)
laml0 <- getTCGA(disease="LAML", data.type="RNASeq", type="RPKM", clinical=TRUE)
lihc0 <- getTCGA(disease="LIHC", data.type="RNASeq", type="RPKM", clinical=TRUE)
luad0 <- getTCGA(disease="LUAD", data.type="RNASeq", type="RPKM", clinical=TRUE)
lusc0 <- getTCGA(disease="LUSC", data.type="RNASeq", type="RPKM", clinical=TRUE)
ov0 <- getTCGA(disease="OV", data.type="RNASeq", type="RPKM", clinical=TRUE)
read0 <- getTCGA(disease="READ", data.type="RNASeq", type="RPKM", clinical=TRUE)
stad0 <- getTCGA(disease="STAD", data.type="RNASeq", type="RPKM", clinical=TRUE)
thca0 <- getTCGA(disease="THCA", data.type="RNASeq", type="RPKM", clinical=TRUE)
ucec0 <- getTCGA(disease="UCEC", data.type="RNASeq", type="RPKM", clinical=TRUE) 

### You will need to change the directory to save all the data files
save.image("/extra/akim127/Project1/TCGA_all.RData")
