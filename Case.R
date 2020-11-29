library(data.table)
gwas=fread("...\\K05_gwas.txt",header=T)
gwas=gwas[which(gwas$pval<0.05),]
gwas$z=as.numeric(gwas$beta)/as.numeric(gwas$se)

eqtl=fread("...\\blood_eqtl.txt",data.table=F)
eqtl$z=eqtl$slope/eqtl$slope_se
threshold1=0.05/length(unique(eqtl$gene_id))   ###阈值
snp1=intersect(eqtl$rs_number,gwas$rs)

Result<-SMR(gwas,eqtl,threshold1,snp1)
