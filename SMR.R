SMR<-function(gwas,eqtl,threshold,snp){
  
  mm<-data.frame()
  for (i in 1:length(snp)){
    index1=which(eqtl[,14]==snp[i])
    index2=which(gwas[,13]==snp[i])
    z_zx=gwas$z[index2]
    z_zy=eqtl$z[index1]
    gene=eqtl$gene_name[index1]
    for (j in 1:length(index1)){
      T=(z_zx^2*z_zy[j]^2)/(z_zx^2+z_zy[j]^2)
      P=pchisq(T,1,lower.tail=F)
    
      result=data.frame(snp[i],P,gene[j])
      mm<-rbind(mm,result)
    
    }
  
  }
  which(mm[,2]<=threshold)->nn
  gene_selected<-mm[nn,]
  gene<-mm
 
  jg<-list(All_genes=gene,Gene_selected=gene_selected) 
  return(jg)
  
}