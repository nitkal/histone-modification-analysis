library("entropy")
load("../../yeast_table.rda")
load(file="../../final_mod_gene_array.rda")
load(file="../../steady_state_headers.rda")
load(file="../../genelist.rda")
no_of_mod = length(steady_state_headers)
mutinfo_gene = array(dim=c(no_of_mod,no_of_mod))

yeast_table1 = read.table('../../data/yeast_db.csv',
                          sep = ',',
                          header = TRUE)

yeast_table1 <- yeast_table1[-1,]

nuc_corr = array(dim=c(no_of_mod,no_of_mod))
  
for ( i in 1: no_of_mod) {
  mod_i = final_mod_gene_array[i,,]
  mod_i[is.na(mod_i)] <- 0
  
  mod_i_name = steady_state_headers[i]
  
  modi = yeast_table1[mod_i_name]
  modi[is.na(modi)] <- 0
  modi = modi[,]
  #pci = (svd(mod_i))$v[,1]
 # npci = (pci - min(pci))/(max(pci)-min(pci))

  for( j in 1:no_of_mod) {
    mod_j_name = steady_state_headers[j]
    
    modj = yeast_table1[mod_j_name]
    modj[is.na(modj)] <- 0
    modj = modj[,]
    nuc_corr[i,j] = cor(modi,modj)
    
    mod_j = final_mod_gene_array[j,,]
    mod_j[is.na(mod_j)] <- 0
  #  pcj = (svd(mod_j))$v[,1]
  #  npcj = (pcj - min(pcj))/(max(pcj)-min(pcj))
    
    d_pcij = discretize2d(mod_i,mod_j,10,10)
    mutinfo_gene[i,j] = mi.plugin(d_pcij)
    
    if( mutinfo_gene[i,j] >=0.2) {
      if ((nuc_corr[i,j]>= -0.5) && (nuc_corr[i,j] <= 0.3)) {
        print(paste(steady_state_headers[i],"-",steady_state_headers[j],"-",round(mutinfo_gene[i,j],3),"-",round(nuc_corr[i,j],3)))
      }
      
    }
  }
  
}
#mutualinfo_gene = mutinfo_gene
#write.csv(mutualinfo_gene,file="mutualinfo_gene.csv")
# testcol<-color.gradient(c(1,0),c(0,0),c(0,1),nslices=5)
 #color2D.matplot(mutualinfo_gene,c(1,0),c(0,0),c(0,1),
  #               main="Gene wise Mutual Information")
# col.labels<-c("High","Moderate","Low")
# color.legend(-1,0,-2,-6,col.labels,testcol,gradient = "Y")


color2D.matplot(bty="L",mutinfo_nuc,extremes = c("lightskyblue","blue","deepskyblue","yellow","pink","red"))
