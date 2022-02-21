library("entropy")
library("plotrix")
load(file="../../steady_state_headers1.rda")
no_of_mod = length(steady_state_headers_order1)
yeast_table1 = read.table('../../data/yeast_db.csv',
                          sep = ',',
                          header = TRUE)

yeast_table1 <- yeast_table1[-1,]
mutinfo_nuc = array(dim=c(no_of_mod,no_of_mod))
nuc_corr = array(dim=c(no_of_mod,no_of_mod))

# 
#     mod_m = 2**yeast_table1['Htz1']
#     plot(mod_m[1:65357,],type="h",xlab="Nucleosomes",ylab="Htz values",main="Htz1 levels across nucleosomes")


for ( i in 1: no_of_mod) {
  mod_i_name = steady_state_headers_order1[i]
  mod_i = yeast_table1[mod_i_name]
  mod_i[is.na(mod_i)] <- 0
  modi = mod_i[,]
  #n_modi = (modi - min(modi)) /(max(modi) - min(modi))
  for( j in 1:(no_of_mod)) {
    mod_j_name = steady_state_headers_order1[j]
    mod_j = yeast_table1[mod_j_name]
    mod_j[is.na(mod_j)] <- 0
    modj = mod_j[,]
    #n_modj = (modj - min(modj)) /(max(modj) - min(modj))
    

    dpij = discretize2d(modi,modj,5,5)
    mutinfo_nuc[i,j] = mi.plugin(dpij)
    
    nuc_corr[i,j] = cor(mod_i,mod_j)
    
    if(steady_state_headers_order1[i] == "H3K36me") {
      print(paste(steady_state_headers_order1[i],"-",steady_state_headers_order1[j]))
      print(paste("mutual information",mutinfo_nuc[i,j]))
      
      
      print(paste("correlation",nuc_corr[i,j]))
      print(" ")
      
    }

  }
  
}

# 
# testcol<-color.gradient(c(1,0),c(0,0),c(0,1),nslices=5)
# color2D.matplot(mutinfo_nuc,c(1,0),c(0,0),c(0,1),
#                 main="Nucleosome wise Mutual Information")
# col.labels<-c("High","Moderate","Low")
# color.legend(-1,0,-2,-6,col.labels,testcol,gradient = "Y")


#color2D.matplot(bty="L",mutinfo_nuc,extremes = c("lightskyblue","blue","deepskyblue","yellow","pink","red"))