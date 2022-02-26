
load("../yeast_table.rda")
load(file="../final_mod_gene_array.rda")
load(file="../steady_state_headers.rda")
load(file="../genelist.rda")
no_of_mod = length(steady_state_headers)
headers = colnames(yeast_table1)
nuc_pos = c(-4:-1,1:8)

no_of_mod = length(steady_state_headers)



for (m in 26:26) {
  mod_name = steady_state_headers[m]
  orig_mod_gene_val = 2**final_mod_gene_array[m,,]
  orig_mod_gene_val=na.omit(orig_mod_gene_val)
  mod_gene_val = colMeans(orig_mod_gene_val)
first_deriv_htz1 = vector(length = round(length(mod_gene_val)/1))
j=1
for(i in seq(from=2, to=length(mod_gene_val), by=1)) {
  # if(i == 3) {
  #   j = 2
  # }
  x_i = mod_gene_val[i]
  x_i1  = mod_gene_val[i-1]
  x_i2 = mod_gene_val[i-2]
  
  #first_deriv_htz1[j] = 0.5 *((3*x_i) - (4*x_i1) + x_i2)
  first_deriv_htz1[j] = 0.5 *((x_i) - (x_i1) )
  j = j+1
}

a = (sort(abs(first_deriv_htz1[abs(first_deriv_htz1)>0])))
a= a[1:1]
par(mfrow=c(1,length(a)))
for (h in 1:length(a)) {
  if(!is.na(a[h])) {
    index = which(first_deriv_htz1 == (1*a[h]))
    print(index)
    if(length(index) ==0) {
      index = which(first_deriv_htz1 == (-1*a[h]))
    }
    if(length(index) > 1) {
      for(m in 1:length(index)) {
        plot(nuc_pos[(index[m]-2):(index[m]+2)],mod_gene_val[(index[m]-2):(index[m]+2)],type="l",xlab="",ylab=steady_state_headers[m])
      }
    }
    else {
      plot(nuc_pos[(index-2):(index+2)],mod_gene_val[(index-2):(index+2)],type="l",xlab="",ylab=steady_state_headers[m])
      
    }
    print(index)
  }
}
}
