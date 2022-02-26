
load("../yeast_table.rda")
load(file="../final_mod_gene_array.rda")
load(file="../steady_state_headers.rda")
load(file="../genelist.rda")
no_of_mod = length(steady_state_headers)
headers = colnames(yeast_table1)

no_of_mod = length(steady_state_headers)
mod_data <- 2**(yeast_table1[, ])
mod_data[is.na(mod_data)] <- 0

mod_doi = mod_data[steady_state_headers]

htz1 = mod_doi[1:66357,"H3K4me3"]

first_deriv_htz1 = vector(length = round(length(htz1)/1))

for(i in seq(from=3, to=length(htz1), by=1)) {
  if(i == 3) {
    j = 2
  }
  x_i = htz1[i]
  x_i1  = htz1[i-1]
  x_i2 = htz1[i-2]
  
  first_deriv_htz1[j] = 0.5 *((3*x_i) - (4*x_i1) + x_i2)
  j = j+1
}
#first_deriv_htz1 = diff(htz1)a = 
#plot(first_deriv_htz1,type="l")
a = (sort(abs(first_deriv_htz1[abs(first_deriv_htz1)>0])))
a = a[1:20]
par(mfrow=c(2,5))
for (h in 1:20) {
  if(!is.na(a[h])) {
  index = which(first_deriv_htz1 == (a[h]))
  if(length(index) ==0) {
    index = which(first_deriv_htz1 == (-1*a[h]))
  }
  if(length(index) > 1) {
    for(m in 1:length(index)) {
    plot(c((index[m]-5):(index[m]+5)),htz1[(index[m]-5):(index[m]+5)],type="l",xlab="",ylab="H3K4me3")
    }
  }
    else {
      plot(c((index-5):(index+5)),htz1[(index-5):(index+5)],type="l",xlab="",ylab="H3K4me3")
      
    }
  
  }
}
