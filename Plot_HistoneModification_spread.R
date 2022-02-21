

#####################
#PLOT FROM 3D MATRIX
#####################


load("../../final_mod_gene_array.rda")
load("../../steady_state_headers.rda")

h3k14ac = final_mod_gene_array[3,,]
h3k14ac= na.omit(h3k14ac)
h3k14ac_gene_avg = colMeans(h3k14ac)

nuc = c(-4:-1,1:8)
plot(nuc,h3k14ac_gene_avg,type="l")


###############################
#PLOT FROM BP (RAW DATA)
###############################

yeast_table1 = read.table('../../data/yeast_db.csv',
                          sep = ',',
                          header = TRUE)
library(sqldf)
library(rapportools)
library(hashmap)
yeast_table1 <- yeast_table1[-1,]
yeast_table_gene <-
  yeast_table1[!(is.empty(as.character(yeast_table1[, "gene"]))), ]
mod_data <- ( yeast_table_gene[, "H3K14ac"])

df  = data.frame(
  position = yeast_table_gene[, "center"],
  chromosome = yeast_table_gene[, "chr"],
  level = mod_data,
  gene = yeast_table_gene[, "gene"]
)
df1 = na.omit(df)

