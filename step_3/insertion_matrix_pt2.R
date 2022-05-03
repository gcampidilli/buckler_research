
wd = "/Users/gcampidilli/Documents/research/mu_tir/step_4"
setwd(wd)

# construct insertion matrix

inbred_list = read.delim("inbred_list.txt", header = F)
total_df = read.csv("prematrix_df/total_prematrix_df.csv", header = T)
total_df = total_df[-c(1),-c(1)]
m_rows = inbred_list$V1

# order data by chromosome and bp position
total_df = total_df[order(total_df$chrom, total_df$pos_thousand),]
m_cols = unique(total_df$id)

insertion_matrix = matrix(0, nrow = length(m_rows), ncol = length(unique(total_df$id)), dimnames = list(m_rows, m_cols))


for(x in 1:length(inbred_list[,1])) {
  inbred_df = total_df[total_df$inbred == inbred_list[x,1],]
  col_index = unlist(lapply(1:length(inbred_df$id), function(x){
    grep(paste("^",inbred_df$id[x],"$",sep=""), colnames(insertion_matrix))
  }))
  row_index = grep(paste("^",inbred_list[x,1],"$",sep=""), row.names(insertion_matrix))
  insertion_matrix[row_index,col_index] = 1
}

write.csv(insertion_matrix, file = "jan21_insertion_matrix.csv")

# exact = read.csv("jan21_insertion_matrix_exact.csv")
# hundred = read.csv("jan21_insertion_matrix_hundred.csv")
# thousand = read.csv("jan21_insertion_matrix_thousand.csv")



insertion_matrix = matrix(0, nrow = length(m_rows), ncol = length(unique(total_df$id)), dimnames = list(m_rows, m_cols))


for(x in 1:length(inbred_list[,1])) {
  inbred_df = total_df[total_df$inbred == inbred_list[x,1],]
  col_index = unlist(lapply(1:length(inbred_df$id), function(x){
    grep(paste("^",inbred_df$id[x],"$",sep=""), colnames(insertion_matrix))
  }))
  row_index = grep(paste("^",inbred_list[x,1],"$",sep=""), row.names(insertion_matrix))
  insertion_matrix[row_index,col_index] = 1
}



