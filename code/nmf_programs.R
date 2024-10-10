# -------------------------------------------------------------------------------------------
# Function for getting heterogeneity programs using nonnegative matrix factorization (NMF)
# ------------------------------------------------------------------------------------------- 

# - cpm = CPM expression matrix (rows = genes, columns = cells) # CPM表达矩阵；列-基因，行-细胞
# - is.log = indicates if the data is log transformed # 指示数据是否经过log转换
# - rank = NMF factorization rank #NMF分解秩
# - method = NMF algorithm # NMF算法

# Returns a list with NMF program scores for genes and cells #返回一个包含基因和细胞NMF程序得分的列表

library(NMF)

nmf_programs <- function(cpm, is.log=T, n=1, rank, method="snmf/r", seed=1) {
    
  if(is.log==F) CP100K_log <- log2((cpm/10) + 1) else CP100K_log <- cpm
  CP100K_log <- CP100K_log[apply(CP100K_log, 1, function(x) length(which(x > n)) > ncol(CP100K_log)*0.02),]
  CP100K_log <- CP100K_log - rowMeans(CP100K_log)
  CP100K_log[CP100K_log < 0] <- 0
 
  nmf_programs <- nmf(CP100K_log, rank=rank, method=method, seed=seed)
  
  nmf_programs_scores <- list(w_basis=basis(nmf_programs), h_coef=t(coef(nmf_programs)))
                                 
  return(nmf_programs_scores)
}
    
