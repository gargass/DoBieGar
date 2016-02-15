pvalue <- read.table('data/p_value_NA.txt')
pvalue$gen <- as.character(pvalue$gen)
class(pvalue$gen)

# pvalue$n_istotnych <- 0
# for(gen in 1:nrow(pvalue)){
# pvalue$n_istotnych[gen] <- sum(as.numeric(ifelse(pvalue[gen,]<0.05,1,0 )), na.rm = TRUE)
# }

pvalue_spis <- NULL
for(i in 1:nrow(pvalue)){
  gen <- pvalue$gen[i]
  for(j in 2:ncol(pvalue)){
    nowotwor <- colnames(pvalue)[j]
    pvalue_spis <- rbind(pvalue_spis, c(nowotwor, gen, pvalue[i,j]))
  }
}
colnames(pvalue_spis) <- c("nowotwor", "gen", "pvalue")
pvalue_spis <- data.frame(pvalue_spis)
pvalue_spis$pvalue <- as.numeric(as.character(pvalue_spis$pvalue))
head(pvalue_spis)
pvalue_spis_istotne <- pvalue_spis[!is.na(pvalue_spis$pvalue) & pvalue_spis$pvalue<0.05, ]

pvalue_spis_istotne$n_cancer <- 0
for(i in 1:nrow(pvalue_spis_istotne)){
  gen <- as.character(pvalue_spis_istotne$gen[i])
  pvalue_spis_istotne$n_cancer[i] <- sum((as.character(pvalue_spis_istotne$gen) == gen))
}

order.n <- order()
write.table(pvalue_spis_istotne, file="p_value_istotne.txt")

pvalue2 <- pvalue
pvalue2$n_cancer <- 0
for(i in 1:nrow(pvalue2)){
  pvalue2$n_cancer[i] <- sum(as.numeric(ifelse(!is.na(pvalue2[i,2:15]) & pvalue2[i,2:15]<0.05, 1, 0)))
}

order.n <- order(pvalue2$n_cancer, decreasing = TRUE)
pvalue2 <- pvalue2[order.n,]
write.table(pvalue2, file='pvalue_NA_n.txt')
