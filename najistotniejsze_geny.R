pvalue <-read.table('p_value/p_value_NA.txt', h=T)

for (i in 1:nrow(pvalue))
{
  for (j in 1:ncol(pvalue))
  {
    if (is.na(pvalue[i, j])==TRUE)
    {
      pvalue[i, j]=1
    }
  }
}


cz = matrix(nrow=nrow(pvalue), ncol=2)
cz[, 1]=as.character(pvalue[, 1])
for (i in 1:nrow(pvalue))
{
  cz[i, 2]=sum(pvalue[i, 2:ncol(pvalue)])
}
colnames(cz)=c('gen', 'pvalue')
cz[,2]= as.numeric(cz[, 2])
cz[order(cz[, 2]),]

cz2 = cz[1:50,]

write.table(cz2[,1], file = "najistotniejsze_geny.txt")