# Training-only EDA keeps held-out protein patterns out of model development.
palette_loc <- c(Cytoplasm='#268B9C',Mitochondrion='#C56B35',Nucleus='#6858A7')
open_plot <- function(name,w=1200,h=800) {
  png(file.path('figures',name),width=w,height=h,res=150)
  par(mar=c(5,5,4,2)+0.1,las=1,cex.axis=0.85)
}
open_plot('01_class_balance.png')
barplot(table(d$localization[train]),col=palette_loc,ylim=c(0,max(table(d$localization[train]))*1.18),
        main='Training proteins by localization',ylab='Number of proteins',border=NA)
dev.off()
open_plot('02_sequence_length.png')
boxplot(log10(length)~localization,data=d[train,],col=palette_loc,outline=FALSE,
        main='Protein length by localization (training)',ylab='log10(sequence length in amino acids)',xlab='')
dev.off()
# Each row is an amino acid and each column is a class. Values are mean percentages.
nt_means <- sapply(classes,function(cl) colMeans(x[train[d$localization[train]==cl],421:440,drop=FALSE])*100)
rownames(nt_means) <- aa
open_plot('03_nterminal_composition.png',1200,1000)
par(mar=c(6,4,4,3))
image(1:3,1:20,t(nt_means),col=colorRampPalette(c('#F2F7F5','#79C8BB','#134C72'))(40),axes=FALSE,
      xlab='',ylab='Amino acid',main='Mean N-terminal amino acid percentage (training)')
axis(1,at=1:3,labels=classes); axis(2,at=1:20,labels=aa,las=1)
for(i in 1:3) for(j in 1:20) text(i,j,sprintf('%.1f',nt_means[j,i]),cex=.72,
                                    col=if(nt_means[j,i]>mean(range(nt_means))) 'white' else '#162331')
box(); dev.off()
length_summary <- do.call(rbind,lapply(classes,function(cl) {
  v <- d$length[train[d$localization[train]==cl]]
  data.frame(Localization=cl,N=length(v),Median=median(v),Q1=unname(quantile(v,.25)),Q3=unname(quantile(v,.75)))
}))
write.csv(length_summary,'results/length_summary.csv',row.names=FALSE)
