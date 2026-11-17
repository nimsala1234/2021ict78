# Fixed comparison: same train/test split and 300 trees for every feature set.
# No tuning on test data. Predictions are probabilities from a ranger forest.
score <- function(truth,pred) {
  cm <- table(factor(truth,levels=classes),factor(pred,levels=classes))
  tp <- diag(cm); actual <- rowSums(cm); predicted <- colSums(cm)
  precision <- ifelse(predicted>0,tp/predicted,0)
  recall <- ifelse(actual>0,tp/actual,0)
  f1 <- ifelse(precision+recall>0,2*precision*recall/(precision+recall),0)
  list(cm=cm, accuracy=sum(tp)/sum(cm), macro_f1=mean(f1), balanced_accuracy=mean(recall),
       per_class=data.frame(Localization=classes,Precision=as.numeric(precision),Recall=as.numeric(recall),F1=as.numeric(f1)))
}
metrics <- data.frame(); per_class <- data.frame(); cms <- list(); predictions <- list()
majority <- names(which.max(table(d$localization[train])))
base_pred <- rep(majority,length(test))
s <- score(d$localization[test],base_pred)
metrics <- rbind(metrics,data.frame(Model='Majority',Accuracy=s$accuracy,Macro_F1=s$macro_f1,Balanced_accuracy=s$balanced_accuracy))
per_class <- rbind(per_class,cbind(Model='Majority',s$per_class)); cms[['Majority']] <- s$cm
for (name in names(feature_sets)) {
  message('Training ',name)
  cols <- feature_sets[[name]]
  fit <- ranger::ranger(x=as.data.frame(x[train,cols,drop=FALSE]),y=d$localization[train],
                        num.trees=300,probability=TRUE,importance=if(name=='AAC_Nterm') 'permutation' else 'none',
                        seed=4133,num.threads=2)
  prob <- predict(fit,data=as.data.frame(x[test,cols,drop=FALSE]),num.threads=2)$predictions
  prob <- prob[,classes,drop=FALSE]
  pred <- classes[max.col(prob,ties.method='first')]
  s <- score(d$localization[test],pred)
  metrics <- rbind(metrics,data.frame(Model=name,Accuracy=s$accuracy,Macro_F1=s$macro_f1,Balanced_accuracy=s$balanced_accuracy))
  per_class <- rbind(per_class,cbind(Model=name,s$per_class)); cms[[name]] <- s$cm
  predictions[[name]] <- data.frame(Model=name,Accession=d$accession[test],Actual=d$localization[test],Predicted=pred,
                                    Cytoplasm=prob[,1],Mitochondrion=prob[,2],Nucleus=prob[,3])
  saveRDS(list(fit=fit,columns=cols,classes=classes),file.path('results',paste0('model_',name,'.rds')))
  if(name=='AAC_Nterm') importance <- sort(fit$variable.importance,decreasing=TRUE)
}
write.csv(metrics,'results/metrics.csv',row.names=FALSE)
write.csv(per_class,'results/per_class.csv',row.names=FALSE)
write.csv(do.call(rbind,predictions),'results/test_predictions.csv',row.names=FALSE)
for(name in names(cms)) write.csv(as.data.frame.matrix(cms[[name]]),file.path('results',paste0('confusion_',name,'.csv')))
write.csv(data.frame(Feature=names(importance),Importance=as.numeric(importance)),'results/feature_importance.csv',row.names=FALSE)
# The primary biological comparison was chosen in advance: AAC_Nterm vs AAC.
open_plot('04_model_comparison.png')
barplot(t(as.matrix(metrics[,c('Accuracy','Macro_F1','Balanced_accuracy')])),beside=TRUE,
        names.arg=metrics$Model,col=c('#268B9C','#6858A7','#C56B35'),ylim=c(0,1.15),
        main='Held-out partition 4: model comparison',ylab='Score (0 to 1)',border=NA,
        legend.text=c('Accuracy','Macro F1','Balanced accuracy'),args.legend=list(x='topright',bty='n',cex=.8))
dev.off()
cm <- cms[['AAC_Nterm']]
open_plot('05_confusion_matrix.png',1100,900)
par(mar=c(6,8,4,2))
image(1:3,1:3,t(cm[3:1,]),col=colorRampPalette(c('#F1F5F9','#77ABD2','#164B77'))(30),axes=FALSE,
      main='AAC + N-terminal: held-out confusion matrix',xlab='Predicted localization',ylab='')
axis(1,at=1:3,labels=classes); axis(2,at=1:3,labels=rev(classes),las=1)
for(i in 1:3) for(j in 1:3) text(i,4-j,cm[j,i],col=if(cm[j,i]>max(cm)*.55) 'white' else '#162331',font=2)
box(); dev.off()
open_plot('06_feature_importance.png',1200,850)
par(mar=c(5,8,4,2)); top <- rev(head(importance,12))
barplot(top,horiz=TRUE,col=ifelse(grepl('^NT_',names(top)),'#6858A7','#268B9C'),las=1,
        main='AAC + N-terminal: permutation importance',xlab='Permutation importance (OOB prediction error)',border=NA)
abline(v=0,col='grey40'); dev.off()
writeLines(capture.output(sessionInfo()),'results/session_info.txt')
