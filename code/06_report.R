# Build the README from actual calculations, never invented results.
# The report is regenerated on each run; edit the template for durable changes.
md_table <- function(df,digits=3) {
  out <- lapply(df,function(x) if(is.numeric(x)) format(round(x,digits),trim=TRUE,scientific=FALSE) else as.character(x))
  out <- as.data.frame(out,stringsAsFactors=FALSE)
  header <- paste0('| ',paste(names(out),collapse=' | '),' |')
  rule <- paste0('| ',paste(rep('---',ncol(out)),collapse=' | '),' |')
  rows <- apply(out,1,function(r) paste0('| ',paste(r,collapse=' | '),' |'))
  paste(c(header,rule,rows),collapse='\n')
}
get_metric <- function(model,metric='Macro_F1') metrics[metrics$Model==model,metric]
base <- get_metric('Majority'); aac <- get_metric('AAC'); combined <- get_metric('AAC_Nterm'); bg <- get_metric('Bigrams')
ml <- metrics[metrics$Model!='Majority',]
best <- paste(ml$Model[ml$Macro_F1==max(ml$Macro_F1)],collapse=' and ')
delta <- combined-aac
fmt <- function(x) sprintf('%.3f',x)
comparison <- if(delta>0) 'higher' else if(delta<0) 'lower' else 'equal'
cm <- cms[['AAC_Nterm']]; off <- cm; diag(off) <- 0
pos <- which(off==max(off),arr.ind=TRUE)
error_text <- if(max(off)==0) 'No off-diagonal errors were observed.' else paste0('The largest off-diagonal count was ',max(off),': ',paste(apply(pos,1,function(p) paste0(rownames(cm)[p[1]],' predicted as ',colnames(cm)[p[2]])),collapse='; '),'.')
nt_span <- apply(nt_means,1,function(v) max(v)-min(v))
nt_residue <- names(which.max(nt_span))
pc <- per_class[per_class$Model=='AAC_Nterm',c('Localization','Precision','Recall','F1')]
replacements <- list(
 ABSTRACT=sprintf('Among %s retained proteins, the highest held-out macro F1 among the forests was %s (%s). The combined representation achieved %s.',nrow(d),fmt(max(ml$Macro_F1)),best,fmt(combined)),
 QC=paste0('**Table 1. Sequential data-quality checks.**\n\n',md_table(qc),'\n\nAfter filtering, ',nrow(d),' proteins remained from ',nrow(raw),' original records.'),
 SPLIT=paste0('**Table 2. Training and test class counts.**\n\n',md_table(split_counts),'\n\nThe training set contained ',length(train),' proteins and the test set contained ',length(test),'.'),
 LENGTH=paste0('**Table 3. Training sequence length summary (amino acids).**\n\n',md_table(length_summary)),
 EDA=sprintf('The largest training class was %s. The residue with the largest range in mean N-terminal percentage across classes was %s (%.2f percentage points). This is descriptive and not a test of statistical significance.',majority,nt_residue,max(nt_span)),
 METRICS=paste0('**Table 4. Held-out classification metrics (0–1).**\n\n',md_table(metrics)),
 PERFORMANCE=sprintf('The highest observed forest macro F1 was %s (%s). AAC scored %s, bigrams %s and AAC + N-terminal %s, compared with the majority baseline of %s.',fmt(max(ml$Macro_F1)),best,fmt(aac),fmt(bg),fmt(combined),fmt(base)),
 ERRORS=paste0(error_text,' The combined model correctly classified ',sum(diag(cm)),' of ',sum(cm),' held-out proteins.'),
 PERCLASS=paste0('**Table 5. Class-specific combined-model performance.**\n\n',md_table(pc)),
 IMPORTANCE=paste0('The highest-ranked feature was `',names(importance)[1],'` (permutation importance ',fmt(importance[1]),'). This ranking is exploratory.'),
 DISCUSSION=sprintf('The combined model had %s macro F1 than AAC alone (difference: %+.3f). %s %s',comparison,delta,if(delta>0) 'The observation is consistent with the primary hypothesis on this split.' else 'The observation does not support the predicted improvement on this split.',if(max(ml$Macro_F1)>base) 'At least one sequence model outperformed the majority baseline in macro F1.' else 'No sequence model outperformed the majority baseline in macro F1.'),
 CONCLUSION=sprintf('In the held-out partition, the best observed forest representation was %s with macro F1 %s. N-terminal augmentation changed macro F1 by %+.3f relative to AAC. The results should be interpreted in the context of the limited three-class dataset.',best,fmt(max(ml$Macro_F1)),delta),
 VERSIONS=paste0('**Run environment:** ',R.version.string,'; ranger ',as.character(packageVersion('ranger')),'.')
)
report <- paste(readLines('report/README_template.md',warn=FALSE,encoding='UTF-8'),collapse='\n')
for(k in names(replacements)) report <- gsub(paste0('@@',k,'@@'),replacements[[k]],report,fixed=TRUE)
if(grepl('@@[A-Z]+@@',report)) stop('Unresolved report placeholder.')
writeLines(enc2utf8(report),'README.md',useBytes=TRUE)
message('Generated README.md from actual data and model results.')
