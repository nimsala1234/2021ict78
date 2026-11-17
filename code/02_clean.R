# Keep unambiguous, single-location proteins with canonical amino acids.
# Multi-label proteins are biologically real, but outside this three-class task.
qc <- data.frame(step = character(), rows = integer())
record <- function(step, n) qc <<- rbind(qc, data.frame(step=step, rows=n))
record('Raw records', nrow(raw))
parse_label <- function(x) {
  x <- tolower(trimws(as.character(x)))
  out <- rep(NA_integer_, length(x))
  out[x %in% c('1','1.0','true','t')] <- 1L
  out[x %in% c('0','0.0','false','f')] <- 0L
  out
}
lab <- vapply(raw[labels], parse_label, integer(nrow(raw)))
valid <- complete.cases(lab)
if (!any(valid)) stop('No valid binary labels. Inspect the source CSV.')
record('Complete binary localization annotations', sum(valid))
single <- valid & rowSums(lab, na.rm=TRUE)==1 & rowSums(lab[,classes,drop=FALSE],na.rm=TRUE)==1
idx <- which(single)
record('Exactly one of the three target locations', length(idx))
d <- data.frame(accession=trimws(as.character(raw$ACC[idx])),
                sequence=toupper(gsub('[[:space:]]','',raw$Sequence[idx])),
                partition=suppressWarnings(as.integer(as.character(raw$Partition[idx]))),
                localization=classes[max.col(lab[idx,classes,drop=FALSE],ties.method='first')], stringsAsFactors=FALSE)
d <- d[!is.na(d$sequence) & !is.na(d$accession) & nzchar(d$accession) &
         grepl('^[ACDEFGHIKLMNPQRSTVWY]+$',d$sequence) & nchar(d$sequence)>40 & d$partition %in% 0:4, ]
record('Canonical sequences, valid accessions and partitions', nrow(d))
# Remove all exact sequence duplicates rather than risk cross-fold leakage or label conflicts.
dup_seq <- duplicated(d$sequence) | duplicated(d$sequence, fromLast=TRUE)
dup_acc <- duplicated(d$accession) | duplicated(d$accession, fromLast=TRUE)
d <- d[!(dup_seq | dup_acc), ]
record('Unique sequences and accessions', nrow(d))
d$localization <- factor(d$localization, levels=classes)
d$length <- nchar(d$sequence)
d <- d[order(d$partition,d$accession), ]
rownames(d) <- NULL
train <- which(d$partition %in% 0:3)
test <- which(d$partition == 4)
stopifnot(length(train)>0, length(test)>0, all(table(d$localization[train])>0), all(table(d$localization[test])>0),
          !any(d$sequence[train] %in% d$sequence[test]))
write.csv(qc,'results/quality_control.csv',row.names=FALSE)
write.csv(d,'data/clean_proteins.csv',row.names=FALSE)
split_counts <- as.data.frame.matrix(table(d$localization, factor(ifelse(d$partition==4,'Test','Train'),levels=c('Train','Test'))))
split_counts <- data.frame(Localization=rownames(split_counts),split_counts,row.names=NULL)
write.csv(split_counts,'results/split_counts.csv',row.names=FALSE)
