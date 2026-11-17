# Convert amino acid sequences into normalized numerical representations.
aa <- strsplit('ACDEFGHIKLMNPQRSTVWY','')[[1]]
# Each bigram has an integer index: (first - 1) * 20 + second.
pair_names <- as.vector(t(outer(aa, aa, paste0)))
feature_vector <- function(sequence) {
  v <- match(strsplit(sequence,'',fixed=TRUE)[[1]],aa)
  n <- length(v)
  stopifnot(n>=2, !anyNA(v))
  full <- tabulate(v,nbins=20)/n
  pairs <- (v[-n]-1L)*20L + v[-1L]
  bigrams <- tabulate(pairs,nbins=400)/(n-1L)
  # A 30-residue window is fixed before model evaluation, not tuned on test data.
  start <- head(v,30)
  nt <- tabulate(start,nbins=20)/length(start)
  c(setNames(full,paste0('AAC_',aa)),setNames(bigrams,paste0('BG_',pair_names)),setNames(nt,paste0('NT_',aa)))
}
# Sanity checks: AAC_A=2/3 and the AC/CA bigrams are equally frequent in ACA.
z <- feature_vector('ACA')
stopifnot(abs(z['AAC_A']-2/3)<1e-10,z['BG_AC']==0.5,z['BG_CA']==0.5,
          abs(sum(z[1:20])-1)<1e-10,abs(sum(z[21:420])-1)<1e-10)
x <- t(vapply(d$sequence,feature_vector,numeric(440)))
colnames(x) <- names(z)
stopifnot(nrow(x)==nrow(d), all(is.finite(x)), max(abs(rowSums(x[,1:20])-1))<1e-8)
feature_sets <- list(AAC=1:20, Bigrams=21:420, AAC_Nterm=c(1:20,421:440))
