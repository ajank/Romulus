#' Annotations of 4,828 NRSF (REST) motif instances.
#'
#' A dataset containing the annotations of 4,828 NRSF (REST) motif instances in the human genome (hg19 assembly).
#'
#' @format A data frame with 4828 rows and 6 variables:
#' \describe{
#'   \item{chrom}{chromosome}
#'   \item{start}{starting base pair (1-based, inclusive)}
#'   \item{end}{ending base pair}
#'   \item{strand}{strand (\code{"-"} or \code{"+"})}
#'   \item{score}{Position Weight Matrix score (log-likelihood)}
#'   \item{signalValue}{ChIP-seq signal value in K562 cells}
#' }
#' @source Motif instances are downloaded from \url{http://homer.salk.edu/homer/} (HOMER Known Motifs track).
#'
#' ChIP-seq signal value was extracted from \url{http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/wgEncodeAwgTfbsHaibK562NrsfV0416102UniPk.narrowPeak.gz}.
"NRSF.anno"

#' DNase I cuts around 4,828 NRSF (REST) motif instances.
#'
#' A dataset containing the exact numbers of DNase I cuts around 4,828 NRSF (REST) motif instances, split into forward and reverse strand cuts.
#'
#' @format An integer matrix with 4828 rows and 838 columns. Columns 1-419 correspond to forward strand cuts (200 bp upstream + 19 bp motif site + 200 bp downstream), while columns 420-838 correspond to reverse strand cuts at the same positions.
#'
#' @source Extracted from \url{http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeOpenChromDnase/wgEncodeOpenChromDnaseK562AlnRep1V2.bam}.
"NRSF.cuts"

#' Number of base pairs of margin for \code{NRSF.cuts}.
#'
#' Number of base pairs of upstream and downstream margin for \code{NRSF.cuts}.
#'
#' @format Integer, equal to 200 (base pairs).
"NRSF.margin"
