% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/motifs_scored_in_bins.R
\name{setup_motifs_scored_in_bins}
\alias{setup_motifs_scored_in_bins}
\title{Setup Motifs Scored in Bins}
\usage{
setup_motifs_scored_in_bins(
  bin_motifs = bin_motifs,
  motifs_scored = motifs_scored,
  contig_bins = contig_bins,
  assembly_stats = assembly_stats
)
}
\arguments{
\item{bin_motifs}{A dataframe of bin motifs.}

\item{motifs_scored}{A dataframe of scored motifs.}

\item{contig_bins}{A dataframe mapping contigs to bins.}

\item{assembly_stats}{A dataframe of assembly statistics.}
}
\value{
A dataframe of motifs scored in bins, with additional information from assembly statistics
and bin-contig mappings. Includes calculated mean methylation for each motif.
}
\description{
This function processes and combines various datasets to calculate the motif methylation degree
for all contigs, considering only methylations found in bins. It integrates motif scores with
assembly statistics and binning information, adjusting for unbinned cases.
}
\examples{
# Example usage:
setup_motifs_scored_in_bins(bin_motifs, motifs_scored, contig_bins, assembly_stats)

}
