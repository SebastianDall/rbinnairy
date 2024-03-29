% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/contig_belonging_score.R
\name{convert_to_binary}
\alias{convert_to_binary}
\title{Convert Methylation Data to Binary Format}
\usage{
convert_to_binary(df, MEAN_METHYLATION_CUTOFF)
}
\arguments{
\item{df}{A dataframe containing at least one column named 'mean', which holds the
mean methylation values.}

\item{MEAN_METHYLATION_CUTOFF}{A numeric threshold for considering a motif as methylated.
Methylation values equal to or above this cutoff are considered as methylated (1),
and values below are considered as unmethylated (0).}
}
\value{
A dataframe similar to the input but with an additional column 'methylation_binary',
indicating the binary methylation status (1 or 0) for each row.
}
\description{
This function processes a dataframe containing methylation data and converts the mean
methylation values into a binary format based on a specified cutoff. Methylation values
equal to or above the cutoff are converted to 1, and values below the cutoff are converted to 0.
}
\examples{
# Assuming 'df' is a dataframe with a 'mean' column containing methylation values:
binary_df <- convert_to_binary(df, MEAN_METHYLATION_CUTOFF)


}
