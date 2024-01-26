


#### Setup 

#' Setup Motifs Scored in Bins
#'
#' This function processes and combines various datasets to calculate the motif methylation degree
#' for all contigs, considering only methylations found in bins. It integrates motif scores with 
#' assembly statistics and binning information, adjusting for unbinned cases.
#'
#' @param bin_motifs A dataframe of bin motifs.
#' @param motifs_scored A dataframe of scored motifs.
#' @param contig_bins A dataframe mapping contigs to bins.
#' @param assembly_stats A dataframe of assembly statistics.
#'
#' @return A dataframe of motifs scored in bins, with additional information from assembly statistics
#' and bin-contig mappings. Includes calculated mean methylation for each motif.
#'
#' @examples
#' # Example usage:
#' setup_motifs_scored_in_bins(bin_motifs, motifs_scored, contig_bins, assembly_stats)
#'
#' @export
setup_motifs_scored_in_bins <- function (
  bin_motifs = bin_motifs,
  motifs_scored = motifs_scored,
  contig_bins = contig_bins,
  assembly_stats = assembly_stats
) {
  
  # Get motifs associated with bins. Only methylations found in bins are considered as they are the basis for scoring contamination and binning
  motifs_in_bins <- bin_motifs %>%
    mutate(motif_mod = paste0(motif,"_",mod_type)) %>% 
    pull(motif_mod) %>%
    unique()
  
  # Calculate motif methylation degree for all contigs.
  motifs_scored_in_bins <- motifs_scored %>%
    mutate(motif_mod = paste0(motif,"_",mod_type)) %>% 
    
    # Only consider motifs scored in bins
    filter(motif_mod %in% motifs_in_bins) %>%
    
    # Add assembly stats
    left_join(contig_bins, by = "contig") %>%
    left_join(assembly_stats %>% select(contig, length), by = "contig") %>%
    
    # If bin is NA it is considered as unbinned
    mutate(
      bin = if_else(is.na(bin), "unbinned", bin),
      bin_contig = paste(bin, contig, length, sep = "_")
    ) %>%
    
    # Calculate mean
    mutate(
      n_motifs = n_mod + n_nomod
    ) %>% 
    mutate(
      mean = n_mod / n_motifs
    )
  
  return(motifs_scored_in_bins)
  
}




