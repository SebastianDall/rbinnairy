

## TODO: Implement more sophisticated scoring system taking n_motifs into account.
## TODO: Deferentiate scoring system based on no methylations. If there is no bins with no methylation patterns, then remove all contigs with no methylation pattern. 



#' Calculate Belonging Score
#'
#' This function performs calculations to determine the belonging scores of motifs 
#' based on binary comparisons of methylation patterns in contigs and bins. It includes 
#' a scoring mechanisms considering the presence or absence of methylation 
#' and handles special cases based on the available data.
#'
#' @param motifs_scored_in_bins A dataframe of motifs scored in bins, containing motif 
#'        modification and bin information.
#' @param MEAN_METHYLATION_CUTOFF A numeric threshold for considering a motif as methylated. 
#'        Default is 0.25.
#' @param N_MOTIFS_CUTOFF A numeric threshold for the minimum number of motifs required for 
#'        consideration in the analysis. Default is 6.
#'
#' @return A dataframe containing the belonging scores for each bin and contig, along with 
#'         additional details such as the number of belonging bins and assembly length.
#'
#' @examples
#' # Assuming you have a properly formatted 'motifs_scored_in_bins' dataframe:
#' calculate_belonging_score(motifs_scored_in_bins)
#'
#' @export
calculate_belonging_score <- function(
  motifs_scored_in_bins = motifs_scored_in_bins,
  MEAN_METHYLATION_CUTOFF = 0.25,
  N_MOTIFS_CUTOFF = 6
) {
  
  # Convert bin motifs to binary
  bin_motif_binary <- motifs_scored_in_bins %>% 
    filter(bin != "unbinned") %>% 
    group_by(bin, motif_mod) %>% 
    convert_to_binary(MEAN_METHYLATION_CUTOFF) %>% 
    ungroup()
  
  # Convert contig methylation to binary
  contig_motif_binary <- motifs_scored_in_bins %>% 
    # Consider only motifs found in bins
    filter(motif_mod %in% bin_motif_binary$motif_mod %>% unique()) %>% 
    group_by(bin_contig, motif_mod) %>% 
    # Place a filter for number of motifs before. Few motifs can by noise have a huge impact
    filter(n_motifs > N_MOTIFS_CUTOFF) %>%
    # Convert to binary
    convert_to_binary(MEAN_METHYLATION_CUTOFF) %>% 
    rename(bin = bin_contig) %>% 
    # Add NA values
    pivot_wider(names_from = motif_mod, values_from = methylation_binary) %>% 
    pivot_longer(-bin, names_to = "motif_mod", values_to = "methylation_binary") %>% 
    # If there is no bin with no methylation pattern, then remove contigs with no methylation pattern
    find_bins_w_no_methylation(bin_motif_binary)
  
  
  # Compute scoring metrics:
  motif_binary_compare <- bin_motif_binary %>%
    left_join(
      contig_motif_binary  %>% rename(bin_compare = bin, methylation_binary_compare = methylation_binary),
      relationship = "many-to-many"
    ) %>%
    mutate(
      motif_comparison_score = case_when(
        methylation_binary == 1 & methylation_binary_compare == 1 ~ 1,    # Methylation on both contig and bin +1
        methylation_binary == 1 & methylation_binary_compare == 0 ~ -1,   # Methylation missing on contig -1
        methylation_binary == 0 & methylation_binary_compare == 1 ~ 0,    # Methylation only found on contig 0 (can be due to noise)
        methylation_binary == 0 & methylation_binary_compare == 0 ~ 0,    # Methylation missing on both contig and bin 0
        methylation_binary == 1 & is.na(methylation_binary_compare) ~ 0,  # Methylation missing on contig due to NA 0 (No penalty means benefit of the doubt)
        methylation_binary == 0 & is.na(methylation_binary_compare) ~ 0   # Methylation missing on contig due to NA 0 (No penalty means benefit of the doubt)
      )
    )
  
  
  
  belonging_score <- motif_binary_compare %>% 
    # Summarize the scores
    group_by(bin, bin_compare) %>% 
    summarise(n = sum(motif_comparison_score), .groups = "drop") %>%
    # Find the max per bin_compare
    group_by(bin_compare) %>%
    filter(n == max(n)) %>%
    mutate(belonging_bins = n()) %>% 
    separate(bin_compare, into = c("bin_id", "contig", "contig_number", "length"), sep = "_", remove = FALSE) %>% 
    mutate(
      contig = paste0("contig_", contig_number),
      length = as.numeric(length)
    ) %>% 
    select(!contig_number)
  
  return(belonging_score)
}



#' Convert Methylation Data to Binary Format
#'
#' This function processes a dataframe containing methylation data and converts the mean 
#' methylation values into a binary format based on a specified cutoff. Methylation values 
#' equal to or above the cutoff are converted to 1, and values below the cutoff are converted to 0.
#'
#' @param df A dataframe containing at least one column named 'mean', which holds the 
#'        mean methylation values.
#' @param MEAN_METHYLATION_CUTOFF A numeric threshold for considering a motif as methylated. 
#'        Methylation values equal to or above this cutoff are considered as methylated (1), 
#'        and values below are considered as unmethylated (0).
#'
#' @return A dataframe similar to the input but with an additional column 'methylation_binary', 
#'         indicating the binary methylation status (1 or 0) for each row.
#'
#' @examples
#' # Assuming 'df' is a dataframe with a 'mean' column containing methylation values:
#' binary_df <- convert_to_binary(df, MEAN_METHYLATION_CUTOFF)
#'
#'
#' @importFrom dplyr summarise mutate select if_else
convert_to_binary <- function(df, MEAN_METHYLATION_CUTOFF){
  df <- df %>% 
    summarise(
      mean_methylation = mean(mean)
    ) %>% 
    mutate(
      methylation_binary = if_else(mean_methylation >= MEAN_METHYLATION_CUTOFF, 1, 0)
    ) %>% 
    select(!mean_methylation)
  
  return(df)
}



#' Find Bins with No Methylation
#'
#' This function identifies bins with no methylation and optionally filters out these bins 
#' from the dataset. It checks for bins where the sum of methylation binary values is zero,
#' indicating no methylation, and then applies this information to filter the input dataframe.
#'
#' @param df A dataframe containing methylation data, including a 'bin' column for bin identifiers 
#'        and a 'methylation_binary' column for methylation status.
#' @param bin_motif_binary A dataframe containing binary methylation data by bin.
#'
#' @return A dataframe either identical to the input or filtered to exclude bins with no 
#'         methylation, depending on whether such bins exist.
#'
#' @examples
#' # Assuming 'df' is a dataframe with methylation data and 'bin_motif_binary' contains binary methylation data:
#' result_df <- find_bins_w_no_methylation(df, bin_motif_binary)
#'
#' @export
#' @importFrom dplyr group_by summarise filter
find_bins_w_no_methylation <- function(df, bin_motif_binary){
  no_methylation_bin_present <- bin_motif_binary %>% 
    group_by(bin) %>% 
    summarise(
      sum_methylation = sum(methylation_binary, na.rm = TRUE), .groups = "drop"
    ) %>% 
    filter(sum_methylation == 0)
  
  if (length(no_methylation_bin_present$bin) == 0) {
    bin_contigs_w_0_or_NA_only <- df %>% 
      group_by(bin) %>% 
      summarise(
        sum_methylation = sum(methylation_binary, na.rm = TRUE), .groups = "drop"
      ) %>% 
      filter(sum_methylation == 0)
    
    df <- df %>% 
      filter(!bin %in% bin_contigs_w_0_or_NA_only$bin)
    
    return(df)
    
  } else {
    return(df)
  }
}