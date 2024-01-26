

#' Determine Contamination Categories
#'
#' This function categorizes data into different contamination groups based on belonging scores.
#' It classifies each item as:
#' - correct: Methylation pattern fits the assigned bin
#' - single_contamination: Contig has a best match based on methylation pattern
#' - multi_contamination: Contig can be assigned to multiple bins based on methylation pattern. 
#'
#'
#' @param belonging_score A dataframe containing the belonging scores for each bin and contig,
#'        as well as additional details such as the number of belonging bins.
#'
#' @return A dataframe similar to the input but with an additional column named 'group', which 
#'         indicates the contamination category for each row.
#'
#' @examples
#' # Assuming 'belonging_score' is a dataframe with the appropriate structure:
#' contamination_df <- determine_contamination(belonging_score)
#'
#' @export
determine_contamination <- function(belonging_score) {
  contamination_df <- belonging_score %>% 
    mutate(
      group = case_when(
        bin == bin_id & belonging_bins == 1 ~ "correct",
        bin != bin_id & belonging_bins == 1 ~ "single_contamination",
        belonging_bins > 1 ~ "multi_contamination"
      )
    )
  
  return(contamination_df)
}

#' Convert Distance Matrix to Long Format
#'
#' This function processes a distance matrix, converting it into a long format dataframe. 
#' It isolates contig distances, rounds values, and summarizes distances by contig, 
#' concatenating them into a single string.
#'
#' @param dist_matrix A distance matrix to be processed.
#'
#' @return A dataframe with two columns: 'contig' and 'matching_bins'. 
#'         The 'matching_bins' column contains concatenated distance values for each contig.
#'
#' @examples
#' # Assuming 'dist_matrix' is a distance matrix:
#' long_format_df <- extract_kmer_dist(dist_matrix)
#'
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr filter select group_by summarize
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_detect
extract_kmer_dist <- function(dist_matrix) {
  dist_matrix %>% 
    # Convert to dataframe
    as.matrix() %>% 
    as.data.frame() %>% 
    # Isolate contig distance
    rownames_to_column("contig") %>% 
    filter(str_detect(contig, "contig")) %>% 
    select(!contains("contig_")) %>% 
    pivot_longer(-contig) %>% 
    # Summarise distance
    mutate(
      value = round(value, 5),
      matching_bins = paste0(name,"|",value)
    ) %>% 
    select(contig,matching_bins) %>% 
    group_by(contig) %>%
    summarize(matching_bins = paste(matching_bins, collapse = ";"), .groups = "drop")
}


#' Calculate K-mer Frequency
#'
#' This function calculates the k-mer frequency in sequences extracted from a FASTA file. 
#' It processes a given dataframe to include k-mer counts for each sequence in the FASTA file 
#' based on a specified k-mer window size.
#'
#' @param df A dataframe with a column 'contig' that corresponds to sequence identifiers in the FASTA file.
#' @param fasta A list or other suitable structure containing sequences from a FASTA file, 
#'        indexed by sequence identifiers.
#' @param kmer_window An integer specifying the k-mer window size for frequency calculation.
#'
#' @return A dataframe similar to the input, but with additional columns for each k-mer 
#'         representing its frequency count in each sequence.
#'
#' @examples
#' # Assuming 'df' is a dataframe and 'fasta' contains sequences:
#' kmer_freq_df <- calculate_kmer_frequency(df, fasta, kmer_window)
#'
#' @importFrom dplyr mutate select
#' @importFrom purrr map
#' @importFrom tidyr unnest pivot_wider
#' @importFrom seqinr count
calculate_kmer_frequency <- function(df, fasta, kmer_window) {
  df <- df %>% 
    mutate(
      contig_seq = map(.x = contig, ~fasta[[paste0(.x)]]),
      kmer_count = map(.x = contig_seq, ~(as.data.frame(seqinr::count(.x, wordsize = kmer_window)) %>% pivot_wider(names_from = Var1, values_from = Freq)))
    ) %>% 
    unnest(kmer_count) %>% 
    select(!contig_seq)
}



# TODO: Implement logic if no contigs could be assigned to multiple bins

#' Calculate distance from Contamination to bin from K-mer Frequency
#'
#' This function assigns contamination categories to sequences based on k-mer frequency analysis. 
#' It involves several steps including calculation of k-mer frequencies for bins and contigs, 
#' determination of contamination categories, and analysis of k-mer distances.
#'
#' @param belonging_score A dataframe containing the belonging scores, used to identify potential 
#'        contamination categories.
#' @param fasta A list or other suitable structure containing sequences from a FASTA file, 
#'        indexed by sequence identifiers.
#' @param contig_bins A dataframe mapping contigs to bins.
#' @param kmer_window An integer specifying the k-mer window size for frequency calculation.
#'
#' @return A dataframe containing the calculated k-mer distances for each contig, which can be used 
#'         to infer contamination.
#'
#' @examples
#' # Assuming you have the necessary 'belonging_score', 'fasta', 'contig_bins', and 'kmer_window':
#' contamination_result <- assign_contamination_from_kmer(belonging_score, fasta, contig_bins, kmer_window)
#'
#' @importFrom dplyr filter mutate select group_by summarize ungroup left_join arrange distinct 
#' @importFrom purrr map map2
#' @importFrom tidyr pivot_longer pivot_wider unnest nest
#' @importFrom tibble rownames_to_column
#' @importFrom stringr str_detect
#' @export
calculate_contamination_dist_from_kmer <- function(
  belonging_score,
  fasta,
  contig_bins,
  kmer_window
) {
  print("Calculating Bin kmer frequency")
  kmer_bin_df <- contig_bins %>% 
    calculate_kmer_frequency(fasta, kmer_window)
  
  # Remove unbinned contigs as kmer frequency does not make sense.
  contamination_multi <- belonging_score %>% 
    filter(!str_detect(bin_compare, "unbinned")) %>% 
    determine_contamination() %>% 
    filter(
      group == "multi_contamination"
    )
  
  if (length(contamination_multi$bin_compare) == 0) {
    print("No contigs could be assigned to multiple bins")
    return(tibble(contig = c(), matching_bins = c()))
  }
  print("Calculating contig kmer frequency")
  kmer_contig_df <- tibble(contig = contamination_multi %>% pull(contig) %>% unique()) %>% 
    # calculate_kmer_frequency(fasta)
    left_join(kmer_bin_df) %>%
    select(!bin)
  
  
  
  cat("Calculating kmer distance to bins\n")
  contig_bins_out_tmp <- contamination_multi %>% 
    ungroup() %>% 
    mutate(
      # Extract kmer bins, where contig could belong (remove contig if found in bin)
      bin_kmer = map2(
        .x = contig, .y = bin, 
        ~(
          kmer_bin_df %>% 
            filter(contig != .x, bin == .y) %>% 
            mutate(contig = .y) %>% 
            select(-bin)
        )
      ),
      # Extract kmer for contig
      bin_contig_kmer = map2(
        .x = bin_kmer, .y = contig, 
        ~(
          kmer_contig_df %>% 
            filter(contig == .y) %>%
            bind_rows(.x)
        )
      ),
      # Summarise kmer counts for the bins and convert to frequency
      bin_contig_kmer = map(
        .x = bin_contig_kmer, 
        ~(
          .x %>% 
            pivot_longer(-contig, names_to = "kmer", values_to = "count") %>% 
            group_by(contig, kmer) %>% 
            summarise(count = sum(count), .groups = "drop") %>% 
            # Calculate total count per contig
            mutate(total_count = sum(count)) %>%
            # Calculate frequency
            mutate(frequency = count / total_count) %>%
            # Remove the total_count column
            select(-c(total_count, count)) %>%
            pivot_wider(names_from = "kmer", values_from = "frequency")
        )
      )
    ) %>% 
    ungroup() %>% 
    select(contig, bin_contig_kmer) %>% 
    arrange(contig)
  
  
  contig_bins_out <- contig_bins_out_tmp %>% 
    rename(unbinned = contig) %>% 
    unnest(bin_contig_kmer) %>%
    group_by(unbinned) %>% 
    # Remove multiple entries for contig
    distinct(contig, .keep_all = TRUE) %>% 
    ungroup() %>% 
    # Print matching bins and distance
    nest(.by = unbinned) %>% 
    mutate(
      df_matrix = map(
        .x = data, 
        ~(.x %>% column_to_rownames("contig"))
      ),
      dist_matrix = map(
        .x = df_matrix, 
        ~(.x %>% dist())
      ),
      kmer_dist = map(.x = dist_matrix, .y = unbinned, ~extract_kmer_dist(.x))
    ) %>% 
    select(kmer_dist) %>% 
    unnest(kmer_dist)
  
  return(contig_bins_out)
}

