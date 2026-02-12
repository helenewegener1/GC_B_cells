# Function to combine follicles within each sample group
combine_follicles <- function(bcr_list) {
  
  # Parse sample names to identify base sample (without follicle number)
  sample_metadata <- data.frame(
    original_name = names(bcr_list),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      # Remove follicle numbers, Negative, and Doublet suffixes
      base_sample = str_remove(original_name, "_Fol-\\d+$"),
      base_sample = str_remove(base_sample, "_Negative$"),
      base_sample = str_remove(base_sample, "_Doublet$")
    )
  
  # Get unique base samples
  unique_samples <- unique(sample_metadata$base_sample)
  
  # Combine data for each base sample
  combined_list <- lapply(unique_samples, function(base) {
    # Find all original samples that belong to this base
    matching_samples <- sample_metadata %>%
      filter(base_sample == base) %>%
      pull(original_name)
    
    # Combine all matching samples
    combined_df <- bind_rows(bcr_list[matching_samples], .id = "original_sample")
    
    return(combined_df)
  })
  
  # Name the list with base sample names
  names(combined_list) <- unique_samples
  
  return(combined_list)
}
