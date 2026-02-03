# Zero point: Point between noise and signal defined on log10 scale.

ADT_zero_point_log10 <- list(
  
  "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH" <- c(
      "Fol-1" = 2.6,
      "Fol-2" = 2.8,
      "Fol-3" = 2.1,
      "Fol-4" = 2.4,
      "Fol-5" = 2.7, 
      "Fol-6" = 2.5,
      "Fol-7" = 2.5, 
      "Fol-8" = 2, # check 
      "Fol-9" = 2.3, 
      "Fol-10" = 2.2,
      "Fol-11" = 2.8,
      "Fol-12" = 2.5, # maybe check this
      "Fol-13" = 2.5, # maybe check this
      "Fol-14" = 2.4,
      "Fol-15" = 1.7,
      "Fol-16" = 2.5,
      "Fol-17" = 2.6,  
      "Fol-18" = 2.4
  ),
  
  # "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH" <- c(
  #   "Fol-1" = 2.6,
  #   "Fol-2" = 2.8,
  #   "Fol-3" = 2.1,
  #   "Fol-4" = 2.4,
  #   "Fol-5" = 2.7, 
  #   "Fol-6" = 2.5,
  #   "Fol-7" = 2.5, 
  #   "Fol-8" = 2, # check 
  #   "Fol-9" = 2.3, 
  #   "Fol-10" = 2.2,
  #   "Fol-11" = 2.8,
  #   "Fol-12" = 2.5, # maybe check this
  #   "Fol-13" = 2.5, # maybe check this
  #   "Fol-14" = 2.4,
  #   "Fol-15" = 1.7,
  #   "Fol-16" = 2.5,
  #   "Fol-17" = 2.6,  
  #   "Fol-18" = 2.4
  # )
  
  
  
  
  
)

# Transform zero points to raw counts 
zero_point <- 10^zero_point_log10 

ADT_zero_point <- lapply(ADT_zero_point_log10, function(x) 10^x)
