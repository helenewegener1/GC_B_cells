# Zero point: Point between noise and signal defined on log10 scale.

ADT_zero_point_log10 <- list(
  
  "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH" = c(
      "Fol-1" = 2.1,
      "Fol-2" = 2.4,
      "Fol-3" = 1.6,
      "Fol-4" = 2,
      "Fol-5" = 2, 
      "Fol-6" = 2.2,
      "Fol-7" = 1.6, 
      "Fol-8" = 1.9,
      "Fol-9" = 1.5, 
      "Fol-10" = 1.5,
      "Fol-11" = 1.9,
      "Fol-12" = 2.1, # maybe check this
      "Fol-13" = 2.45, # maybe check this
      "Fol-14" = 1.6,
      "Fol-15" = 0.8,
      "Fol-16" = 1.2,
      "Fol-17" = 1.7,  
      "Fol-18" = 1.9
  ),
  
  "HH119-SI-PP-CD19-Pool1" = c(
    "Fol-1" = 2.1,
    "Fol-2" = 2.8,
    "Fol-3" = 1.8,
    "Fol-4" = 2.2,
    "Fol-5" = 1.5,
    "Fol-6" = 2,
    "Fol-7" = 2.1,
    "Fol-8" = 2,
    "Fol-9" = 2,
    "Fol-10" = 2,
    "Fol-11" = 2.1,
    "Fol-12" = 2.2,
    "Fol-13" = 1.7,
    "Fol-14" = 2.4,
    "Fol-15" = 2.5,
    "Fol-16" = 2.2,
    "Fol-17" = 2
  ),

  "HH119-SI-PP-CD19-Pool2" = c(
    "Fol-18" = 1.9,
    "Fol-19" = 2.3,
    "Fol-20" = 1.5,
    "Fol-21" = 1.5,
    "Fol-22" = 2,
    "Fol-23" = 2,
    "Fol-24" = 2.1,
    "Fol-25" = 2.2,
    "Fol-26" = 1.7,
    "Fol-27" = 2,
    "Fol-28" = 1.8,
    "Fol-29" = 2,
    "Fol-30" = 1.4,
    "Fol-31" = 2,
    "Fol-32" = 2,
    "Fol-33" = 2,
    "Fol-34" = 2
  ),
  
  "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1" = c(
    "Fol-1" = 2.2,
    "Fol-2" = 3,
    "Fol-3" = 2,
    "Fol-4" = 2.4,
    "Fol-5" = 1.7,
    "Fol-6" = 2,
    "Fol-7" = 2.3,
    "Fol-8" = 2,
    "Fol-9" = 2.3,
    "Fol-10" = 2,
    "Fol-11" = 2,
    "Fol-12" = 1.8,
    "Fol-13" = 1.9,
    "Fol-14" = 2.4,
    "Fol-15" = 2.6,
    "Fol-16" = 2.3,
    "Fol-17" = 2
  ),
  
  "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2" = c(
    "Fol-18" = 2.8,
    "Fol-19" = 2.2,
    "Fol-20" = 1.5,
    "Fol-21" = 1.5,
    "Fol-22" = 2.1,
    "Fol-23" = 1.9,
    "Fol-24" = 2.3,
    "Fol-25" = 2.2,
    "Fol-26" = 1.7,
    "Fol-27" = 2,
    "Fol-28" = 1.8,
    "Fol-29" = 2,
    "Fol-30" = 1.7,
    "Fol-31" = 1.9,
    "Fol-32" = 1.9,
    "Fol-33" = 1.9,
    "Fol-34" = 2
  )
  
  # "" = c(
  #   "Fol-1" = ,
  #   "Fol-2" = ,
  #   "Fol-3" = ,
  #   "Fol-4" = ,
  #   "Fol-5" = ,
  #   "Fol-6" = ,
  #   "Fol-7" = ,
  #   "Fol-8" = , 
  #   "Fol-9" = ,
  #   "Fol-10" = ,
  #   "Fol-11" = ,
  #   "Fol-12" = , 
  #   "Fol-13" = , 
  #   "Fol-14" = ,
  #   "Fol-15" = ,
  #   "Fol-16" = ,
  #   "Fol-17" = ,
  #   "Fol-18" = 
  # )
  
  
  
)

# Transform zero points to raw counts 
# zero_point <- 10^zero_point_log10 

ADT_zero_point <- lapply(ADT_zero_point_log10, function(x) 10^x)

