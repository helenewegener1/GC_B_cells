# Zero point: Point between noise and signal defined on log10 scale.

ADT_zero_point_log10 <- list(
  
  "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH" = c(
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
  
  "HH119-SI-PP-CD19-Pool1" = c(
    "Fol-1" = 2.7,
    "Fol-2" = 3.5,
    "Fol-3" = 2.3,
    "Fol-4" = 3.5,
    "Fol-5" = 2,
    "Fol-6" = 2.5,
    "Fol-7" = 2.5,
    "Fol-8" = 2.3,
    "Fol-9" = 2.5,
    "Fol-10" = 2.5,
    "Fol-11" = 2.4,
    "Fol-12" = 2.7,
    "Fol-13" = 2.5,
    "Fol-14" = 2.7,
    "Fol-15" = 2.9,
    "Fol-16" = 2.9,
    "Fol-17" = 2.5
  ),

  "HH119-SI-PP-CD19-Pool2" = c(
    "Fol-18" = 2.5,
    "Fol-19" = 2.3,
    "Fol-20" = 2.1,
    "Fol-21" = 2,
    "Fol-22" = 2.2,
    "Fol-23" = 2.2,
    "Fol-24" = 2.5,
    "Fol-25" = 2.4,
    "Fol-26" = 2.2,
    "Fol-27" = 2.3,
    "Fol-28" = 2.4,
    "Fol-29" = 2.5,
    "Fol-30" = 2,
    "Fol-31" = 2.4,
    "Fol-32" = 2.5,
    "Fol-33" = 2.5,
    "Fol-34" = 2.35
  ),
  
  "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1" = c(
    "Fol-1" = 2.7,
    "Fol-2" = 3.3,
    "Fol-3" = 2.5,
    "Fol-4" = 2.8,
    "Fol-5" = 2.5,
    "Fol-6" = 2.5,
    "Fol-7" = 2.8,
    "Fol-8" = 2.65,
    "Fol-9" = 2.6,
    "Fol-10" = 2.6,
    "Fol-11" = 2.6,
    "Fol-12" = 2.6,
    "Fol-13" = 2.4,
    "Fol-14" = 2.7,
    "Fol-15" = 2.85,
    "Fol-16" = 2.8,
    "Fol-17" = 2.8
  ),
  
  "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2" = c(
    "Fol-18" = 2.5,
    "Fol-19" = 2.7,
    "Fol-20" = 2.5,
    "Fol-21" = 2.1,
    "Fol-22" = 2.6,
    "Fol-23" = 2.5,
    "Fol-24" = 2.55,
    "Fol-25" = 2.8,
    "Fol-26" = 2,
    "Fol-27" = 2.5,
    "Fol-28" = 2.5,
    "Fol-29" = 2.6,
    "Fol-30" = 2,
    "Fol-31" = 2.5,
    "Fol-32" = 2.5,
    "Fol-33" = 2.6,
    "Fol-34" = 2.5
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

