##### Library dependencies
library(flowCore)

#####Functions to run flowCore metadata search
flpaths_FCS_list <- function(fld_path) {list.files(path=fld_path, pattern = ".fcs", recursive = T, full.names = T)}

make_sampleName <- function(string) {
  basename(string) |> gsub(pattern = ".fcs", replacement = "") |> gsub(pattern = "[abcdefghABCDEFGH]{1}[0123456789]{1,2}", replacement = "") |> gsub(pattern = "_", replacement = "")
}

std_readFCS <- function(fls, transf = "linearize") {
  flowCore::read.FCS(fls,
                     alter.names = TRUE,
                     transformation = transf)
}

loadFF <- function(fls) {
  frames<-std_readFCS(fls)
  return(frames)
}

get_metaData <- function(fF) {
  pe<-function(x) paste0("fF@description$`$P",x,"E`")
  evaluate <- function(x) eval(parse(text = x))
  
  instrument = list(
    OS = fF@description$`$SYS`,
    TOT_PARAMETERS = fF@description$`$PAR`,
    DEV_SOFTWARE = fF@description$CREATOR,
    DEV = fF@description$`$CYT`,
    DEV_SERIAL = fF@description$`CYTNUM`,
    LASERS = c(fF@description$`LASER1NAME`,fF@description$`LASER2NAME`,fF@description$`LASER3NAME`,fF@description$`LASER4NAME`)
  )
  
  PNE <- sapply(1:length(instrument$PARAMETER_NAMES), pe)
  
  experiment = list(
    exp_name = fF@description[["EXPERIMENT NAME"]]
    ,file_name = fF@description[["FILENAME"]]
    ,esoteric_id = fF@description[["FILENAME"]] |> make_sampleName()
    ,comp.matrix = fF@description$SPILL
    ,PNEvalues = sapply(PNE, evaluate)
    ,colors = fF@parameters@data[["name"]][grep("FSC.*|SSC.*|Time",fF@parameters@data[["name"]], invert = T)]
    ,markers = markernames(fF)
  )
  
  ls <- list(instrument = instrument, experiment = experiment)
  return(ls)
}

#####Functions to run flowCore workflow

get_hypsin_tfx_2 <- function(fF, fluo) {
  lgcl <- estimateLogicle(fF, fluo)
  tfx_data <- transform(fF, lgcl)
  return(tfx_data)
}

make_lymph_pop_data <- function(file) {
  #### step 1 Load file into a flowFrame object
  sampleName <- make_sampleName(file)
  print(paste("Loading sample:", sampleName, "into flowFrame"))
  exp_fF <- loadFF(file)
  
  #### step 2 Get experimental meta data parameters
  meta_data <- get_metaData(exp_fF)
  
  #### step 3 Compensate and transform cytometry parameter data
  stat_raw <- summary(exp_fF@exprs[,c(11,10,14)])
  qc1_fF <- flowCore::compensate(exp_fF, exp_fF@description[["SPILL"]])
  stat_comp <-summary(qc1_fF@exprs[,c(11,10,14)])
  print(paste("Compensation check:", ifelse(!identical(stat_raw, stat_comp),"PASS, altered signal data", "FAIL; see compensation matrix - signal data is exactly the same")))
  comp_check <- ifelse(!identical(stat_raw, stat_comp), "PASS", "FAIL")
  
  tryCatch({
    skip <- F
    error_df <- matrix(data = c(sampleName, rep(NA, 1), basename(file), rep(NA, 8)),nrow = 1, ncol = 11) |> as.data.frame()
    qcFinal_fF <- get_hypsin_tfx_2(qc1_fF, unname(exp_fF@parameters@data[["name"]]))}
    ,warning = function(w) {
      print("Warning: Transform did not work and returning empty matrix")
      skip <<- T}
    ,error = function(e) {
      print("Error: Transform did not work and returning empty matrix")
      skip <<- T
      })
  
  if(skip == T) {return(error_df)} else {
    stat_trsfm <-summary(qcFinal_fF)
    print(paste("Transformation check:", ifelse(!identical(stat_raw, stat_trsfm),"PASS, altered signal data", "FAIL; see get_hypsin_tfx wrapper - signal data is exactly the same")))
    trsfm_check <- ifelse(!identical(stat_raw, stat_trsfm), "PASS", "FAIL")
    
    #### step 4 filter data for positive populations using kmeansFilter function
    kmf <- kmeansFilter("FSC.A" = c("debris","lymph","macrophages"), filterId = "kfilter" )
    kmf_results <- filter(qcFinal_fF, kmf)
    debris_filtered <- split(qcFinal_fF, kmf_results)
    
    
    kmf_cd3 <- kmeansFilter("Am.Cyan.A" = c("CD3neg","CD3pos"), filterId = "CD3" )
    cd3_filter <- filter(debris_filtered$lymph, kmf_cd3)
    CD3_pops <- split(debris_filtered$lymph, cd3_filter)
    
    kmf_cd8 <- kmeansFilter("APC.Cy7.A" = c("CD8neg","CD8pos"), filterId = "CD8" )
    flt_CD8 <- filter(CD3_pops$CD3pos, kmf_cd8)
    CD3posCD8 <- split(x = CD3_pops$CD3pos, flt_CD8)
    
    kmf_cd4 <- kmeansFilter("Pacific.Blue.A" = c("CD4neg","CD4pos"), filterId = "CD4" )
    flt_CD8negCD4pos <- filter(CD3posCD8$CD8neg, kmf_cd4)
    CD3posCD8negCD4 <- split(x = CD3posCD8$CD8neg, flt_CD8negCD4pos)
    flt_CD8posCD4 <- filter(CD3posCD8$CD8pos, kmf_cd4)
    CD3posCD8posCD4 <- split(x = CD3posCD8$CD8pos, flt_CD8posCD4)
    
    ls_qG_params <- list("APC.Cy7.A" = 2.5,  "Pacific.Blue.A" =  1.85)
    qG_CD8_CD4 <- quadGate(.gate = ls_qG_params)
    flt_CD8_CD4 <- filter(CD3_pops$CD3pos, qG_CD8_CD4)
    CD8_CD4_pop <- split(x = CD3_pops$CD3pos, flt_CD8_CD4)
    
    #### step 5 generate table
    
    # Tbl 1 Measured classified counts
    
    pop_names <- names(CD8_CD4_pop)
    mtx_meta_colnames <- c("Sample", "Experiment", "Filename", "Device", "Device Serial", "Target(s)", "Colors", pop_names)
    df <- matrix(data = replicate(11, "NA"), nrow = 1, ncol = 11) |> `colnames<-`(mtx_meta_colnames) |> as.data.frame()
    for (names in names(CD8_CD4_pop)) {
      df[names] <- nrow(CD8_CD4_pop[[names]]) |> as.numeric()
    }
    df["Filename"] <- meta_data$experiment$file_name |> basename()
    df["Experiment"] <- meta_data$experiment$exp_name
    df["Device"] <- meta_data$instrument$DEV
    df["Target(s)"] <- meta_data$experiment$markers |> unname() |> paste(collapse = " ")
    df["Colors"] <- meta_data$experiment$markers |> names() |> paste(collapse = " ")
    df["Device Serial"] <- meta_data$instrument$DEV_SERIAL
    df["Sample"] <- basename(file) |> gsub(pattern = ".fcs", replacement = "") |> gsub(pattern = "[abcdefghABCDEFGH]{1}[0123456789]{1,2}", replacement = "") |> gsub(pattern = "_", replacement = "")
    
    print(df)
    return(df)}
}

#####Functions to perform data organization and basic statistical analysis

make_lymphDF <- function(lst) {
  df <- matrix(nrow = 0, ncol = 11) |> as.data.frame()
  for(j in 1:length(lst)) {
    row_being_added <- lst[[j]]
    colNames <- names(row_being_added)
    names(df) <- colNames
    df <- rbind(df, row_being_added)
  }
  colnames(df) <- c("Sample", "Experiment", "Filename", "Device", "Device.Serial",
                     "Target(s)", "Colors", "CD8+CD4+", "CD8-CD4+", "CD8+CD4-", 
                     "CD8-CD4-")
  return(df)
}

sort_by_exp <- function(df) {
  exp <- unique(df$Experiment)
}

grab_df <- function(exp, df) {
  rows <- is.na(df$`CD8+CD4+`)
  filtered_df <- df[!rows,]
  filtered_df <- filtered_df[filtered_df$Experiment == exp,]
  return(filtered_df)
}

make_5_num_summary_per_col <- function(df) {
  colname <- c("Min", "Lower Hinge", "Median", "Upper Hinge", "Max")
  rowname <- names(df) |> grep(pattern = "CD", value = T, invert = F)
  ignore <- names(df) |> grep(pattern = "CD", value = T, invert = T)
  dim_names <- list(rowname, colname)
  fivenum_df <- matrix(nrow = length(rowname), ncol = length(colname), dimnames = dim_names) |> as.data.frame()
  counter = 1
  for(col in names(df)) {
    if(any(col == ignore)) {next}
    values <- df[[col]] |> as.numeric()
    fivenum_df[counter,] <- fivenum(values, na.rm = T)
    counter = counter + 1
  }
  return(as.data.frame(fivenum_df))
}

calc_perc <- function(x,y) (x/y *100) |> round(digits = 1)

make_percentages_across_rows <- function(df) {
  colnames <- c(names(df), 'Total')
  rownames <- df[['Sample']]
  dim_names <- list(rownames, colnames)
  perc_df <- matrix(nrow = nrow(df), ncol = length(colnames), dimnames = dim_names)
  lymphs <- names(df) |> grep(pattern = "CD", value = T, invert = F)
  for(row in 1:nrow(df)) {
    total_cells <- sum(df[row,lymphs] |> as.numeric(), na.rm = T) 
    for(col in 1:length(names(df))) {
      col_name <- names(df)[col]
      if(col_name %in% lymphs) {
        vals <- df[row,col] |> as.numeric() |> calc_perc(y = total_cells)
        perc_df[row,col] <- vals
        next
      } else {
        perc_df[row,col] <- df[row,col]
        next}
    }
    total <- sum(perc_df[row, 8:11]|>as.numeric())
    if(is.na(total)){
      perc_df[row,12] <- NA
    } else { perc_df[row,12] <- total |> round(digits = 0) }
  }
  perc_df <- perc_df |> as.data.frame()
  return(perc_df)
}

convert_char_num <- function(df, ignore) {
  for(col in names(df)) {
    if(col != ignore) df[[col]] <- as.numeric(df[[col]])
  }
  return(df)
}

#Main
dir <- "~/Documents/Proj/DS/FCM/Flow_studies/"
exp_lymph_raw_df <-flpaths_FCS_list(dir) |> lapply(FUN = make_lymph_pop_data) |> make_lymphDF()
exp_lymph_perc_df <- make_percentages_across_rows(exp_lymph_raw_df)
exp_names <- unique(exp_lymph_raw_df$Experiment)
vec <- unique(exp_lymph_raw_df$Experiment) |> is.na()
exp_names<- exp_names[!vec]
grouped_exp_raw_lst <- lapply(X = exp_names, FUN = grab_df, df = exp_lymph_raw_df) |> `names<-`(value = exp_names)
grouped_exp_perc_lst <- lapply(X = exp_names, FUN = grab_df, df = exp_lymph_perc_df) |> `names<-`(value = exp_names)
raw_five_num_lst <- lapply(X = grouped_exp_raw_lst, FUN = make_5_num_summary_per_col) |> `names<-`(value = exp_names)
perc_five_num_lst <- lapply(X = grouped_exp_perc_lst, FUN = make_5_num_summary_per_col) |> `names<-`(value = exp_names)
