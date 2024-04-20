record_stats <- function(statistics, title = "SimulationReport", 
                         outdir = directory, fileform = "csv", delim = ";"){
  filename <- paste(outdir, "/", title, ".", fileform, sep = "")
  
  # Check if the file exists
  if (!file.exists(filename)) {
    # If the file doesn't exist, create it with column names from statistics
    if (is.null(dim(statistics)[1])) {
      report <- data.frame(matrix(nrow = 0, ncol = length(statistics)))
      names(report) <- names(statistics)
    } else {
      report <- data.frame(matrix(nrow = 0, ncol = dim(statistics)[2]))
      names(report) <- colnames(statistics)
    }
    
    write.table(report, 
                filename,
                append = FALSE,
                sep = delim,
                row.names = FALSE,
                col.names = !file.exists(filename))
    
    # cat("Report file was created. \n")
  } 
  
  # Read the existing report file
  report <- read.csv(filename, header = FALSE, sep = delim, encoding="UTF-8")
  
  if (is.null(dim(statistics)[1])) {
    colnames(report) <- names(statistics)
  } else {
    colnames(report) <- colnames(statistics)
  }
  
  # Append statistics to the report dataframe
  report <- rbind(report, statistics)
  
  # Write the updated report to the file
  write.table(report,
              filename,
              append = FALSE,
              sep = delim,
              row.names = FALSE,
              col.names = !file.exists(filename))
  
  # cat("Report file was updated. \n")
}