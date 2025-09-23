# ---------- #
# Function to read in PopCluster results
# ---------- #

readPopCluster <- function(path, file, genind) {
  
  # Read the lines of the text file
  lines <- readLines(paste0(path, file))
  
  # Find the index of the marker that indicates the start of the table
  start_index <- grep("Inferred ancestry of individuals", lines)
  
  # Find the index of the first empty line after the start of the table
  end_index <- start_index + which(lines[start_index:length(lines)] == "")[1] - 1
  
  # Skip lines before the start of the table and read the table into a data frame
  table_data <- read.table(text = lines[(start_index + 2):end_index], header = FALSE, fill = TRUE)
  
  # Subset columns
  admix_data <- subset(table_data, select = c(3,7:ncol(table_data)))
  colnames(admix_data) <- c("Individual", paste(rep("Cluster"), 1:(ncol(admix_data)-1)))
  
  # Add a column for site labels
  admix_data$Site <- as.character(genind$pop)
  
  admix_data <- dplyr::select(admix_data, Site, Individual, everything())
  
  # Return table
  return(admix_data)
}