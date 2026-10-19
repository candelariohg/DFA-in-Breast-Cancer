################################################################################
#                                                                              #
# This program calculates the mutual information of the entire chromosome 8    #
# with the entire chromosome 8 for expression values of cancer variants.       #
#                                                                              #
# The result is saved in an external array compatible with R.                  #
#                                                                              #
#                                                                              #
# What the program receives: A gene expression matrix.                         #
# What this program gives: An array of mutual information values.              #
#                                                                              #
# Made by:                                                                     #
# Last modified: September 27, 2025                                            #
################################################################################

library(infotheo)

################################################################################
# 1) Arrays are read from an external file.                                    #
################################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Please provide a filename, e.g.: Rscript myscript.R datos.tsv")
}

filename <- args[1]

# Read expression matrix
x <- read.table(filename, header = TRUE, sep = "\t", row.names = 1)

R <- nrow(x)   # Number of rows
C <- ncol(x)   # Number of columns

cat("Rows:", R, "\n")
cat("Columns:", C, "\n")

################################################################################
# 2) The Mutual Information (MI) matrix is created.                            #
################################################################################

IM <- matrix(0, nrow = R, ncol = R)   # square matrix initialized with zeros

################################################################################
# 3) The MI is calculated.                                                     #
################################################################################

for (i in 1:R) {
  fila1 <- unname(unlist(x[i,]))    # i-th row
  fila1d <- discretize(fila1)       # discretize values for MI calculation
  
  j <- i + 1
  while (j <= R) {
    fila2 <- unname(unlist(x[j,]))  # j-th row
    fila2d <- discretize(fila2)
    
    infMut <- mutinformation(fila1d, fila2d)
    IM[i, j] <- infMut
    IM[j, i] <- infMut   # matrix is symmetric
    
    j <- j + 1
  }
}

################################################################################
# 4) Assign row/column names and save to file.                                 #
################################################################################

rownames(IM) <- rownames(x)
colnames(IM) <- rownames(x)

# Create output filename: input base name + "MI.tsv"
base_name <- tools::file_path_sans_ext(basename(filename))
output_file <- paste0(base_name, "MI.tsv")

write.table(IM, file = output_file, sep = "\t", dec = ".", 
            row.names = TRUE, col.names = NA, quote = FALSE)

cat("Mutual Information matrix saved to", output_file, "\n")
################################################################################

