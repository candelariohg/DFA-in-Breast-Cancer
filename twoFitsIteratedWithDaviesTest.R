# ---------------------------------------------------------------
# Script: run_dfa_segmented.R (Versión corregida para CIs)
# ---------------------------------------------------------------

suppressPackageStartupMessages(library(segmented, lib.loc="~/R/libs"))

# Define chromosomes and subtypes
chromosomes <- c(1:22, "X")
subtypes <- c("Basal", "Healthy", "Her2", "LumB", "LumA", "Normal")

# Data frame para almacenar todos los resultados y armar la tabla final
master_results <- data.frame(
  File = character(),
  Breakpoint = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  Y_Break = numeric(),
  Slope1 = numeric(),
  Slope2 = numeric(),
  Davies_PValue = numeric(),
  stringsAsFactors = FALSE
)

# Loop through all chromosome/subtype combinations
for (chr in chromosomes) {
  for (subtype in subtypes) {
    infile <- sprintf("chr%s%sMISortedByDistance.dat", chr, subtype)
    
    if (!file.exists(infile)) {
      message(sprintf("⚠️ File not found: %s", infile))
      next
    }
    
    message(sprintf("🔹 Processing %s ...", infile))
    
    # Run external DFA program and capture its output
    dfa_output <- tryCatch({
      system(paste("./dfa < ", infile), intern = TRUE)
    }, error = function(e) NULL)
    
    if (is.null(dfa_output) || length(dfa_output) == 0) {
      message(sprintf("❌ DFA failed or produced no output for %s", infile))
      next
    }
    
    # Convert DFA output (two columns) to numeric data frame
    con <- textConnection(dfa_output)
    data <- tryCatch(read.table(con, header = FALSE), error = function(e) NULL)
    close(con)
    
    if (is.null(data) || ncol(data) < 2) {
      message(sprintf("⚠️ Invalid DFA output for %s", infile))
      next
    }
    
    x <- data[[1]]
    y <- data[[2]]
    
    # Fit simple linear model
    fit <- lm(y ~ x)
    
    # Davies test (Score test)
    davies_res <- tryCatch(davies.test(fit), error = function(e) NULL)
    p_val_davies <- if (!is.null(davies_res)) davies_res$p.value else NA
    
    # Fit segmented model
    segfit <- tryCatch(segmented(fit), error = function(e) NULL)
    
    if (is.null(segfit)) {
      message(sprintf("⚠️ Segmented fit failed for %s", infile))
      next
    }
    
    # Extract results
    psi_x <- segfit$psi[2]     # Breakpoint estimate
    coef_seg <- coef(segfit)
    slopes <- slope(segfit)$x
    y_break <- coef_seg[1] + slopes[1] * psi_x
    
    # Confidence intervals for the breakpoint (CORREGIDO: usando confint(segfit) directo)
    ci_psi <- tryCatch(confint(segfit, parm = "x"), error = function(e) NULL)
    ci_lower <- if (!is.null(ci_psi) && nrow(ci_psi) > 0) ci_psi[1, 1] else NA
    ci_upper <- if (!is.null(ci_psi) && nrow(ci_psi) > 0) ci_psi[1, 2] else NA
    
    # Guardar en el data frame global
    master_results <- rbind(master_results, data.frame(
      File = infile,
      Breakpoint = psi_x,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      Y_Break = y_break,
      Slope1 = slopes[1],
      Slope2 = slopes[2],
      Davies_PValue = p_val_davies
    ))
    
    # Save individual result to text file
    result_line <- sprintf("%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f",
                           infile, psi_x, ci_lower, ci_upper, y_break, slopes[1], slopes[2], p_val_davies)
    result_file <- sub(".dat$", "_segmented_results.txt", infile)
    writeLines(result_line, result_file)
    
    # Save plot
    plot_file <- sub(".dat$", "_segmented_plot.png", infile)
    png(plot_file, width = 800, height = 600)
    plot(x, y, main = paste("Segmented Regression for", infile),
         xlab = "x", ylab = "y", pch = 16, col = "darkgray")
    plot(segfit, add = TRUE, col = "red", lwd = 2)
    points(psi_x, y_break, col = "blue", pch = 19, cex = 1.2)
    legend("topleft",
           legend = c("Data (DFA output)", "Segmented fit", "Breakpoint"),
           col = c("darkgray", "red", "blue"),
           pch = c(16, NA, 19),
           lty = c(NA, 1, NA),
           lwd = c(NA, 2, NA),
           bty = "n")
    dev.off()
  }
}

# Guardar la tabla maestra consolidada
write.table(master_results, file = "resumen_global_segmented.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

message("\n🏁 All files processed successfully! Master table saved as 'resumen_global_segmented.tsv'.\n")

# Imprimir la tabla limpia y bonita en consola al terminar
print("=== TABLA RESUMEN DE RESULTADOS SEGMENTADOS ===")
print(master_results)
