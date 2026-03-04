# ====================================================
# False Discovery Rate Case Study: 
# Comparing Treatments Across Multiple Populations
# ====================================================
install.packages("DiscreteFDR")
# Load required packages
library(DiscreteFDR)
library(DiscreteTests)
library(ggplot2)
library(dplyr)

# ----------------------------------------------------
# 1. Create the dataset
# ----------------------------------------------------
# Data represents 9 populations with responder/non-responder counts
# for Treatment 1 (X1, Y1) and Treatment 2 (X2, Y2) [citation:8]

df <- data.frame(
  Population = 1:9,
  X1 = c(4, 2, 2, 14, 6, 9, 4, 0, 1),   # Responders - Treatment 1
  Y1 = c(144, 146, 146, 134, 142, 139, 144, 148, 147), # Non-responders - Treatment 1
  X2 = c(0, 0, 1, 3, 2, 1, 2, 2, 2),    # Responders - Treatment 2
  Y2 = c(132, 132, 131, 129, 130, 131, 130, 130, 130)  # Non-responders - Treatment 2
)

# Calculate total sample sizes
df$N1 <- df$X1 + df$Y1  # Total Treatment 1
df$N2 <- df$X2 + df$Y2  # Total Treatment 2
df$Total <- df$N1 + df$N2

print("Dataset: Responder counts by population and treatment")
print(df[, c("Population", "X1", "Y1", "X2", "Y2", "N1", "N2")])

# ----------------------------------------------------
# 2. Perform Fisher's exact tests for each population
# ----------------------------------------------------
# Test H0: Response rates are equal between treatments

p_values <- numeric(9)
test_results <- list()

for(i in 1:9) {
  # Create 2x2 contingency table
  table_i <- matrix(c(df$X1[i], df$Y1[i], df$X2[i], df$Y2[i]), 
                    nrow = 2, byrow = TRUE)
  colnames(table_i) <- c("Responders", "Non-responders")
  rownames(table_i) <- c("Treatment 1", "Treatment 2")
  
  # Perform Fisher's exact test
  test <- fisher.test(table_i, alternative = "two.sided")
  p_values[i] <- test$p.value
  test_results[[i]] <- test
  
  # Print for first few populations
  if(i <= 3) {
    cat("\nPopulation", i, "Contingency Table:\n")
    print(table_i)
    cat("P-value:", round(test$p.value, 4), "\n")
  }
}

# ----------------------------------------------------
# 3. Apply multiple testing corrections
# ----------------------------------------------------
alpha <- 0.05

# Bonferroni correction (FWER control)
bonferroni_p <- p.adjust(p_values, method = "bonferroni")
bonferroni_reject <- bonferroni_p < alpha

# Benjamini-Hochberg (FDR control)
bh_p <- p.adjust(p_values, method = "BH")
bh_reject <- bh_p < alpha

# Also try Holm correction (step-down FWER)
holm_p <- p.adjust(p_values, method = "holm")
holm_reject <- holm_p < alpha

# ----------------------------------------------------
# 4. Create results table
# ----------------------------------------------------
results <- data.frame(
  Population = 1:9,
  Raw_P = round(p_values, 5),
  BH_Adjusted = round(bh_p, 5),
  Bonferroni = round(bonferroni_p, 5),
  Holm = round(holm_p, 5),
  BH_Reject = bh_reject,
  Bonferroni_Reject = bonferroni_reject,
  Holm_Reject = holm_reject
)

# Sort by raw p-value for easier interpretation
results <- results[order(results$Raw_P), ]
row.names(results) <- NULL

print("\n=== Multiple Testing Results ===")
print(paste("FDR Control Level:", alpha))
print(results)

# Summary statistics
print("\n=== Summary ===")
print(paste("Total hypotheses tested:", 9))
print(paste("Rejections (FDR/BH):", sum(bh_reject)))
print(paste("Rejections (Bonferroni):", sum(bonferroni_reject)))
print(paste("Rejections (Holm):", sum(holm_reject)))

# ----------------------------------------------------
# 5. Visualize results
# ----------------------------------------------------
# Create plotting data
plot_df <- data.frame(
  Population = factor(1:9),
  P_Value = p_values,
  BH_Threshold = (rank(p_values)/9) * alpha,
  Rejected_BH = bh_reject,
  Rejected_Bonf = bonferroni_reject
)
plot_df <- plot_df[order(plot_df$P_Value), ]
plot_df$Rank <- 1:9

# Plot 1: P-values with BH threshold line
p1 <- ggplot(plot_df, aes(x = Rank, y = P_Value)) +
  geom_point(aes(color = Rejected_BH), size = 3) +
  geom_line(aes(y = BH_Threshold), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = alpha/9, linetype = "dotted", color = "red") +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "green")) +
  labs(title = "FDR Control: Benjamini-Hochberg Procedure",
       subtitle = "Blue dashed line: BH threshold (α × rank/m)\nRed dotted line: Bonferroni threshold (α/m)",
       x = "Rank (sorted by p-value)",
       y = "P-value",
       color = "Rejected") +
  theme_minimal()

print(p1)

# Plot 2: Comparison of adjusted p-values
comparison_df <- data.frame(
  Population = rep(1:9, 3),
  Method = factor(c(rep("BH", 9), rep("Bonferroni", 9), rep("Holm", 9))),
  Adjusted_P = c(bh_p, bonferroni_p, holm_p)
)

p2 <- ggplot(comparison_df, aes(x = Population, y = Adjusted_P, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = alpha, linetype = "dashed") +
  labs(title = "Comparison of Multiple Testing Corrections",
       y = "Adjusted P-value",
       x = "Population") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")

print(p2)

# ----------------------------------------------------
# 6. Interpretation and discussion
# ----------------------------------------------------
cat("\n\n=== INTERPRETATION ===\n")
cat("The Benjamini-Hochberg procedure controls the expected proportion of\n")
cat("false discoveries among all rejected hypotheses. At FDR = 0.05, we expect\n")
cat("that no more than 5% of our discoveries are false positives.\n\n")

if(sum(bh_reject) > 0) {
  cat("Populations identified as significant by FDR:",
      paste(which(bh_reject), collapse = ", "), "\n")
  cat("Expected false discoveries among these: ~",
      round(sum(bh_reject) * alpha, 1), "\n")
}

if(sum(bonferroni_reject) < sum(bh_reject)) {
  cat("\nBonferroni is more conservative - it controls the probability of ANY\n")
  cat("false positive, which is why it identifies fewer populations.\n")
}

# ----------------------------------------------------
# 7. Additional: Discrete FDR (more powerful for discrete tests)
# ----------------------------------------------------
# Fisher's exact test produces discrete p-values, so we can use
# the DiscreteFDR package for potentially more power [citation:8]

if(require(DiscreteFDR, quietly = TRUE)) {
  # This would require additional steps to get the full distribution
  # But conceptually, it demonstrates advanced methods
  cat("\nNote: For discrete p-values (like from Fisher's exact test),\n")
  cat("specialized procedures like DiscreteFDR can provide more power\n")
  cat("by accounting for the discrete nature of the test statistics.\n")
}
