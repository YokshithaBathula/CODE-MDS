install.packages("data.table")
library(data.table)
WT <- fread("~/Desktop/TP53_analysis.tab", header = FALSE)
Variant <- fread("~/Desktop/TP53_R158Q_em_analysis.tab", header = FALSE)
colnames(WT) <- c("Time", "Energy_kJmol", "Bond", "Angle", "Dihedral", "Planarity", "Coulomb", "VdW", "RMSD", "Backbone", "HeavyAtoms", "RMSD_Avg", "Backbone_Avg", "HeavyAtoms_Avg")
colnames(Variant) <- c("Time", "Energy_kJmol", "Bond", "Angle", "Dihedral", "Planarity", "Coulomb", "VdW", "RMSD", "Backbone", "HeavyAtoms")
TP53 <- WT[c(1:801),c("Time","RMSD")]
TP53_R158Q <- Variant[c(1:804), c("Time","RMSD")]
plot(TP53$Time, TP53$RMSD, type = "l", main = "RMSD", xlab = "Time (ps)", ylab = "RMSD", col = "red")
lines(TP53_R158Q$Time, TP53_R158Q$RMSD, type = "l")
labels <- c("TP53", "TP53_R158Q")
legend("bottomright", legend = labels, col = c("red", "black"), lty = c(1, 1), cex = 0.6, inset = 0.1)

library(ggplot2)
RMSF_WT <- fread("~/Desktop/TP53_analysisres.tab", header = FALSE)
RMSF_Variant <- fread("~/Desktop/TP53_R158Q_em_analysisres.tab", header = FALSE)
colnames(RMSF_WT) <- c("Residue", "Position no", "Letter", "RMSDs", "Backbone", "HeavyAtoms", "RMSF")
colnames(RMSF_Variant) <- c("Residue", "Position no", "Letter", "RMSDs", "Backbone", "HeavyAtoms", "RMSF")
#RMSF_TP53 <- RMSF_WT[c(43:236),c("Residue","RMSF")]
RMSF_TP53 <- RMSF_WT[c(43:236),c("Position no","RMSF")]
RMSF_TP53$Fluctuation <- "RMSF_WT"
#RMSF_TP53_R158Q <- RMSF_Variant[c(1:194),c("Residue","RMSF")]
RMSF_TP53_R158Q <- RMSF_Variant[c(1:194),c("Position no","RMSF")]
RMSF_TP53_R158Q$Fluctuation <- "RMSF_Variant"
combined_data <- rbind(RMSF_TP53, RMSF_TP53_R158Q)
combined_data$dataset <- rep(c("RMSF_TP53", "RMSF_TP53_R158Q"), each = nrow(RMSF_TP53))
ggplot(combined_data, aes(x = `Position no`, y = RMSF, color = Fluctuation, group = Fluctuation)) +
  geom_line() +
  labs(x = "Residue", y = "RMSF", title = "RMSF") +
  theme_minimal()
