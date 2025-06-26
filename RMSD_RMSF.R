#RMSD

#Install and load package
install.packages("data.table")
library(data.table)

#Import .tab files as tables without a header row
WT <- fread("~/Desktop/TP53_analysis.tab", header = FALSE)
Variant <- fread("~/Desktop/TP53_R158Q_em_analysis.tab", header = FALSE)

#Assign column names
colnames(WT) <- c("Time", "Energy_kJmol", "Bond", "Angle", "Dihedral", "Planarity", "Coulomb", "VdW", "RMSD", "Backbone", "HeavyAtoms", "RMSD_Avg", "Backbone_Avg", "HeavyAtoms_Avg")
colnames(Variant) <- c("Time", "Energy_kJmol", "Bond", "Angle", "Dihedral", "Planarity", "Coulomb", "VdW", "RMSD", "Backbone", "HeavyAtoms")

#Subset Time and RMSD columns
TP53 <- WT[c(1:801),c("Time","RMSD")]
TP53_R158Q <- Variant[c(1:804), c("Time","RMSD")]

#Plot a line graph of the wild-type subset
plot(TP53$Time, TP53$RMSD, type = "l", main = "RMSD", xlab = "Time (ps)", ylab = "RMSD", col = "red")

#Add variant subset graph to the above plot
lines(TP53_R158Q$Time, TP53_R158Q$RMSD, type = "l")

#Add a legend to identify the two lines
labels <- c("TP53", "TP53_R158Q")
legend("bottomright", legend = labels, col = c("red", "black"), lty = c(1, 1), cex = 0.6, inset = 0.1)


#RMSF

#Install and load ggplot2
install.packages("ggplot2")
library(ggplot2)

#Import .tab files as tables without a header row
RMSF_WT <- fread("~/Desktop/TP53_analysisres.tab", header = FALSE)
RMSF_Variant <- fread("~/Desktop/TP53_R158Q_em_analysisres.tab", header = FALSE)

#Assign column names
colnames(RMSF_WT) <- c("Residue", "Position no", "Letter", "RMSDs", "Backbone", "HeavyAtoms", "RMSF")
colnames(RMSF_Variant) <- c("Residue", "Position no", "Letter", "RMSDs", "Backbone", "HeavyAtoms", "RMSF")

#Isolate amino acid position and RMSF columns 
RMSF_TP53 <- RMSF_WT[c(43:236),c("Position no","RMSF")]
RMSF_TP53_R158Q <- RMSF_Variant[c(1:194),c("Position no","RMSF")]

#Combine both subsets into one table. 
combined_data <- rbind(RMSF_TP53, RMSF_TP53_R158Q)

#Add a third column called "dataset" which identifies and separates the two sets
combined_data$dataset <- rep(c("RMSF_TP53", "RMSF_TP53_R158Q"), each = nrow(RMSF_TP53))

#Plot a line graph with legend
ggplot(combined_data, aes(x = `Position no`, y = RMSF, color = dataset, group = dataset)) +
  geom_line() +
  labs(x = "Residue", y = "RMSF", title = "RMSF") +
  theme_minimal()
