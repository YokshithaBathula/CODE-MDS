# CODE-MDS
HudsonAlpha Characterizing Our DNA Exceptions (CODE) project - R script for plotting Molecular Dynamic Simulations (MDS) of wild-type and variant protein

Molecular dynamics is a simulation methodology often used for studying the conformational rearrangements of molecules and their interactions with other molecular species in a range of environments. 

YASARA has built in macros to simulate how the protein will move inside the cell. Some areas of a protein will move quite a bit (the ends of the protein are often unstructured and “floppy”). Some areas are very stable. The stable areas are important because they are more likely to be the location of interaction with other molecules.

If you ran a MDS macro through YASARA or through HudsonAlpha, you should receive 2 folders, one containing the files for the wild-type model and one for the variant model. Each folder will contain 800 simulation files, images of the protein, and summary reports. We are interested in the "analysis.tab" - that contains the Root Mean Square Deviations (RMSD) of the protein over 20 nanoseconds, and "analysisres.tab" files - that contain the Root Mean Square Fluctuations (RMSF) of the amino acids in the protein. 

RMSD is a metric that assesses the degree of similarity between two structures. RMSD is a measure of the deviation of the positions of atoms in a molecule over the course of a molecular dynamics simulation relative to a reference structure. The RMSD calculation involves aligning the coordinates of a molecule at each time step to a reference structure and calculating the average distance between the aligned atoms. The resulting value provides information about how much the structure has deviated from the reference structure over the simulation time. In particular, the RMSD data can be used to assess the conformational changes of the molecule over time, the stability of the simulation, and the degree of convergence of the simulation. For example, a low RMSD value indicates that the structure of the molecule has not changed significantly during the simulation, while a high RMSD value suggests that the molecule has undergone significant conformational changes.

RMSF is a measurement similar to RMSD, but rather than measuring whole structural differences over time, RMSF calculates the movement of each individual residue during a simulation. Therefore, RMSF provides an understanding of the flexibilty of a residue in a protein. Higher RMSF values indicate a flexible region, while lower values correspond to a more rigid area in the protein. 

Here's a guide to creating RMSD and RMSF plots of your protein's Molecular Dynamic Simulations using RStudio. 

---
![RMSD 1](RStudio_images/RMSD%201.png)

Open RStudio. <br>
Install "data.table" package and load the package.
```{r}
install.packages("data.table")
library(data.table)
```

---
![RMSD 2](RStudio_images/RMSD%202.png)

Import your wild-type and variant "analysis.tab" files. <br>
Tab files are separated by a tab-delimited text files. We need to convert these into a table format. To do that, we will import them without a header and assign the column names in the next step.
```{r}
WT <- fread("path_to_the_file", header = FALSE)
Variant <- fread("path_to_the_file", header = FALSE)
```

Click on WT or Variant in the top right panel to view the table.

---
![RMSD 3](RStudio_images/RMSD%203.png)

---
The wild-type analysis has 14 columns and the variant analysis has 11 columns (which is missing the last 3 columns - average RMSD, average backbone, average heavy atoms).

Assign each column its name:
```{r}
colnames(WT) <- c("Time", "Energy_kJmol", "Bond", "Angle", "Dihedral", "Planarity", "Coulomb", "VdW", "RMSD", "Backbone", "HeavyAtoms", "RMSD_Avg", "Backbone_Avg", "HeavyAtoms_Avg")
colnames(Variant) <- c("Time", "Energy_kJmol", "Bond", "Angle", "Dihedral", "Planarity", "Coulomb", "VdW", "RMSD", "Backbone", "HeavyAtoms")
```

![RMSD 4](RStudio_images/RMSD%204.png)

---

Since we are only interested in RMSD over time, let's subset only those two columns
```{r}
TP53 <- WT[c(1:801),c(1,9)]
TP53_R158Q <- Variant[c(1:804), c(1,9)]
```

![RMSD 7](RStudio_images/RMSD%207.png)

---

Plot a line graph of the wild-type subset data titled "RMSD" with Time (in picoseconds) as the x-axis, RMSD as the y-axis, and the line in red
```{r}
plot(TP53$Time, TP53$RMSD, type = "l", main = "RMSD", xlab = "Time (ps)", ylab = "RMSD", col = "red")
```

![RMSD 9](RStudio_images/RMSD%209.png)

---

Add the variant subset data to the previous graph
```{r}
lines(TP53_R158Q$Time, TP53_R158Q$RMSD, type = "l")
```

![RMSD 10](RStudio_images/RMSD%2010.png)

---

Add a legend in the bottom-right corner to identify the two lines
```{r}
labels <- c("TP53", "TP53_R158Q")
legend("bottomright", legend = labels, col = c("red", "black"), lty = c(1, 1), cex = 0.6, inset = 0.1)
```

![RMSD 11](RStudio_images/RMSD%2011.png)

---

For RMSF analysis, we will use *ggplot2* package to graph the movement of the amino acids. If not already installed, install *ggplot2* and run the package. 
```{r}
install.packages("ggplot2")
library(ggplot2)
```

Import your wild-type and variant "analysisres.tab" files as a table
```{r}
RMSF_WT <- fread("~/Desktop/TP53_analysisres.tab", header = FALSE)
RMSF_Variant <- fread("~/Desktop/TP53_R158Q_em_analysisres.tab", header = FALSE)
```

![RMSF 2](RStudio_images/RMSF%202.png)

---

Assign column names to each column
```{r}
colnames(RMSF_WT) <- c("Residue", "Position no", "Letter", "RMSDs", "Backbone", "HeavyAtoms", "RMSF")
colnames(RMSF_Variant) <- c("Residue", "Position no", "Letter", "RMSDs", "Backbone", "HeavyAtoms", "RMSF")
```

![RMSF 3](RStudio_images/RMSF%203.png)

---

Since we are plotting RMSF's for each amino acid, let's isolate the amino acid position and RMSF columns. 

[!NOTE]  
Amino acids in the wild-type dataset starts on row 43

```{r}
RMSF_TP53 <- RMSF_WT[c(43:236),c("Position no","RMSF")]
RMSF_TP53_R158Q <- RMSF_Variant[c(1:194),c("Position no","RMSF")]
```

![RMSF 5](RStudio_images/RMSF%205.png)

---

Combine both subsets into one table. Add a new column called *dataset* in which all wild-type values are assigned "RMSF_TP53" and all variant values are assigned "RMSF_TP53_R158Q"
```{r}
combined_data <- rbind(RMSF_TP53, RMSF_TP53_R158Q)
combined_data$dataset <- rep(c("RMSF_TP53", "RMSF_TP53_R158Q"), each = nrow(RMSF_TP53))
```

![RMSF 6](RStudio_images/RMSF%206.png)

---

Plot a line graph of the combined data titled "RMSF" with *position no* in the x-axis, *RMSF* values in the y-axis, and the color of the lines differentiated based on the grouping in the *dataset* column
```{r}
ggplot(combined_data, aes(x = `Position no`, y = RMSF, color = dataset, group = dataset)) +
  geom_line() +
  labs(x = "Residue", y = "RMSF", title = "RMSF") +
  theme_minimal()
  ```
  
  ![RMSF 7](RStudio_images/RMSF%207.png)
  
  ---