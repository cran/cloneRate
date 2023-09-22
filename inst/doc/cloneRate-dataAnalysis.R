## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  fig.align = "center",
  collapse = TRUE,
  comment = "#>",
  out.width = "80%",
  dpi = 150
)

## ----setup, results=FALSE, message = FALSE------------------------------------
# Load and attach our package cloneRate
library(cloneRate)

# Load and attach ape
library(ape)

# Install from CRAN if necessary, then load and attach with library()
library(ggplot2)
library(survival)
library(ggsurvfit)
library(car)

## ----setColors----------------------------------------------------------------
colorPal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## ----dataPreview--------------------------------------------------------------
summary(cloneRate::realCloneData)

## ----dataSplit----------------------------------------------------------------
fullTrees.list <- cloneRate::realCloneData$fullTrees
cloneTrees.list <- cloneRate::realCloneData$cloneTrees

## ----classPhylo---------------------------------------------------------------
# Let's see the class of the tree
PD9478 <- fullTrees.list$PD9478_1
class(PD9478)

## ----firstTreePlot, fig.asp = 0.8, fig.width = 5.5----------------------------
# Plot the full tree from individual PD9478 at timepoint 1
PD9478 <- fullTrees.list$PD9478_1

# Plot, then add scale and title
plot.phylo(PD9478,
  direction = "downwards",
  show.tip.label = FALSE, edge.width = 2,
  edge.color = c(rep("black", 2), "red", "#0070C0", rep("black", (200)))
)
axisPhylo(side = 2, backward = FALSE, las = 1)
title(main = "PD9478 full tree", ylab = "Time (years)")

## ----plotSubclone, fig.asp = 0.8, fig.width = 5.5-----------------------------
# Load the clone tree from our cloneTrees.list
PD9478_subClone <- cloneTrees.list[["PD9478_1_clone1"]]

# Plot, then add scale and title
plot.phylo(PD9478_subClone,
  direction = "downwards",
  show.tip.label = FALSE, edge.width = 2,
)
axisPhylo(side = 2, backward = FALSE, las = 1)
title(main = "PD9478 JAK2 & DNMT3A clone tree", ylab = "Time (years)")

## ----applyMethods-------------------------------------------------------------
# Get maximum likelihood and internal lengths estimates
print(round(maxLikelihood(PD9478_subClone)$estimate, 3))
print(round(internalLengths(PD9478_subClone)$estimate, 3))

## ----applyLengths-------------------------------------------------------------
# Apply each of our estimates
resultsLengths <- internalLengths(cloneTrees.list)

## ----plotNonBinary, fig.asp = 0.8, fig.width = 5.5----------------------------
# Plot, highlighting the 3 descendants of node 76. Then add scale and title
plot.phylo(cloneTrees.list$PD5847_1_clone1,
  edge.width = 2,
  direction = "downwards", show.tip.label = FALSE,
  edge.color = c(
    rep(1, 4), "darkorange", rep(1, 34),
    "darkorange", rep(1, 6), "darkorange", rep(1, 1000)
  )
)
axisPhylo(side = 2, backward = FALSE, las = 1)
title(main = "PD5847 clone 1 tree (multi-branching)", ylab = "Time (years)")

## ----applyMaxLike-------------------------------------------------------------
resultsMaxLike <- suppressWarnings(maxLikelihood(cloneTrees.list))

# Preview the output
print(head(resultsMaxLike[, c(1:3)]))

# Print correlation coefficient of the estimates from the two methods
print(stats::cor.test(resultsLengths$estimate, resultsMaxLike$estimate)$estimate)

## ----applyBDMCMC, eval=FALSE--------------------------------------------------
#  # Longer running code. Not run as part of vignette.
#  resultsBDMCMC <- suppressWarnings(birthDeathMCMC(cloneTrees.list))
#  
#  # Print correlation coefficient of the estimates from the two methods
#  print(stats::cor.test(resultsBDMCMC$estimate, resultsMaxLike$estimate)$estimate)
#  # >  cor
#  # > 0.9964975

## ----columnsInfo--------------------------------------------------------------
# Check that column names are the same
stopifnot(all(colnames(resultsMaxLike) == colnames(resultsLengths)))

# Print column names
colnames(resultsMaxLike)

## ----viewMetdata--------------------------------------------------------------
# Show the metadata of an individual without MPN
print(cloneTrees.list$PD34493_clone1$metadata)

# And the metadata of an individual with MPN
print(cloneTrees.list$PD9478_1_clone1$metadata)

## ----combineMetdata-----------------------------------------------------------
# Combine all metadata into a single data.frame
metadataAll <- do.call(rbind, lapply(cloneTrees.list, function(x) {
  return(x$metadata)
}))

## ----joinMetaResults----------------------------------------------------------
# Combine metadata with estimates using cbind. Check if cloneNames match
resultsLengthsMeta <- cbind(resultsLengths, metadataAll)
stopifnot(resultsLengthsMeta$cloneName_result == resultsLengthsMeta$cloneName_meta)

resultsMaxLikeMeta <- cbind(resultsMaxLike, metadataAll)
stopifnot(resultsLengthsMeta$cloneName_result == resultsLengthsMeta$cloneName_meta)

# Because max. likelihood performs slightly better, use that going forward
results <- resultsMaxLikeMeta

## ----splitMalNorm-------------------------------------------------------------
# If patient ID and clone number are the same, even if timepoint
# differs, then we have a duplicate.
# Make a new column for patient, removing anything after "_" from cloneName_result
results$patient <- gsub("_.*", "", results$cloneName_result)
results$cloneNumber <- gsub(".*_clone", "", results$cloneName_result)

# Combining patient ID (without timepoint) and clone number
# will give us a unique ID for the clone regardless of sampling time
results$uniqueCloneID <- paste0(results$patient, "_", results$cloneNumber)

# Find which clone IDs appear twice
tmp <- table(results$uniqueCloneID)
repeatsVec <- names(tmp)[tmp == 2]

# Remove the duplicate with fewer number of cells sequenced, n
rowsRemove <- c()
for (cloneID in repeatsVec) {
  duplicateRows <- which(results$uniqueCloneID == cloneID)
  removeIndex <- duplicateRows[which.min(results$n[duplicateRows])]
  rowsRemove <- c(rowsRemove, removeIndex)
}

uniqueResults <- results[!c(1:nrow(results)) %in% rowsRemove, ]
# Check that each unique clone now only appears once (no duplicates)
stopifnot(all(table(uniqueResults$uniqueCloneID) == 1))

## ----meanInd, out.width = '50%'-----------------------------------------------
# Get the unique individuals, and initialize vectors
uniqueIndividuals <- unique(uniqueResults$patient)
individualMeans <- c()
malNorm <- c()

# Fill vectors with mean growth rate and MPN status of each individual
for (ind in uniqueIndividuals) {
  individualMeans <- c(individualMeans, mean(uniqueResults$estimate[uniqueResults$patient == ind]))
  malNorm <- c(malNorm, uniqueResults$malnorm[uniqueResults$patient == ind][1])
}

# Combine results into a data.frame
mal_vs_norm.df <- data.frame(
  "Patient" = uniqueIndividuals,
  "meanEstimate" = individualMeans,
  "malNorm" = malNorm
)

# Run a non-parametric Mann-whitney test to see if they're significantly different
nonparamTest <- wilcox.test(meanEstimate ~ malNorm, data = mal_vs_norm.df)
print(nonparamTest)

# Set factor ordering for plot and plot using ggplot (Fig. 5E)
mal_vs_norm.df$malNorm <- factor(mal_vs_norm.df$malNorm, levels = c("Normal", "Malignant"))
ggplot(mal_vs_norm.df, aes(x = malNorm, y = meanEstimate)) +
  geom_label(label = paste0("p = ", round(nonparamTest$p.value, 3)), x = 1.5, y = 1.2) +
  geom_boxplot(width = 0.3, aes(color = malNorm), outlier.shape = NA) +
  geom_jitter(aes(x = malNorm, y = meanEstimate, color = malNorm), width = .1, size = 2) +
  scale_color_manual(values = c("black", "red"), labels = c("Normal", "MPN")) +
  theme_bw() +
  ylab("Mean net growth rate (r)") +
  xlab("") +
  theme(
    legend.position = "none",
    axis.ticks.x = element_blank(), axis.title.x = element_blank()
  )

## ----multiDriver, out.width = '50%'-------------------------------------------
# Set the multi driver variable according to whether multiple drivers are annotated
uniqueResults$multiDriver <- "Single or unknown"
uniqueResults$multiDriver[grepl("AND", uniqueResults$cloneDriver)] <- "Multiple"

# Run a non-parametric test to see if the differences are significant
nonParamTest <- wilcox.test(estimate ~ multiDriver, data = uniqueResults)
print(nonParamTest)

# Set factor levels to control plot order
uniqueResults$multiDriver <- factor(uniqueResults$multiDriver, levels = c("Single or unknown", "Multiple"))

# Plot (Fig. 5F)
ggplot(uniqueResults, aes(x = multiDriver, y = estimate)) +
  geom_boxplot(width = 0.3, aes(color = multiDriver), outlier.shape = NA) +
  geom_jitter(aes(x = multiDriver, y = estimate, color = multiDriver), width = .1) +
  scale_color_manual(values = c(colorPal[4], colorPal[8])) +
  theme_bw() +
  geom_label(label = paste0("p=", round(nonParamTest$p.value, 7)), x = 1.5, y = 2) +
  ylab("Net growth rate (r)") +
  xlab("") +
  ggtitle("Driver mutations") +
  theme(
    legend.position = "none", axis.ticks.x = element_blank(),
    axis.title.x = element_blank(), plot.title = element_text(hjust = .5)
  )

## ----getMax-------------------------------------------------------------------
# Get the clone IDs of the fittest clone in each patient
cloneIDs_max <- sapply(unique(uniqueResults$patient), function(x) {
  patient.df <- uniqueResults[uniqueResults$patient == x, ]
  cloneID <- patient.df$cloneName_result[which.max(patient.df$estimate)]
  cloneID
})

# Get the clone IDs of the youngest clone in each patient
cloneIDs_youngest <- sapply(unique(uniqueResults$patient), function(x) {
  patient.df <- uniqueResults[uniqueResults$patient == x, ]
  cloneID <- patient.df$cloneName_result[which.min(patient.df$cloneAgeEstimate)]
  cloneID
})

# See how much overlap there is between the most recent clone and the highest growth rate
table(cloneIDs_max == cloneIDs_youngest)

## ----maxStratify, out.width = '50%'-------------------------------------------
# Subset the results to only have the fittest clone for each patient
maxResults <- uniqueResults[uniqueResults$cloneName_result %in% cloneIDs_max, ]

# Subset to the MPN patients
malMax <- maxResults[maxResults$malnorm == "Malignant", ]

# Determine time from clone initiation to diagnosis
malMax$latency_to_dx <- malMax$cloneAgeEstimate - (malMax$age - malMax$diagnosis.age)

# Stratify the population by the mean growth rate, adding a column "aboveMean"
malMax$aboveMean <- malMax$estimate > mean(malMax$estimate)

# Use the survival package to see if differences are significant
malMax$status <- 1
survivalTest <- survival::survdiff(survival::Surv(time = latency_to_dx, event = status) ~ aboveMean, data = malMax)
print(survivalTest$pvalue)

# Plot survival curves
survfit2(Surv(time = latency_to_dx, event = status) ~ aboveMean, data = malMax) %>%
  ggsurvfit() +
  labs(
    y = "Probability diagnosis-free",
    x = "Time from clone initiation (yrs.)"
  ) + scale_color_manual(values = c(colorPal[7], colorPal[3])) +
  theme(legend.position = "none") +
  geom_label(label = paste0("p = ", round(survivalTest$pvalue, 4)), x = 45, y = .85, color = "black", fill = "white") +
  geom_label(label = "r < mean", x = 40, y = .375, color = colorPal[7], fill = "white") +
  geom_label(label = "r > mean", x = 12, y = .35, color = colorPal[3], fill = "white")

## ----youngestStratify---------------------------------------------------------
youngestResults <- uniqueResults[uniqueResults$cloneName_result %in% cloneIDs_youngest, ]

# Subset to the MPN patients
malYoungest <- youngestResults[youngestResults$malnorm == "Malignant", ]

# Check how much overlap there is between youngest clone and fittest clone in the MPN dataset
table(malYoungest$cloneName_result == malMax$cloneName_result)

## ----consistency, out.width = '50%'-------------------------------------------
# Like before, we can find which clone IDs appear twice
tmp <- table(results$uniqueCloneID)
repeatsVec <- names(tmp)[tmp == 2]

# Subset the results to include only repeated clones
repeatResults <- results[results$uniqueCloneID %in% repeatsVec, ]

# Plot the repeated clones
ggplot(repeatResults) +
  geom_pointrange(data = repeatResults, aes(
    x = cloneName_result, y = estimate, ymin = lowerBound,
    ymax = upperBound, color = uniqueCloneID
  )) +
  scale_color_manual(values = colorPal[c(1, 6, 7)]) +
  theme_bw() +
  ylab("Net growth rate estimate (r)") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.title = element_blank())

## ----loadLong, out.width = '50%'----------------------------------------------
# Get data from PD9478
PD9478_long <- longitudinalData[longitudinalData$Sample.ID == "PD9478", ]

# Get a rough estimate for starting params
startingParams <- suppressWarnings(coef(lm(logit(VAF / .5) ~ Age, data = PD9478_long)))

# Fit to a three parameter logistic curve
fit <- nls(VAF ~ K / (1 + exp(-(phi2 + r * Age))),
  start = list(K = 0.2, phi2 = min(c(-0.000001, startingParams[1])), r = max(c(0.00001, startingParams[2]))),
  data = PD9478_long, trace = FALSE, algorithm = "port", lower = c(0, -500, 0.00001), upper = c(0.5, 0, 5)
)

# Set output equal to summary of longitudinal models
output <- summary(fit)

# Assign fitted params
r <- output$coefficients["r", 1]
t_m <- -output$coefficients["phi2", 1] / r # Get midpoint time by dividing by -1/r
K <- output$coefficients["K", 1]

# Get CI for growth rate of longitudinal/logistic model
stdError <- output$coefficients["r", "Std. Error"]
lb <- r - 1.96 * stdError
ub <- r + 1.96 * stdError

# Prepare a df for plotting the fit line along with the data
x <- c((min(PD9478_long$Age) - 10):(max(PD9478_long$Age) + 10)) # construct a range of x values bounded by the data
y <- K / (1 + exp(-(x - t_m) * r)) # curve VAF (3 param model)
predict.df <- data.frame("x" = x, "y" = y)

# Set color for longitudinal fit info
fitColor <- colorPal[6]

# Plot
ggplot(PD9478_long, aes(x = Age, y = VAF)) +
  theme_bw() +
  coord_cartesian(xlim = c(min(x), max(x)), ylim = c(-0.01, 0.52), expand = 0) +
  labs(x = "Person Age (yr)", y = "Variant Allele Frequency (VAF)") +
  ggtitle(gsub("_", " & ", paste0(PD9478_long$Sample.ID[1], " ", PD9478_long$Gene[1]))) +
  geom_line(data = predict.df, aes(x = x, y = y), color = fitColor, linewidth = 1, show.legend = TRUE) +
  geom_point(color = "#808080", shape = 18, size = 1.5) +
  geom_vline(xintercept = t_m, linetype = "dotted", color = fitColor, linewidth = .6) +
  geom_hline(yintercept = K, linetype = "dotted", color = fitColor, linewidth = .6)

## ----validateLong, out.width = '50%'------------------------------------------
# Get the single cell estimate from PD9478 clone 1
scPD9478 <- results[results$cloneName_result == "PD9478_1_clone1", ]

# Combine orthogonal estimates into one df for plotting
combine.df <- data.frame(
  "Clone" = "PD9478_1_clone1",
  "r" = c(r, scPD9478$estimate),
  "lowerBound" = c(lb, scPD9478$lowerBound),
  "upperBound" = c(ub, scPD9478$upperBound),
  "method" = c("longitudinal", "max. likelihood")
)

# Plot
ggplot(combine.df) +
  geom_pointrange(aes(x = method, y = r, ymin = lowerBound, ymax = upperBound)) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  ylab("Net growth rate (r)")

