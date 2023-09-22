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

# Load and attach ape, which will be installed if you've installed cloneRate
library(ape)

# Install ggplot2 from CRAN if necessary, then load and attach it with library()
library(ggplot2)

## ----setColors----------------------------------------------------------------
colorPal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## ----simOne-------------------------------------------------------------------
# Generate the tree
tree <- simUltra(a = 1, b = 0, cloneAge = 20, n = 100, addStem = TRUE)

# We see that the tree is of class phylo
class(tree)

# Preview the tree
print(tree)

## ----previewMeta--------------------------------------------------------------
print(tree$metadata)

## ----plotOne, fig.asp = 0.8, fig.width = 5.5----------------------------------
plot.phylo(tree, direction = "downwards", show.tip.label = FALSE)
axisPhylo(side = 2, backward = FALSE, las = 1)
title(main = "Simulated tree", ylab = "Time")

## ----estimateOne--------------------------------------------------------------
maxLikelihood(tree)$estimate

internalLengths(tree)$estimate

## ----sim100-------------------------------------------------------------------
ptm <- proc.time()
tree.list <- simUltra(a = 1, b = 0, cloneAge = 20, n = 100, nTrees = 100, nCores = 1, addStem = TRUE)
print(proc.time()[["elapsed"]] - ptm[["elapsed"]])

## ----applyAll-----------------------------------------------------------------
resultsMaxLike <- maxLikelihood(tree.list)

resultsLengths <- internalLengths(tree.list)

## ----plotEstimates, fig.asp = 0.8, fig.width = 5.5----------------------------
# Combine for ggplot formatting
resultsCombined <- rbind(resultsLengths, resultsMaxLike)

# Plot, adding a vertical line at r=1 because that's the true growth rate
ggplot(resultsCombined) +
  geom_density(aes(x = estimate, color = method), linewidth = 1.5) +
  geom_vline(xintercept = exampleUltraTrees[[1]]$metadata$r) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    legend.title = element_blank()
  ) +
  xlab("Net growth rate estimate (r)") +
  ylab("Density") +
  scale_color_manual(labels = c("Internal lengths", "Max. likelihood"), values = c("black", "#009E73"))

## ----ultra2mut----------------------------------------------------------------
# Set mutation rate equal to 10 muts/year for all trees
mutTree.list <- ultra2mut(tree.list, nu = 10)

## ----simMut, eval = FALSE-----------------------------------------------------
#  # Set params for ultra tree + a mutation rate
#  mutTree.list2 <- simMut(a = 1, b = 0, cloneAge = 20, n = 100, nTrees = 100, nu = 100, nCores = 1)

## ----applyShared--------------------------------------------------------------
resultsShared <- sharedMuts(mutTree.list, nu = 10)

## ----plotAll, fig.asp = 0.8, fig.width = 5.5----------------------------------
# Combine the columns with the estimates
colsUse <- c("lowerBound", "estimate", "upperBound", "method")
resultsAll <- rbind(resultsShared[, colsUse], resultsCombined[, colsUse])

# Plot, adding a vertical line at r=1 because that's the true growth rate
ggplot(resultsAll) +
  geom_density(aes(x = estimate, color = method), linewidth = 1.5) +
  geom_vline(xintercept = exampleUltraTrees[[1]]$metadata$r) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    legend.title = element_blank()
  ) +
  xlab("Net growth rate estimate (r)") +
  ylab("Density") +
  scale_color_manual(labels = c("Internal lengths", "Max. likelihood", "Shared Muts."), values = c(colorPal[1], colorPal[4], colorPal[6]))

## ----simRMSE------------------------------------------------------------------
# Calculate the RMSE
groundTruth <- 1
rmse <- unlist(lapply(
  split(resultsAll, resultsAll$method),
  function(x) {
    sqrt(sum((x$estimate - groundTruth)^2) / length(x$estimate))
  }
))

print(rmse)

## ----simMean------------------------------------------------------------------
# Calculate the mean
simMean <- unlist(lapply(
  split(resultsAll, resultsAll$method),
  function(x) {
    mean(x$estimate)
  }
))

print(simMean)

## ----simSD--------------------------------------------------------------------
# Calculate the standard deviation of our estimates
simSD <- unlist(lapply(
  split(resultsAll, resultsAll$method),
  function(x) {
    sd(x$estimate)
  }
))

print(simSD)

## ----simCoverage--------------------------------------------------------------
# Calculate the coverage probability of our estimates
groundTruth <- 1
simCoverage <- unlist(lapply(
  split(resultsAll, resultsAll$method),
  function(x) {
    # Set insideInterval to TRUE if inside interval and FALSE if outside
    insideInterval <- (x$lowerBound < groundTruth & x$upperBound > groundTruth)
    # Count TRUE as 1 and FALSE as 0. See the fraction inside interval
    sum(insideInterval) / length(insideInterval)
  }
))

print(simCoverage)

## ----gen_n_vec----------------------------------------------------------------
n_vec <- c(rep(10, 100), rep(30, 100))

## ----a_vec--------------------------------------------------------------------
a_vec <- stats::runif(n = 200, min = 1, max = 4)
b_vec <- a_vec - 1

## ----vary_n-------------------------------------------------------------------
n10_n30_trees <- simUltra(
  a = a_vec, b = b_vec, cloneAge = 20, n = n_vec,
  nTrees = length(n_vec)
)

## ----mergeOutput--------------------------------------------------------------
# Apply our estimates for ultrametric trees
n10_n30_maxLike <- maxLikelihood(n10_n30_trees)
n10_n30_lengths <- internalLengths(n10_n30_trees)

# Combine the estimates
n10_n30_both <- rbind(n10_n30_maxLike, n10_n30_lengths)

## ----mergeAll-----------------------------------------------------------------
# Merge the results data.frames
results_vary_n <- rbind(n10_n30_both, resultsCombined)

## ----plot_vary_n, fig.asp = 0.8, fig.width = 5.5------------------------------
ggplot(results_vary_n) +
  geom_density(aes(x = estimate, color = method), linewidth = 1.5) +
  geom_vline(xintercept = exampleUltraTrees[[1]]$metadata$r) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    legend.title = element_blank()
  ) +
  xlim(0, 3) +
  xlab("Net growth rate estimate (r)") +
  ylab("Density") +
  scale_color_manual(
    labels = c("Internal lengths", "Max. likelihood"),
    values = c("black", "#009E73")
  ) +
  facet_wrap(~ factor(paste0("n = ", n), levels = c("n = 10", "n = 30", "n = 100")),
    ncol = 1, strip.position = "bottom", scales = "free", dir = "v"
  )

## ----sim_vary_r---------------------------------------------------------------
# First tree, r = a - b = 2
tree1 <- simUltra(a = 2.5, b = .5, cloneAge = 30, n = 50)

# Second tree, with r = a - b = 0.5
tree2 <- simUltra(a = 1, b = .5, cloneAge = 30, n = 50)

# Third tree, with r = 0.5 but cloneAge = 120
tree3 <- simUltra(a = 1, b = .5, cloneAge = 120, n = 50)

## ----plot_vary_r, fig.asp = 0.4, fig.width = 5.5------------------------------
oldpar <- graphics::par(mfrow = c(1, 3))
ape::plot.phylo(tree1, direction = "downwards", show.tip.label = F, main = "r = 2, cloneAge = 30")
ape::plot.phylo(tree2, direction = "downwards", show.tip.label = F, main = "r = 0.5, cloneAge = 30")
ape::plot.phylo(tree3, direction = "downwards", show.tip.label = F, main = "r = 0.5, cloneAge = 120")
# reset par settings
graphics::par(oldpar)

## ----vary_r, eval = FALSE-----------------------------------------------------
#  # Uniform ditribution of r used to generate b_vec
#  r_vec <- stats::runif(n = 1000, min = 0.1, max = 1)
#  a_vec <- stats::runif(n = 1000, min = 1, max = 3)
#  b_vec <- a_vec - r_vec
#  
#  # Input to simUltra()
#  vary_r_trees <- simUltra(a = a_vec, b = b_vec, cloneAge = 40, n = 50, nTrees = length(a_vec))

## ----slowClone, fig.asp = 0.8, fig.width = 5.5--------------------------------
# Let's call it slowClone
slowClone <- simUltra(a = .15, b = .05, n = 50, cloneAge = 40)

# Plot the tree
ape::plot.phylo(slowClone, direction = "downwards", show.tip.label = F)

## ----throwWarning-------------------------------------------------------------
# Apply our estimates
maxLikelihood(slowClone)
internalLengths(slowClone)

