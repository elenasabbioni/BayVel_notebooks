# ---------------------------------------------
# Script to produce Table 5 of BayVel paper, with the values of WAIC for all the different BayVel models
#
# OUTPUT: .txt file with the LaTeX code to reproduce table 5
#
# INPUTS:
# - pathToYourDirectory: path to the working directory
#
# DEPENDECIES
# - External functions loaded from "functions.R" (e.g. u(), s(), u0(), s0())
# - packages: xtable
#             data.table
# ---------------------------------------------

rm(list = ls())
seed = 1234
set.seed(seed)

# -----------------------------
#  PATH 
# -----------------------------
# Set working directory and load file with auxiliary functions
pathToYourDirectory <- "pathToYourDirectory"
pathToBayVel_results <- paste0(pathToYourDirectory, "/real data")
setwd(paste0(pathToYourDirectory))
source(paste0(pathToYourDirectory, "/functions.R"))

# Set the output path where you will save the table
pathOutput <- paste0(pathToYourDirectory, "/tablesPaper/")
# -----------------------------
#  PACKAGES 
# -----------------------------
library(xtable)
library(data.table)

# -----------------------------
# TABLE 5
# Combine results of WAIC in the different models (obtained with the file "computeWAIC.r")
# -----------------------------
SW_vec <- c("SW1", "SW2")
T_vec <- c("T1", "T2", "T3")
g <- expand.grid(SW_vec, T_vec)
g$Var3 <- "D4"
combinations <- paste(g$Var1, g$Var2, g$Var3, sep = "-")

tableWAIC <- data.frame(mod = combinations, WAIC = rep(NA, length(combinations)))

for(i in combinations){
    tmp <- fread(paste0(pathToBayVel_results, "/", i, "/output/WAIC.csv"))
    tableWAIC[tableWAIC$mod == i, ]$WAIC = tmp
    tableWAIC[tableWAIC$mod == i,]$mod = nameReal_Pancreas(i)
}

# save the code to generate the Table in TeX
sink(paste0(pathOutput, "/Tab5.txt"))
print(xtable(tableWAIC))
sink()

