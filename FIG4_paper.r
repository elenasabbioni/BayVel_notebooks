# ---------------------------------------------
# Script to produce Fig. 4 of BayVel paper
#
# OUTPUT: 6 .pdf files with the plots in the 3 columns of Figure 4
#
#
# INPUTS (set at the beginning or controlled in loops):
# - pathToYourDirectory: path to the working directory
#
# DEPENDECIES
# - External functions loaded from "functions.R" (e.g. u(), s(), u0(), s0())
# - packages: ggplot2
#             LaplacesDemon
#             cowplot
# ---------------------------------------------

rm(list = ls())
seed = 1234
set.seed(seed)

# -----------------------------
#  PATH 
# -----------------------------
# Set working directory and load file with auxiliary functions
pathToYourDirectory <- "pathToYourDirectory"
setwd(paste0(pathToYourDirectory))
pathToResults <- paste0(pathToYourDirectory, "/simulations/")
source(paste0(pathToYourDirectory, "/functions.R"))
# Set the output path where you will save the images
pathOutput <- paste0(pathToYourDirectory, "/figuresPaper/")


# -----------------------------
#  PACKAGES 
# -----------------------------
library(ggplot2)
library(LaplacesDemon)
library(cowplot)


# -----------------------------
# TYPE OF SIMULATION WE ARE CONSIDERING
# -----------------------------
n_genes <- 2000                   # number of genes
mcmc <- 150                       # number of iterations of the mcmc we have run
typeSIM <- "sim"                  # specify that we are considering simulated data
typeSW <- "SW1"                   # common or cluster-specific switching points
typeT <- "T1"                     # number of subgroups
typeD <- "D1"                     # type of data distribution
genesToPlot <- c(3, 12)           # indexes of the genes that we want to plot        

addBayVel <- TRUE                 # specify if you want to plot the results of BayVel


# -----------------------------
# LOAD THE DATA 
# -----------------------------
nameSim <- paste0(typeSW, "-", typeT, "-", typeD)

# environment where we will load the results from the different methods
bv <- new.env(parent = baseenv()) 

# --- BayVel results
maxLog <- FALSE
source(paste0(pathToYourDirectory, "/loadBayVel_results.r"))
iterToPlot <- which.max(bv$loglikTOT[,1])


# --- real parameters
real <- new.env(parent = baseenv())
pathReal <- paste0(pathToYourDirectory, "/simulations/", nameSim, "/", nameSim, "_", n_genes, ".RData")
load(paste0(pathReal), envir = real)


# upload the functions in order to check the version of the functions is the last one and they were not overwritten when the data were loaded.
source(paste0(pathToYourDirectory, "/functions.R"))


# -----------------------------
# Fig. 4, columns 1: variability of the steady states coordinates and of the switching point coordinates
# -----------------------------
for(g in genesToPlot){
    for(tyT0_off in 1:max(bv$typeCellT0_off)){

        df <- data.frame(
            s_SS_off = bv$sSS_off_chain[g,], 
            u_SS_off = bv$uSS_off_chain[g,], 
            s_SS_on = bv$sSS_on_chain[g,], 
            u_SS_on = bv$uSS_on_chain[g,], 
            s0_off = bv$s0_off_chain[tyT0_off, g, ], 
            u0_off = bv$u0_off_chain[tyT0_off, g, ]
        )

        # variability steady state coordinates and switching point
        gg <- ggplot() +
              geom_point(data = df, aes(x = s_SS_off, y = u_SS_off), col = "#FF9900", alpha = 0.2, size = 1) +
              geom_point(data = df, aes(x = s_SS_on, y = u_SS_on),   col = "#666666", alpha = 0.2, size = 1) +
              geom_point(data = df, aes(x = s0_off, y = u0_off),     col = "#3399FF", alpha = 0.2, size = 1) + 
              labs(x = "s", y = "u", title = paste0("Gene ", g))


        # real gene-dynamic
        gg <- plot_sVSu(t0_off = real$t0_off_real[tyT0_off, g], t0_on = 0, alpha = real$alpha_real[g,], beta = real$beta_real[g], gamma = real$gamma_real[g], pos_u = NA, pos_s = NA, g = g, subGrLabels = bv$subtypeCell, add = TRUE, colCell = NA, xlim = NA, ylim = NA, axisTitle.size = 20, axisText.size = 10, title.size = 30, colDyn = "red", lineSize = 1, shapePoint = 22, sizePoint = 3, gg = gg) 
        
        # BayVel gene-dynamic
        gg <- plot_sVSu(t0_off = bv$T0_off_chain[tyT0_off, g, iterToPlot], t0_on = 0, alpha = bv$alpha_chain[g, , iterToPlot], beta = bv$beta_chain[g, iterToPlot], gamma = bv$gamma_chain[g, iterToPlot], pos_u = NA, pos_s = NA, g = g, subGrLabels = bv$subtypeCell, add = TRUE, colCell = NA, , xlim = 10, ylim = 10, axisTitle.size = 20, axisText.size = 20, title.size = 30, colDyn = "darkgreen", lineSize = 1, shapePoint = 21, sizePoint = 3, gg = gg)

        if(g == 12){
            # uniform graphics parameters
            gg <- gg + scale_y_continuous(breaks= c(2, 5, 7, 10), labels = c("2", "5", "7", "10"))
        }           

        # personalize text dimension and font in the plot
        gg <- gg + 
            theme(plot.title = element_text(family = "serif", size=50,  hjust = 0.5), axis.text = element_text(family = "serif", size = 35), axis.title = element_text(family = "serif", size = 45), legend.text = element_text(family = "serif", size = 38), legend.title = element_text(family = "serif", size = 38), plot.margin = unit(c(0.2,0.5,0,0), "lines")) 

        # extract legend 
        gg2 <- gg + 
            geom_point(data = data.frame(x = bv$alpha_chain[g, 1, iterToPlot]/bv$gamma_chain[g, iterToPlot], y = bv$alpha_chain[g, 1, iterToPlot]/bv$beta_chain[g, iterToPlot], group = c("Lower steady state", "Upper steady state", "Switching point")), aes(x = x, y = y, fill = group), size = 0.1) + 
            scale_fill_manual(values = c("Lower steady state" = "#FF9900", "Upper steady state" = "#666666", "Switching point" = "#3399FF"), guide = guide_legend(order = 1, override.aes = list(size = 8, shape = 21), nrow = 3)) + 
            labs(fill = "Posterior samples of") + theme(legend.position = c(0.05, 0.3), legend.text = element_text(family = "serif", size = 38), legend.title = element_text(family = "serif", size = 40), legend.key = element_rect(colour = "transparent", fill = NA), legend.background=element_blank())

        leg2 <- cowplot::get_plot_component(gg2, "guide-box-inside")

        # Combine the plots using ggdraw and draw_plot
        final_plot <- ggdraw() +
                      draw_plot(gg, 0, 0, 1, 1) # Main plot takes full space

        if(g == 3){
            final_plot <- final_plot + 
            draw_plot(leg2, 0.41, 0.73, 0.2, 0.3)     # Inset 2 at bottom left
        }else if(g == 12){
            final_plot <- final_plot + 
                draw_plot(leg2, 0.45, 0.73, 0.2, 0.3)
        }

        pdf(paste0(pathOutput, "/Fig4_col1_gene", g, "_tyT0_off", tyT0, ".pdf"), width = 7*1.2, height = 7*1.5)
            plot(final_plot)
        dev.off() 
    }   
}


# -----------------------------
# Fig. 4, columns 1: estimated and real position of the subgroups coordinates
# -----------------------------
for(g in genesToPlot){
    for(tyT0_off in 1:max(bv$typeCellT0_off)){

        # graphical parameters that are different between the two plots
        if(g == 3){
            scaling_y <- c(1.01, 0.97, 1, 0.97, 1, 1, 1, 1, 0.98, 1.14)
            scaling_x <- c(1.12, 1.12, 0.8,0.75, 1.15, 1.05, 1.06, 1.1, 1.05, 0.8)
            scaling_y_bv <- c(1, 1, 1, 1.08, 1, 1, 1, 1.02, 1, 0.98)
            scaling_x_bv <- c(0.86, 0.88, 1.27,1.35, 0.84, 0.93, 0.92, 0.92, 0.95, 1.3)
            margin <- c(0.2,0.5,0,0)
        }else if(g == 12){
            scaling_y <- c(1, 0.97, 0.95, 1, 1, 1, 0.96, 1, 1.05, 1)
            scaling_x <- c(1.03, 1.04, 0.99, 1.06, 1.06, 1.04, 1.03,  1.05, 1.06, 1.03)
            scaling_y_bv <- c(1.01,    1,    1.01, 1.03, 1,    1.035, 1.035, 1,    1,    1.01)
            scaling_x_bv <- c(0.97, 0.94, 0.96, 0.975, 0.94, 1.02,  0.99, 0.93, 0.95, 0.97)
            margin <- c(0.1,0.5,0,0)
        }
        gg <- ggplot() + labs(x = "s", y = "u", title = paste0("Gene ", g))

        # real gene-dynamic
        gg <- plot_sVSu(t0_off = real$t0_off_real[tyT0_off, g], t0_on = 0, alpha = real$alpha_real[g,], beta = real$beta_real[g], gamma = real$gamma_real[g], pos_u = real$pos_u_real[,g], pos_s = real$pos_s_real[,g], g = g, subGrLabels = real$subtypeCell, add = TRUE, colCell = "red", xlim = NA, ylim = NA, axisTitle.size = 20, axisText.size = 10, title.size = 30, colDyn = "red", lineSize = 0.5, shapePoint = 22, sizePoint = 3.5, gg = gg)  

        #  BayVel
        gg <- plot_sVSu(t0_off = bv$T0_off_chain[tyT0_off, g, iterToPlot], t0_on = 0, alpha = bv$alpha_chain[g, , iterToPlot], beta = bv$beta_chain[g, iterToPlot], gamma = bv$gamma_chain[g, iterToPlot], pos_u = bv$u_chain[,g,iterToPlot], pos_s = bv$s_chain[,g,iterToPlot], g = g, subGrLabels = bv$subtypeCell, add = TRUE, colCell = "darkgreen", , xlim = 10, ylim = 10, axisTitle.size = 20, axisText.size = 20, title.size = 30, colDyn = "darkgreen", lineSize = 0.5, shapePoint = 21, sizePoint = 4, gg = gg)

        gg <- gg + 
            annotate(geom = "text", 
                y = unique(real$pos_u_real[,g])*scaling_y, 
                x = unique(real$pos_s_real[,g])*scaling_x, 
                label = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "L"), size = 14, parse = TRUE,  family = "serif")  + 
            annotate(geom = "text", 
            y = unique(bv$u_chain[,g, iterToPlot])*scaling_y_bv, 
            x = unique(bv$s_chain[,g, iterToPlot])*scaling_x_bv, 
            label =  c("paste(A^'~')", "paste(B^'~')", "paste(C^'~')", "paste(D^'~')", "paste(E^'~')", "paste(F^'~')","paste(G^'~')","paste(H^'~')","paste(I^'~')","paste(L^'~')"), size = 14, parse = TRUE,  family = "serif")

        gg <- gg + 
            geom_point(data = data.frame(x = bv$alpha_chain[g, 1, iterToPlot]/bv$gamma_chain[g, iterToPlot], y = bv$alpha_chain[g, 1, iterToPlot]/bv$beta_chain[g, iterToPlot], group = c("Real", "BayVel")), aes(x = x, y = y, color = group), size = 0.1) + 
            scale_color_manual(values = c("Real" = "red", "BayVel" = "darkgreen")) + 
            labs(color = "") + 
            guides(color = guide_legend(override.aes = list(size = 10, shape = c(16, 15)), nrow = 2)) + 
            theme(plot.title = element_text(family = "serif", size=50,  hjust = 0.5), axis.text = element_text(family = "serif", size = 35), axis.title = element_text(family = "serif", size = 45), legend.position = c(0.82, 0.08), legend.text = element_text(family = "serif", size = 38), legend.title = element_blank(), legend.key = element_rect(colour = "transparent", fill = NA), legend.background=element_blank(), plot.margin = unit(margin, "lines")) 

        pdf(paste0(pathOutput, "/Fig4_col2_gene", g, "_tyT0_off", tyT0, ".pdf"), width = 7*1.2, height = 7*1.5)
            plot(gg)
        dev.off() 

    }
}

# -----------------------------
# Fig. 4, columns 3: variability of the subgroups coordinates
# -----------------------------
for(g in genesToPlot){
    for(tyT0_off in 1:max(bv$typeCellT0_off)){
        # graphical parameters that are different between the two plots
        if(g == 3){
            subgroups <- c(4, 5, 7, 8)
            scaling_y <- c(0.97, 1, 1, 1)
            scaling_x <- c(0.75, 1.15, 1.06, 1.1)
            labels_real <- c("D", "E", "G", "H")
            scaling_y_bv <- c(1.08, 1, 1, 1.02) 
            scaling_x_bv <- c(1.35, 0.84,0.92,  0.92)
            labels_bv <- c("paste(D^'~')", "paste(E^'~')","paste(G^'~')", "paste(H^'~')")
            legend_y <- 5.5
            legend_x <- 2.5
        }else if(g == 12){
            subgroups <- c(1, 6, 8, 9)
            scaling_y <- c(1, 1, 1, 1.05)
            scaling_x <- c(1.03, 1.04,  1.05, 1.06)
            labels_real <- c("A", "F", "H", "I")
            scaling_y_bv <- c(1.01,  1.035, 1,    1)
            scaling_x_bv <- c(0.97, 1.02, 0.93, 0.95)
            labels_bv <- c("paste(A^'~')", "paste(F^'~')","paste(H^'~')", "paste(I^'~')")
            legend_y <- 5.2
            legend_x <- 2.8
        }

        df <- data.frame(
            s1 = bv$s_chain[subgroups[1],g,], 
            s2 = bv$s_chain[subgroups[2],g,], 
            s3 = bv$s_chain[subgroups[3],g,], 
            s4 = bv$s_chain[subgroups[4],g,], 
            u1 = bv$u_chain[subgroups[1],g,], 
            u2 = bv$u_chain[subgroups[2],g,], 
            u3 = bv$u_chain[subgroups[3],g,], 
            u4 = bv$u_chain[subgroups[4],g,]
        )

        # plot variability of the positions of 4 subgroups estimated by BayVel (the number of subgroups has been decided for graphical reasons)
        gg <- ggplot() +
        geom_point(data = df, aes(x = s1, y = u1), col = "#9966CC", alpha = 0.2, size = 1) +
        geom_point(data = df, aes(x = s2, y = u2), col = "#66CCCC", alpha = 0.2, size = 1) +
        geom_point(data = df, aes(x = s3, y = u3), col = "#FF66CC", alpha = 0.2, size = 1) +
        geom_point(data = df, aes(x = s4, y = u4), col = "#CC9900", alpha = 0.2, size = 1) + 
        labs(x = "s", y = "u")

        sty <- c(which(bv$subtypeCell == subgroups[1])[1], which(bv$subtypeCell == subgroups[2])[1], which(bv$subtypeCell == subgroups[3])[1], which(bv$subtypeCell == subgroups[4])[1])

        # real gene-dynamic
        gg <- plot_sVSu(t0_off = real$t0_off_real[tyT0_off, g], t0_on = 0, alpha = real$alpha_real[g,], beta = real$beta_real[g], gamma = real$gamma_real[g], pos_u = real$pos_u_real[sty,g], pos_s = real$pos_s_real[sty,g], g = g, subGrLabels = bv$subtypeCell[sty], add = TRUE, colCell = "red", xlim = NA, ylim = NA, axisTitle.size = 20, axisText.size = 10, title.size = 30, colDyn = "red", lineSize = 1, shapePoint = 22, sizePoint = 5, gg = gg) + labs(title = "")

        # BayVel gene-dynamic
        gg <- plot_sVSu(t0_off = bv$T0_off_chain[tyT0_off, g, iterToPlot], t0_on = 0, alpha = bv$alpha_chain[g, , iterToPlot], beta = bv$beta_chain[g, iterToPlot], gamma = bv$gamma_chain[g, iterToPlot], bv$u_chain[subgroups,g,iterToPlot], pos_s = bv$s_chain[subgroups,g,iterToPlot], g = g, subGrLabels = bv$subtypeCell[sty], add = TRUE, colCell = "darkgreen", , xlim = 10, ylim = 10, axisTitle.size = 20, axisText.size = 20, title.size = 30, colDyn = "darkgreen", lineSize = 1, shapePoint = 21, sizePoint = 5, gg = gg)

        # add labels of the estimated and real positions
        gg <- gg + 
            annotate(geom = "text", 
                y = unique(real$pos_u_real[sty,g])*scaling_y, 
                x = unique(real$pos_s_real[sty,g])*scaling_x,
                label = labels_real, size = 14, parse = TRUE,  family = "serif")  + 
            annotate(geom = "text", 
                y = unique(bv$u_chain[subtypes, g, iterToPlot])*scaling_y_bv, 
                x = unique(bv$s_chain[subtypes, g, iterToPlot])*scaling_x_bv, 
                label =  labels_bv, size = 14, parse = TRUE,  family = "serif") + 
            annotate(geom = "text", y =  legend_y, x = legend_x,  label = "paste('Posterior samples of ','(', s^'~', ', ', u^'~', ')')", size = 14, parse = TRUE,  family = "serif") 

        # graphical parameters
        gg <- gg + 
            theme(plot.title = element_text(family = "serif", size=50,  hjust = 0.5), axis.text = element_text(family = "serif", size = 35), axis.title = element_text(family = "serif", size = 45), legend.text = element_text(family = "serif", size = 38), legend.title = element_blank(), legend.key = element_rect(colour = "transparent", fill = NA), legend.background=element_blank(), plot.margin = unit(c(0.2,0.5,0,0), "lines")) 
        
        if(g == 12)
        {
            gg <- gg + scale_x_continuous(breaks= c(2, 2.5, 3, 3.5, 4))
        }

        # save the output
        pdf(paste0(pathOutput, "/Fig4_col3_gene", g, "_tyT0_off", tyT0, ".pdf"), width = 7*1.2, height = 7*1.5)
            plot(gg)
        dev.off() 
    }
}

