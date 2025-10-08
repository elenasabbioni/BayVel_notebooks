# ---------------------------------------------
# Script to produce Fig. 1a, 1b and 1c of BayVel paper
#
# OUTPUT
# - 3 PDF files, each for each panel in Figure 1
#
# INPUTS (set at the beginning or controlled in loops):
# - pathToYourDirectory: path to the working directory
#
# DEPENDECIES
# - External functions loaded from "functions.R" (e.g. u(), s(), u0(), s0())
# - packages: ggplot2
#             latex2exp
#             RColorBrewer
#             grid
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
source(paste0(pathToYourDirectory, "/functions.R"))
# Set the output path where you will save the images
pathOutput <- paste0(pathToYourDirectory, "/figuresPaper/")


# -----------------------------
#  PACKAGES 
# -----------------------------
library(ggplot2)
library(latex2exp)
library(RColorBrewer)
library(grid)

# ----------------------------------------------
# FIG. 1A
# Example of time evolution of alpha(t), u(t) and s(t)
# ----------------------------------------------

# decide values for transcription, splicing and degradation rate and for off-switching time
df <- data.frame(alphaON = 3, alphaOFF = 1, gamma = 0.75, beta = 1, t0_off = 6, t0_on = 6)

# time limit fot the plots
xlim <- 4*df$t0_off + 1.9

#################################################################
#################################################################
themeGGPLOT <- theme(
  axis.text.x = element_text(face = "bold", size = 25),
  axis.text.y = element_text(face = "bold", size = 25),
  axis.title.x = element_text(face = "bold", size = 25),
  axis.title.y = element_text(face = "bold", size = 25),
  legend.text = element_text(face = "bold", size = 25),
  legend.title = element_text(face = "bold", size = 25)
)

themeLEGEND <- theme(legend.key.size = unit(1, 'cm'), legend.title = element_text(size = 15), legend.text = element_text(size = 15))

colGGPLOT <- c("#ca0020", "#f4a582", "#f7f7f7", "#92c5de", "#0571b0")[-3]

#####################################################################


# compute the switching coordinates
df$u0 <- u0(t0_off = df$t0_off, t0_on = 0, u0_off = NA, u0_on = NA, k = 0, alpha = c(df$alphaOFF, df$alphaON), beta = df$beta)
df$s0 <- s0(t0_off = df$t0_off, t0_on = 0, u0_off = NA, u0_on = NA, s0_off = NA, s0_on = NA, k = 0, alpha = c(df$alphaOFF, df$alphaON), beta = df$beta, gamma = df$gamma)

# compute u(t)
uON <- u(t = seq(0, df$t0_off, 0.01), t0_off = df$t0_off, t0_on = 0, u0_off = df$u0, u0_on = NA, k = rep(2, length(seq(0, df$t0_off, 0.01))), alpha = c(df$alphaOFF, df$alphaON), beta = df$beta)
uOFF <- u(t = seq(df$t0_off, xlim, 0.01), t0_off = df$t0_off, t0_on = 0, u0_off = df$u0, u0_on = NA, k = rep(0, length(seq(df$t0_off, xlim, 0.01))), alpha = c(df$alphaOFF, df$alphaON), beta = df$beta)

df_uON  <- data.frame(uON = uON,   x = seq(df$t0_on, df$t0_on + df$t0_off, 0.01))
df_uOFF <- data.frame(uOFF = uOFF, x = seq(df$t0_on + df$t0_off, df$t0_on + xlim, 0.01))

# compute s(t)
sON <- s(t = seq(0, df$t0_off, 0.01), t0_off = df$t0_off, t0_on = 0, u0_off = df$u0, u0_on = NA, s0_off = df$s0, s0_on = NA, k = rep(2, length(seq(0, df$t0_off, 0.01))), alpha = c(df$alphaOFF, df$alphaON), beta = df$beta, gamma = df$gamma)
sOFF <- s(t =  seq(df$t0_off, xlim, 0.01), t0_off = df$t0_off, t0_on = 0, u0_off = df$u0, u0_on = NA, s0_off = df$s0, s0_on = NA, k = rep(0, length(seq(df$t0_off, xlim, 0.01))), alpha = c(df$alphaOFF, df$alphaON), beta = df$beta, gamma = df$gamma)

df_sON  <- data.frame(sON  = sON,  x = seq(df$t0_on, df$t0_on + df$t0_off, 0.01))
df_sOFF <- data.frame(sOFF = sOFF, x = seq(df$t0_on + df$t0_off, df$t0_on + xlim, 0.01))


# --- Plot alpha
gg_alpha <- ggplot(df) + 
  geom_segment(aes(x = 0, xend = t0_on, y = alphaOFF, yend = alphaOFF, color = "blue"), linewidth = 3) + 
  geom_segment(aes(x = t0_on, xend = t0_on + t0_off, y = alphaON, yend = alphaON, color = "red"), linewidth = 3)  + 
  geom_segment(aes(x = t0_on + t0_off, xend = xlim, y = alphaOFF, yend = alphaOFF, color = "blue"), linewidth = 3) + 
  geom_text(y = df$alphaOFF*1.3,  x = df$t0_off/4 + 1.5, label = "alpha^{off}", size = 35, parse = TRUE, family = "serif") + 
  geom_text(y = df$alphaOFF*1.3,  x = 2*df$t0_off + 5,   label = "alpha^{off}", size = 35, parse = TRUE, family = "serif") + 
  geom_text(y = df$alphaON*0.9,   x = df$t0_off + 3,     label = "alpha^{on}",  size = 35, parse = TRUE, family = "serif") + 
  geom_text(y = df$alphaON - 0.2, x = df$t0_off*3.8,     label = "alpha(t)",    size = 35, parse = TRUE, family = "serif") + 
  coord_cartesian(ylim = c(df$alphaOFF, 3.1), xlim = c(1.1, xlim - 1), clip = "off") +
  labs(x = "", y = "") + 
  theme(plot.margin = unit(c(1,1,1,1), "lines"), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
       ) #+ themeGGPLOT


# --- Plot u(t)
df_uON <- df_uON[which(df_uON$x < xlim),]
df_uOFF <- df_uOFF[which(df_uOFF$x < xlim),]

gg_u <- ggplot(df_uON, aes(x = x, y = uON)) + 
        geom_segment(x = df$t0_on, xend = df$t0_on, y = df$alphaOFF - 0.08, yend = 2*df$alphaON - 0.38, linetype = "dotted", color = "grey", linewidth = 2) + 
        geom_segment(x = df$t0_on + df$t0_off, xend =  df$t0_on + df$t0_off, y = df$alphaOFF - 0.08, yend = 2*df$alphaON - 0.38, linetype = "dotted", color = "grey", linewidth = 2) +
        geom_segment(aes(x = 0, xend =  df$t0_on, y = df$alphaOFF/df$beta, yend = df$alphaOFF/df$beta, color = "blue"), linewidth = 3) + 
        geom_path(aes(color = "red"), linewidth = 3) + 
        geom_path(data = df_uOFF, aes(x = x, y = uOFF, color = "blue"), size = 3) + 
        annotate(geom='text', y = df$alphaON - 0.2, x = df$t0_off*3.8,  label = TeX("$u(t)$", output='character'), parse = TRUE, size = 35, color = "black", family = "serif") + 
        coord_cartesian(ylim=c(df$alphaOFF,3.1), xlim = c(0+1.1, xlim-1),clip="off") + 
        labs(x = "", y = "") + 
        theme(
            plot.margin = unit(c(1,1,1,1), "lines"),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            legend.position = "none",
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()
        )


# --- Plot u(t)
df_sON  <- df_sON[ which(df_sON$x  < xlim),]
df_sOFF <- df_sOFF[which(df_sOFF$x < xlim),]

scaleFUN <- function(x) sprintf("%.1f", x)


gg_s <- ggplot(df_sON, aes(x = x, y = sON)) + 
        geom_segment(x = df$t0_on, xend = df$t0_on, y = df$alphaOFF/df$gamma - 0.15, yend = df$alphaON/df$gamma + 1, linetype = "dotted", color = "grey", linewidth = 2) + 
        geom_segment(x = df$t0_on + df$t0_off, xend =  df$t0_on + df$t0_off, y = df$alphaOFF/df$gamma - 0.15, yend = df$alphaON/df$gamma + 1, linetype = "dotted", color = "grey", linewidth = 2) + 
        geom_segment(aes(x = 0, xend = df$t0_on, y = df$alphaOFF/df$gamma, yend = df$alphaOFF/df$gamma, colour = "blue", linetype = "b-g"), size = 3) + 
        geom_path(aes(colour = "red"), linewidth = 3) + 
        geom_path(data = df_sOFF, aes(x = x, y = sOFF, colour = "blue"), linewidth = 3) +
        annotate(geom='text', y = df$alphaON/df$gamma - 0.47, x = df$t0_off*3.8 ,  label = TeX("$s(t)$", output='character'), parse = TRUE, size = 35, color = "black", family = "serif") +
        annotation_custom(textGrob(TeX("$t_{0}^{on}$"), gp = gpar(fontsize=90, fontfamily = "serif")), xmin = df$t0_on - 0.5, xmax = df$t0_on - 0.4, ymin = df$alphaOFF/df$gamma - 0.7, ymax = df$alphaOFF/df$gamma - 0.7) + 
        annotation_custom(textGrob(TeX("$t_{0}^{on} + omega$"), gp=gpar(fontsize=90, fontfamily = "serif")), xmin = df$t0_on + df$t0_off + 0.6, xmax = df$t0_on + df$t0_off + 0.5, ymin = df$alphaOFF/df$gamma - 0.7, ymax = df$alphaOFF/df$gamma - 0.7) + 
        geom_segment(aes(x = df$t0_on, xend = df$t0_on, y = df$alphaOFF/df$gamma - 0.15, yend = df$alphaOFF/df$gamma - 0.25), color = "black", linetype = "solid") + 
        geom_segment(aes(x = df$t0_on + df$t0_off, xend = df$t0_on + df$t0_off, y = df$alphaOFF/df$gamma - 0.15, yend = df$alphaOFF/df$gamma - 0.25), color = "black", linetype = "solid") + 
        coord_cartesian(ylim=c(df$alphaOFF/df$gamma,4), xlim = c(0+1.1, xlim-1),clip="off") + 
        labs(x = "", y = "") + 
        theme(
            plot.margin = unit(c(1,1,7.5,1), "lines"),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            legend.position = "none",
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()
        )

pdf(paste0(pathOutput, "./Fig_1a.pdf"), width = 7*3, height = 7*3 )
    grid.arrange(gg_alpha, gg_u, gg_s)
dev.off()


# ----------------------------------------------
# FIG. 1B 
# Example of gene dynamic with a single switch and with "inductive" and "repression" phases highlighted
# ----------------------------------------------

# decide values for transcription, splicing and degradation rate and for off-switching time
df <- data.frame(alphaON = 3, alphaOFF = 1, gamma = 0.75, beta = 1, t0_off = 2)

# compute the switching coordinates
df$u0 <- u0(t0_off = df$t0_off, t0_on = 0, u0_off = NA, u0_on = NA, k = 0, alpha = c(df$alphaOFF, df$alphaON), beta = df$beta)
df$s0 <- s0(t0_off = df$t0_off, t0_on = 0, u0_off = NA, u0_on = NA, s0_off = NA, s0_on = NA, k = 0, alpha = c(df$alphaOFF, df$alphaON), beta = df$beta, gamma = df$gamma)

# plot the dynamic
pl <- plot_GeneDynamic_withNotes(t0_off = df$t0_off, u0_off = df$u0, s0_off = df$s0, t0_on = 0, alpha = c(df$alphaOFF, df$alphaON), beta = df$beta, gamma = df$gamma, r = NA, g = NA) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(plot.margin = unit(c(1,1,6.5,1), "lines"))

pdf(paste0(pathOutput, "./Fig_1b.pdf"), width = 7*3, height = 7*3 )
    print(pl)
dev.off()


# ----------------------------------------------
# FIG. 1C
# Example of gene dynamic with diffferent switching points and subgroup's positions
# ----------------------------------------------

# decide values for transcription, splicing and degradation rate and for the first off-switching time
df <- data.frame(alphaON = 3, alphaOFF = 1, gamma = 0.75, beta = 1, t0_off = 2)
# compute the associated switching coordinates
df$u0 <- u0(t0_off = df$t0_off, t0_on = 0, u0_off = NA, u0_on = NA, k = 0, alpha = c(df$alphaOFF, df$alphaON), beta = df$beta)
df$s0 <- s0(t0_off = df$t0_off, t0_on = 0, u0_off = NA, u0_on = NA, s0_off = NA, s0_on = NA, k = 0, alpha = c(df$alphaOFF, df$alphaON), beta = df$beta, gamma = df$gamma)

# add another switching time and compute the associated switching points
t0_off_v = c(df$t0_off, 1.05)
u0_off_v = rep(NA, length(t0_off_v))
s0_off_v = rep(NA, length(t0_off_v))
for(i in 1:length(t0_off_v)){
  u0_off_v[i] = u0(t0_off = t0_off_v[i], t0_on = 0, u0_off = NA, u0_on = NA, k = 0, alpha = c(df$alphaOFF, df$alphaON), beta = df$beta)
  s0_off_v[i] = s0(t0_off = t0_off_v[i], t0_on = 0, u0_off = NA, u0_on = NA, s0_off = NA, s0_on = NA, k = 0, alpha = c(df$alphaOFF, df$alphaON), beta = df$beta, gamma = df$gamma)
}

# graphic parameters
palette = brewer.pal(length(t0_off_v), "Set1")
sizeLines = c(3,3)

# add the position of te different subgroups
t_v <- list()
u_v <- list()
s_v <- list()
k_v <- list()

t_v[[1]] <- c(0.35, 1.5, 2.5)    # subtype1
t_v[[2]] <-  c(0.1, 1.65, 2.5)   # subtype2
for(i in 1:length(t0_off_v)){
  k_v[[i]] <- ifelse(t_v[[i]] < t0_off_v[i], 2, 0)
  u_v[[i]] <- u(t_v[[i]], t0_off = t0_off_v[i], t0_on = 0, u0_off_v[i], u0_on = NA, k_v[[i]], alpha = c(df$alphaOFF, df$alphaON), beta = df$beta)
  s_v[[i]] <- s(t_v[[i]], t0_off = t0_off_v[i], t0_on = 0, u0_off_v[i], u0_on = NA, s0_off_v[i], s0_on = NA, k_v[[i]],  alpha = c(df$alphaOFF, df$alphaON), beta = df$beta, gamma = df$gamma)
}
df_points <- data.frame(u = unlist(u_v), s = unlist(s_v),fill = brewer.pal(length(t0_off_v)*length(t_v[[1]]), "Set2"), color = rep(palette[1:length(t0_off_v)], each = length(t_v[[1]])))
pl <- pl + geom_point(data = df_points, aes(x = s, y = u), color = df_points$color, fill = df_points$fill, pch = 16, size = 10, stroke = 3) 

# position of the steady states
u_SS_off <- df$alphaOFF/df$beta
u_SS_on  <- df$alphaON/df$beta
s_SS_off <- df$alphaOFF/df$gamma        
s_SS_on  <- df$alphaON/df$gamma        

# compute gene dynamic on the on branch (from lower to upper steady state)
u_plotON <- seq(u_SS_off, u_SS_on, 0.01)
s_plotON <- s_di_u(u_plotON, u_SS_on, s_SS_on, u_SS_off, s_SS_off, k = 2, c(df$alphaOFF, df$alphaON), df$beta, df$gamma)
# compute gene dynamic off the on branch (from upper to lower steady state)
u_plotOFF <- seq(u_SS_off, u_SS_on, 0.01)
s_plotOFF <- s_di_u(u_plotOFF, u_SS_on, s_SS_on, u_SS_off, s_SS_off, k = 0, c(df$alphaOFF, df$alphaON), df$beta, df$gamma)

df_switching <- data.frame(uON_SS = as.vector(u_plotON), sON_SS = as.vector(s_plotON), uOFF_SS = as.vector(u_plotOFF), sOFF_SS = as.vector(s_plotOFF))

pl <- ggplot(df_switching)  +
    geom_path(aes(x = sON_SS, y = uON_SS), color = "grey", linetype = "dotted", linewidth = 3) + 
    geom_path(aes(x = sOFF_SS, y = uOFF_SS), color = "grey", linetype = "dotted", linewidth = 3)

# add the switching points
for(i in 1:length(t0_off_v)){
    if(t0_off_v[i] != Inf){ # the dynamic switches before arriving to the upper steady state
        # compute gene dynamic on the on branch (from lower steady state to off-switching point)
        u_plotON_switch <- seq(u_SS_off, u0_off_v[i], 0.001)
        s_plotON_switch <- s_di_u(u_plotON_switch, u0_off_v[i], s0_off_v[i], u_SS_off, s_SS_off, k = 2, c(df$alphaOFF, df$alphaON), df$beta, df$gamma)
        dfSwitchON <- data.frame(uON_switch = u_plotON_switch, sON_switch = s_plotON_switch)

        # compute gene dynamic on the off branch (from the off-switching point to lower steady state)
        u_plotOFF_switch <- seq(u_SS_off, u0_off_v[i], 0.001)
        s_plotOFF_switch <- s_di_u(u_plotOFF_switch, u0_off_v[i], s0_off_v[i], u_SS_off, s_SS_off, k = 0, c(df$alphaOFF, df$alphaON), df$beta, df$gamma)
        dfSwitchOFF <- data.frame(uOFF_switch = u_plotOFF_switch, sOFF_switch = s_plotOFF_switch)
        maxSoff <- which.max(s_plotOFF_switch)

        pl <- pl +
        geom_path(data = dfSwitchON, aes(x = sON_switch, y = uON_switch), color = palette[i], size = sizeLines[i]) +     
        geom_path(data = dfSwitchOFF, aes(x = sOFF_switch, y = uOFF_switch), color = palette[i], size = sizeLines[i])   

    }else{
        pl <- pl +
            geom_path(aes(x = sON_SS, y = uON_SS, color = "red"), linewidth = sizeLines[i]) +
            geom_path(aes(x = sOFF_SS, y = uOFF_SS, color = "blue"), linewidth = sizeLines[i])
    }
}

# add labels for switching points
l1 <- TeX(r"($(s_{1g}^{omega}, u_{1g}^{omega})$)", output = "character")
l2 <- TeX(r"($(s_{2g}^{omega}, u_{2g}^{omega})$)", output = "character")
labels <- c(l1, l2)
for(i in 1:length(t0_off_v)){
    pl <- pl + annotate(geom = "text", y = u0_off_v[i] + 0.1, x = s0_off_v[i] - 0.4, label = labels[i], size = 35, parse = TRUE, family = "serif") 
}
# add labels for subgroups'positions
l11 <- "paste('(', s['11g']^'~', ', ', u['11g']^'~', ')')"
l12 <- "paste('(', s['12g']^'~', ', ', u['12g']^'~', ')')"
l13 <- "paste('(', s['13g']^'~', ', ', u['13g']^'~', ')')"
l21 <- "paste('(', s['24g']^'~', ', ', u['24g']^'~', ')')"
l22 <- "paste('(', s['25g']^'~', ', ', u['25g']^'~', ')')"
l23 <- "paste('(', s['26g']^'~', ', ', u['26g']^'~', ')')"
labels_point <- c(l11, l12, l13, l21, l22, l23)
for(i in 1:length(unlist(u_v))){
    pl <- pl + annotate(geom = "text", y = unlist(u_v)[i] + 0.09, x = unlist(s_v)[i] - 0.2, label = labels_point[i], size = 27, parse = TRUE, family = "serif") 
}


# add graphical parameters
pl <- pl + 
    coord_cartesian(xlim =  c(min(rbind(df_switching$sON_SS, df_switching$sOFF_SS)) - 0.5, max(rbind(df_switching$sON_SS, df_switching$sOFF_SS))), ylim = c(min(rbind(df_switching$uON_SS, df_switching$uOFF_SS)), max(rbind(df_switching$uON_SS, df_switching$uOFF_SS)))) +
    labs(x = "", y = "") + 
    theme(legend.position = "none", plot.title = element_text(family = "serif", size=70,  hjust = 0.5), 
          axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), 
          axis.text.x = element_blank(), plot.margin = unit(c(1,1,6.5,1), "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()
          ) 

pdf(paste0(pathOutput, "./Fig_1c.pdf"), width = 7*3, height = 7*3 )
    print(pl)
dev.off()















