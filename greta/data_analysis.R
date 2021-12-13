library(ggplot2)
library(patchwork)
library(RColorBrewer)
#library(hrbrthemes)
library(viridis)

rm(list=ls()) #reset working space
graphics.off() #closing current or all graphical windows
setwd('C:/Users/Gretka/Desktop/combined')

data = read.table(file = 'data_whole_combined.tsv', sep = '\t', header = TRUE, comment.char = '@')
tempura = read.table(file = 'matched_temp_multi.tsv', sep = '\t', header = TRUE, comment.char = '@')

data = data[data$temp > -21, ] 
data = data[data$temp < 110, ] 
data = data[!grepl('Eukaryota',data$taxonomy),]

data=data[data$OTU_ID %in% names(which(table(data$OTU_ID) > 10)), ]
#n <- 60
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))

otu_temp = intersect(unique(data$OTU_ID),unique(tempura$OTU_id_5))
col_vector = brewer.pal(n = 8, name = "Dark2")
for (nr in 1:20){
  otu_sample = sample(otu_temp, 5, replace=FALSE)
  counter = 1
  for (otu in otu_sample){
    plot_name = paste("plot_", counter, sep = "")
    plot = ggplot(data[data$OTU_ID == otu, ], aes(x = "", y = temp), ylim = c(tempura[tempura$OTU_id_5 == otu,]$Tmin -2, tempura[tempura$OTU_id_5 == otu,]$Tmax + 2)) + 
      geom_boxplot() +
      geom_jitter(width = 0.2,colour = sample(col_vector, 1), size = 3) +
      ggtitle(otu) +
      geom_hline(yintercept = tempura[tempura$OTU_id_5 == otu,]$Tmin, linetype="dashed", size = 1, colour="blue") +
      geom_hline(yintercept = tempura[tempura$OTU_id_5 == otu,]$Topt_ave, linetype="dashed", size = 1, colour="black") +
      geom_hline(yintercept = tempura[tempura$OTU_id_5 == otu,]$Tmax, linetype="dashed", size = 1, colour="red") +
      theme_bw() +
      theme(axis.title.y=element_blank()) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    assign(plot_name, plot)
    counter = counter + 1
  }
  png(file=paste("C:/Users/Gretka/Desktop/combined/data_distrb", nr,".png", sep = ""),
      width=1100, height=1100 * 0.5833333)
  print(plot_1 + plot_2 + plot_3 + plot_4 + plot_5 + plot_layout(ncol = 5))
  dev.off()
}

################################################################################
min_opt = read.table(file = 'min_opt_z.tsv', sep = '\t', header = TRUE, comment.char = '@')
max_opt = read.table(file = 'max_opt_z.tsv', sep = '\t', header = TRUE, comment.char = '@')
opt_opt = read.table(file = 'opt_opt_z.tsv', sep = '\t', header = TRUE, comment.char = '@')

min_opt = read.table(file = 'big_regression_table/min_opt_z.tsv', sep = '\t', header = TRUE, comment.char = '@')
max_opt = read.table(file = 'big_regression_table/max_opt_z.tsv', sep = '\t', header = TRUE, comment.char = '@')
opt_opt = read.table(file = 'big_regression_table/opt_opt_z.tsv', sep = '\t', header = TRUE, comment.char = '@')

# https://www.r-graph-gallery.com/79-levelplot-with-ggplot2.html

# Heatmap 
ggplot(opt_opt, aes(min_reads, min_sampl_nr, fill= r_sq)) + 
  geom_tile() +
  ggtitle("Topt") +
  scale_fill_viridis(discrete=FALSE)

ggplot(min_opt, aes(min_reads, min_sampl_nr, fill= r_sq)) + 
  geom_tile() +
  ggtitle("Tmin") +
  scale_fill_viridis(discrete=FALSE)

ggplot(max_opt, aes(min_reads, min_sampl_nr, fill= r_sq)) + 
  geom_tile() +
  ggtitle("Tmax") +
  scale_fill_viridis(discrete=FALSE)
##############################################################################
ggplot(max_opt, aes(min_reads, min_sampl_nr, fill= TEMPURA_matches)) + 
  geom_tile() +
  ggtitle("TEMPURA matches") +
  scale_fill_viridis(discrete=FALSE)

new <- opt_opt
new$r_sq[new$TEMPURA_matches < 40] <- 0
ggplot(new, aes(min_reads, min_sampl_nr, fill= r_sq)) + 
  geom_tile() +
  ggtitle("Topt") +
  scale_fill_viridis(discrete=FALSE)
  scale_fill_viridis(discrete=FALSE)

new <- max_opt
new$r_sq[new$TEMPURA_matches < 40] <- 0
ggplot(new, aes(min_reads, min_sampl_nr, fill= r_sq)) + 
  geom_tile() +
  ggtitle("Tmax") +
  scale_fill_viridis(discrete=FALSE)
  scale_fill_viridis(discrete=FALSE)
  
new <- min_opt
new$r_sq[new$TEMPURA_matches < 40] <- 0
ggplot(new, aes(min_reads, min_sampl_nr, fill= r_sq)) + 
  geom_tile() +
  ggtitle("Tmin") +
  scale_fill_viridis(discrete=FALSE)
  scale_fill_viridis(discrete=FALSE)

#https://stackoverflow.com/questions/21119095/overlaying-an-image-plot-on-a-perspective-plot-in-r
################################################################################
library(plotly)
library(DescTools)
opt_data = as.matrix(xtabs(r_sq ~ min_reads + min_sampl_nr, data = opt_opt))
plot_ly(z = ~opt_data, type = "surface") %>% 
  layout(scene = list(xaxis=list(title ="min_sampl_nr"),yaxis=list(title ="min_reads"),zaxis=list(title ="R²")), title = "Topt")

min_data = as.matrix(xtabs(r_sq ~ min_reads + min_sampl_nr, data = min_opt))
plot_ly(z = ~min_data, type = "surface") %>% 
  layout(scene = list(xaxis=list(title ="min_sampl_nr"),yaxis=list(title ="min_reads"),zaxis=list(title ="R²")), title = "Tmin")

max_data = as.matrix(xtabs(r_sq ~ min_reads + min_sampl_nr, data = max_opt))
plot_ly(z = ~max_data, type = "surface") %>% 
  layout(scene = list(xaxis=list(title ="min_sampl_nr"),yaxis=list(title ="min_reads"),zaxis=list(title ="R²")), title = "Tmax")

################################################################################

opt_opt_match_filt=opt_opt[opt_opt$TEMPURA_matches>390,]
min_opt_match_filt=min_opt[min_opt$TEMPURA_matches>60,]
max_opt_match_filt=max_opt[max_opt$TEMPURA_matches>60,]


df <- data.frame(x1 = numeric(), x2 = numeric())
for (i in 1:487){
  max_r2 <- c(i, max(min_opt[opt_opt$TEMPURA_matches>i,]$r_sq)) 
  df[i, ] <- max_r2
}
df <- df[is.finite(rowSums(df)),]
plot(df$x1, df$x2, type = "S",
     xlab = "minimum matches", ylab = "R²", main = "min")

##############
df <- data.frame(x1 = numeric(), x2 = numeric())
for (i in 2:488){
  max_r2 <- c(i, max(opt_opt[opt_opt$TEMPURA_matches  == i,]$r_sq))
  df[i, ] <- max_r2
}
df <- df[is.finite(rowSums(df)),]
plot(df$x1, df$x2, type = "S", xaxp = c(0, 500, 10),
     xlab = "matches", ylab = "maximum R²", main = "opt")
######
df <- data.frame(x1 = numeric(), x2 = numeric())
for (i in 2:488){
  max_r2 <- c(i, max(max_opt[max_opt$TEMPURA_matches  == i,]$r_sq))
  df[i, ] <- max_r2
}
df <- df[is.finite(rowSums(df)),]
plot(df$x1, df$x2, type = "S", xaxp = c(0, 500, 10),
     xlab = "matches", ylab = "maximum R²", main = "max")
######
df <- data.frame(x1 = numeric(), x2 = numeric())
for (i in 2:488){
  max_r2 <- c(i, max(min_opt[min_opt$TEMPURA_matches  == i,]$r_sq))
  df[i, ] <- max_r2
}
df <- df[is.finite(rowSums(df)),]
plot(df$x1, df$x2, type = "S", xaxp = c(0, 500, 10),
     xlab = "matches", ylab = "maximum R²", main = "min")

max_opt[opt_opt$TEMPURA_matches==91,]
opt_opt[opt_opt$TEMPURA_matches==110,][which.max(opt_opt[opt_opt$TEMPURA_matches==110,]$r_sq),]
min_opt[min_opt$TEMPURA_matches==110,][which.max(min_opt[min_opt$TEMPURA_matches==110,]$r_sq),]
max_opt[max_opt$TEMPURA_matches==110,][which.max(max_opt[max_opt$TEMPURA_matches==110,]$r_sq),]


## Create a simple surface  f(x,y) = -x^2 - y^2
## Colour the surface according to x^2 only

x = max_opt$min_reads
y = max_opt$min_sampl_nr
z = max_opt$r_sq
## Fourth dim
z_col = max_opt$TEMPURA_matches
nx = 480; ny = 480
## Average the values at the corner of each facet
## and scale to a value in [0, 1].  We will use this
## to select a gray for colouring the facet. 
hgt = 0.25 * (z_col[-nx,-ny] + z_col[-1,-ny] + z_col[-nx,-1] + z_col[-1,-1])
hgt = (hgt - min(hgt))/ (max(hgt) - min(hgt))

##  Plot the surface with the specified facet colours.
persp(x, y, z, col = gray(1 - hgt))
persp(x, y, z, col=cm.colors(32)[floor(31*hgt+1)], theta=-35, phi=10)


