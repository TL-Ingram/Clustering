library <- c("tidyverse", "data.table", "hrbrthemes", "GGally", "factoextra", "mlr", "clue", "cowplot", "ggrepel", "ggforce", "here")
lapply(library, require, character.only = TRUE)

############################################################
### Data cleaning /wrangling
############################################################
data <- do.call(rbind, 
                lapply(list.files(path = here("raw_data"), full.names = TRUE), read.csv)) %>%
      select(17, 7, 23, 26, 28, 31, 33, 36) %>%
      separate(Well, c("Y", "X"), 1, remove = F)

data$Dilution[data$Y <= "P" & data$X == 1] <- "120"
data$Dilution[data$Y <= "P" & data$X == 3] <- "60"
data$Dilution[data$Y <= "P" & data$X == 5] <- "30"
data$Dilution[data$Y <= "P" & data$X == 7] <- "15"
data$Dilution[data$Y <= "P" & data$X == 9] <- "7.5"
data$Dilution[data$Y <= "P" & data$X == 11] <- "3.5"
data$Dilution[data$Y <= "P" & data$X == 13] <- "1.8"
data$Dilution[data$Y <= "P" & data$X == 15] <- "0.9"
data$Dilution[data$Y <= "P" & data$X == 17] <- "0.45"
data$Dilution[data$Y <= "P" & data$X == 19] <- "0.225"
data$Dilution[data$Y <= "P" & data$X == 21] <- "0.113"

gathered_data <- data %>%
     filter(!is.na(Dilution)) %>%
     filter(!(Y == "B" | Y == "D" | Y == "F" | Y == "H" | Y == "J" | Y == "L" | Y == "N" | Y == "P")) %>%
     mutate(Dilution = factor(Dilution)) %>%
     gather(key = "Gene", value = "Fluorescence", 6,8,10) %>%
     select(4,8,9,10) %>%
     rename("lot" = 1) %>%
     slice(-(c(1,1057,2113,1834,778,2890,128,2285,2219,1382,326,3075,2022,2000)))

gathered_data$Gene[gathered_data$Gene %like% "Channel_1_EndFluoBC"] <- "ORF1ab"
gathered_data$Gene[gathered_data$Gene %like% "Channel_2_EndFluoBC"] <- "N"
gathered_data$Gene[gathered_data$Gene %like% "Channel_3_EndFluoBC"] <- "S"
gathered_data$lot[gathered_data$lot %like% "QC20211011-QNostics-Panel.eds"] <- "8185"
gathered_data$lot[gathered_data$lot %like% "QC20211014-Qnostics-Panel-6171.eds"] <- "6171"
gathered_data$lot[gathered_data$lot %like% "QC20211020-Qnostics-Panel-5158.eds"] <- "5158"
gathered_data$lot[gathered_data$lot %like% "QC20211104-Qnostics-Panel-1145.eds"] <- "1145"
gathered_data$lot[gathered_data$lot %like% "QC20211104-Qnostics-Panel-1142.eds"] <- "1142"
gathered_data$lot[gathered_data$lot %like% "QC20211104-Qnostics-Panel-4155.eds"] <- "4155"
gathered_data$lot[gathered_data$lot %like% "QC20211104-Qnostics-Panel-6166.eds"] <- "6166"
gathered_data$lot[gathered_data$lot %like% "QC20211014-Qnostics-Panel-7177.eds"] <- "7177"
gathered_data$lot[gathered_data$lot %like% "QC20211014-Qnostics-Panel-7179.eds"] <- "7179"
gathered_data$lot[gathered_data$lot %like% "QC20211020-Qnostics-Panel-6162.eds"] <- "6162"
gathered_data$lot[gathered_data$lot %like% "QC20211124-Qnostics-Panel-8189.eds"] <- "8189"
gathered_data$lot[gathered_data$lot %like% "QC20211124-Qnostics-Panel-8184.eds"] <- "8184"

# boxplots for outliers
# slice is removed outlier rows
box_outlier_graph <- gathered_data %>%
     slice(-(c(1,1057,2113,1834,778,2890,128,2285,2219,1382,326,3075,2022,2000))) %>%
     ggplot(aes(fct_reorder(Dilution, Fluorescence), Fluorescence, col = Gene)) +
     geom_boxplot(fill = "grey95", alpha = 0.7, size = 0.6) +
     facet_wrap(.~lot) +
     theme_ipsum(axis_title_just = "cc",
                 axis_title_face = "bold", 
                 axis_text_size = 14, 
                 axis_title_size = 18) +
     theme(panel.grid.major = element_line(colour = "grey90"),
           panel.grid.minor = element_line(colour = "grey90"),
           legend.title = element_text(size = 16, face = "bold"),
           legend.text = element_text(size = 16)) +
     guides(x = guide_axis(angle = 45)) +
     labs(title = "",
          x = "Concentration (copies/reaction)",
          y = "Fluorescence",
          colour = "Cluster")
box_outlier_graph
ggsave(here("results/graphs", "boxplot_outlier_removal.tiff"), dpi = 300, width = 16, height = 10)

### PCA
PCA_ready_data <- gathered_data %>%
      slice(-(c(1,1057,2113,1834,778,2890,128,2285,2219,1382,326,3075,2022,2000))) %>%
     unite(Dilution_Gene, Dilution, Gene) %>%
     lapply(function(x) ifelse(x<0, 0, x)) %>%
     as_tibble(.) %>%
     group_by(Dilution_Gene, lot) %>%
     summarise(mean = mean(Fluorescence)) %>%
     filter(!grepl("MS2", Dilution_Gene)) %>%
     ungroup(.) %>%
     spread(Dilution_Gene, mean)
summary(PCA_ready_data)

pca <- prcomp(PCA_ready_data[,2:34], center = TRUE, scale = TRUE)
summary(pca)
dimnames(pca$rotation)[[1]]

pca_biplot <- fviz_pca_biplot(pca, col.var = "contrib", label = "var",
             gradient.cols = c("white", "blue", "red"), repel = T, labelsize = 4, pointsize = 4, geom = c("point"), alpha.var = 0.2) +
     theme_ipsum(axis_title_just = "cc",
                 axis_title_face = "bold",
                 axis_text_size = 16,
                 axis_title_size = 18) +
     theme(panel.grid.major = element_line(colour = "grey90"),
           panel.grid.minor = element_line(colour = "grey90")) +
     labs(title = "",
          col = "Contribution",
            x = "PC1 (80.0%)",
            y = "PC2 (10.2%)")
pca_biplot

eigen_scree <- fviz_screeplot(pca, addlabels = TRUE, choice = "eigenvalue") +
     theme_ipsum(axis_title_just = "cc",
                 axis_title_face = "bold", 
                 axis_text_size = 16, 
                 axis_title_size = 18) +
     theme(panel.grid.major = element_line(colour = "grey90"),
           panel.grid.minor = element_line(colour = "grey90")) +
     labs(title = "",
          x = "Principal component",
          y = "Eigenvalue")
variance_scree <- fviz_screeplot(pca, addlabels = TRUE, choice = "variance") +
     theme_ipsum(axis_title_just = "cc",
                 axis_title_face = "bold", 
                 axis_text_size = 16, 
                 axis_title_size = 18) +
     theme(panel.grid.major = element_line(colour = "grey90"),
           panel.grid.minor = element_line(colour = "grey90")) +
     labs(title = "",
          x = "Principal component",
          y = "% of variance")
PCA_complete <- PCA_ready_data %>%
     mutate(PCA1 = pca$x[, 1], PCA2 = pca$x[, 2]) %>%
     select(1,35,36)

# pca_dots <- ggplot(PCA_complete, aes(PCA1, PCA2, col = lot)) +
#      geom_point(size = 3) +
#      theme_ipsum(axis_title_just = "cc",
#                  axis_title_face = "bold", 
#                  axis_text_size = 16, 
#                  axis_title_size = 18) +
#      theme(panel.grid.major = element_line(colour = "grey90"),
#            panel.grid.minor = element_line(colour = "grey90")) +
#      labs(title = "",
#           col = "Lot",
#           x = "PC1 (eigenvalue)",
#           y = "PC2 (eigenvalue)")
scree_biplot <- plot_grid(variance_scree, pca_biplot, ncol = 2, labels = c("A","B"))
scree_biplot
save_plot(here("results/graphs", "scree_pca_plots.tiff") , plot = scree_biplot, dpi = 300, base_width = 20, base_height = 10)


##PCA DONE FOR NOW


### clustering
cluster_task <- PCA_complete %>%
     select(2,3) %>%
     as_tibble()
cluster_task_two <- PCA_complete
cluster_task <- makeClusterTask(data = as.data.frame(cluster_task))


kMeans <- makeLearner("cluster.kmeans",
                      par.vals = list(iter.max = 500, nstart = 25))
########
kMeansParamSpace <- makeParamSet(
     makeDiscreteParam("centers", values = c(3)),
     makeDiscreteParam("algorithm",
                       values = c("Hartigan-Wong", "Lloyd", "MacQueen", "Forgy")))
gridSearch <- makeTuneControlGrid()
kFold <- makeResampleDesc("CV", iters = 3)
########
tunedK <- tuneParams(kMeans, task = cluster_task,
                     resampling = kFold,
                     par.set = kMeansParamSpace,
                     control = gridSearch,
                     measures = list(db, G1, G2))
tunedK
#######
kMeansTuningData <- generateHyperParsEffectData(tunedK)
kMeansTuningData$data

gatheredTuningData <- gather(kMeansTuningData$data,
                             key = "Metric",
                             value = "Value",
                             c(-centers, -iteration, -algorithm))

ggplot(gatheredTuningData, aes(centers, Value, col = algorithm)) +
     facet_wrap(~ Metric, scales = "free_y") +
     geom_line() +
     geom_point() +
     theme_bw()
#######
tunedKMeans <- setHyperPars(kMeans, par.vals = tunedK$x)
tunedKMeansModel <- train(tunedKMeans, cluster_task)

kMeansModelData <- getLearnerModel(tunedKMeansModel)
tunedKMeansModel$learner.model
kMeansModelData$cluster
kMeansModelData$centers


#######
cluster_task_two <- cluster_task_two %>%
      mutate(kMeansCluster = as.factor(kMeansModelData$cluster)) #%>%
      # mutate(kMeansCluster = ifelse(kMeansCluster == 1, 2,
      #                               ifelse(kMeansCluster == 3, 1, 3)))
write.csv(cluster_task_two, here("results/csv", "lots_corresponding_cluster.csv"), row.names = F)

pairs_plot <- ggpairs(cluster_task_two, aes(col = kMeansCluster),
        upper = list(continuous = "density")) +
     theme_bw()
pairs_plot
ggsave("pairs_plot_three_clusters.tiff", plot = pairs_plot, dpi = 300, width = 16, height = 12)
cluster_task_two %>%
     ggplot(aes(PCA1, PCA2, colour = kMeansCluster)) +
     geom_point(size = 4, alpha = 0.7) +
     geom_mark_ellipse(size = 1) +
     xlim(-8,10) +
     theme_ipsum(axis_title_just = "cc",
                 axis_title_face = "bold", 
                 axis_text_size = 16, 
                 axis_title_size = 18) +
     theme(panel.grid.major = element_line(colour = "grey90"),
           panel.grid.minor = element_line(colour = "grey90"),
           legend.title = element_text(size = 16, face = "bold"),
           legend.text = element_text(size = 16))+
     labs(title = "",
          x = "PC1 (eigenvalue)",
          y = "PC2 (eigenvalue)",
          colour = "Cluster")
ggsave(here("results/graphs", "Kmeans_three_cluster.tiff"), dpi = 300, height = 10, width = 10)




# looks like three clusters. One has all three targets doing badly; one has just S doing bad; one has all three targets doing well

# show sig differences between PCA 1 clusters and between PCA 2 clusters
# Then take all data for each of these clusters. Group the values together and do statistics to show differences between S.
# Say something like "this method could be extended to predict new Kits quality upon entering the laboratory. 

## take lots by group


