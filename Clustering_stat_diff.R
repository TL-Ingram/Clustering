library <- c("tidyverse", "data.table", "hrbrthemes", "GGally", "factoextra", "mlr", "clue", "cowplot", "ggrepel", "ggforce", "here", "lubridate")
lapply(library, require, character.only = TRUE)

next_step <- gathered_data %>%
     dplyr::inner_join(cluster_task_two, by="lot") %>%
     select(-PCA1, -PCA2)
clustered_boxplot <- next_step %>%
     group_by(kMeansCluster, Dilution) %>%
     ggplot(aes(fct_reorder(Dilution, Fluorescence), Fluorescence, colour = kMeansCluster)) +
     geom_boxplot(fill = "grey95", alpha = 0.8, size = 0.6) +
     facet_wrap(Gene~.) +
     theme_ipsum(axis_title_just = "cc",
                 axis_title_face = "bold", 
                 axis_text_size = 16, 
                 axis_title_size = 18) +
     theme(panel.grid.major = element_line(colour = "grey90"),
           panel.grid.minor = element_line(colour = "grey90"),
           strip.text = element_text(size = 18, face = "bold"),
           legend.title = element_text(size = 16, face = "bold"),
           legend.text = element_text(size = 16)) +
     guides(x = guide_axis(angle = 45)) +
     labs(x = "Concentration (copies/reaction)",
          colour = "Cluster")
ggsave(here("results/graphs", "clusters_boxplot.tiff"), dpi = 300, height = 12, width = 14)
clustered_boxplot

next_steps <- next_step %>%
     mutate(Id= row_number()) %>%
     spread(Dilution, Fluorescence) %>%
     select(-1,-4) %>%
     group_by(Gene) %>%
     # mutate_if(is.numeric, scale1) %>%
     mutate_at(vars(contains(c("0", "1", "5"))), scale) %>%
     gather(dilution, value, -Gene, -kMeansCluster)
     # group_by(Gene, kMeansCluster) %>%

ggplot(next_steps, aes(Gene, value, colour = kMeansCluster)) +
      geom_boxplot(fill = "grey95", alpha = 0.8, size = 0.6) +
      theme_ipsum(axis_title_just = "cc",
                  axis_title_face = "bold", 
                  axis_text_size = 16, 
                  axis_title_size = 18) +
      theme(panel.grid.major = element_line(colour = "grey90"),
            panel.grid.minor = element_line(colour = "grey90"),
            strip.text = element_text(size = 18, face = "bold"),
            legend.title = element_text(size = 16, face = "bold"),
            legend.text = element_text(size = 16)) +
      labs(x = "Gene",
           y = "Fluorescence (scaled)",
           colour = "Cluster")
ggsave(here("results/graphs", "_overall_differences_scaled_genes.tiff"), dpi = 300, height = 10, width = 10)


#### Next thing to do is work out statistical differences of this using one-way ANOVA and do it three times for each gene.
#### Also make above graph look better

stats_comparison <- next_steps %>%
     select(1,2,4) %>%
     drop_na()
# normality test qq plot
# homoscedasity test
# if both okay then do one way ANOVA for each gene individually

ggplot(stats_comparison, aes(sample = value, colour = kMeansCluster)) +
      facet_grid(Gene ~ .) +
      stat_qq(alpha = 0.5) +
      stat_qq_line() +
      theme_ipsum(axis_title_just = "cc",
                  axis_title_face = "bold", 
                  axis_text_size = 16, 
                  axis_title_size = 18) +
      theme(panel.grid.major = element_line(colour = "grey90"),
            panel.grid.minor = element_line(colour = "grey90"),
            strip.text = element_text(size = 18, face = "bold"),
            legend.title = element_text(size = 16, face = "bold"),
            legend.text = element_text(size = 16)) +
      labs(x = "z-score",
           y = "Fluorescence (scaled)",
           colour = "Cluster")
ggsave(here("results/graphs", "_qq_normality.tiff"), dpi = 300, height = 10, width = 10)
      
# Clusters generally normally distributed.

res.aov <- stats_comparison
tidy_tukey <- res.aov %>% 
      nest(data = c(kMeansCluster,value)) %>% 
      mutate(model = map(data, ~TukeyHSD(aov(value ~ kMeansCluster, .))),
             tidy = map(model, broom::tidy))%>%
      select(Gene, tidy) %>% 
      unnest(tidy)
write.csv(tidy_tukey, here("results/csv", "tukey_sig.csv"), row.names = F)


tidy_resid <- res.aov %>% 
      group_by(Gene) %>%
      nest(data = c(kMeansCluster,value)) %>% 
      mutate(model = map(data, ~lm(value ~ kMeansCluster, .)),
             # tidy = map(model, broom::tidy),
             resid = map(model, residuals),
             fitted = map(model, fitted)) %>%
      select(Gene, resid, fitted) %>% 
      unnest(c(resid, fitted)) %>%
      ungroup() %>%
      mutate(Cluster = ifelse(fitted > 0.5, "3",
                              ifelse((Gene == "N" | Gene == "ORF1ab") & fitted < -0.25, "1",
                                     ifelse(Gene == "S" & fitted > -0.5 & fitted < -0.25, "1", "2"))))
###
tidy_resid %>%
      ggplot(aes(fitted, resid, col = Cluster)) +
      geom_point(alpha = 0.5) +
      facet_grid(Gene~.) +
      theme_ipsum(axis_title_just = "cc",
                  axis_title_face = "bold", 
                  axis_text_size = 16, 
                  axis_title_size = 18) +
      theme(panel.grid.major = element_line(colour = "grey90"),
            panel.grid.minor = element_line(colour = "grey90"),
            strip.text = element_text(size = 18, face = "bold"),
            legend.title = element_text(size = 16, face = "bold"),
            legend.text = element_text(size = 16)) +
      labs(x = "Fitted values",
           y = "Residuals",
           colour = "Cluster")
ggsave(here("results/graphs", "residuals_fitted.tiff"), dpi = 300, height = 12, width = 14)

var_ratio <- tidy_resid %>%
      group_by(Gene, Cluster) %>%
      summarise((round(var(resid),2))) %>%
      rename(Variance = 3) %>%
      mutate(Standardised_variance = ifelse(Gene == "N", round(Variance/0.14,2), 
                              ifelse(Gene == "ORF1ab", Variance/0.2, Variance/0.1)))
write.csv(var_ratio, here("results/csv", "var_ratio.csv"), row.names = F)