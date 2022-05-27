###
library(tidyverse)
library(hrbrthemes)
library(scales)
library(rstatix)
library(data.table)
library(lubridate)
library(ggrepel)
library(reshape2)
library(viridis)
library(purrr)
library(tidyr)
library(here)
library(ggridges)
###


### Point to directory .csv is contained within
dir <- read_csv(here("pc_data", "PosCon data DEC2020-OCT2021.csv"))
dir$Experiment.created.at..UTC. <- ymd_hms(dir$Experiment.created.at..UTC.)
dir$Experiment.created.at..UTC. <- date(dir$Experiment.created.at..UTC.)


### Removing non-normal or misnamed (therefore duplicated) plates
dir <- dir[!grepl("AD20", dir$File.name),]
dir <- dir[!grepl("Training", dir$File.name),]
dir <- dir[!grepl("COMPETENCY", dir$File.name),]
dir <- dir[!grepl("36H3U0E0FF", dir$File.name),]
dir <- dir[!grepl("38H3N03LE9", dir$File.name),]
dir <- dir[!grepl("2021-06-10 212656", dir$File.name),]
dir <- dir[!grepl("2021-09-19 113545", dir$File.name),]
dir <- dir[!grepl("ML-", dir$File.name),]
dir <- dir[!grepl("QC2", dir$File.name),]
dir <- dir[!grepl("QC1", dir$File.name),]
dir <- dir[!grepl("QC_", dir$File.name),]
dir <- dir[!grepl("QC-", dir$File.name),]
dir <- dir[!grepl("VP2", dir$File.name),]
dir <- dir[!grepl("VP3", dir$File.name),]
dir <- dir[!grepl("VP5", dir$File.name),]
dir <- dir[!grepl("VT_H3U0", dir$File.name),]
dir <- dir[!grepl("PCR_75_test_1", dir$File.name),]
dir <- dir[!grepl("WS3_", dir$File.name),]
dir <- dir[!grepl("Competency", dir$File.name),]
dir <- dir[!grepl("AUTO", dir$File.name),]
dir <- dir[!grepl("AP-pcr-", dir$File.name),]
dir <- dir[!grepl("AT_", dir$File.name),]
dir <- dir[!grepl("concordance", dir$File.name),]

### Cleaning and pre-processing data 
raw <- dir %>%
        filter(Experiment.created.at..UTC. >= "2020-12-01") %>%
        filter(Result == "Positive") %>%
        filter(Well == "O23" | Well == "O24" | Well == "P23" | Well == "P24") %>%
        distinct(,.keep_all = TRUE) %>%
        select(4,9,10,17,20,23,27,32,37,42) %>%
        rename("ORF1ab" = `Channel_1_EndFluoBC`, "N" = `Channel_2_EndFluoBC`, "S" = `Channel_3_EndFluoBC`, "MS2" = `Channel_4_EndFluoBC`) %>%
        gather(key = "Gene", value = "Fluorescence", 7,8,9,10) %>%
        filter(Fluorescence >= 15000)
raw$Instrument.ID <- as.character(raw$Instrument.ID)

### Extreme outlier removal: values above Q3 + 3xIQR or below Q1 - 3xIQR are considered as extreme points (or extreme outliers).
output <- raw %>%
        group_by(Instrument.ID, Gene) %>%
        identify_outliers(Fluorescence)
output$is.extreme[output$is.extreme %like% "TRUE"] <- "yes"
output$is.extreme[output$is.extreme %like% "FALSE"] <- "no"
combine_df <- left_join(raw, output)
combine_df$is.extreme[is.na(combine_df$is.extreme)] <- "no"

output %>%
      group_by(`is.extreme`) %>%
      summarise(n = n())
### Development of thresholds using AP MU determination (March 2021)
instrument_summary <- combine_df %>%
        filter(`is.extreme` != "yes") %>%
        group_by(Instrument.ID, Gene) %>%
        summarise(Count = n(),
                  mean = mean(Fluorescence, na.rm = T),
                  `1SD` = sd(Fluorescence, na.rm = T),
                  Threshold = `mean`-`1SD`,
                  First_quantile = quantile(Fluorescence, 0.25, na.rm = T),
                  Third_quantile = quantile(Fluorescence, 0.75, na.rm = T),
                  max = max(Fluorescence),
                  min = min(Fluorescence)) %>%
        transform(Threshold_exp_uncert=ifelse(Gene=="ORF1ab", Threshold*0.9649, 
                                              ifelse(Gene=="N", Threshold*0.9398,
                                                     ifelse(Gene=="S", Threshold*0.9404, Threshold*1))))
#write.csv(instrument_summary, "Summary_BC.csv", row.names = FALSE)
combine_df_2 <- left_join(combine_df, instrument_summary)
combine_df_2$Experiment.created.at..UTC. <- month(combine_df_2$Experiment.created.at..UTC.)
combine_df_2 <- combine_df_2 %>%
        rename("Month" = Experiment.created.at..UTC.)

n_month_table <- combine_df_2 %>%
      filter(!(`is.extreme` == TRUE)) %>%
      group_by(Month) %>%
      summarise(n = n())
#write.csv(n_month_table, "n_month_table.csv", row.names = FALSE)

dist_data <- combine_df_2 %>%
        select(3,7,8,9,10)

with(dist_data, month.abb[Month])
dist_data <- transform(dist_data, MonthAbb = month.abb[Month])
x <- c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")
dist_data <- dist_data %>%
      mutate(MonthAbb = factor(MonthAbb, levels = x)) %>%
      arrange(MonthAbb)


#####################################DONE####################################################
dist_graph <- dist_data %>%
        filter(is.extreme != "yes") %>%
     filter(Gene != "MS2") %>%
     mutate(MonthAbb = fct_reorder(MonthAbb, desc(MonthAbb))) %>%
        ggplot(aes(x=Fluorescence, y=MonthAbb, fill = Gene)) +
     geom_density_ridges(scale = 2, rel_min_height = 0.01, alpha = 0.5, size = 0.5, colour = "grey50") +
        facet_grid(.~ Gene, scales = "free") +
      scale_x_continuous(labels = scales::comma) +
        scale_fill_manual(values = c("yellow3", "#52B119", "orange")) +
     theme_ipsum(axis_title_just = "cc",
                 axis_title_face = "bold", 
                 axis_text_size = 16, 
                 axis_title_size = 18) +
     theme(panel.grid.major = element_line(colour = "grey90"),
           panel.grid.minor = element_line(colour = "grey90"),
           strip.text = element_text(size = 18, face = "bold"),
           legend.position = "none") +
     guides(x = guide_axis(angle = 45)) +
        ggtitle("") +
        xlab("Fluorescence (baseline corrected)") +
        ylab("Month (Dec 2020 : Oct 2021")
ggsave(here("results/graphs", "1_PC_fluor_dist_months.tiff"), dpi = 300, height = 12, width = 14)
dist_graph





##########
##########
##########
################################################################################

S <- dist_data %>%
     filter(Gene == "S")
one.way <- aov(Fluorescence ~ MonthAbb, data = S)
summary(one.way)
Tukey.one.way_S <- TukeyHSD(one.way)
modify_S <- as.data.frame(Tukey.one.way_S[1])
#write.csv(modify_S, "anova_Cq_between_months_S.csv", row.names = TRUE)

subsequent_months_S <- modify_S[row.names(modify_S) %in% c('Jan-Dec','Feb-Jan','Mar-Feb','Apr-Mar','May-Apr','Jul-Jun','Aug-Jul','Sep-Aug','Oct-Sep'), ]
#write.csv(subsequent_months_S, "stats_subsequent_months_S.csv", row.names = TRUE)
subsequent_months_S_stats <- subsequent_months_S %>%
     select(1) %>%
     abs %>%
     summarise(mean = mean(MonthAbb.diff),
               sd = sd(MonthAbb.diff))
#write.csv(subsequent_months_S_stats, "stats_subsequent_months_MEAN+SD_S.csv", row.names = TRUE)

ORF1ab <- dist_data %>%
     filter(Gene == "ORF1ab")
one.way_ORF1ab <- aov(Fluorescence ~ MonthAbb, data = ORF1ab)
Tukey.one.way_ORF1ab <- TukeyHSD(one.way_ORF1ab)
modify_ORF1ab <- as.data.frame(Tukey.one.way_ORF1ab[1])
#write.csv(modify_ORF1ab, "anova_Cq_between_months_ORF1ab.csv", row.names = TRUE)
subsequent_months_ORF1ab <- modify_ORF1ab[row.names(modify_ORF1ab) %in% c('Jan-Dec','Feb-Jan','Mar-Feb','Apr-Mar','May-Apr','Jul-Jun','Aug-Jul','Sep-Aug','Oct-Sep'), ]
#write.csv(subsequent_months_ORF1ab, "stats_subsequent_months_ORF1ab.csv", row.names = TRUE)
subsequent_months_ORF1ab_stats <- subsequent_months_ORF1ab %>%
     select(1) %>%
     abs %>%
     summarise(mean = mean(MonthAbb.diff),
               sd = sd(MonthAbb.diff))
#write.csv(subsequent_months_ORF1ab_stats, "stats_subsequent_months_MEAN+SD_ORF1ab.csv", row.names = TRUE)

N <- dist_data %>%
     filter(Gene == "N")
one.way_N <- aov(Fluorescence ~ MonthAbb, data = N)
Tukey.one.way_N <- TukeyHSD(one.way_N)
modify_N <- as.data.frame(Tukey.one.way_N[1])
#write.csv(modify_N, "anova_Cq_between_months_N.csv", row.names = TRUE)
subsequent_months_N <- modify_N[row.names(modify_N) %in% c('Jan-Dec','Feb-Jan','Mar-Feb','Apr-Mar','May-Apr','Jul-Jun','Aug-Jul','Sep-Aug','Oct-Sep'), ]
#write.csv(subsequent_months_N, "stats_subsequent_months_N.csv", row.names = TRUE)
subsequent_months_N_stats <- subsequent_months_N %>%
     select(1) %>%
     abs %>%
     summarise(mean = mean(MonthAbb.diff),
               sd = sd(MonthAbb.diff))
#write.csv(subsequent_months_N_stats, "stats_subsequent_months_MEAN+SD_N.csv", row.names = TRUE)


# summary_plot <- rbind(subsequent_months_S, subsequent_months_ORF1ab, subsequent_months_N)
# summary_plot$Months <- rownames(summary_plot)
# summary_plot$Gene[summary_plot$Months %like% "1"] <- "ORF1ab"
# summary_plot$Gene[summary_plot$Months %like% "2"] <- "N"
# summary_plot$Gene[is.na(summary_plot$Gene)] <- "S"
# summary_plot$Months <- gsub('[0-9.]', '', summary_plot$Months)
# summary_plot$Months <-factor(summary_plot$Months, levels=c("Jan-Dec", "Feb-Jan", "Mar-Feb", "Apr-Mar", "May-Apr", "Jun-May", "Jul-Jun", "Aug-Jul", "Sep-Aug", "Oct-Sep"))
# 
# ggplot(summary_plot, aes(x=Months, y=MonthAbb.diff, colour = Gene)) + 
#      geom_errorbar(aes(ymin=MonthAbb.lwr, ymax=MonthAbb.upr), width=.2) +
#      geom_line(size = 4) +
#      geom_point(size = 4) +
#      facet_grid(Gene ~ . , scales = "free_y") +
#      scale_y_continuous(labels = scales::comma) +
#      geom_hline(yintercept = 0, size = 1, linetype = "dashed", colour = "#2D708EFF") +
#      scale_colour_manual(values = c("yellow3", "#52B119", "orange")) +
#      theme_ipsum(axis_title_just = "cc", 
#                  axis_title_face = "bold", 
#                  axis_text_size = 18, 
#                  axis_title_size = 30) +
#      theme(legend.position="none",
#            legend.title = element_text(size=30),
#            legend.text = element_text(size=24), 
#            plot.title = element_text(size=), 
#            axis.text.x=element_text(size = 20, angle=45, vjust=+0.5),
#            strip.text = element_text(face="bold", 
#                                      size=24,
#                                      lineheight=5.0),) +
#      guides(x = guide_axis(angle = 45)) +
#      ggtitle("") +
#      xlab("Month comparison") +
#      ylab("Delta Flourescence")
# #ggsave(width = 20, height = 14, dpi = 600, file="Summary_stats_figure_monthly_changes_Fluorescence.PNG")