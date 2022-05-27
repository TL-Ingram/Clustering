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
###


### Point to directory .csv is contained within
setwd("C:/Users/tom.ingram/OneDrive - Medicines Discovery Catapult/Desktop/R scripts/Kit threshold script/Full project_1/Kit changes per month")
dir <- read_csv("Kit_change_per_month_data_edited.csv")
dir$`Creation Date` <- ymd_hms(dir$`Creation Date`)
dir$`Creation Date` <- month(dir$`Creation Date`)

### Cleaning and pre-processing data 
raw <- dir %>%
        select(34,37) %>%
        rename("Kit" = 1, "Month" = 2) %>%
      distinct(,.keep_all = TRUE)
with(raw, month.abb[Month])
dist_data <- transform(raw, MonthAbb = month.abb[Month])
x <- c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")
dist_data <- dist_data %>%
      mutate(MonthAbb = factor(MonthAbb, levels = x)) %>%
      arrange(MonthAbb)
summary <- dist_data %>%
      group_by(MonthAbb) %>%
      summarise(count = n())



bar_graph <- summary %>%
        ggplot(aes(x=MonthAbb, y=count)) +
        geom_col(color = "white", fill = "black", alpha = .5, width = 0.75) +
        theme_ipsum(axis_title_just = "cc", 
                    axis_title_face = "bold", 
                    axis_text_size = 30, 
                    axis_title_size = 38) +
        theme(legend.position="none",
              legend.title = element_text(size=30),
              legend.text = element_text(size=24), 
              plot.title = element_text(size=), 
              axis.text.x=element_text(size = 26, angle=45, vjust=+0.5)) +
        ggtitle("") +
        xlab("Month") +
        ylab("Kit Lots used")
ggsave(width = 20, height = 20, dpi = 600, file="Lots_per_month_Paper.PNG")
bar_graph
