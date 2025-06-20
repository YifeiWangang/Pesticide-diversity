library(ggplot2)
library(ggalluvial)
library(dplyr)

# Load data
data <- read.table(file.choose(), header = TRUE, sep = "\t", row.names = 1)

# Convert data to alluvial format
convertion <- to_lodes_form(data,                        
                            axes = 1:ncol(data),                       
                            key = "x",                       
                            value = "stratum",                       
                            id = "alluvium")

# Create alluvial plot
ggplot(convertion,
       aes(x = x, stratum = stratum, alluvium = alluvium,
           fill = stratum, label = stratum)) +
  scale_fill_manual(values = color_palette) +
  geom_flow(width = 0.4, aes.flow = "forward") +
  # Customize strata: add white borders
  geom_stratum(alpha = 0.8, width = 0.4, color = "white", size = 0.3) +
  # Customize labels: smaller font and adjusted position
  geom_text(stat = "stratum", 
            size = 3,  # reduced font size
            color = "black", 
            family = "serif",
            position = position_nudge(y = -0.05)) +  # slight downward adjustment
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(size = 15, family = "serif", colour = "black")
  ) +
  xlab("") +
  ylab("") +
  ggtitle("") +
  guides(fill = "none")  # Remove legend