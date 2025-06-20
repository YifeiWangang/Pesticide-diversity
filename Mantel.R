library(linkET)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Load data files
arg_data <- read.table(file.choose(), header = TRUE, sep = "\t", row.names = 1)
mge_data <- read.table(file.choose(), header = TRUE, sep = "\t", row.names = 1)

# Select specific MGE columns
selected_mge <- mge_data[, c(1, 2, 4, 5)] 
colnames(selected_mge) <- c("tnpA27", "tnpA3", "tnpA", "tnpA2")

# Perform Mantel test (each selected MGE vs all ARGs)
mantel_result <- mantel_test(selected_mge, 
                             arg_data,
                             spec_select = list(
                               tnpA27 = 1,
                               tnpA3 = 2,
                               tnpA = 3,
                               tnpA2 = 4
                             )) %>% 
  mutate(
    rd = cut(r, breaks = c(-Inf, 0.1, 0.3, Inf),
             labels = c("< 0.1","0.1 - 0.3", ">= 0.3")),
    pd = cut(p, breaks = c(-Inf,0.05, Inf),
             labels = c("< 0.05", ">= 0.05")),
    r.sign = case_when(
      r > 0  ~ "Positive",
      r <= 0 ~ "Negative"),
    p.level = ifelse(p < 0.05, "sig", "ns")
  ) %>%
  mutate(
    r.sign = factor(r.sign, levels = c("Positive", "Negative")),
    p.level = factor(p.level, levels = c("sig", "ns"))
  )

# Create visualization
qcorrplot(cor(arg_data), 
          type = "upper", 
          grid.col = "#BFBFBF",
          diag = FALSE) +
  geom_square() +
  geom_couple(
    aes(colour = r.sign, size = rd,linetype = p.level), 
    data = mantel_result, 
    curvature = 0.1
  ) +
  scale_fill_gradientn(
    colours = RColorBrewer::brewer.pal(11, "RdBu")
  ) +
  scale_size_manual(
    values = c(0.5, 1, 2)
  ) +
  scale_colour_manual(
    values = color_pal(3)
  ) +
  guides(
    size = guide_legend(
      title = "Mantel's r",
      override.aes = list(colour = "grey35"),
      order = 2
    ),
    colour = guide_legend(
      title = "Mantel's p",
      override.aes = list(size = 3),
      order = 1
    ),
    fill = guide_colorbar(
      title = "Pearson's r", 
      order = 3
    )
  )