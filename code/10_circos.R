library(data.table)
library(dplyr)
library(tidyr)
library(circlize)
library(RColorBrewer)

plotting_data <- fread("data/compensatory_pairs.csv")

# Filter relevant rows (same population group + SNV types)
data <- plotting_data %>%
  filter(grpmax_joint_MIR == grpmax_joint_TAR,
         snv_type %in% c("snp_both", "mixed"))

# Format chromosome names
data$Chr_MIR <- paste0("chr", gsub("chr", "", data$Chr_MIR))
data$Chr_TAR <- paste0("chr", gsub("chr", "", data$Chr_TAR))

# Initialize Circos
circos.par(start.degree = 106)
circos.initializeWithIdeogram(plotType = NULL)

# Outer track with colored chromosome labels
circos.track(
  ylim = c(0, 1), track.height = 0.10, bg.border = NA,
  panel.fun = function(x, y) {
    chr <- gsub("chr", "", CELL_META$sector.index)
    xlim <- CELL_META$xlim
    ylim <- CELL_META$ylim
    circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))
    circos.text(mean(xlim), mean(ylim), chr, cex = 0.6, col = "white", facing = "inside", niceFacing = TRUE)
  }
)

# Add ideogram
circos.genomicIdeogram(track.height = 0.08)

# OCS value → color mapping
color_levels <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
color_palette <- c("#66C2FF", "#3399FF", "#9966CC", "#CC3399", "#FF6633", "#FF0000")

map_color <- function(value) {
  index <- findInterval(unlist(value), color_levels, rightmost.closed = TRUE)
  return(color_palette[index])
}

# Prepare point data (miRNA positions)
point_data <- data.frame(
  chr = data$Chr_MIR,
  start = data$Start_mirna,
  end = data$Start_mirna + 1,
  value = pmin(pmax(data$OCS, 0.4), 1.0)
) %>% distinct()

# Plot OCS points
circos.genomicTrack(
  point_data, numeric.column = 4, jitter = 2,
  ylim = c(0, 1), track.height = 0.1,
  panel.fun = function(region, value, ...) {
    circos.genomicPoints(region, value, col = map_color(value), pch = 16, cex = 1.2)
  }
)

# Plot miRNA labels
mirna_data <- data %>%
  select(Chr_MIR, Start_mirna, mirna_name) %>%
  distinct() %>%
  mutate(mirna_name = gsub("hsa-", "", mirna_name))

if (nrow(mirna_data) > 0) {
  circos.track(
    factors = mirna_data$Chr_MIR,
    x = mirna_data$Start_mirna,
    ylim = c(0, 1),
    track.height = 0.1,
    panel.fun = function(x, y) {
      chr <- CELL_META$sector.index
      subset <- mirna_data[mirna_data$Chr_MIR == chr, ] %>% arrange(Start_mirna)
      last_pos <- -Inf
      shift_bp <- 1e7
      
      for (i in seq_len(nrow(subset))) {
        pos <- subset$Start_mirna[i]
        if (abs(pos - last_pos) < shift_bp) {
          pos <- last_pos + shift_bp
        }
        last_pos <- pos
        circos.text(pos, 0.5, subset$mirna_name[i], sector.index = chr,
                    facing = "clockwise", niceFacing = TRUE, col = "black", cex = 0.3)
      }
    }
  )
}

# Color setup for population groups
group_levels <- c("afr", "amr", "eas", "mid", "nfe", "sas")
group_colors <- setNames(c("darkgrey", "darkblue", "forestgreen", "orange", "red", "violet"), group_levels)

# Filter relevant groups and compute transparency
filtered_data <- data %>%
  filter(grpmax_joint_MIR %in% group_levels) %>%
  mutate(
    max_allele_freq = pmax(Allele_freq_MIR, Allele_freq_TAR, na.rm = TRUE),
    alpha = 0.2 + (max_allele_freq - min(max_allele_freq, na.rm = TRUE)) /
      (max(max_allele_freq, na.rm = TRUE) - min(max_allele_freq, na.rm = TRUE)) * 0.8
  )

# Draw links between miRNAs and targets
link_data_1 <- data.frame(chr = filtered_data$Chr_MIR,
                          start = filtered_data$Start_mirna,
                          end = filtered_data$End_mirna)

link_data_2 <- data.frame(chr = filtered_data$Chr_TAR,
                          start = filtered_data$target_start,
                          end = filtered_data$target_end)

for (i in seq_len(nrow(filtered_data))) {
  circos.genomicLink(
    link_data_1[i, ], link_data_2[i, ],
    col = adjustcolor(group_colors[filtered_data$grpmax_joint_MIR[i]], alpha.f = filtered_data$alpha[i]),
    lwd = 1
  )
}

# Clear plot
circos.clear()
