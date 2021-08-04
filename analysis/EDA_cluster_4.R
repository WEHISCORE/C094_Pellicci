# EDA to look for technical factors driving cluster 4.
# Peter Hickey
# 2021-08-04

library(here)
library(scater)

sce <- readRDS(here("data", "SCEs", "C094_Pellicci.cells_selected.SCE.rds"))

# There's perhaps a bias towards cells being in cluster 4 when they are near
# the edges of the plate, specifically towards columns 1-3 and in particular
# on plate LCE514.
p <- lapply(levels(sce$plate_number), function(p) {
  z <- sce[, sce$plate_number == p]
  plotPlatePosition(
    z,
    as.character(z$well_position),
    point_size = 2,
    point_alpha = 1,
    theme_size = 5,
    colour_by = "cluster") +
    ggtitle(p) +
    theme(
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 6)) +
    guides(
      colour = guide_legend(override.aes = list(size = 4)),
      shape = guide_legend(override.aes = list(size = 4)))
})

wrap_plots(p, ncol = 3) + plot_layout(guides = "collect")

# LCE508 and LCE514 have the most cluster 4 cells.
table(sce$cluster, sce$plate_number)
round(proportions(table(sce$cluster, sce$plate_number), 2), 2)
round(proportions(table(sce$cluster, sce$plate_number), 1), 2)
table(sce$cluster == 4, sce$plate_number)
round(proportions(table(sce$cluster == 4, sce$plate_number), 2), 2)
round(proportions(table(sce$cluster == 4, sce$plate_number), 1), 2)

# Columns 1-3 have the most cluster 4 cells ...
table(sce$cluster, factor(substr(sce$well_position, 2, 3), 1:24))
table(sce$cluster == 4, factor(substr(sce$well_position, 2, 3), 1:24))
p1 <- proportions(table(sce$cluster == 4, factor(substr(sce$well_position, 2, 3), 1:24)), 2)
p1
par(mfrow = c(1, 1))
plot(p1[2, ], type = "b", ylim = c(0, 1))
abline(h = median(p1[2, ]), col = "blue")
# But it depends on plate
table(sce$cluster, factor(substr(sce$well_position, 2, 3), 1:24), sce$plate_number)
table(sce$cluster == 4, factor(substr(sce$well_position, 2, 3), 1:24), sce$plate_number)
p2 <- proportions(table(sce$cluster, factor(substr(sce$well_position, 2, 3), 1:24), sce$plate_number), 2:3)
p2
par(mfrow = c(2, 3))
lapply(levels(sce$plate_number), function(k) {
  plot(p2[4, , k], ylim = c(0, 1), type = "b", main = k, ylab = "prop. cluster 4", xlab = "column")
  abline(h = median(p2[4, , ]), col = "blue")
  abline(h = median(p2[4, , k]), col = "red")
})

# No strong association with row.
table(sce$cluster, substr(sce$well_position, 1, 1))
table(sce$cluster, substr(sce$well_position, 1, 1), sce$plate_number)
table(sce$cluster == 4, factor(substr(sce$well_position, 2, 3), 1:24))
q1 <- proportions(table(sce$cluster == 4, substr(sce$well_position, 1, 1)), 2)
q1
par(mfrow = c(1, 1))
plot(x = q1[2, ], type = "b", ylim = c(0, 1))
abline(h = median(q1[2, ]), col = "blue")
# But it depends on plate
table(sce$cluster, substr(sce$well_position, 1, 1), sce$plate_number)
table(sce$cluster == 4, substr(sce$well_position, 1, 1), sce$plate_number)
q2 <- proportions(table(sce$cluster, substr(sce$well_position, 1, 1), sce$plate_number), 2:3)
q2
par(mfrow = c(2, 3))
lapply(levels(sce$plate_number), function(k) {
  plot(q2[4, , k], ylim = c(0, 1), type = "b", main = k, ylab = "prop. cluster 4", ylab = "row")
  abline(h = median(q2[4, , ], na.rm = TRUE), col = "blue")
  abline(h = median(q2[4, , k], na.rm = TRUE), col = "red")
})


