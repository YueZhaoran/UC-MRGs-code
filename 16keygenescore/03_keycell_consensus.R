library(png)
library(grid)

if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
project_dir <- dirname(getwd())

a <- readPNG("02_Plots/05_Method_Comparison_Plate.png")
b <- readPNG("02_Plots/04_Violin_Combined_FinalScore_by_Group.png")

w_cm <- 15.2
h_a <- w_cm * dim(a)[1] / dim(a)[2]
h_b <- w_cm * dim(b)[1] / dim(b)[2]
h_cm <- h_a + h_b + 0.6

pdf("02_Plots/06_KeyCell_Consensus_Plate.pdf", width = w_cm / 2.54, height = h_cm / 2.54)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(h_a, h_b), "cm"))))
pushViewport(viewport(layout.pos.row = 1))
grid.raster(a, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = FALSE)
grid.text("A", x = unit(1, "mm"), y = unit(1, "npc") + unit(1, "mm"), just = c("left", "bottom"),
          gp = gpar(fontface = "bold", fontsize = 11))
popViewport()
pushViewport(viewport(layout.pos.row = 2))
grid.raster(b, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = FALSE)
popViewport(2)
dev.off()

png("02_Plots/06_KeyCell_Consensus_Plate.png", width = w_cm, height = h_cm, units = "cm", res = 600)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(h_a, h_b), "cm"))))
pushViewport(viewport(layout.pos.row = 1))
grid.raster(a, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = FALSE)
grid.text("A", x = unit(1, "mm"), y = unit(1, "npc") + unit(1, "mm"), just = c("left", "bottom"),
          gp = gpar(fontface = "bold", fontsize = 11))
popViewport()
pushViewport(viewport(layout.pos.row = 2))
grid.raster(b, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = FALSE)
popViewport(2)
dev.off()
