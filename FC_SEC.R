library(dplyr)
library(janitor)
library(purrr)
library(tidyr)
library(VennDiagram)
library(RColorBrewer)

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# Max value
SEC_df_HHPV_max <- SEC_df_HHPV %>% rowwise() %>% 
  mutate(nc_max = max(c_across(starts_with("nc")))) %>% 
  mutate(x1_max = max(c_across(starts_with("x1")))) %>%
  mutate(x2_max = max(c_across(starts_with("x2")))) %>%
  mutate(x3_max = max(c_across(starts_with("x3")))) %>%
  select(-(nc_sec8:x3_sec10)) %>% 
  ungroup()
  
# Create Name list
name_list <- SEC_df_HHPV_max %>% 
      select(nc_max:x3_max) %>% 
      names()

dim(SEC_df_HHPV_max)
(unique(SEC_df_HHPV_max$entry_name))
(unique(SEC_df_HHPV_max$entry))

# Extract protein list

mapped_entry <- SEC_df_HHPV_max %>% select(entry, nc_max:x3_max) %>% 
 pivot_longer(-entry, names_to = "category", values_to = "intensity") %>% 
  group_by(category) %>% 
  filter(intensity != 0) %>% 
  nest() 

mapped_entry$data
SEC_protein_list <- setNames(map(mapped_entry$data, ~.x %>% pull(entry)), name_list)
SEC_protein_list

# Venn Diagram ------------------------------------------------------------

## Palette
Col4 <- brewer.pal(4, "Dark2")
Col3 <- brewer.pal(3, "Dark2")
Col2 <- Col3[1:2]
display.brewer.all()
names(SEC_protein_list)

## Drawing Venn

# All stage
draw_4venn(SEC_protein_list, file = "All_stage.jpg", 
           category.names = c("NC", "Stage 1", "Stage 2", "Stage 3"))

# Generate pair list
pair_list <- expand.grid(name_list, name_list) %>% 
  filter(Var1 != Var2) %>% 
  filter(as.character(Var1) < as.character(Var2)) %>% # Keep only unique pair 
  t %>% as.data.frame %>% 
  map(. , ~.x[])
names(pair_list) <- unlist(map(pair_list, ~paste(.x, collapse = "_vs_")))

# Draw pair venn
map(pair_list, ~draw_2venn(SEC_protein_list[.x], 
                           file = paste(.x[1], "vs", .x[2], ".jpg", sep = "_"),
                           category.names = c(.x[1],.x[2])))
      
SEC_diff_list <- map(pair_list, ~SEC_protein_list[.x]) %>% 
  map(~gplots::venn(.x, show.plot = FALSE) %>% 
        attributes)
SEC_diff_list

# CSV
SEC_diff_list
SEC_write <- map(SEC_diff_list, ~.x[][[3]][]) %>% 
  map(~map(., ~.x %>% tibble(entry = .))) %>% 
  unlist(recursive = FALSE)

names(SEC_write) <- gsub("\\.", "_", names(SEC_write))

iwalk(SEC_write, ~write_csv(.x, file = paste0(.y, ".csv")))


SEC_intersect_write <- map(SEC_list_intersect, ~tibble(entry = .))

iwalk(SEC_intersect_write, ~write_csv(.x, file = paste0(.y, ".csv")))

# DEG ---------------------------------------------------------------------
# Create intersection list

SEC_list_intersect <- map(pair_list, ~SEC_protein_list[.x]) %>% 
                            map(~gplots::venn(.x, show.plot = FALSE) %>% 
                            attributes %>% 
                            .[["intersections"]] %>% .[[1]])
names(SEC_list_intersect) <- names(pair_list)

iwalk(SEC_list_intersect, ~write_csv(.x, file = paste0(.y, ".csv")))

# Equal to VennDiagram
sapply(SEC_list_intersect, length) 
pair_list
# 

library(ggplot2)
library(readr)

walk2(pair_list, SEC_list_intersect,
      ~fc_calculator(.x, .y, data = SEC_df_HHPV_max, 50))

# Heatmap
library(ComplexHeatmap)

SEC_max_matrix<- SEC_df_HHPV_max %>% select(nc_max:x3_max) %>% 
                  as.matrix 
rownames(SEC_max_matrix) <-  SEC_df_HHPV_max$entry

SEC_nctox3 <- SEC_df_HHPV_max %>% select(entry, nc_max:x3_max) %>% 
                        filter((nc_max > x1_max) & (x1_max > x2_max) & (x2_max > x3_max))
SEC_nctox3_matrix <- SEC_nctox3 %>% select(nc_max:x3_max) %>% as.matrix 
rownames(SEC_nctox3_matrix) <- SEC_nctox3$entry

SEC_x3tonc <- SEC_df_HHPV_max %>% select(entry, nc_max:x3_max) %>% 
                        filter((nc_max < x1_max) & (x1_max < x2_max) & (x2_max < x3_max))
SEC_x3tonc_matrix <- SEC_x3tonc %>% select(x3_max:nc_max) %>% as.matrix 
rownames(SEC_x3tonc_matrix) <- SEC_x3tonc$entry

jpeg("SEC_Heatmap.jpeg",width=2000,height=10000, units="px", pointsize = 1, res=300)
Heatmap(SEC_max_matrix, name = "Intensity", row_dend_width = unit(4, "cm"),
        column_dend_height = unit(5, "cm"),
        row_names_gp = gpar(fontsize = 5), km = 20)
dev.off()

jpeg("SEC_nctox3.jpeg",width=2000,height=10000, units="px", pointsize = 1, res=300)
Heatmap(SEC_nctox3_matrix, name = "Intensity", row_dend_width = unit(4, "cm"),
             column_dend_height = unit(5, "cm"),
        row_names_gp = gpar(fontsize = 5), km = 20)
dev.off()

jpeg("SEC_x3tonc.jpeg",width=2000,height=3000, units="px", pointsize = 1, res=300)
Heatmap(SEC_x3tonc_matrix, name = "Intensity", row_dend_width = unit(4, "cm"),
        column_dend_height = unit(5, "cm"), colorRamp2(c(0, 12, 25)),
        row_names_gp = gpar(fontsize = 5), km = 20)
dev.off()

# PCA ---------------------------------------------------------------------

## PCAtools
library(PCAtools)
library(ggplot2)
x <- SEC_df_HHPV %>% select(nc_sec8:x3_sec10) %>% t
colnames(x) <- SEC_df_HHPV$entry
x
metadata <- data.frame(row.names = rownames(x))
rownames(x)
metadata <- metadata %>% mutate(Group = rep(c("nc", "x1", "x2", "x3"), each = 3))
metadata

pc <- pca(x, metadata, center = TRUE, scale = FALSE, transposed = TRUE, removeVar = 0.1)

screeplot(pc, axisLabSize = 18, components = getComponents(pc,1:12))

biplot(pc, colby = "Group", 
       colkey = c("forestgreen", "purple", "darkblue", "red"), 
       encircle = TRUE,
       hline = 0, vline = 0,
       showLoadings = FALSE,
)
ggsave("Result/PCA.jpg", device = "jpeg")

plotloadings(pc, drawConnectors=TRUE, returnPlot = TRUE)

pairsplot(pc,
          components = getComponents(pc, c(1:6)),
          triangle = TRUE, trianglelabSize = 20,
          hline = 0, vline = 0,
          pointSize = 1,
          gridlines.major = FALSE, gridlines.minor = FALSE,
          colby = 'Group',
          title = 'Pairs plot', plotaxes = FALSE,
          margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm'))
ggsave("Result/Pairsplot.jpg", device = "jpeg")


elbow <- findElbowPoint(pc$variance)
elbow
horn <- parallelPCA(t(x))
horn

library(ggplot2)

screeplot(pc,
          components = getComponents(pc, 1:12),
          vline = c(horn$n, elbow, which(cumsum(pc$variance) > 80)[1]), hline = 80) +
  geom_label(aes(x = horn$n+0.5, y = 50,
                 label = 'Horn\'s', vjust = -1, size = 8)) +
  geom_label(aes(x = elbow, y = 50,
                 label = 'Elbow method', vjust = -1, size = 8)) 

which(cumsum(pc$variance) > 80)[1] # 8

ggsave("Result/Screeplot.jpg", device = "jpeg", dpi = 300)

# Kmeans ------------------------------------------------------------------
library(cluster)
library(fpc)
library(dplyr)
library(tribble)
library(ggplot2)
library(factoextra)

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:6

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, ~wss(x, .x))

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

fviz_nbclust(x, kmeans, method = "silhouette")
ggsave("Result/Optimal_cluster.jpg", device = "jpeg")

km <- kmeans(x, 4, iter.max = 1000)
km$cluster
View(x)
cluster <- fviz_cluster(km, x) + theme_bw() + ggtitle("Kmean-clustering") + 
  geom_text_repel(label = rownames(x), aes(col = factor(km$cluster))) 
cluster$layers[[4]] <- NULL
cluster

ggsave("Result/Kmeans_4.jpg", device = "jpeg")


