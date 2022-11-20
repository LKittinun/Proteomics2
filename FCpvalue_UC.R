library(janitor)

# DEP --------------------------------------------------------------------------

# BiocManager::install("DEP")
library(DEP)
library(SummarizedExperiment)
library(dplyr)

# The raw data from Aj. Aim is log2 transformed, so we need to transform back first
UC_dfpwr2 <- UC_df_HHPV %>% mutate(across(.cols = starts_with(c("nc","x")), 
                                          ~ifelse(.x == 0, NA, .))) %>% 
                  mutate(across(.cols = starts_with(c("nc","x")), ~2^(.x)))
glimpse(UC_dfpwr2)

unique_UC_df <- make_unique(UC_dfpwr2, "gene_names", "entry")
glimpse(unique_UC_df)
colnames(UC_dfpwr2)
UC_dfpwr$gene_names

data_col <- UC_dfpwr2 %>% select(nc_ev1:x3_ev3) %>% colnames
data_col
which(colnames(unique_UC_df) %in% data_col) # 7:18

metadata_UC <- data.frame(label = data_col, 
                        condition = sub("([a-z0-9]+_[a-z]+).*", "\\1", data_col),
                        replicate = rep(1:3,4)
                        )
metadata_UC

# Make_SE performs log2 transformation to the data
df_SE <- make_se(unique_UC_df, 
                 which(colnames(unique_UC_df) %in% data_col), 
                 metadata_UC)
get_df_wide(df_SE) %>% View()

View(assay(df_filt))
plot_numbers(df_SE)
plot_coverage(df_SE)
plot_missval(df_SE)
plot_detect(df_SE)

df_filt <- filter_missval(df_SE, thr = 1)
plot_missval(df_filt)
get_df_wide(df_filt) %>% View

df_impute <- impute(df_filt, fun = "knn") #impute data
plot_missval(df_impute)
get_df_wide(df_impute) %>% View()

df_norm <- normalize_vsn(df_impute) #normalize the data
get_df_wide(df_norm) %>% View

plot_normalization(df_impute, df_norm)
meanSdPlot(df_impute)
meanSdPlot(df_norm)

data_diff_all_contrasts_norm <- test_diff(df_norm, type = "all") #normalized data

data_diff_all_contrasts <- test_diff(df_impute, type = "all") #un-normalized data

data_diff_mini <- analyze_dep(df_impute, type = "all")

dep <- add_rejections(data_diff_all_contrasts_norm, alpha = 0.05, lfc = 2) #labeled data
get_df_wide(dep)
### Visualization
library(ggplot2)
library(dplyr)

## Plot p-value
plot_p_hist(dep, adjust = FALSE, wrap = TRUE) + theme_bw() # un-adjusted
plot_p_hist(dep, adjust = TRUE, wrap = TRUE) + theme_bw() # adjusted
ggsave("Result/p_value_adjusted_hist.jpg", device = "jpeg", 
       width = 1920/72, height = 1080/72)

## Heatmap
jpeg("Heatmap.jpeg",width=2000,height=3000,units="px", pointsize = 1, res=300)
plot_heatmap(dep, kmeans = FALSE, clustering_distance = "euclidean") 
dev.off()

##
UC_matrix <- UC_df_HHPV %>% select(nc_ev1:x3_ev3) %>% 
  as.matrix 
rownames(UC_matrix) <-  UC_df_HHPV$entry

jpeg("Heatmap.jpeg",width=2000,height=10000, units="px", pointsize = 1, res=300)
Heatmap(UC_matrix, name = "Intensity", row_dend_width = unit(4, "cm"),
        column_dend_height = unit(5, "cm"),
        row_names_gp = gpar(fontsize = 5), km = 20)
dev.off()
# Plot all graph
plot_all(dep, plots = "volcano")

# Get list of sig gene
dep_df <- get_results(dep) 
dep_list <- unique(colnames(dep_df)[grepl("_p.val" , colnames(dep_df), fixed = TRUE)]) %>% 
            gsub("(.)_p.val", "\\1", .)
dep_list

# Write a list of significant proteins
library(purrr)
dep_sig <- setNames(map(dep_list, ~plot_volcano(dep, contrast = .x, adjusted = TRUE, plot = FALSE) %>% 
              arrange(desc(significant))), dep_list)
dep_sig
walk2(dep_sig, dep_list, ~write.csv(.x, file = paste0(.y, ".csv")))

walk(dep_list, function(x) {
  plot_volcano(dep, contrast = x, adjusted = FALSE, plot = TRUE) +
  geom_point(color = "darkred") 
  ggsave(filename = paste0(x, ".png"))} ) 


# Test 

test <- plot_volcano(dep, contrast = "nc_ev_vs_x3_ev", adjusted = TRUE, plot = FALSE) %>% 
  filter(significant == TRUE)
test

cat(paste0(unlist(strsplit(test$protein, split = " ")),"\n"))
# Corplot
colData(dep)

plot_cor(dep, pal = "RdBu", significant = TRUE, indicate = "condition")
?plot_cor

# EnrichR ---------------------------------------------------------------------------------------------------------
library(enrichR)
vignette("enrichR")
db <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021")

enriched_GO <- enrichr(test$protein, databases = db)
dblist <- listEnrichrDbs()
dblist
write.csv(dblist, file = "dblist.csv")
listEnrichrSites()
# test
plotEnrich(enriched_GO[[1]], y = "Ratio", orderBy = "Combined.Score")
plotEnrich(enriched_GO[[2]], y = "Ratio", orderBy = "Combined.Score")
plotEnrich(enriched_GO[[3]], y = "Ratio", orderBy = "Combined.Score")

# make enrichR list
dep_sig_protein <- map(dep_sig, ~.x %>% filter(significant == TRUE) %>% 
                         select(protein) %>% pull() %>% strsplit(split = " ") %>% unlist)
dep_sig_protein

enrich_list_GO <- map(dep_sig_protein, ~enrichr(.x, databases = db)) %>% 
                          unlist(., recursive = FALSE)
enrich_list_KEGG <- map(dep_sig_protein, ~enrichr(.x, databases = "KEGG_2021_Human")) %>% 
                          unlist(., recursive = FALSE)

library(ggplot2)
# raw data
iwalk(enrich_list_GO, ~write.csv(.x, file = paste0(.y, ".csv")))
iwalk(enrich_list_KEGG, ~write.csv(.x, file = paste0(.y, ".csv")))

# by combine score
iwalk(enrich_list_GO, function(x, z){
  plotEnrich(x, y = "Ratio", orderBy = "Combined.Score", title = z) 
  ggsave(filename = paste0(z, ".png"))
})

iwalk(enrich_list_KEGG, function(x, z){
  plotEnrich(x, y = "Ratio", orderBy = "Combined.Score", title = z) 
  ggsave(filename = paste0(z, ".png"))
})

# by P-value
iwalk(enrich_list_GO, function(x, z){
  plotEnrich(x, y = "Ratio", orderBy = "P.value", title = z) 
  ggsave(filename = paste0(z, ".png"))
})

iwalk(enrich_list_KEGG, function(x, z){
  plotEnrich(x, y = "Ratio", orderBy = "P.value", title = z) 
  ggsave(filename = paste0(z, ".png"))
})


# STRINGdb ----------------------------------------------------------------
library(STRINGdb)
library(purrr)
test
test_protein <- data.frame(test = unlist(strsplit(test$protein, " ")))
cat(test_protein$test, sep = "\n")

string_db <- STRINGdb$new(version = "11", 
                          species = 9606, score_threshold = 200, 
                          input_directory = "")

# test
mapped <- string_db$map(test_protein, "test", removeUnmappedRows = TRUE )
string_db$plot_network(mapped)
string_db$get_png(mapped, file = "test.png")

# PCA -------------------------------------------------------------------------------------------------------------
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/factoextra")
library(corrplot)
library(factoextra)

forPCA_norm <- get_df_wide(df_norm) %>% select(nc_ev_1:x3_ev_3) %>% as.matrix %>% t
colnames(forPCA_norm) <- get_df_wide(df_norm)$name

## PCAtools
library(PCAtools)
library(ggplot2)
x <- forPCA_norm
View(x)
metadata <- data.frame(row.names = rownames(x))
rownames(x)
metadata <- metadata %>% mutate(Group = rep(c("nc", "x1", "x2", "x3"), each = 3))
metadata

pc <- pca(x, metadata, center = TRUE, scale = FALSE, transposed = TRUE, removeVar = 0.1)
str(pc)
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
View(df)
library(ggplot2)
df %>%filter(gene_names == "HYLS1 HLS") %>%  select(gene_names)
screeplot(pc,
          components = getComponents(pc, 1:12),
          vline = c(horn$n, elbow, which(cumsum(pc$variance) > 80)[1]), hline = 80) +
  geom_label(aes(x = horn$n, y = 50,
                 label = 'Horn\'s', vjust = -1, size = 8)) +
  geom_label(aes(x = elbow, y = 50,
                 label = 'Elbow method', vjust = -1, size = 8)) 

which(cumsum(pc$variance) > 80)[1] # 6

ggsave("Result/Screeplot.jpg", device = "jpeg", dpi = 300)


# Kmeans
library(cluster)
library(fpc)
library(dplyr)
library(tribble)
library(ggplot2)
library(factoextra)

km <- kmeans(forPCA_norm, 4, iter.max = 1000)
fviz_cluster(km, forPCA_norm) + theme_bw() + ggtitle("Kmean-clustering")

set.seed(123)

# function to compute total within-cluster sum of square 
library(purrr)
wss <- function(k) {
  kmeans(forPCA_norm, k, nstart = 10 )$tot.withinss
}
View(forPCA_norm)
# Compute and plot wss for k = 1 to k = 15
k.values <- 1:6

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

fviz_nbclust(forPCA_norm, kmeans, method = "silhouette")
ggsave("Result/Optimal_cluster.jpg", device = "jpeg")

