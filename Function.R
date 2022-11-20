## Function
library(VennDiagram)
library(RColorBrewer)

# Color palette --------------------------------------------------------------------------------------------------

Col4 <- brewer.pal(4, "Dark2")
Col3 <- brewer.pal(3, "Dark2")
Col2 <- Col3[1:2]
display.brewer.all()

# Draw function ---------------------------------------------------------------------------------------------------
draw_4venn <- function(lists, filename = NULL, font = 0.25, 
                       fontfamily = "sans", category.names = names(lists)){
  venn.diagram(lists, 
               filename = filename, imagetype = "tiff",
               cex = font, cat.cex = font,
               fontfamily = fontfamily, cat.fontfamilty = fontfamily,
               cat.fontface = "bold",
               lwd = 1.5, fill = Col4, lty = 'blank',
               category.names = category.names,
               height = 700, width = 700, resolution = 500,
               print.mode = c("raw","percent")
  )
}

draw_3venn <- function(lists, filename = NULL, font = 0.25, 
                       fontfamily = "sans", category.names = names(lists)){
  venn.diagram(lists, 
               filename = filename, imagetype = "tiff",
               cex = font, cat.cex = font,
               fontfamily = fontfamily, cat.fontfamilty = fontfamily,
               cat.fontface = "bold",
               lwd = 1.5, fill = Col3, lty = 'blank',
               category.names = category.names,
               height = 500, width = 500, resolution = 500,
               print.mode = c("raw","percent")
  )
}

draw_2venn <- function(lists, filename = NULL, font = 0.25, 
                       fontfamily = "sans", category.names = names(lists),
                       cat.pos = c(-30,30), cat.dist = 0.04, rotation.degree = 0, scaled = FALSE){
  venn.diagram(lists, 
               filename = filename, scaled = scaled, imagetype = "tiff",
               cex = font, cat.cex = font,
               fontfamily = fontfamily, cat.fontfamilty = fontfamily,
               cat.fontface = "bold",
               lwd = 1, fill = Col2, lty = 'blank',
               category.names = category.names,
               cat.pos = cat.pos,
               cat.dist = cat.dist, cat.default.pos = "outer",
               height = 500, width = 500, resolution = 500, 
               rotation.degree = rotation.degree,
               print.mode = c("raw","percent")
  ) 
}

draw_2venn_prop <-function(lists, filename = NULL, cex = 0.25, cat.cex = 0.3,
                           fontfamily = "sans", category.names = names(lists),
                           cat.dist = 0.04, rotation.degree = 0, scaled = FALSE,
                           height = 600, width = 600, margin = 0.05){
  venn.diagram(lists, 
               filename = filename, scaled = scaled, imagetype = "tiff",
               cex = cex, cat.cex = cat.cex,
               fontfamily = fontfamily, cat.fontfamilty = fontfamily,
               cat.fontface = "bold",
               lwd = 1, fill = Col2, lty = 'blank',
               category.names = category.names,
               cat.dist = cat.dist, cat.default.pos = "outer",
               height = height, width = width, resolution = 500, 
               rotation.degree = rotation.degree,
               margin = margin,
               print.mode = c("raw","percent"))
}

fc_calculator <- function(columns, filters = NULL, data, n){
  
  list <- data %>% 
    filter(if (!is.null(filters)) entry_name %in% filters else entry_name %in% entry_name) %>% 
    select(entry_name, contains(columns[1]), contains(columns[2]), entry, protein_name) %>%  
    print() %>% 
    mutate(log2fc = .[[2]] - .[[3]], abs_log2fc = abs(log2fc)) %>%   # column 2 and column 3
#           log2fc = log2(fc), abs_log2fc = abs(log2fc)) %>% # already log2fc
    mutate(Change = ifelse(log2fc < 0, "Downregulated", "Upregulated")) %>%
    group_by(Change)
  
  message("Make list: pass")
  
  top_list <- list %>%
    slice_max(order_by = abs_log2fc, n = n) %>%
    arrange(desc(Change), desc(abs_log2fc)) %>%
    mutate(entry_name = forcats::fct_reorder(entry_name, log2fc)) %>% 
    relocate(c("entry", "protein_name"), .after = entry_name)

  message("Top list production: pass")
  
  half_top_list <- list %>%
    slice_max(order_by = abs_log2fc, n = floor(n/2)) %>%
    arrange(desc(Change), desc(abs_log2fc)) %>%
    mutate(entry_name = forcats::fct_reorder(entry_name, log2fc))
  
  message("Half top list production: pass")
  
  max_value <- max(half_top_list$abs_log2fc)
  p <- ggplot(half_top_list, aes(x = entry_name, y = log2fc, fill = Change)) +
    geom_col() +
    coord_flip() +
    scale_fill_discrete(breaks=c("Upregulated","Downregulated")) +
    scale_y_continuous(limits = c(-max_value, max_value)) +
    labs(x = "Entry Name", y = "Log2 Fold Change") +
    ggtitle(paste0(columns[1], " VS ", columns[2])) +
    theme(text = element_text(size= 8))

  ggsave(filename = paste0(columns[1], "_", columns[2], ".png"))      
  
  message("Plot graph: pass")
  
  write_csv(top_list, paste0(columns[1], "_", columns[2], ".csv"))
  
  message("Write top_n protein: pass")
  
}

wss <- function(data, k) {
  kmeans(data, k, nstart = 10 )$tot.withinss
}
