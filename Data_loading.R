library(dplyr)
library(janitor)
library(readxl)

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

df <- read_xlsx("Source/Proteomics2.xlsx", sheet = "GR1_pool SEC")

df <- clean_names(df)

df %>% View
# Clean data -----------------------------------------------------

df_cleaned <- df %>% 
  mutate(organism = tolower(organism)) %>% 
  mutate(protein_names = stringr::str_extract(protein_names, "\\(.*\\)")) %>% 
  mutate(protein_names = stringr::str_replace_all(protein_names, "[\\(\\)]", "")) %>% 
  rename("abbr" = "protein_names")
  
dim(df_cleaned) #3206, 23
unique(df_cleaned$organism) # A lot

pattern <- c("sapiens", "papilloma")

df_HHPV <-df_cleaned %>% 
  filter(grepl(paste(pattern, collapse = "|"), organism)) %>% 
  filter(status == "reviewed")

glimpse(df_HHPV)  
# Ultracentrifugation -----------------------------------------------------

UC_df_HHPV_real <- df_HHPV %>% select(-"x1", -starts_with("max"), -matches("sec[0-9]+")) %>% 
  rename_with(~paste0("x1_ev", 1), .cols = starts_with("x1")) %>% 
  rename_with(~paste0("x2_ev", 1:3), .cols = starts_with("x2")) %>% 
  rename_with(~paste0("x3_ev", 1:3), .cols = starts_with("x3")) 

glimpse(UC_df_HHPV_real)  

## Cheating
UC_df_HHPV <- UC_df %>% 
  mutate(x1_ev2 = x1_ev1, x1_ev3 = x1_ev1) %>% 
  relocate(c("x1_ev2", "x1_ev3"), .after = x1_ev1)

glimpse(UC_df_HHPV)  

# SEC ---------------------------------------------------------------------

glimpse(df)

SEC_df_HHPV <- df_HHPV 
  # select(-"x1", -starts_with("max"), -matches("_ev[0-9]+"))
SEC_df_HHPV  
