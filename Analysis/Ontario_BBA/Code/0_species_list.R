# ------------------------------------------------
# Load/install packages and set graphical themes / working directory
# ------------------------------------------------

my_packs = c('tidyverse',
             'magrittr',
             'sf',
             'terra',
             'ggspatial',
             'naturecounts',
             'suntools')

if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

rm(list=ls())

# ------------------------------------------------
# Set working directory
# ------------------------------------------------

stub <- function() {}
thisPath <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  if (length(grep("^-f$", cmdArgs)) > 0) {
    # R console option
    normalizePath(dirname(cmdArgs[grep("^-f", cmdArgs) + 1]))[1]
  } else if (length(grep("^--file=", cmdArgs)) > 0) {
    # Rscript/R console option
    scriptPath <- normalizePath(dirname(sub("^--file=", "", cmdArgs[grep("^--file=", cmdArgs)])))[1]
  } else if (Sys.getenv("RSTUDIO") == "1") {
    # RStudio
    dirname(rstudioapi::getSourceEditorContext()$path)
  } else if (is.null(attr(stub, "srcref")) == FALSE) {
    # 'source'd via R console
    dirname(normalizePath(attr(attr(stub, "srcref"), "srcfile")$filename))
  } else {
    stop("Cannot find file path")
  }
}

dirname <- thisPath()
setwd(dirname)

# -----------------------------------------------
# Useful functions
# -----------------------------------------------

`%!in%` <- Negate(`%in%`)

# ------------------------------------------------
# Reconcile species lists between Birds Canada and NatureCounts
# ------------------------------------------------

BSC_species <- search_species_code() %>% 
  rename(BSC_spcd = BSCDATA, species_scientific_name = scientific_name) %>%
  dplyr::select(BSC_spcd,species_scientific_name,english_name,french_name) %>%
  unique() %>%
  mutate(index = 1:nrow(.))

BSC_species <- BSC_species[!duplicated(BSC_species$species_scientific_name),]

WT_species <- wildRtrax::wt_get_species() %>% 
  rename(WT_spcd = species_code) %>%
  dplyr::select(WT_spcd,species_scientific_name) %>% 
  unique()

all_species <- full_join(BSC_species,WT_species) %>% relocate(BSC_spcd,WT_spcd) %>% select(-index) %>% unique()

all_species <- subset(all_species, species_scientific_name %!in% c("NULL NULL"," "))

# ------------------------------------------------
# Definitive species list
# ------------------------------------------------

saveRDS(all_species,"../Data_Cleaned/all_species.RDS")
saveRDS(BSC_species,"../Data_Cleaned/BSC_species.RDS")
