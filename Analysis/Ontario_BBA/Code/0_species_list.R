# ------------------------------------------------
# Load/install packages and set graphical themes / working directory
# ------------------------------------------------

my_packs = c('tidyverse',
             'magrittr',
             'sf',
             'terra',
             'wildrtrax',
             'ggspatial',
             'naturecounts',
             'suntools')

if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

rm(list=ls())

# -----------------------------------------------
# Useful functions
# -----------------------------------------------

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/Analysis/Ontario_BBA/Code")

`%!in%` <- Negate(`%in%`)

# -----------------------------------------------
# WildTrax credentials
# -----------------------------------------------

wt_auth() # Need your wildtrax username

# ------------------------------------------------
# Reconcile species lists between Birds Canada and NatureCounts
# ------------------------------------------------

BSC_species <- search_species_code() %>% 
  rename(BSC_spcd = BSCDATA, species_scientific_name = scientific_name)

WT_species <- wildrtrax::wt_get_species() %>% 
  rename(WT_spcd = species_code) %>%
  dplyr::select(WT_spcd,species_common_name, species_class)




# Fix several species with incorrect english names
WT_species$species_common_name[WT_species$species_common_name == "Bald eagle"] <- "Bald Eagle"
WT_species$species_common_name[WT_species$species_common_name == "Black-crowned Night-Heron"] <- "Black-crowned Night Heron"
WT_species$species_common_name[WT_species$species_common_name == "Rock Pigeon"] <- "Rock Pigeon (Feral Pigeon)"
WT_species$species_common_name[WT_species$species_common_name == "Eastern Screech-owl"] <- "Eastern Screech-Owl"
WT_species$species_common_name[WT_species$species_common_name == "Northern Goshawk"] <- "American Goshawk"





# subset(WT_species, species_common_name == "American Pipit")
# 
# "American Pipit"
# "Arctic Tern"
# "Black-crowned Night-Heron"
# "Black Scoter"
# "Brant"
# "Chimney Swift"
# "Double-crested Cormorant"
# "Eastern Screech-owl"
# "Great Egret"
# "Northern Goshawk"
# "Rock Pigeon"


all_species <- full_join(BSC_species,WT_species, by = c("english_name" = "species_common_name")) %>% 
  relocate(BSC_spcd,WT_spcd) %>%
  subset(!is.na(species_id))

sp = "NOGO"
subset(all_species, BSC_spcd == sp)
subset(BSC_species,BSC_spcd == sp)
subset(WT_species,WT_spcd == sp)


# ------------------------------------------------
# Definitive species list
# ------------------------------------------------
saveRDS(all_species,"../Data_Cleaned/all_species.RDS")