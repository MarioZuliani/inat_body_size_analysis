# Read in all iNaturalist data

library(tidyverse)
library(httr)
library(jsonlite)
library(purrr)
library(dplyr)
library(ratelimitr)

data <- read_tsv("Data/GBIF/0000420-260623161305970/obs.csv")

View(data[1:50,])

count_of_obs <- data %>%
  group_by(recordedBy) %>%
  summarise(count=n_distinct(species))

# let's only examine users with at least 20 species observed
obs_20 <- count_of_obs %>%
  filter(count >= 20)

users <- unique(obs_20$recordedBy)

# Users lookup --------------------------------------------------------------

# give that the recordedBy column contains a mix of usernames and names, we have to ensure that the recordedBY
# users are truly unique. To do this, we will look up users given the username/name in the data to see how many
# individuals that recordedBy could be associated with. If more then one user ID is found, then we will scrap that user

# WARNING: This part of the script will take a while to run (~12-14 hours for ~40,000 users)

get_user_match <- function(username){
  
  tryCatch({
    
    response <- GET(
      "https://api.inaturalist.org/v2/users/autocomplete",
      query = list(q = username)
    )
    
    stop_for_status(response)
    
    dat <- fromJSON(
      content(response, as = "text", encoding = "UTF-8"),
      flatten = TRUE
    )
    
    if(nrow(dat$results) == 0){
      return(
        data.frame(
          queried_name = username,
          id = NA
        )
      )
    }
    
    dat$results %>%
      mutate(queried_name = username)
    
  }, error = function(e){
    
    data.frame(
      queried_name = username,
      id = NA
    )
    
  })
}

safe_get_user <- limit_rate(
  get_user_match,
  rate(n = 1, period = 1)
)

chunk_size <- 1000

user_chunks <- split(
  users,
  ceiling(seq_along(users) / chunk_size)
)

for(i in seq_along(user_chunks)){
  
  cat("Starting chunk", i, "of", length(user_chunks), "\n")
  
  user_info <- map_dfr(
    user_chunks[[i]],
    safe_get_user
  )
  
  saveRDS(
    user_info,
    paste0("user_info", i, ".rds")
  )
  
  rm(user_info)
  gc()
  
  cat("Finished chunk", i, "\n")
}



# Clean data --------------------------------------------------------------

# read in all user information 
user_table <- bind_rows(lapply(list.files("Data/user_info", full.names=TRUE), readRDS))

# now, determine how many iNaturalist IDs are associted with each recordedBy user
user_freq <- user_table %>%
  group_by(queried_name) %>%
  summarise(count=n_distinct(id))

# distinct users
user_distinct <- user_freq %>%
  filter(count == 1)

# what percentage of users are distinct
nrow(user_distinct)/nrow(user_freq) * 100
# nearly 70%

# now filter the gbif data
gbif_filtered <- data %>%
  filter(recordedBy %in% user_distinct$queried_name)

# select relevant columns
gbif_filtered <- gbif_filtered %>%
  dplyr::select(gbifID, occurrenceID, kingdom, phylum, class, order, family, genus, species, taxonRank,
                scientificName, verbatimScientificName, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters,
                eventDate, day, month, year, recordedBy) %>%
  mutate(
    verbatimScientificName = ifelse(
      str_count(verbatimScientificName, "\\S+") < 2,
      NA,
      str_extract(verbatimScientificName, "^(\\S+\\s+\\S+)")
    )
  )

# now let's see how many users we have and their observation count for birds and butterflies
gbif_filtered <- gbif_filtered %>%
  mutate(group=ifelse(class=="Aves", "bird", "butterfly"))

gbif_filtered_count <- gbif_filtered %>%
  group_by(group, recordedBy) %>%
  summarise(number_of_species = n_distinct(verbatimScientificName)) %>%
  filter(number_of_species >= 20)

summary(gbif_filtered_count[gbif_filtered_count$group=="bird",]$number_of_species)
summary(gbif_filtered_count[gbif_filtered_count$group=="butterfly",]$number_of_species)

length(unique(gbif_filtered_count[gbif_filtered_count$group=="bird",]$recordedBy))
length(unique(gbif_filtered_count[gbif_filtered_count$group=="butterfly",]$recordedBy))

# how many have over 100 species
nrow(gbif_filtered_count[gbif_filtered_count$group=="bird" & gbif_filtered_count$number_of_species>99,])
nrow(gbif_filtered_count[gbif_filtered_count$group=="butterfly" & gbif_filtered_count$number_of_species>99,])


# Add body size data ------------------------------------------------------


## Birds -------------------------------------------------------------------

# get just bird data
birds <- gbif_filtered %>%
  filter(group=="bird")

# bird body size data
bird_body_size_avonet <- read_csv("Data/body_size_data/AVONET_bird_body_size_data.csv")
bird_body_size_amniote <- read_csv("Data/body_size_data/Amniote_Database_Aug_2015.csv")

# start with AVONET
bird_body_size_avonet <- bird_body_size_avonet %>%
  dplyr::select(Species1, Mass) %>%
  rename(species=Species1, bird_mass=Mass) %>%
  mutate(source="AVONET")

# now amniote
bird_body_size_amniote <- bird_body_size_amniote %>%
  # they use -999 to indicate no date, so let's replace this with NA
  dplyr::mutate(adult_body_mass_g = ifelse(adult_body_mass_g == -999, NA, adult_body_mass_g)) %>%
  dplyr::filter(complete.cases(adult_body_mass_g)) %>%
  dplyr::mutate(species=paste(genus, species), source="Amniote") %>%
  dplyr::select(species, adult_body_mass_g, source) %>%
  dplyr::rename(bird_mass=adult_body_mass_g) 

# combine body size across the two sources
bird_body_size <- rbind(bird_body_size_avonet, bird_body_size_amniote)
bird_body_size <- bird_body_size %>% group_by(species) %>%
  dplyr::summarise(body_size=mean(bird_mass), source=toString(unique(source))) %>% 
  dplyr::mutate(group="Aves")

birds <- left_join(birds, bird_body_size, by=c("verbatimScientificName"="species"))

# what percentage of observations have body size
nrow(birds[complete.cases(birds$body_size),])/nrow(birds[complete.cases(birds$species),])*100
# 97.4%

birds_clean <- birds %>%
  select(-scientificName, -species, -kingdom, -phylum, -order, -class, -group.y, -source,
         -taxonRank, -family, -genus, -occurrenceID, -day, -month, -year,
         -decimalLatitude, -decimalLongitude, -coordinateUncertaintyInMeters, 
         -group.x, -gbifID) %>%
  rename(species_name=verbatimScientificName) %>%
  # add user ID
  left_join(user_table, by=c("recordedBy"="queried_name")) %>%
  rename(user_id=id) %>%
  select(-recordedBy)

# save the data
saveRDS(birds_clean, "Data/body_size_birds.RDS")


## Butterflies -------------------------------------------------------------

# get just butterfly data
butterflies <- gbif_filtered %>%
  filter(group=="butterfly")

# butterfly body size data
butterfly_size <- read_csv("Data/body_size_data/LepTraist_body_size_data.csv")

butterfly_size <- as.data.frame(butterfly_size) %>% 
  dplyr::group_by(Species) %>%
  dplyr::summarise(WingSpanUpper_Female=mean(WingSpanUpper_Female, na.rm=TRUE),
                   WingSpanUpper_Male=mean(WingSpanUpper_Male, na.rm=TRUE),
                   WingSpanUpper_Unspecified=mean(WingSpanUpper_Unspecified, na.rm=TRUE),
                   Genus=first(Genus)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(WingSpanUpper = rowMeans(dplyr::select(., starts_with("WingSpanUpper")), na.rm=TRUE),
                group="Lepidoptera", source="LepTraits") %>%
  dplyr::filter(complete.cases(WingSpanUpper)) %>%
  dplyr::select(Species, WingSpanUpper, group, source) %>% 
  dplyr::rename(species=Species, body_size=WingSpanUpper)

butterflies <- left_join(butterflies, butterfly_size, by=c("verbatimScientificName"="species"))

# what percentage of observations have body size
nrow(butterflies[complete.cases(butterflies$body_size),])/nrow(butterflies[complete.cases(butterflies$species),])*100
# 80.03%

# clean up final dataset
butterflies_clean <- butterflies %>%
  select(-scientificName, -species, -kingdom, -phylum, -order, -class, -group.y, -source,
         -taxonRank, -family, -genus, -occurrenceID, -day, -month, -year,
         -decimalLatitude, -decimalLongitude, -coordinateUncertaintyInMeters, 
         -group.x, -gbifID) %>%
  rename(species_name=verbatimScientificName) %>%
  # add user ID
  left_join(user_table, by=c("recordedBy"="queried_name")) %>%
  rename(user_id=id) %>%
  select(-recordedBy)

# save the data
saveRDS(butterflies_clean, "Data/body_size_butterflies.RDS")


