# get user's start date

# Check and install required packages
required_packages <- c("httr", "jsonlite", "dplyr")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing missing package: ", pkg)
    install.packages(pkg)
  }
}

# Load relevant packages
library(httr)
library(jsonlite)
library(dplyr)

# Get the base URL
base_url <- "https://api.inaturalist.org/v2/observations/observers"

# read in body size data 
birds <- readRDS("Data/body_size_birds.RDS")
butterflies <- readRDS("Data/body_size_butterflies.RDS")

users <- unique(c(birds$user_id, butterflies$user_id))

# Before creating a customized request, the v2 API requires that users specify desired data fields
# The code below will provide an example of all data available from iNaturalist along with field names
fields <- fromJSON(content(GET("https://api.inaturalist.org/v2/observations/observers?fields=all"), as = "text", encoding = "UTF-8"))
colnames(fields$results)

# Use this code to customize the request
# See https://api.inaturalist.org/v2/docs/#!/Observations/get_observations_observers for all possible parameters
params <- list(
  
  # Research Grade observations only
  quality_grade = "research",
  
  # Specify desired data fields (see lines 24-27 above for all available fields)
  fields = paste("observation_count", "species_count", "user.id", "user.login", 
                 "user.created_at", "user.name", "user.observations_count",
                 "user.identifications_count", sep=","),
  
  # Maximum results for this endpoint
  per_page = 500                         
)

response <- GET(base_url, query = params)
stop_for_status(response)

data_parsed <- fromJSON(content(response, as = "text", encoding = "UTF-8"), flatten = TRUE)

res <- data_parsed$results
output_data <- if (is.null(res)) tibble() else as_tibble(res)


# Examine data
head(output_data)
cat("Total rows fetched:", nrow(output_data), "\n")







library(httr)
library(jsonlite)
library(dplyr)
library(purrr)

# Split users into groups of 50
user_chunks <- split(users, ceiling(seq_along(users) / 50))

all_results <- vector("list", length(user_chunks))

for(i in seq_along(user_chunks)) {
  
  params <- list(
    
    # Request only these users
    user_id = paste(user_chunks[[i]], collapse = ","),
    
    quality_grade = "research",
    
    fields = paste(
      "observation_count",
      "species_count",
      "user.id",
      "user.login",
      "user.created_at",
      "user.name",
      "user.observations_count",
      "user.identifications_count",
      sep = ","
    ),
    
    per_page = 500
  )
  
  response <- GET(base_url, query = params)
  stop_for_status(response)
  
  data_parsed <- fromJSON(
    content(response, as = "text", encoding = "UTF-8"),
    flatten = TRUE
  )
  
  res <- data_parsed$results
  
  all_results[[i]] <-
    if (is.null(res)) tibble() else as_tibble(res)
  
  message(sprintf("Finished %d of %d", i, length(user_chunks)))
  
  Sys.sleep(1)   # be polite to the API
}

output_data <- bind_rows(all_results)

# save the data
saveRDS(output_data, "Data/inat_user_info.RDS")
