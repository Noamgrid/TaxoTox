#title: "taxotox"
#author: "Yair&Noam"
#date: "2025-03-11"
#output: html_document

#Libraries

library(openxlsx)
library(tidyr)
library(dplyr)
library(ggplot2)
library(scales)
library(ggthemes)
library(tidyverse)
library(ggpattern)
library(ggpubr)
library(ggpmisc)
library(hrbrthemes)
library(remotes)
library(purrr)  
library(readxl)
library(stringr)
library(RSQLite)
library(data.table)
library(stringdist)
library(fst)
library(tcltk)
library(ECOTOXr)


# Step 1 - Data loading and arrangment
## example data - the green basine project, paths will have to become relative

GB <- read.xlsx("C:/Users/owner/OneDrive - huji.ac.il/PhD/reports/GreenBasine/CW_nabulus_raw.xlsx")

gb_p <- GB %>% 
  select(1, 4, 21:ncol(.)) %>% 
  mutate(across(3:ncol(.), ~as.numeric(.) %>% replace_na(0))) %>% 
  filter(rowSums(across(3:ncol(.))) > 0) %>% 
  pivot_longer(cols = 3:length(.), names_to = "PREFERRED_NAME", values_to = "Concentration") %>% 
  pivot_wider(names_from = c(1,2), values_from = Concentration)

# NUKA data - 2500 relevant organic pollutants in the EU and their CASRN

NUKA <- read.fst("C:/Users/owner/Desktop/TaxoTox/Data/NUKA.fst") %>% 
  rename("CASRN" = CAS)

# DSSTox - database of over million pollutants, their names and IUPAC names (updating DB, needs to create a new table each time)
DSSTox <- read.fst("C:/Users/owner/Desktop/DSSTox.fst")

Internal_data <- gb_p # The user's data

# Convert to data.table for more efficient operations
setDT(NUKA)
setDT(DSSTox)
setDT(Internal_data)





# Step 2 - CASRN search

p_vector <- Internal_data[[1]] # Create a vector of compounds from the user's data for CAS search

## First search - NUKA (exact matches)
internal_list <- NUKA[NUKA$PREFERRED_NAME %in% p_vector, ]
p_vector_found <- internal_list$PREFERRED_NAME # Creating a vector of all pollutants that were found
unfound <- p_vector[!p_vector %in% p_vector_found] # Saving the pollutants that weren't found in NUKA for further search in DSSTox

# Create a data.table for unfounded pollutants for easier manipulation
unfound_dt <- data.table(PREFERRED_NAME = unfound)


# Second search - DSSTox (fuzzy)
# Interactive fuzzy matching function 
fuzzy_match_interactive <- function(source_names, target_dt, match_col, threshold = 0.05, confirm_threshold = 0.01) {
  result <- data.table(
    source_name = source_names,
    matched_name = NA_character_,
    distance = NA_real_,
    confirmed = FALSE
  )
  
  for (i in 1:length(source_names)) {
    name <- source_names[i]
    
    if (is.na(name)) {
      next
    }
    
    valid_targets <- target_dt[[match_col]][!is.na(target_dt[[match_col]])] # Calculate distances to all target names (excluding NA targets)
    
    if (length(valid_targets) == 0) { # Handle case where there are no valid targets
      next
    }
    
    distances <- stringdist(name, valid_targets, method = "jw")
    
    if (length(distances) > 0 && !all(is.na(distances))) { # Find the best match if any valid distances exist
      min_dist <- min(distances, na.rm = TRUE)
      best_idx <- which.min(distances)
      best_match <- valid_targets[best_idx]
      match_quality <- round((1 - min_dist) * 100, 1)
      
      if (!is.na(min_dist) && min_dist <= threshold) {  # Get CASRN number for the matched compound
        casrn_number <- target_dt[get(match_col) == best_match, CASRN]
        casrn_display <- if(length(casrn_number) > 0 && !is.na(casrn_number[1])) casrn_number[1] else "Not available"
        
        if (match_quality == 100) { # Auto-accept perfect matches
          confirmed <- TRUE
        } else if (match_quality >= threshold) { # Show message box for high confidence matches (95% and above, but not 100%)
          message_text <- paste0(
            "Match found:\n\n",
            "Source compound:  ", name, "\n",
            "Matched compound: ", best_match, "\n",
            "CASRN number: ", casrn_display, "\n",
            "Match confidence: ", match_quality, "%\n\n",
            "Do you want to accept this match?"
          )
          
          answer <- tcltk::tkmessageBox( # Display message box with consistent format
            title = paste0("Compound Match (", match_quality, "% confidence)"),
            message = message_text,
            icon = "question", 
            type = "yesno"
          )
          
          confirmed <- as.character(answer) == "yes"
        } else {
          confirmed <- FALSE
        }
        
        if (confirmed) {
          result[i, `:=`(
            matched_name = best_match,
            distance = min_dist,
            confirmed = TRUE
          )]
        }
      }
    }
  }
  
  return(result)
}


# Run the interactive matching
matches <- fuzzy_match_interactive(
  unfound,
  DSSTox, 
  "PREFERRED_NAME",
  threshold = 0.1,      # Maximum allowed distance - needs to be defind by the user?
  confirm_threshold = 0  # Auto-accept 
)


# Create final result data.table with CASRN numbers
# First, get the CASRN numbers for exact matches from NUKA
exact_matches <- NUKA[PREFERRED_NAME %in% p_vector_found, .(PREFERRED_NAME, CASRN)]

# Then get CASRN numbers for fuzzy matches from DSSTox
fuzzy_confirmed <- matches[confirmed == TRUE]
fuzzy_matches <- DSSTox[DSSTox$PREFERRED_NAME %in% fuzzy_confirmed$matched_name, 
                        .(PREFERRED_NAME, CASRN)]

# Map the original names to the matched DSSTox names to get CASRN numbers
fuzzy_results <- merge(
  fuzzy_confirmed[, .(source_name, matched_name)],
  fuzzy_matches,
  by.x = "matched_name",
  by.y = "PREFERRED_NAME",
  all.x = TRUE
)

# Rename columns for consistency
setnames(fuzzy_results, "source_name", "PREFERRED_NAME")

# Combine exact and fuzzy matches
final_search_results <- rbindlist(
  list(
    exact_matches,
    fuzzy_results[, .(PREFERRED_NAME, CASRN)]
  ),
  use.names = TRUE,
  fill = TRUE
)

#search CARSN result:
percentage_found <- round((nrow(final_search_results) * 100 / length(p_vector)), 1)

#information about first stage- how many found and not found in exact search from NUKA
print(paste("Number of compounds found in NUKA:", length(p_vector_found), "out of", length(p_vector)))
# Count confirmed matches
match_count <- sum(matches$confirmed, na.rm = TRUE)
print(paste("Number of confirmed fuzzy matches (exact match from DSSTox):", match_count))
# Count unexact matches (matches accepted through message box, excluding 100% matches)
unexact_matches <- sum(matches$confirmed == TRUE & matches$distance > 0, na.rm = TRUE)
print(paste("Number of unexact matches accepted by user:", unexact_matches))
print(paste("Total pollutants with CASRN numbers:", nrow(final_search_results), "out of", length(p_vector), "(",percentage_found,"%)"))




# Step 3 - Auto fill of the rest of the pollutnat's CASRN

all_p_df <- data.frame(
  PREFERRED_NAME = p_vector,
  CASRN = NA,
  stringsAsFactors = FALSE
)

# Combine and remove duplicates, keeping original CASRN values when they exist
manual_fill_df <- bind_rows(final_search_results, all_p_df) %>%
  group_by(PREFERRED_NAME) %>%
  summarise(CASRN = ifelse(any(!is.na(CASRN) & CASRN != ""), 
                           first(CASRN[!is.na(CASRN) & CASRN != ""]), 
                           ""), 
            .groups = 'drop')

#write.csv(manual_fill_df, "manual_fill.csv") #fill
all_cas <- read_csv("../Data/manual_fill.csv") # the ".." - route the path to one directory level higher 




# Step 4 - Ecotox toxicity value search

#R TO SQL
conn <- dbConnect(RSQLite::SQLite(), "C:/Users/owner/AppData/Local/ECOTOXr/ECOTOXr/Cache/ecotox_ascii_03_13_2025.sqlite")

#for field examintaion:
#dbListFields(conn, "tests")

#A list of Ecotox organisms groups:
mysearch <- paste0("SELECT DISTINCT ecotox_group
FROM species
ORDER BY ecotox_group;")

unique_ecotox_group <- dbGetQuery(conn, mysearch)


#toxicity data search

p_final <- as.vector(all_cas$CASRN) %>%  #save as vector- without "-" (as appear in ecotox)
 gsub("-", "", .)

# Convert the vector to a comma-separated string for the SQL IN clause
cas_list <- paste0("'", p_final, "'", collapse = ", ")





#CREATING A FILTERD DATA FOR TOX EVALUETION OF AQUATIC ORGANISEMS

#combining all relevant columns from result, test and chemical df

mysearch <- paste0("WITH enriched_results AS (
    SELECT 
        r.result_id,
        t.test_id,
        t.reference_number,
        t.exposure_type,
        r.endpoint,
        r.trend,
        r.effect,
        r.conc1_unit,
        r.measurement,
        r.measurement_comments,
        r.conc1_mean,
        r.conc2_mean,
        r.conc3_mean,
        r.conc1_min,
        r.conc1_max,
        s.species,
        s.ecotox_group,
        c.chemical_name,
        -- Step 1: compute adjusted conc1_mean
        CASE 
            WHEN (r.conc1_mean IS NULL OR r.conc1_mean <= 0) 
                 AND r.conc1_min > 0 AND r.conc1_max > 0 
            THEN (r.conc1_min + r.conc1_max) / 2.0
            ELSE r.conc1_mean
        END AS calc_conc1_mean
    FROM results r
    JOIN tests t ON r.test_id = t.test_id
    JOIN species s ON t.species_number = s.species_number
    JOIN chemicals c ON t.test_cas = c.cas_number
)

SELECT 
    er.result_id,

    -- Step 2: calculate min_concentration using calc_conc1_mean
    CASE 
        WHEN er.calc_conc1_mean > 0 AND er.conc2_mean > 0 AND er.conc3_mean > 0 THEN
            CASE 
                WHEN er.calc_conc1_mean <= er.conc2_mean AND er.calc_conc1_mean <= er.conc3_mean THEN er.calc_conc1_mean
                WHEN er.conc2_mean <= er.conc3_mean THEN er.conc2_mean
                ELSE er.conc3_mean
            END
        WHEN er.calc_conc1_mean > 0 AND er.conc2_mean > 0 THEN
            CASE WHEN er.calc_conc1_mean <= er.conc2_mean THEN er.calc_conc1_mean ELSE er.conc2_mean END
        WHEN er.calc_conc1_mean > 0 AND er.conc3_mean > 0 THEN
            CASE WHEN er.calc_conc1_mean <= er.conc3_mean THEN er.calc_conc1_mean ELSE er.conc3_mean END
        WHEN er.conc2_mean > 0 AND er.conc3_mean > 0 THEN
            CASE WHEN er.conc2_mean <= er.conc3_mean THEN er.conc2_mean ELSE er.conc3_mean END
        WHEN er.calc_conc1_mean > 0 THEN er.calc_conc1_mean
        WHEN er.conc2_mean > 0 THEN er.conc2_mean
        WHEN er.conc3_mean > 0 THEN er.conc3_mean
        ELSE NULL
    END AS min_concentration,
    
    er.endpoint,
    er.trend,
    er.effect,
    er.exposure_type,
    er.measurement,
    er.measurement_comments,
    er.test_id,
    er.reference_number,
    er.conc1_unit,
    t.test_cas,
    er.species,
    
    CASE 
        WHEN LOWER(er.ecotox_group) LIKE '%fish%' THEN 'fish'
        WHEN LOWER(er.ecotox_group) LIKE '%algae%' THEN 'algae'
        WHEN LOWER(er.ecotox_group) LIKE '%crustacean%' THEN 'crustacean'
        ELSE er.ecotox_group
    END AS ecotox_group,
    
    er.chemical_name

FROM enriched_results er
JOIN tests t ON er.test_id = t.test_id

WHERE LOWER(er.ecotox_group) LIKE '%fish%' 
   OR LOWER(er.ecotox_group) LIKE '%algae%' 
   OR LOWER(er.ecotox_group) LIKE '%crustacean%';
")


filterd_ecotox_data <-dbGetQuery(conn, mysearch)


#unit analysis- what is relevent
relevent_units <- c("ng/L", "ug/L", "mg/L", "g/L")  
  
  
filterd_ecotox_data_conc_unit <- filterd_ecotox_data %>%
  mutate(
    conc1_unit = str_remove(conc1_unit, c("^AI\\s*")), # we assume that the lab mesures for active indgredient and not product
         conc1_unit = str_replace_all(conc1_unit, "\\bdm3\\b", "L"),
         conc1_unit = str_replace_all(conc1_unit, "ppt", "ng/L"), 
# This SQL will help you find publication to understand wht ppm mean
#  SELECT 
#  t.test_cas,
#  r.conc1_unit,
#  ref.title
#FROM results r
#JOIN tests t ON r.test_id = t.test_id
#JOIN [references] ref ON t.reference_number = ref.reference_number
#WHERE r.conc1_unit LIKE '%ppm%';
         conc1_unit = str_replace_all(conc1_unit, "ppm", "mg/L"),
         conc1_unit = str_replace_all(conc1_unit, "ppb", "ug/L"),
         conc1_unit = str_replace_all(conc1_unit, "0/00", "g/L")
         ) %>%
  filter(conc1_unit %in% relevent_units)



conversion_df <- tibble(
  conc1_unit = relevent_units,
  factor_to_ng_L = c(1, 1e3, 1e6, 1e9)
)


irelevant_endpoints = c("NOEC","NR","LOEC","BCF","NOEL","NR-LETH","NR-ZERO","LOEL","LT50")



final_ecotox_data <- filterd_ecotox_data_conc_unit %>% 
  left_join(conversion_df, by = "conc1_unit") %>%
  mutate(min_concentration = as.numeric(min_concentration)) %>% 
  mutate(conc_ng_L = min_concentration * factor_to_ng_L) %>% 
  mutate(endpoint = str_replace_all(endpoint, "[*/]", "")) %>% 
  filter(!endpoint %in% irelevant_endpoints) %>% 
  mutate(effect = str_replace_all(effect, "[~/]", "")) %>% 
  filter(effect %in% ef)







#filter by ecotox groups

algae <- ec50 %>% 
  filter(ecotox_group == "AlgaeStandard Test Species")

length(unique(algae$species))



crustaceans <- filterd_ecotox_data %>% 
  filter(ecotox_group == "CrustaceansStandard Test Species")

length(unique(crustaceans$chemical_name))
length(unique(crustaceans$reference_number))
length(unique(crustaceans$species))


fish <- filterd_ecotox_data %>% 
  filter(ecotox_group== "FishStandard Test Species")

length(unique(fish$chemical_name))
length(unique(fish$reference_number))
length(unique(fish$species))



  
  algae_by_count <- algae %>% 
    group_by(chemical_name) %>% 
    filter(n() > 3) %>% 
    ungroup()
    
unique(algae_by_count$conc1_unit)
    
    
    
    ggplot(aes(x=chemical_name, y=conc1_mean))+
    geom_boxplot()



    
###################################old search
    
    mysearch <- paste0("SELECT
    -- Select the relevant endpoint concentration, unit and measurment columns
    results.conc1_mean AS Endpoint_Concentration,
    results.conc1_unit AS Endpoint_Unit,
    results.endpoint AS Endpoint,
    results.measurement AS Measurement,
    results.effect AS Effect,
    tests.test_type AS TestType,
    
    -- Select details about the chemical
    chemicals.chemical_name AS ChemicalName,
    chemicals.cas_number AS CAS_Number,
    
    -- Select details about the organism
    species.common_name AS OrganismCommonName,
    species.latin_name AS OrganismLatinName,
    species.phylum_division AS OrganismPhylum,
    species.ecotox_group AS SpeciesGroup,
    
    -- Select details about the test
    tests.test_id AS TestID,
    results.result_id AS ResultID
    
FROM
    results -- Start with the results table where ED50 data is found
JOIN
    tests ON results.test_id = tests.test_id -- Join to get test details, chemical CAS, and species number
JOIN
    chemicals ON tests.test_cas = chemicals.cas_number -- Join to get chemical name from CAS number
JOIN
    species ON tests.species_number = species.species_number -- Join to get species phylum/group
WHERE
    chemicals.cas_number IN (", cas_list, ")
    AND species.ecotox_group LIKE '%Crustacea%'
    AND (results.endpoint LIKE '%EC50%' OR results.endpoint LIKE '%LC50%')
    AND (results.measurement LIKE '%IMBL%' OR '%MORT%' OR '%SURV%')
    AND (results.effect LIKE '%ITX%' OR '%MOR%')
ORDER BY
    results.conc1_mean -- Order results for easier review")
    
    output <- dbGetQuery(conn, mysearch) 
    
    output <- output %>%  
      group_by(Endpoint, Effect, Measurement) %>%
      summarise(n = n(), .groups = "drop") 
    
    unique(output$ChemicalName)
    
    
    
    
    
    
    ###  Filtering the unit that are used in test of aqueous organisms
    mysearch <- paste0("SELECT
    results.conc1_unit AS Endpoint_Unit
FROM
    results
JOIN
    tests ON results.test_id = tests.test_id
JOIN
    species ON tests.species_number = species.species_number
WHERE
    species.ecotox_group LIKE '%Crustacea%' 
    OR species.ecotox_group LIKE '%Algae%' 
    OR species.ecotox_group LIKE '%Fish%'")
    
    
    conc_unit <- dbGetQuery(conn, mysearch)
    
    conc_unit <- conc_unit %>%
      group_by(Endpoint_Unit) %>%
      summarise(n = n(), .groups = "drop") 
    
    