#title: "taxotox_install"
#author: "Yair&Noam"
#date: "2025-03-11"
#output: html_document



#installation - different script?

#######################################################################

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

#######################################################################

#Filtering and treating ECOTOX data for aquatic organisms, relevant units and saving the filtered data as csv for quicker search - RUN ONCE IN 3 MONTHS

#R TO SQL
conn <- dbConnect(RSQLite::SQLite(), "C:/Users/owner/AppData/Local/ECOTOXr/ECOTOXr/Cache/ecotox_ascii_03_13_2025.sqlite")

#combining all relevant columns from result, test and chemical df

mysearch <- paste0("
-- ============================================================
-- TaxoTox ECOTOX extraction query
-- Purpose: pull all aquatic toxicity results (fish, algae,
--   crustacean) with concentration values converted to a
--   single numeric field (min_concentration) ready for unit
--   conversion in R.
-- Source tables (ECOTOXr SQLite cache):
--   results   - one row per measured endpoint result
--   tests     - experimental conditions, links to species/chemical
--   species   - taxonomic group (ecotox_group)
--   chemicals - chemical name and CAS number
-- ============================================================

-- CTE: enrich each result row with species, chemical, and
-- test metadata, and compute calc_conc1_mean.
-- calc_conc1_mean: use conc1_mean when valid (>0);
--   fall back to (conc1_min + conc1_max)/2 when mean is
--   missing or zero but range bounds are available.
WITH enriched_results AS (
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
        c.cas_number,
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

    -- Step 2: min_concentration = lowest available concentration
    -- across the up to three concentration columns (conc1/2/3).
    -- Using the minimum follows the conservative approach for
    -- mixture toxicity: we prefer the most sensitive value when
    -- multiple concentrations are reported for a single result.
    -- Only positive (non-zero) values are considered.
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
    er.cas_number,
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

    -- Normalise ecotox_group labels: ECOTOX uses varied strings
    -- (e.g. 'Fish', 'Saltwater Fish'); map all to lowercase tokens
    -- so downstream R code can filter with simple equality checks.
    CASE
        WHEN LOWER(er.ecotox_group) LIKE '%fish%' THEN 'fish'
        WHEN LOWER(er.ecotox_group) LIKE '%algae%' THEN 'algae'
        WHEN LOWER(er.ecotox_group) LIKE '%crustacean%' THEN 'crustacean'
        ELSE er.ecotox_group
    END AS ecotox_group,

    er.chemical_name

FROM enriched_results er
JOIN tests t ON er.test_id = t.test_id

-- Keep only the three aquatic organism groups relevant to TaxoTox.
-- All other ecotox_groups (birds, mammals, etc.) are excluded.
WHERE LOWER(er.ecotox_group) LIKE '%fish%'
   OR LOWER(er.ecotox_group) LIKE '%algae%'
   OR LOWER(er.ecotox_group) LIKE '%crustacean%';
")


filterd_ecotox_data <-dbGetQuery(conn, mysearch) #creats a table of all endpoints for fish, algae and crustacean

relevent_units <- c("ng/L", "ug/L", "mg/L", "g/L")  

filterd_ecotox_data_conc_unit <- filterd_ecotox_data %>% #filtering for relevant concentration units (dissolved)
  mutate(
    conc1_unit = str_remove(conc1_unit, c("^AI\\s*")), # we assume that the lab measure for active ingredient and not product
    conc1_unit = str_replace_all(conc1_unit, "\\bdm3\\b", "L"), #decimeter^3, equivalent to L
    conc1_unit = str_replace_all(conc1_unit, "ppt", "ng/L"), 
    conc1_unit = str_replace_all(conc1_unit, "ppm", "mg/L"),
    conc1_unit = str_replace_all(conc1_unit, "ppb", "ug/L"),
    conc1_unit = str_replace_all(conc1_unit, "0/00", "g/L")
  ) %>%
  filter(conc1_unit %in% relevent_units)

conversion_df <- tibble(
  conc1_unit = relevent_units,
  factor_to_ng_L = c(1, 1e3, 1e6, 1e9)
)


endpoint_count <- filterd_ecotox_data_conc_unit %>% 
  group_by(endpoint) %>% 
  count(endpoint)

final_ecotox_data <- filterd_ecotox_data_conc_unit %>%
  left_join(conversion_df, by = "conc1_unit") %>%
  mutate(min_concentration = as.numeric(min_concentration)) %>%
  # Back-transform (log)-reported values BEFORE unit conversion.
  # "(log)EC50" / "(log)LC50" mean the concentration is stored as log10(conc).
  mutate(min_concentration = if_else(
    str_detect(endpoint, "^\\(log\\)"),
    10 ^ min_concentration,
    min_concentration
  )) %>%
  mutate(conc_ng_L = min_concentration * factor_to_ng_L) %>%
  # Normalise endpoint labels:
  #   (log)EC50 -> EC50,  (log)LC50 -> LC50,  LC50* -> LC50
  mutate(endpoint = str_remove(endpoint, "^\\(log\\)"),
         endpoint = str_replace_all(endpoint, "[*/]", "")) %>%
  # Fix C: keep only true acute benchmarks
  filter(endpoint %in% c("LC50", "EC50")) %>%
  mutate(effect = str_replace_all(effect, "[~/]", ""))


taxotox_data <- final_ecotox_data %>% 
  group_by(cas_number, ecotox_group) %>% 
  summarize(median_min_conc = median(min_concentration, na.rm = TRUE),
            observation_count = n()) 



write_fst(final_ecotox_data, "../Data/final_ecotox_data.fst", compress = 50)

#filter dsstox data to casrn that appear in ecotox

cas_ecotox <- final_ecotox_data$cas_number
cas_ecotox <- as.character(cas_ecotox)

#######################################################################

#combining the DSSTox df into one and filtering the relevent CAS according to ecotox- also include in the updates

# Set the folder path containing your Excel files
folder_path <- "../../Data/DSSTox_Feb_2024" #should be on the desktop of the virtual computer


# Get all Excel file names in the folder
excel_files <- list.files(path = folder_path, 
                          pattern = "\\.(xlsx|xls)$", 
                          full.names = TRUE)

# Read all Excel files and store them in a list
all_data <- map(excel_files, ~{
  file_name <- basename(.x)
  data <- read_excel(.x)
  # Add a column to identify the source file if needed  - remove?
  # data$source_file <- file_name - remove?
  return(data)
})

# Name each dataframe in the list according to the file name
names(all_data) <- basename(excel_files)

# combine all data into a single dataframe and add a CAS column
combined_df <- bind_rows(all_data) %>% 
  select(-c(1,4,6:length(.))) 

#add IUPAC names to the name column
DSSTox <- rbind(
  combined_df[, c("PREFERRED_NAME", "CASRN")],              
  data.frame(PREFERRED_NAME = combined_df$IUPAC_NAME,       
             CASRN = combined_df$CASRN)) %>% 
  mutate("cas" = gsub("-", "", CASRN)) 

cas_dsstox <- DSSTox$cas

common <- intersect(cas_dsstox, cas_ecotox)

DSSTox <- DSSTox %>% 
  filter(cas %in% common)

write.fst(DSSTox, "../Data/DSSTox.fst") #save filtered and updated fst file with names and CAS of pollutants from DSSTox database


#######################################################################


