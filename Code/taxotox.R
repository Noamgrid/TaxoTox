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

#functions:
#######################################################################
##load_file_as_df =  the user way to path to its data, various formats (xlsx, txt, delim)

load_file_as_df <- function() {
  # Required packages
  if (!requireNamespace("readxl", quietly = TRUE)) {
    install.packages("readxl")
  }
  
  file_path <- choose.files()
  
  # Extract file extension
  ext <- tools::file_ext(file_path)
  
  # Load based on file type
  df <- switch(tolower(ext),
               "csv" = read.csv(file_path, stringsAsFactors = FALSE),
               "txt" = read.delim(file_path, stringsAsFactors = FALSE),
               "tsv" = read.delim(file_path, sep = "\t", stringsAsFactors = FALSE),
               "xls" = readxl::read_excel(file_path),
               "xlsx" = readxl::read_excel(file_path),
               stop("Unsupported file type: ", ext)
  )
  
  return(df)
}

#######################################################################

# fuzzy_match_interactive = Interactive fuzzy matching function (for DSSTox search)
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

#######################################################################

#add_casrn_interactive = window popup massage for manually insert the unfound contaminants CASRN for further ECOTOX search

add_casrn_interactive <- function(compounds_df, casrn_column = "CASRN", compound_column = "PREFERRED_NAME") {
  result_df <- compounds_df
  
  # Find rows where CASRN is missing
  missing_casrn <- is.na(result_df[[casrn_column]]) | result_df[[casrn_column]] == ""
  missing_indices <- which(missing_casrn)
  
  if (length(missing_indices) == 0) {
    tcltk::tkmessageBox(
      title = "Complete",
      message = "All compounds already have CASRN numbers!",
      icon = "info",
      type = "ok"
    )
    return(result_df)
  }
  
  for (i in missing_indices) {
    compound_name <- result_df[i, compound_column]
    
    # Create input dialog for CASRN
    tt <- tcltk::tktoplevel()
    tcltk::tkwm.title(tt, "Enter CASRN")
    
    tcltk::tklabel(tt, text = paste("Compound:", compound_name)) %>%
      tcltk::tkpack(pady = 5)
    
    tcltk::tklabel(tt, text = "Enter CASRN number:") %>%
      tcltk::tkpack()
    
    casrn_var <- tcltk::tclVar("")
    entry_widget <- tcltk::tkentry(tt, textvariable = casrn_var, width = 30)
    tcltk::tkpack(entry_widget, pady = 5)
    
    # Buttons
    button_frame <- tcltk::tkframe(tt)
    tcltk::tkpack(button_frame, pady = 10)
    
    saved <- FALSE
    
    ok_button <- tcltk::tkbutton(button_frame, text = "Save", command = function() {
      casrn_value <- tcltk::tclvalue(casrn_var)
      if (casrn_value != "") {
        result_df[i, casrn_column] <<- casrn_value
        saved <<- TRUE
      }
      tcltk::tkdestroy(tt)
    })
    
    skip_button <- tcltk::tkbutton(button_frame, text = "Skip", command = function() {
      tcltk::tkdestroy(tt)
    })
    
    tcltk::tkpack(ok_button, side = "left", padx = 5)
    tcltk::tkpack(skip_button, side = "right", padx = 5)
    
    tcltk::tkwait.window(tt)
    
    if (saved) {
      cat(paste("CASRN added for", compound_name, "\n"))
    }
  }
  
  return(result_df)
}

#######################################################################

#Filtering and treating ECOTOX data for aquatic organisms, relevant units and saving the filtered data as csv for quicker search - RUN ONCE IN 3 MONTHS

#R TO SQL
conn <- dbConnect(RSQLite::SQLite(), "C:/Users/owner/AppData/Local/ECOTOXr/ECOTOXr/Cache/ecotox_ascii_03_13_2025.sqlite")

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


filterd_ecotox_data <-dbGetQuery(conn, mysearch) #creats a table of all endpoints for fish, algae and crustacean

relevent_units <- c("ng/L", "ug/L", "mg/L", "g/L")  

filterd_ecotox_data_conc_unit <- filterd_ecotox_data %>% #filtering for relevant concentration units (dissolved)
  mutate(
    conc1_unit = str_remove(conc1_unit, c("^AI\\s*")), # we assume that the lab measure for active ingredient and not product
    conc1_unit = str_replace_all(conc1_unit, "\\bdm3\\b", "L"),
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

irelevant_endpoints = c("NOEC","NR","LOEC","BCF","NOEL","NR-LETH","NR-ZERO","LOEL","LT50") #filtering the relevant endpoints for toxicity assessment

final_ecotox_data <- filterd_ecotox_data_conc_unit %>% 
  left_join(conversion_df, by = "conc1_unit") %>%
  mutate(min_concentration = as.numeric(min_concentration)) %>% 
  mutate(conc_ng_L = min_concentration * factor_to_ng_L) %>% 
  mutate(endpoint = str_replace_all(endpoint, "[*/]", "")) %>% 
  filter(!endpoint %in% irelevant_endpoints) %>% 
  mutate(effect = str_replace_all(effect, "[~/]", ""))


write.csv(final_ecotox_data, "../Data/final_ecotox_data.csv")

###############################################################################





# T A X O T O X




# Step 1 - Data loading and arrangement

NUKA <- read.fst("../Data/NUKA.fst") %>% # NUKA data - 2500 relevant organic pollutants in the EU and their CASRN
  rename("CASRN" = CAS)

DSSTox <- read.fst("C:/Users/owner/Desktop/DSSTox.fst") # DSSTox - database of over million pollutants, their names and IUPAC names (updating DB, needs to create a new table each time)- outside the taxotox file so full path is needed


Internal_data <- load_file_as_df() #user's data

setDT(NUKA) # Convert to data.table for more efficient operations
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
# Run the interactive matching
matches <- fuzzy_match_interactive(
  unfound,
  DSSTox, 
  "PREFERRED_NAME",
  threshold = 0.1,      # Maximum allowed distance - needs to be defined by the user?
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

#search CASRN result:
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






# Step 3 - Interactive fill the rest of pollutant's CASRN
all_p_df <- data.frame(
  PREFERRED_NAME = p_vector,
  CASRN = NA,
  stringsAsFactors = FALSE
)

# Combine and remove duplicates
all_cas <- bind_rows(final_search_results, all_p_df) %>%
  group_by(PREFERRED_NAME) %>%
  summarise(CASRN = ifelse(any(!is.na(CASRN) & CASRN != ""), 
                           first(CASRN[!is.na(CASRN) & CASRN != ""]), 
                           ""), 
            .groups = 'drop')

# Interactive CASRN addition
all_cas <- add_casrn_interactive(all_cas) %>% 
  mutate(
    CASRN = gsub("-", "", CASRN),                   
    cas_number = as.numeric(CASRN)                  
  ) %>% 
  select(-CASRN)
  
all_cas <- read_csv("../Data/manual_fill.csv") %>%  ##### TEMPORARY- manually filled data to skip the pop-up messages
  mutate(
    CASRN = gsub("-", "", CASRN),                   
    cas_number = as.numeric(CASRN)                  
  ) %>% 
  select(-CASRN)

p_final <- as.vector(all_cas$CASRN) %>%  #save as vector- without "-" (as appear in ecotox)
  gsub("-", "", .)







# Step 4 - ECOTOX endpoint value search

final_ecotox_data <- read.csv("../Data/final_ecotox_data.csv")


endpoint_data <- final_ecotox_data %>% #filtering ECOTOX data according to the user's CASRN
  filter(cas_number %in% p_final)

endpoint_data_algae <- endpoint_data %>% 
  filter(ecotox_group == "algae") %>% 
  left_join(all_cas) %>% 
  group_by(PREFERRED_NAME, cas_number) %>%  #, ecotox_group) %>% 
  summarise(
    median_conc = median(conc_ng_L, na.rm = TRUE),
    #count = sum(!is.na(conc_ng_L)),
    .groups = "drop"
  )

endpoint_data_crustacean <- endpoint_data %>% 
  filter(ecotox_group == "crustacean") %>% 
  left_join(all_cas) %>% 
  group_by(PREFERRED_NAME, cas_number) %>% #, ecotox_group) %>% 
  summarise(
    median_conc = median(conc_ng_L, na.rm = TRUE),
    #count = sum(!is.na(conc_ng_L)),
    .groups = "drop"
  )

endpoint_data_fish <- endpoint_data %>% 
  filter(ecotox_group == "fish") %>% 
  left_join(all_cas) %>% 
  group_by(PREFERRED_NAME, cas_number) %>% #, ecotox_group) %>% 
  summarise(
    median_conc = median(conc_ng_L, na.rm = TRUE),
   # count = sum(!is.na(conc_ng_L)),
    .groups = "drop"
  )






# Step 5 - Toxicity calculations:


tox_cal_algae <- endpoint_data_algae %>% 
  left_join(Internal_data, by = "PREFERRED_NAME") %>% 
  mutate(across(4:length(.), ~ .x / median_conc)) %>%  # Divide each column by conc_ng_L median value
  select(-c(cas_number, median_conc)) %>% 
  pivot_longer(cols = 2:length(.), names_to = "Sample", values_to = "TU") %>% 
  pivot_wider(names_from = PREFERRED_NAME, values_from = TU) %>% 
  rowwise() %>%
  mutate(RQtg = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
  ungroup()
 
tox_cal_crustacean <- endpoint_data_crustacean %>% 
  left_join(Internal_data, by = "PREFERRED_NAME") %>% 
  mutate(across(4:length(.), ~ .x / median_conc)) %>%  # Divide each column by conc_ng_L median value
  select(-c(cas_number, median_conc)) %>% 
  pivot_longer(cols = 2:length(.), names_to = "Sample", values_to = "TU") %>% 
  pivot_wider(names_from = PREFERRED_NAME, values_from = TU) %>% 
  rowwise() %>%
  mutate(RQtg = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
  ungroup()

tox_cal_fish <- endpoint_data_fish %>% 
  left_join(Internal_data, by = "PREFERRED_NAME") %>% 
  mutate(across(4:length(.), ~ .x / median_conc)) %>%  # Divide each column by conc_ng_L median value
  select(-c(cas_number, median_conc)) %>% 
  pivot_longer(cols = 2:length(.), names_to = "Sample", values_to = "TU") %>% 
  pivot_wider(names_from = PREFERRED_NAME, values_from = TU) %>% 
  rowwise() %>%
  mutate(RQtg = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
  ungroup()







#data visalization:

algea_box <- endpoint_data_algae %>% 
  ggplot(aes(x = PREFERRED_NAME, y = conc_ng_L)) +
  geom_boxplot()


crustacean_box <- endpoint_data_crustacean %>% 
  ggplot(aes(x = PREFERRED_NAME, y = conc_ng_L)) +
  geom_boxplot()


fish_box <- endpoint_data_fish %>% 
  ggplot(aes(x = PREFERRED_NAME, y = conc_ng_L)) +
  geom_boxplot()

##################################
