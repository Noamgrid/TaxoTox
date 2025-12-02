#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
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
library(ECOTOXr)
library(DT)


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("TaxoTox"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            fileInput("user_file", "Choose your data file",
                      accept = c(
                          "text/csv",
                          "text/comma-separated-values,text/plain",
                          ".csv",
                          ".xlsx",
                          ".xls",
                          ".txt",
                          ".tsv"
                      )),
            actionButton("start_processing", "1. Load Data & Find CASRN"),
            hr(),
            uiOutput("interactive_casrn_ui"),
            hr(),
            actionButton("start_toxicity_calc", "3. Calculate Toxicity"),
            hr(),
            downloadButton("download_results", "4. Download Results")
        ),

        # Show a plot of the generated distribution
        mainPanel(
           tabsetPanel(id = "main_tabs",
               tabPanel("Instructions",
                        h3("Workflow"),
                        p("1. Upload your data file using the 'Choose your data file' button."),
                        p("2. Click 'Load Data & Find CASRN' to start the process. The app will find exact matches from the Identified CASRNs database."),
                        p("3. 'Interactive CASRN Matching' tab will enabble you to approve identification of polutant names similar to the ones in the the input resolved using fuzzy matches in the dsstox database of pollutants names and CASRNs."), 
                        p("4. For any remaining compounds, use the 'Manual Entry' section in the sidebar to add CASRNs one by one."),
                        p("5. Once all CASRNs are resolved, click 'Calculate Toxicity'."),
                        p("6. View the results in the 'Toxicity Plots' and 'Toxicity Tables' tabs."),
                        p("7. Click 'Download Results' to get the final Excel file.")
                        ),
               tabPanel("CASRN Search Summary", verbatimTextOutput("cas_search_summary")),
               tabPanel("Interactive CASRN Matching", 
                        h4("Fuzzy Matches"),
                        p("Review the suggested matches below and confirm or reject them by selecting rows and clicking the button."),
                        DTOutput("fuzzy_match_table"),
                        actionButton("submit_fuzzy_matches", "Submit Fuzzy Match Selections"),
                        hr()
                        ),
               tabPanel("Toxicity Plots", 
                        plotOutput("algae_plot"),
                        plotOutput("crustacean_plot"),
                        plotOutput("fish_plot")
                        ),
               tabPanel("Toxicity Tables",
                        h4("Algae"),
                        DTOutput("tox_table_algae"),
                        h4("Crustacean"),
                        DTOutput("tox_table_crustacean"),
                        h4("Fish"),
                        DTOutput("tox_table_fish")
                        )
           )
        )
    )
)

# Define server logic
server <- function(input, output, session) {
    
    v <- reactiveValues(
        user_data = NULL,
        p_vector = NULL,
        summary_log = c("Welcome to TaxoTox! Please upload a file to begin."),
        casrn_results = NULL,
        fuzzy_to_review = NULL,
        manual_to_fill = NULL,
        final_search_results = data.table(),
        tox_results = NULL,
        plots = NULL,
        manual_additions = data.table(PREFERRED_NAME=character(), CASRN=character())
    )

    # Load static data once
    Known_CAS <- read.fst("../Data/Known_CAS.fst", as.data.table = TRUE)
    DSSTox <- read.fst("../Data/DSSTox.fst", as.data.table = TRUE)
    final_ecotox_data <- read.fst("../Data/final_ecotox_data.fst", as.data.table = TRUE)

    # Functions adapted from original script
    
    load_user_file <- function(file_path) {
        ext <- tools::file_ext(file_path)
        df <- switch(tolower(ext),
                     "csv" = read.csv(file_path, stringsAsFactors = FALSE),
                     "txt" = read.delim(file_path, stringsAsFactors = FALSE),
                     "tsv" = read.delim(file_path, sep = "\t", stringsAsFactors = FALSE),
                     "xls" = readxl::read_excel(file_path),
                     "xlsx" = readxl::read_excel(file_path),
                     stop("Unsupported file type: ", ext)
        )
        return(as.data.table(df))
    }

    fuzzy_match_non_interactive <- function(source_names, target_dt, match_col, threshold = 0.05) {
        results <- list()
        for (i in 1:length(source_names)) {
            name <- source_names[i]
            if (is.na(name)) next

            valid_targets <- target_dt[[match_col]][!is.na(target_dt[[match_col]])]
            if (length(valid_targets) == 0) next

            distances <- stringdist(name, valid_targets, method = "jw")
            if (length(distances) > 0 && !all(is.na(distances))) {
                min_dist <- min(distances, na.rm = TRUE)
                if (!is.na(min_dist) && min_dist <= threshold) {
                    best_idx <- which.min(distances)
                    best_match <- valid_targets[best_idx]
                    casrn_number <- target_dt[get(match_col) == best_match, CASRN][1]
                    
                    results[[length(results) + 1]] <- data.table(
                        source_name = name,
                        matched_name = best_match,
                        distance = min_dist,
                        CASRN = casrn_number
                    )
                }
            }
        }
        if (length(results) > 0) return(rbindlist(results))
        return(data.table())
    }

    observeEvent(input$start_processing, {
        req(input$user_file)
        
        v$summary_log <- c("Processing started...")
        
        tryCatch({
            v$user_data <- load_user_file(input$user_file$datapath)
            v$p_vector <- v$user_data[[1]]
            
            v$summary_log <- c(v$summary_log, paste("Loaded", length(v$p_vector), "compounds from user file."))
            
            # Step 2 - CASRN search
            internal_list <- Known_CAS[Known_CAS$PREFERRED_NAME %in% v$p_vector, ]
            p_vector_found <- internal_list$PREFERRED_NAME
            unfound <- v$p_vector[!v$p_vector %in% p_vector_found]
            
            v$final_search_results <- internal_list[, .(PREFERRED_NAME, CASRN)]
            
            v$summary_log <- c(v$summary_log, paste("Found", length(p_vector_found), "exact matches in Known_CAS database."))
            
            if (length(unfound) > 0) {
                v$summary_log <- c(v$summary_log, paste("Searching for", length(unfound), "compounds in DSSTox (fuzzy match)..."))
                fuzzy_matches <- fuzzy_match_non_interactive(
                    unfound,
                    DSSTox,
                    "PREFERRED_NAME",
                    threshold = 0.1
                )
                
                if (nrow(fuzzy_matches) > 0) {
                    v$fuzzy_to_review <- fuzzy_matches
                    v$summary_log <- c(v$summary_log, paste("Found", nrow(v$fuzzy_to_review), "potential fuzzy matches for review."))
                } else {
                    v$summary_log <- c(v$summary_log, "No fuzzy matches found.")
                    v$fuzzy_to_review <- NULL
                }
                
                names_in_fuzzy_review <- v$fuzzy_to_review$source_name
                v$manual_to_fill <- data.table(PREFERRED_NAME = unfound[!unfound %in% names_in_fuzzy_review])
                
            } else {
                v$summary_log <- c(v$summary_log, "All compounds found in Known_CAS database.", no fuzzy search needed.")
                v$fuzzy_to_review <- NULL
                v$manual_to_fill <- data.table(PREFERRED_NAME=character())
            }

            updateTabsetPanel(session, "main_tabs", selected = "Interactive CASRN Matching")
            
        }, error = function(e) {
            v$summary_log <- c(v$summary_log, "An error occurred:", e$message)
        })
    })

    output$cas_search_summary <- renderPrint({
        cat(v$summary_log, sep = "\n")
    })

    output$fuzzy_match_table <- renderDT({
        req(v$fuzzy_to_review)
        v$fuzzy_to_review
    }, selection = 'multiple', options = list(pageLength = 10))

    observeEvent(input$submit_fuzzy_matches, {
        req(v$fuzzy_to_review)
        
        # Assume if no rows are selected, none are confirmed.
        if (!is.null(input$fuzzy_match_table_rows_selected) && length(input$fuzzy_match_table_rows_selected) > 0) {
            confirmed_matches <- v$fuzzy_to_review[input$fuzzy_match_table_rows_selected, ]
            newly_confirmed <- confirmed_matches[, .(PREFERRED_NAME = source_name, CASRN)]
            v$final_search_results <- rbindlist(list(v$final_search_results, newly_confirmed), use.names = TRUE, fill = TRUE)
            v$summary_log <- c(v$summary_log, paste("User confirmed", nrow(newly_confirmed), "fuzzy matches."))
            
            # Update the list of compounds needing manual entry with the rejected ones
            rejected_names <- v$fuzzy_to_review[-input$fuzzy_match_table_rows_selected, .(PREFERRED_NAME = source_name)]
            v$manual_to_fill <- rbind(v$manual_to_fill, rejected_names)
        } else {
             v$summary_log <- c(v$summary_log, "No fuzzy matches were selected/confirmed.")
             # All fuzzy matches are now considered rejected and need manual entry
             v$manual_to_fill <- rbind(v$manual_to_fill, v$fuzzy_to_review[, .(PREFERRED_NAME = source_name)])
        }
        
        v$fuzzy_to_review <- NULL # Clear the review table
    })

    observeEvent(input$add_manual_casrn, {
        req(input$manual_name, input$manual_casrn)
        name <- trimws(input$manual_name)
        casrn <- trimws(input$manual_casrn)
        
        if (name != "" && casrn != "") {
            # Add directly to the final results
            v$final_search_results <- rbind(v$final_search_results, data.table(PREFERRED_NAME = name, CASRN = casrn))
            # Track for display
            v$manual_additions <- rbind(v$manual_additions, data.table(PREFERRED_NAME = name, CASRN = casrn))
            
            # Remove from the list of those to fill
            v$manual_to_fill <- v$manual_to_fill[PREFERRED_NAME != name]
            
            v$summary_log <- c(v$summary_log, paste("Manually added CASRN for", name))

            # Clear input fields for next entry
            updateTextInput(session, "manual_casrn", value = "")
        }
    })

    output$manual_added_list <- renderDT({
        v$manual_additions
    })
    
    observeEvent(input$start_toxicity_calc, {
        req(nrow(v$final_search_results) > 0)
        
        v$summary_log <- c(v$summary_log, "Starting toxicity calculations...")

        tryCatch({
            
            p_final <- as.data.table(v$final_search_results) %>%
              na.omit() %>%
              distinct(PREFERRED_NAME, .keep_all = TRUE) %>% # Ensure unique compounds
              mutate(cas_number = gsub("-", "", CASRN)) %>%
              select(PREFERRED_NAME, cas_number)
            
            cas_search_list <- as.vector(p_final$cas_number)

            endpoint_data <- final_ecotox_data %>%
              mutate(cas_number = as.character(cas_number)) %>%
              filter(cas_number %in% cas_search_list)
            
            # Algae
            endpoint_data_algae <- endpoint_data %>% 
              filter(ecotox_group == "algae") %>% 
              left_join(p_final, by = "cas_number", relationship = "many-to-many") %>% 
              group_by(PREFERRED_NAME, cas_number) %>%
              summarise(median_conc = median(conc_ng_L, na.rm = TRUE), .groups = "drop")
            
            tox_cal_algae <- endpoint_data_algae %>% 
              left_join(v$user_data, by = "PREFERRED_NAME") %>% 
              mutate(across(where(is.numeric) & !c(cas_number, median_conc), ~ .x / median_conc)) %>%
              select(-c(cas_number, median_conc)) %>% 
              pivot_longer(cols = 2:length(.), names_to = "Sample", values_to = "TU") %>% 
              pivot_wider(names_from = PREFERRED_NAME, values_from = TU) %>% 
              rowwise() %>%
              mutate(RQtg = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
              ungroup()
            
            # Crustacean
            endpoint_data_crustacean <- endpoint_data %>% 
              filter(ecotox_group == "crustacean") %>% 
              left_join(p_final, by = "cas_number", relationship = "many-to-many") %>% 
              group_by(PREFERRED_NAME, cas_number) %>%
              summarise(median_conc = median(conc_ng_L, na.rm = TRUE), .groups = "drop")
              
            tox_cal_crustacean <- endpoint_data_crustacean %>% 
              left_join(v$user_data, by = "PREFERRED_NAME") %>% 
              mutate(across(where(is.numeric) & !c(cas_number, median_conc), ~ .x / median_conc)) %>%
              select(-c(cas_number, median_conc)) %>% 
              pivot_longer(cols = 2:length(.), names_to = "Sample", values_to = "TU") %>% 
              pivot_wider(names_from = PREFERRED_NAME, values_from = TU) %>% 
              rowwise() %>%
              mutate(RQtg = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
              ungroup()

            # Fish
            endpoint_data_fish <- endpoint_data %>% 
              filter(ecotox_group == "fish") %>% 
              left_join(p_final, by = "cas_number", relationship = "many-to-many") %>% 
              group_by(PREFERRED_NAME, cas_number) %>%
              summarise(median_conc = median(conc_ng_L, na.rm = TRUE), .groups = "drop")

            tox_cal_fish <- endpoint_data_fish %>% 
              left_join(v$user_data, by = "PREFERRED_NAME") %>% 
              mutate(across(where(is.numeric) & !c(cas_number, median_conc), ~ .x / median_conc)) %>%
              select(-c(cas_number, median_conc)) %>% 
              pivot_longer(cols = 2:length(.), names_to = "Sample", values_to = "TU") %>% 
              pivot_wider(names_from = PREFERRED_NAME, values_from = TU) %>% 
              rowwise() %>%
              mutate(RQtg = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
              ungroup()
            
            v$tox_results <- list(
                "Original data" = v$user_data,
                "Algae toxicity" = tox_cal_algae,
                "Crustacean toxicity" = tox_cal_crustacean,
                "Fish toxicity" = tox_cal_fish
            )

            v$plots <- list(
                algae = ggplot(endpoint_data_algae, aes(x = reorder(PREFERRED_NAME, median_conc, FUN=median), y = median_conc)) + geom_boxplot() + theme_minimal() + labs(title="Algae Endpoint Concentrations (EC50, ng/L)", x="Compound", y="Median EC50") + coord_flip(),
                crustacean = ggplot(endpoint_data_crustacean, aes(x = reorder(PREFERRED_NAME, median_conc, FUN=median), y = median_conc)) + geom_boxplot() + theme_minimal() + labs(title="Crustacean Endpoint Concentrations (EC50, ng/L)", x="Compound", y="Median EC50") + coord_flip(),
                fish = ggplot(endpoint_data_fish, aes(x = reorder(PREFERRED_NAME, median_conc, FUN=median), y = median_conc)) + geom_boxplot() + theme_minimal() + labs(title="Fish Endpoint Concentrations (EC50, ng/L)", x="Compound", y="Median EC50") + coord_flip()
            )
            
            v$summary_log <- c(v$summary_log, "Toxicity calculations complete. See tabs for results.")
            updateTabsetPanel(session, "main_tabs", selected = "Toxicity Plots")

        }, error = function(e) {
            v$summary_log <- c(v$summary_log, "An error occurred during toxicity calculation:", e$message)
        })
    })
    
    output$tox_table_algae <- renderDT({ req(v$tox_results); v$tox_results[["Algae toxicity"]] })
    output$tox_table_crustacean <- renderDT({ req(v$tox_results); v$tox_results[["Crustacean toxicity"]] })
    output$tox_table_fish <- renderDT({ req(v$tox_results); v$tox_results[["Fish toxicity"]] })
    
    output$algae_plot <- renderPlot({ req(v$plots); v$plots$algae })
    output$crustacean_plot <- renderPlot({ req(v$plots); v$plots$crustacean })
    output$fish_plot <- renderPlot({ req(v$plots); v$plots$fish })

    output$download_results <- downloadHandler(
        filename = function() {
            paste("TaxoTox_results_", Sys.Date(), ".xlsx", sep = "")
        },
        content = function(file) {
            req(v$tox_results)
            write.xlsx(v$tox_results, file)
        }
    )
    
    output$interactive_casrn_ui <- renderUI({
      if (!is.null(v$manual_to_fill) && nrow(v$manual_to_fill) > 0) {
        tagList(
          hr(),
          h4("Manual Entry"),
          p("Add CASRNs for remaining compounds."),
          selectInput("manual_name", "Select Compound", choices = v$manual_to_fill$PREFERRED_NAME),
          textInput("manual_casrn", "Enter CASRN (e.g., 123-45-6)"),
          actionButton("add_manual_casrn", "Add CASRN"),
          h5("Manually Added:"),
          DTOutput("manual_added_list")
        )
      }
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
