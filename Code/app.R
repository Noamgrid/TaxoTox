# app.R
# This is a single-file Shiny application that allows a user to upload an Excel file
# and then displays the first row of the uploaded file.

# Load necessary libraries
library(shiny)
library(readxl) # For reading Excel files (.xlsx and .xls)
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
# Define the User Interface (UI) of the Shiny application
ui <- fluidPage(
  
  # Application title
  titlePanel("Testing taxotox on R-shiny"),
  
  # Sidebar layout with a file input and a main panel for output
  sidebarLayout(
    
    # Sidebar panel for file input
    sidebarPanel(
      # fileInput: Widget for users to upload a file
      # 'file1' is the input ID, which will be used to access the uploaded file in the server
      # 'Choose Excel File' is the label displayed to the user
      # 'multiple = FALSE' ensures only one file can be uploaded at a time
      # 'accept' specifies the accepted file types, making it easier for users to select correct files
      fileInput("file1", "Choose Excel File",
                multiple = FALSE,
                accept = c(".xlsx", ".xls", ".csv")) # Accepts both newer (.xlsx) older (.xls) Excel formats and csv
    ),
    
    # Main panel for displaying output
    mainPanel(
      # verbatimTextOutput: Displays preformatted text, useful for printing R objects directly
      # 'first_row_output' is the output ID, linked to the server logic
      h4("First Row of Uploaded Excel Data:"),
      verbatimTextOutput("first_row_output")
    )
  )
)

# Define the Server-side logic of the Shiny application
server <- function(input, output) {
  
  # Reactive expression to read the uploaded file
  # This code will only run when 'input$file1' changes (i.e., when a file is uploaded)
  uploaded_data <- reactive({
    
    # req(input$file1) ensures that this reactive expression only proceeds
    # if a file has actually been uploaded. If no file, it stops execution gracefully.
    req(input$file1)
    
    # Get the file path of the uploaded file
    # 'datapath' is the temporary location where Shiny stores the uploaded file
    filepath <- input$file1$datapath
    
    # Determine the file extension to choose the appropriate read function
    # Although read_excel handles both, this is good practice for more complex scenarios
    file_ext <- tools::file_ext(filepath)
    
    # Read the Excel file using read_excel from the readxl package
    # read_excel can automatically detect the sheet and format
    # For simplicity, we assume the data is in the first sheet.
    # If your Excel files have multiple sheets, you might need to add an input for sheet selection.
    
    # Use a try-catch block for robust error handling in case of an invalid Excel file
    tryCatch({
      data <- readxl::read_excel(filepath)
      return(data)
    }, error = function(e) {
      # Display an error message if reading fails
      showNotification(
        paste("Error reading Excel file:", e$message),
        type = "error",
        duration = NULL
      )
      return(NULL) # Return NULL to indicate failure
    })
  })
  #### Section 1 from Noams code

  Sys.setlocale("LC_TIME", "English")
  getwd()
  #load internal data first
  # Loading data
  # paths will have to become relative
  NUKA <- read.fst("C:/Users/owner/Desktop/TaxoTox/Data/NUKA.fst") %>% 
    rename("CASRN" = CAS)# CAS from NUKA
  DSSTox <- read.fst("C:/Users/owner/Desktop/DSSTox.fst") # CAS from DSSTox - needs to figure out how to pull
  Internal_data <- gb_p # The user's data (in this case- the GreenBasine project)
  p_vector <- Internal_data[[1]] # Create a vector of compounds from the user's data for CAS search
  
  # Convert to data.table for more efficient operations
  setDT(NUKA)
  setDT(DSSTox)
  setDT(Internal_data)
  #CASRN SEARCH
  # First search - NUKA (exact search)
  internal_list1 <- NUKA[NUKA$PREFERRED_NAME %in% p_vector, ]
  p_vector_found <- internal_list1$PREFERRED_NAME # Creating a vector of all p found
  unfound <- p_vector[!p_vector %in% p_vector_found] # Saving the p that weren't found in NUKA for further search in DSSTox
  
  # Create a data.table for unfounded pollutants for easier manipulation
  unfound_dt <- data.table(PREFERRED_NAME = unfound)
  # Second search - DSSTox (fuzzy)
  # Interactive fuzzy matching function with confirmation for uncertain matches
  fuzzy_match_interactive <- function(source_names, target_dt, match_col, threshold = 0.05, confirm_threshold = 0.01) {
    # Create a result data table
    result <- data.table(
      source_name = source_names,
      matched_name = NA_character_,
      distance = NA_real_,
      confirmed = FALSE
    )
    # Process each source name row by row
    for (i in 1:length(source_names)) {
      name <- source_names[i]
      
      # Skip NA values in source
      if (is.na(name)) {
        next
      }
      
      # Calculate distances to all target names (excluding NA targets)
      valid_targets <- target_dt[[match_col]][!is.na(target_dt[[match_col]])]
      
      # Handle case where there are no valid targets
      if (length(valid_targets) == 0) {
        next
      }
      
      distances <- stringdist(name, valid_targets, method = "jw")
      
      # Find the best match if any valid distances exist
      if (length(distances) > 0 && !all(is.na(distances))) {
        min_dist <- min(distances, na.rm = TRUE)
        best_idx <- which.min(distances)
        best_match <- valid_targets[best_idx]
        match_quality <- round((1 - min_dist) * 100, 1)
        
        if (!is.na(min_dist) && min_dist <= threshold) {
          # Get CASRN number for the matched compound
          casrn_number <- target_dt[get(match_col) == best_match, CASRN]
          casrn_display <- if(length(casrn_number) > 0 && !is.na(casrn_number[1])) casrn_number[1] else "Not available"
          
          # Auto-accept 100% matches, show message box for matches >= 95%
          if (match_quality == 100) {
            # Auto-accept perfect matches
            confirmed <- TRUE
          } else if (match_quality >= threshold) {
            # Show message box for high confidence matches (95% and above, but not 100%)
            message_text <- paste0(
              "Match found for compound matching:\n\n",
              "Source compound: ", name, "\n",
              "Matched compound: ", best_match, "\n",
              "CASRN number: ", casrn_display, "\n",
              "Match confidence: ", match_quality, "%\n\n",
              "Do you want to accept this match?"
            )
            
            # Display message box with consistent format
            answer <- tcltk::tkmessageBox(
              title = paste0("Compound Match (", match_quality, "% confidence)"),
              message = message_text,
              icon = "question", 
              type = "yesno"
            )
            
            confirmed <- as.character(answer) == "yes"
          } else {
            # For matches below 95%, don't show message box (reject automatically)
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
  
  
  
  #### End of Noams section 1
  # Render the output for 'first_row_output' in the UI
  # renderPrint is used because it formats R objects (like data frames/tibbles) nicely
  output$first_row_output <- renderPrint({
    
    # Get the data from the reactive expression
    data <- uploaded_data()
    
    # req(data) ensures that data is not NULL (i.e., file was successfully read)
    req(data)
    
    # Check if the data has any rows before trying to access the first one
    if (nrow(data) > 0) {
      # Select the first row of the data frame
      first_row <- data[1, ]
      
      # Print the first row. renderPrint will capture this output.
      print(first_row)
    } else {
      # Message if the file is empty or has no rows
      "The uploaded Excel file is empty or contains no rows of data."
    }
  })
  
}

# Run the Shiny application
shinyApp(ui = ui, server = server)