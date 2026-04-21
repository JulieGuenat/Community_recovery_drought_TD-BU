################################################################################
##################  Cleaning environmental metadata  ###########################
##################       Julie Morgane Guenat       ############################
#################              Mars 2025            ############################             
################################################################################

#1.Load Packages ####

library(readODS)
library(readxl)
library(tidyverse) 
library(lubridate) 

#2. set path ####

setwd("C:/Users/jguenat/OneDrive - Université de Lausanne/PhD_Analyses/")
multiparameter_file_2022 <- "Environmental_data/Multiparameter_2022_JG.ods"
multiparameter_file_2023 <- "Environmental_data/Multiparameter_2023_JG.ods"
samples_info_file <- "Cleaning_data/Samples_info_Planaqua.csv"


#3. Define month mapping between French and English ####

month_mapping <- list( 
  "Janvier" = "December",
  "Janvier_(2)" = "January",
  "Février" = "February", 
  "Mars" = "March",
  "Avril" = "April", 
  "Mai" = "May",
  "Juin" = "June",
  "Juillet" = "July",
  "Aout" = "August",
  "Septembre" = "September",
  "Octobre" = "October",
  "Novembre" = "November"
)

#4. Analysis of data structure #### 

#Function created by Claude: 
examine_sheet_structure <- function(file_path) {
  sheets <- ods_sheets(file_path)
  sheet_structures <- list()
  
  for (sheet in sheets) {
    cat("\n===== EXAMINING SHEET:", sheet, "=====\n")
    data <- read_ods(file_path, sheet = sheet)
    
    # Find rows containing "Lac"
    lac_rows <- which(apply(data, 1, function(x) any(grepl("Lac", x))))
    
    if (length(lac_rows) == 0) {
      cat("No 'Lac' row found\n")
      next
    }
    
    header_row <- lac_rows[1]
    header_values <- as.character(unlist(data[header_row, ]))
    header_values <- header_values[!is.na(header_values)]
    
    # Print headers
    cat("Headers (first 10):", paste(head(header_values, 10), collapse = ", "), "...\n")
    
    # Check for specific parameters
    temp_cols <- grep("°C", header_values, value = TRUE)
    ph_cols <- grep("pH", header_values, value = TRUE)
    do_cols <- grep("DO", header_values, value = TRUE)
    chl_cols <- grep("Chl", header_values, value = TRUE)
    
    cat("Temperature columns:", paste(temp_cols, collapse = ", "), "\n")
    cat("pH columns:", paste(ph_cols, collapse = ", "), "\n")
    cat("DO columns:", paste(do_cols, collapse = ", "), "\n")
    cat("Chlorophyll columns:", paste(chl_cols, collapse = ", "), "\n")
    
    # Save structure info
    sheet_structures[[sheet]] <- list(
      header_row = header_row,
      header_values = header_values,
      temp_cols = temp_cols,
      ph_cols = ph_cols,
      do_cols = do_cols,
      chl_cols = chl_cols
    )
  }
  
  return(sheet_structures)
}

# Examine all sheets to understand structure differences
sheet_structures_2022 <- examine_sheet_structure(multiparameter_file_2022)
sheet_structures_2023 <- examine_sheet_structure(multiparameter_file_2023)

#5. Mappings of parameters Names ####

#Function created by Claude
# Now create a more comprehensive parameter mapping based on what we found
create_detailed_parameter_mapping <- function(sheet_structures) {
  cat("\n===== CREATING PARAMETER MAPPING =====\n")
  
  # Initialize mapping
  mapping <- list(
    "temperature" = character(0),
    "do_percent" = character(0),
    "do_mg_l" = character(0),
    "ph" = character(0),
    "conductivity" = character(0),
    "turbidity" = character(0),
    "chlorophyll_rfu" = character(0),
    "chlorophyll_ug_l" = character(0)
  )
  
  # Collect all possible parameter names across sheets
  for (sheet_name in names(sheet_structures)) {
    structure <- sheet_structures[[sheet_name]]
    
    # Add temperature columns
    mapping$temperature <- unique(c(mapping$temperature, structure$temp_cols))
    
    # Add DO columns - separate percent and mg/L
    do_percent_cols <- grep("DO %", structure$do_cols, value = TRUE)
    do_percent_cols <- do_percent_cols[!grepl("DO %L", do_percent_cols)]  # Exclude DO %L
    mapping$do_percent <- unique(c(mapping$do_percent, do_percent_cols))
    
    do_mg_cols <- grep("DO mg/L", structure$do_cols, value = TRUE)
    mapping$do_mg_l <- unique(c(mapping$do_mg_l, do_mg_cols))
    
    # Add pH columns (exclude pH mV)
    ph_cols <- structure$ph_cols
    ph_cols <- ph_cols[!grepl("pH mV", ph_cols)]
    mapping$ph <- unique(c(mapping$ph, ph_cols))
    
    # Add conductivity columns
    conductivity_cols <- grep("SPC-uS/cm|SPC-uS/cm", structure$header_values, value = TRUE)
    mapping$conductivity <- unique(c(mapping$conductivity, conductivity_cols))
    
    # Add turbidity columns
    turbidity_cols <- grep("FNU", structure$header_values, value = TRUE)
    mapping$turbidity <- unique(c(mapping$turbidity, turbidity_cols))
    
    # Add chlorophyll columns
    chl_rfu_cols <- grep("Chl RFU", structure$chl_cols, value = TRUE)
    mapping$chlorophyll_rfu <- unique(c(mapping$chlorophyll_rfu, chl_rfu_cols))
    
    chl_ug_cols <- grep("Chl ug/L", structure$chl_cols, value = TRUE)
    mapping$chlorophyll_ug_l <- unique(c(mapping$chlorophyll_ug_l, chl_ug_cols))
  }
  
  # Print the mapping
  for (param in names(mapping)) {
    cat(param, ":", paste(mapping[[param]], collapse = ", "), "\n")
  }
  
  return(mapping)
}

# Create detailed parameter mapping
param_mapping_2022 <- create_detailed_parameter_mapping(sheet_structures_2022)
param_mapping_2023 <- create_detailed_parameter_mapping(sheet_structures_2023)

#6. Process the data formatting sheet by sheet ####

#This first function is used in the next function to process all sheet at once. 
#Function to process a sheet with better handling of different formats
process_sheet <- function(sheet_name, data, param_mapping) {
  cat("\nProcessing sheet:", sheet_name, "\n")
  
  # Find header row containing "Lac"
  header_rows <- which(apply(data, 1, function(x) any(grepl("Lac", x))))
  
  if (length(header_rows) == 0) {
    warning("No header row found in sheet ", sheet_name)
    return(NULL)
  }
  
  header_row <- header_rows[1]
  headers <- as.character(unlist(data[header_row, ]))
  
  # Get data rows (skip header)
  data_rows <- data[(header_row+1):nrow(data), ]
  
  # Create data frame with proper column names
  df <- as.data.frame(data_rows, stringsAsFactors = FALSE)
  colnames(df) <- headers
  
  # Remove rows with missing Lac values
  df <- df[!is.na(df$Lac), ]
  
  if (nrow(df) == 0) {
    warning("No valid data rows in sheet ", sheet_name)
    return(NULL)
  }
  
  # Add month
  df$Month <- sheet_name
  
  # Convert specific columns to numeric based on patterns in column names
  numeric_patterns <- c("°C", "DO", "mg/L", "pH", "mV", "uS/cm", "FNU", "TSS", "RFU", "ug/L")
  
  for (pattern in numeric_patterns) {
    cols <- grep(pattern, names(df), value = TRUE)
    for (col in cols) {
      if (col %in% names(df)) {
        df[[col]] <- as.numeric(as.character(df[[col]]))
      }
    }
  }
  
  return(df)
}

# Function to read and process all sheets
read_all_sheets <- function(file_path, param_mapping) {
  sheets <- ods_sheets(file_path)
  all_data <- list()
  
  for (sheet in sheets) {
    tryCatch({
      data <- read_ods(file_path, sheet = sheet)
      processed <- process_sheet(sheet, data, param_mapping)
      
      if (!is.null(processed)) {
        all_data[[sheet]] <- processed
      }
    }, error = function(e) {
      warning("Error processing sheet ", sheet, ": ", e$message)
    })
  }
  
  if (length(all_data) == 0) {
    stop("No data successfully processed from any sheet")
  }
  
  return(all_data)
}

# Read and process all sheets
monthly_data_2022 <- read_all_sheets(multiparameter_file_2022, param_mapping_2022)
monthly_data_2023 <- read_all_sheets(multiparameter_file_2023, param_mapping_2023)

#7. Mean computation for each parameters ####

# Calculate monthly averages for each lake
extract_lake_monthly_averages <- function(monthly_data, param_mapping) {
  cat("\n===== EXTRACTING MONTHLY AVERAGES FOR EACH LAKE =====\n")
  
  results <- list()
  
  for (month in names(monthly_data)) {
    cat("Processing month:", month, "\n")
    month_df <- monthly_data[[month]]
    
    # Add Lake_ID column
    month_df$Lake_ID <- paste0("L", month_df$Lac)
    
    # Find which parameters exist in this month's data
    available_params <- list()
    for (param_type in names(param_mapping)) {
      param_cols <- param_mapping[[param_type]]
      avail_cols <- param_cols[param_cols %in% names(month_df)]
      
      if (length(avail_cols) > 0) {
        available_params[[param_type]] <- avail_cols[1]  # Use first available column
      }
    }
    
    # Process each lake in this month
    for (lake in unique(month_df$Lake_ID)) {
      lake_data <- month_df[month_df$Lake_ID == lake, ]
      
      if (nrow(lake_data) == 0) next
      
      # Calculate averages for each parameter
      # Parse the real environmental measurement date, handling the three
      # inconsistent formats present across sheets:
      #   MM/DD/YYYY  --> most sheets (e.g. "04/21/2022")
      #   DD.MM.YY    --> Juillet 2022 (e.g. "20.07.22")
      #   DD.MM.YYYY  --> Janvier/Mars/Mai 2023 (e.g. "03.06.2023")
      raw_date <- as.character(lake_data$Date[1])
      
      parsed_date <- tryCatch({
        if (grepl("^\\d{2}/\\d{2}/\\d{4}$", raw_date)) {
          # Format MM/DD/YYYY
          as.Date(raw_date, format = "%m/%d/%Y")
        } else if (grepl("^\\d{2}\\.\\d{2}\\.\\d{2}$", raw_date)) {
          # Format DD.MM.YY (Juillet 2022 only)
          as.Date(raw_date, format = "%d.%m.%y")
        } else if (grepl("^\\d{2}\\.\\d{2}\\.\\d{4}$", raw_date)) {
          # Format DD.MM.YYYY (some 2023 sheets)
          as.Date(raw_date, format = "%d.%m.%Y")
        } else {
          # Fallback: let lubridate try
          lubridate::parse_date_time(raw_date, orders = c("mdy", "dmy", "ymd")) %>%
            as.Date()
        }
      }, error = function(e) {
        warning("Could not parse date '", raw_date, "' for lake ", lake,
                " month ", month, ": ", e$message)
        NA_Date_
      })
      
      result <- list(
        Month    = month,
        Lake_ID  = lake,
        month_date = parsed_date  # Real env measurement date, properly parsed
      )
      
      for (param_type in names(available_params)) {
        col_name <- available_params[[param_type]]
        values <- lake_data[[col_name]]
        
        if (!is.null(values) && length(values) > 0) {
          result[[param_type]] <- mean(values, na.rm = TRUE)
        } else {
          result[[param_type]] <- NA
        }
      }
      
      results[[paste(month, lake, sep = "_")]] <- result
    }
  }
  
  # Convert results to data frame
  result_df <- do.call(bind_rows, results)
  
  # Add English month names
  result_df$Month_English <- sapply(result_df$Month, function(m) month_mapping[[m]])
  
  # Check summary
  cat("\nSummary of monthly averages:\n")
  cat("Total rows:", nrow(result_df), "\n")
  for (param in names(param_mapping)) {
    if (param %in% names(result_df)) {
      cat(param, "- Non-NA values:", sum(!is.na(result_df[[param]])), "\n")
    } else {
      cat(param, "- Not available in results\n")
    }
  }
  
  return(result_df)
}

# Extract monthly averages
monthly_averages_2022 <- extract_lake_monthly_averages(monthly_data_2022, param_mapping_2022)
monthly_averages_2023 <- extract_lake_monthly_averages(monthly_data_2023, param_mapping_2023)

#8. Clean the processed parameter datasets ####
# Clean up the final dataset
env_clean_2022 <- monthly_averages_2022 %>%
  select(Lake_ID, Month_English, month_date, everything(), -Month) %>%
  rename(Month = Month_English)

env_clean_2023 <- monthly_averages_2023 %>%
  select(Lake_ID, Month_English, month_date, everything(), -Month) %>%
  rename(Month = Month_English)

# Handle any remaining NaN values
env_clean_2022 <- env_clean_2022 %>%
  mutate(across(where(is.numeric), ~replace(., is.nan(.), NA)))

env_clean_2023 <- env_clean_2023 %>%
  mutate(across(where(is.numeric), ~replace(., is.nan(.), NA)))

#Since we have twice the month of April, we need to label it April22 and April23
#The sample taken in April23 were take on the 3rd of May 2023, as well as the environmental data. 
env_clean_2022$Month[env_clean_2022$Month == "April"] <- "April22"
env_clean_2022$Month[env_clean_2022$Month == "March"] <- "March22"
env_clean_2023$Month[env_clean_2023$Month == "May"] <- "April23"

#Then we'll bind both datasets together
combined_env<-bind_rows(env_clean_2022, env_clean_2023)

combined_env <- combined_env %>%
  rename(date_env_data = month_date)

# Save the processed data
write.csv(combined_env, "Environmental_data/environmental_monthly_lake_data.csv", row.names = FALSE)

#9. bind with the water samples metadata #### 
# Read sample data
samples_data <- read.csv(samples_info_file, header=T, sep = ";")

#I need to change the names for the March for the year 2022 to distinguish between
#both March samples. 

sample_pattern <- "T\\d+_\\d+\\.\\d+_(blue|red|orange)_March22"
# Update the 'Collection month' column for matching Sample_IDs
samples_data$Collection_month <- ifelse(
  grepl(sample_pattern, samples_data$Samples_ID), 
  "March22",
  samples_data$Collection_month
)


# Join with samples
combined_data <- samples_data %>%
  left_join(combined_env, by = c("Lake_ID", "Collection_month" = "Month"))

# Save combined data
write.csv(combined_data, "Chapter_2_ecological_successions/data/samples_with_environmental_data.csv", row.names = FALSE)

