################################################################################
############## GRAPHIC SETTINGS FOR COHERENT DISPLAY ###########################
##################       Julie Morgane Guenat        ###########################
#################              Avril 2025             ###########################             
################################################################################

setwd("C:/Users/jguenat/OneDrive - Université de Lausanne/PhD_Analyses/")


# Nutrients 
nutrient_colors <- c("0" = "#3F72B3", "1" = "#9BA447")  

# Components, i.e. turnover and Nestedness 
component_colors <- c("Turnover" = "#E76F51", "Nestedness" = "#7B2CBF")

# labels
predator_labels <- c("1" = "Pred+", "0" = "Pred-")
nutrient_labels <- c("0" = "Nutr-", "1" = "Nutr+")


# Create a custom month vector in English
english_months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

custom_colors <- c(
  "Pred+/Nutr-" = "#1A4785",  
  "Pred-/Nutr-" = "#84a0e8",  
  "Pred+/Nutr+" = "#4E6C1E",  
  "Pred-/Nutr+" = "#adcb78")

theme_jg_bottom <- theme_minimal() + 
  theme(
    legend.position = "bottom",
    axis.title.x = element_text(size = 14, margin = margin(t = 15, b = 5)), # Add top margin
    axis.title.y = element_text(size = 14, margin = margin(r = 15, l = 5)), # Add right margin
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 0),
    plot.subtitle = element_text(size = 12, face = "italic"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    title = element_text(size = 14),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    # Additional spacing elements
    plot.margin = margin(10, 10, 10, 10),  # Add overall plot margins (top, right, bottom, left)
    plot.title = element_text(margin = margin(b = 15)),  # Add bottom margin to title
    legend.margin = margin(t = 15)  # Add top margin to legend
  )

theme_jg_right <- theme_minimal() + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size = 14, margin = margin(t = 15, b = 5)), # Add top margin
    axis.title.y = element_text(size = 14, margin = margin(r = 15, l = 5)), # Add right margin
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 0),
    plot.subtitle = element_text(size = 12, face = "italic"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    title = element_text(size = 14),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    # Additional spacing elements
    plot.margin = margin(10, 10, 10, 10),  # Add overall plot margins (top, right, bottom, left)
    plot.title = element_text(margin = margin(b = 15)),  # Add bottom margin to title
    legend.margin = margin(t = 15)  # Add top margin to legend
  )

month_labels=c("Apr22", "May22", "Jun22", "Jul22", "Aug22", "Sep22",
               "Oct22", "Nov22", "Dec22", "Jan23", "Feb23", "March23",
               "Apr23")

palette_colorful <- # Define colors for effect types
  effect_colors <- c(
    "Predator only" = "#FFC100",
    "Nutrients only" = "#86c75f",
    "Interaction only" = "#FF9800", #locked
    "Predator + Nutrients" = "#0091a8", #locked
    "Predator + Interaction" = "#D35400",
    "Nutrients + Interaction" = "#9B59B6",
    "Predator + Nutrients + Interaction" = "#c20089",# locked
    "Not highlighted" = "gray60"
  )

save.image("C:/Users/jguenat/OneDrive - Université de Lausanne/PhD_Analyses/Graph_settings.RData")

