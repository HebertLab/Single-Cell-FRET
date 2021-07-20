## Pre-processing code for analysis of high-content microscopy images processed using Columbus and exported as text files.

## Libraries required
library(tidyverse)
library(here)

# Set working directory

setwd("/Users/jacejones-tabah/OneDrive - McGill University/Documents/SingleCell/2021-07-19 Annotated code for Github/Data")


## Function required
## x is the input list 
## d is the distance allowed between coordinate location of a single object between subsequent time points.
## delta is the allowed threshold for relative change in fluorescence intensity allowed for a single object between timepoints

dfret <- function(x, distance, delta) {
  # Select columns of interest from dataframe
  y <- x %>%
    select('PlateID', 'WellName', 'Field', 'Timepoint', 'X', 'Y', 'CellsSelected...FRET1', 'CellsSelected...CFP.Mean',
           'CellsSelected...YFP.Mean', 'CellsSelected...Cell.Roundness', 'CellsSelected...Cell.Area..µm..') %>%
    # Simplify column names
    rename('FRET1' = 'CellsSelected...FRET1',
           'CFP.Mean' = 'CellsSelected...CFP.Mean',
           'YFP.Mean' = 'CellsSelected...YFP.Mean', 
           'Cell.Roundness' = 'CellsSelected...Cell.Roundness',
           'Cell.Area.µm' = 'CellsSelected...Cell.Area..µm..')
  # Split dataframe by timepoint, now have list of dataframes, one for each time point
  y <- split(y, f = y$Timepoint) 
  
  # Loop through list of dataframes
  for (i in seq_along(y)){
    # Start with assigning a new object for the first time point
    if (i == 1){
      full <- y[[i]]
    }
    #For each time point, full_join by Field to match all objects between each time point
    else{
      filter <- full_join(full, y[[i]], by = c('Field'), suffix = c('.x', '.y')) %>%
        # Calculate the distance between the two points
        mutate(dist.y = sqrt((X.y - X.x)^2 + (Y.y - Y.x)^2),
               delta_CFP = abs(CFP.Mean.y - CFP.Mean.x) / CFP.Mean.x) %>%
        # Filter based off distance parameter and change in CFP < 30%
        filter(dist.y < distance,
               delta_CFP < delta) %>%
        # Group_by cell coordinates of preceding time point
        group_by(X.x, Y.x) %>%
        # select pairing with closest pairings
        top_n(-1, dist.y) %>%
        ungroup() %>%
        # Group_by cell coordinates of matched time point
        group_by(X.y, Y.y) %>%
        # select pairing with closest pairings, these top_n make sure that each cell is only kept once
        # in case it matches to two cells within the d parameter
        top_n(-1, dist.y) %>%
        # Calculate dfret and df_fret
        mutate(dfret.y = FRET1.y - FRET1.x,
               df.fret.y = 100*(dfret.y / FRET1.x))
      
      # Select columns from the matched time point
      matched <- filter[,str_detect(colnames(filter), ".y")]
      # Add back the coordinates from the previous time point
      matched$X <- filter$X.x
      matched$Y <- filter$Y.x
      
      # Left join by the coordinates from the previous time point
      # This will add NAs for cells that have not been matched
      full <- left_join(full, matched, by = c('X', 'Y'))
      # rename previous coordinates with the timepoint suffix
      names(full)[names(full) == 'X'] <- paste('X', i - 1, sep = '_')
      names(full)[names(full) == 'Y'] <- paste('Y', i - 1, sep = '_')
      # replace coordinates used for full_join in the next iteration
      names(full)[names(full) == 'X.y'] <- 'X'
      names(full)[names(full) == 'Y.y'] <- 'Y'
      
      # drop unnecessary columns
      full$Timepoint.y <- NULL
      full$WellName.y <- NULL
      full$PlateID.y <- NULL
      
      # Replace '.y' suffix of matched columns with the timepoint
      full <- full %>%
        rename_at(vars(ends_with('.y')), ~str_replace(., '.y', paste0("_", i))) 
    }
  }
  # The last time point matched will not have the timepoint suffix added, this adds it
  len <- length(y)
  names(full)[names(full) == 'X'] <- paste('X', len, sep = "_")
  names(full)[names(full) == 'Y'] <- paste('Y', len, sep = "_")
  
  # Add the timepoint suffix to the first time point columns
  names(full)[names(full) == 'FRET1'] <- paste('FRET1', 1, sep = "_")
  names(full)[names(full) == 'CFP.Mean'] <- paste('CFP.Mean', 1, sep = "_")
  names(full)[names(full) == 'YFP.Mean'] <- paste('YFP.Mean', 1, sep = "_")
  names(full)[names(full) == 'Cell.Roundness'] <- paste('Cell.Roundness', 1, sep = "_")
  names(full)[names(full) == 'Cell.Area.µm'] <- paste('Cell.Area.µm', 1, sep = "_")
  
  # Pivot table to tidy format 
  full <- full %>%
    rownames_to_column() %>%
    pivot_longer(cols = 6:ncol(.)) %>%
    separate(name, into = c("axis", "timepoint"), sep = "_") %>% 
    pivot_wider(names_from = axis, values_from = value) %>%
    mutate(df.fret = replace_na(df.fret, 0),
           dfret = replace_na(dfret, 0),
           dist = replace_na(dist, 0))  %>%
    drop_na(c(FRET1, CFP.Mean, YFP.Mean, Cell.Roundness, Cell.Area.µm)) 
  
  # Return the table
  return(full)
}

# List all files in directory for analysis
files <- list.files(path = here("2021-07-19 Annotated code for Github/Data"), pattern = "*].txt", full.names = T, recursive = T)

# Remove non-single cell files
akar_files <- files[grep("Population", files)]

# Read in files, skip if there is an error
akar <- map(akar_files, read.delim)

# Iterate through all files and apply the dfret function
# Alter the distance parameter
akar_fret <- map_dfr(akar, dfret, distance = 100, delta = 0.7)


# Quick plots to inspect data

# Plots the movement of matched cells in each well across time points
akar_fret %>%
  ggplot(aes(x = X, y = Y, color = rowname)) +
  geom_line(size = 1) +
  facet_wrap(~ WellName) +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        legend.position = "none")

# Plot fret responses over time for each cell
akar_fret_grouped <- akar_fret %>% group_by(WellName, rowname)

ggplot(data=akar_fret_grouped, aes(y=df.fret, x= timepoint, group=rowname))+
  geom_path()+
  facet_wrap(~ WellName) +
  scale_color_gradient(low = '#f7fefe', high = '#1d366a') +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))


# Assign additional variables (e.g. Treatment conditions, time points)

akar1 <- akar_fret %>%
  mutate(WellName = as.character(WellName),
         Well2 = substr(WellName, 2, 3),
         Treatment = ifelse(Well2 == 2, "DMSO", 
                            ifelse(Well2 == 3, "SKF",
                                   ifelse(Well2 == 4, "FSK",
                                          ifelse(Well2 == 5, "DMSO",
                                                 ifelse(Well2 == 6, "SKF",
                                                        ifelse(Well2 == 7, "FSK", 
                                                               ifelse(Well2 == 8, "DMSO",
                                                                      ifelse(Well2 == 9, "SKF",
                                                                             ifelse(Well2 == 10, "FSK",
                                                                                    ifelse(Well2 == 11, NA,NA))))))))))) 



akar1 <- akar1 %>%
  mutate(time = ifelse(timepoint == 1, 0,
                       ifelse(timepoint == 2, 2,
                              ifelse(timepoint == 3, 6,
                                     ifelse(timepoint == 4, 10, NA)))))


# Write resulting dataframe to csv for further analysis

write.csv(akar1, "/Users/jacejones-tabah/OneDrive - McGill University/Documents/SingleCell/2021-07-19 Annotated code for Github/AKAR-NLS_Sample_Data_Github.csv")
