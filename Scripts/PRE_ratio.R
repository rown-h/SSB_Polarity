# PARAMAGNETIC RELAXATION ENHANCEMENT NMR RATIOS ===============================
# R Heggen, 14 May 2026
# Takes input of POKY list files and exports ratio of peak height differences
# between oxidised and reduced PRE-DNA substrates.

rm(list = ls())

# Packages ----
library(here)
library(tidyverse)

# INPUTS =======================================================================
# Selected protein folder
protein <- "YB-1"

# Selected PRE-DNA (5-PRE or 3-PRE)
PREDNA <- "3-PRE"

# Save the plot as a pdf?
save_ratio_graph <- TRUE

# Scale so that the greatest peak height ratio is set to 1? 
# Required for PyMOL colouring. Recommended ON, but off will not overwrite the
# "_scaled" files.
scale_ratio <- FALSE

# Set the y-axis range to [0,1]
standardise_y_axis <- FALSE

# Show lines of mean and mean - SD on graph?
plot_cutoffs <- FALSE

# Show the distance between DNA ends and each residue's backbone on the graph?
include_distance_curves <- TRUE

# Show ribbon graph for distances (distance ± SD)
graph_ribbons <- TRUE

# Show both the distance to the 5' and 3' ends of DNA
graph_both_distances <- TRUE

# Save PyMOL text file?
export_pymol_txt <- FALSE

# Save preratio csv file?
export_preratio_csv <- FALSE

# LOADING FILES ================================================================
# Import CSV with amino acid sequences and position shifts (i.e., a conversion
# number for when the residue numbering differs in the POKY list file because of
# a non-native N-terminal sequence)
seqshift <- read.csv(here('Data', 'protein_sequences_position_shifts.csv'))
positionshift <- seqshift$PositionShift[seqshift$Protein == protein]
aseq <- seqshift$Sequence[seqshift$Protein == protein]

# Finding filepaths of all .list files in the appropriate folder.
folderfilepath <- here('Data', 'List_Files', protein)

files = list.files(path = folderfilepath,
                   pattern = "\\.list$",
                   full.names = TRUE)

# Creating a function to read a .list file (skipping the blank line POKY inserts)
read_list_file <- function(file_path) {
  file_name <- basename(file_path)
  file_name <- gsub("\\.list", "", file_name)
  data <- read.table(file_path, header = FALSE, skip = 1)
}

# Importing the data
datalist <- lapply(files, read_list_file)

## Adding data names ----
#Define the new column names
new_colnames <- c("peak", "N", "H", "Data.Height")

#Use lapply to apply the column name change to each dataframe in the list
datalist <- lapply(datalist, function(df) {
  colnames(df) <- new_colnames
  return(df)
})

#Naming each dataframe within the list
extractsubstring <- function(file_path) {
  file_name <- basename(file_path)
  file_name <- gsub("\\.list$", "", file_name) # Remove ".list" extension
  return(file_name)
}

names(datalist) <- sapply(files, extractsubstring)

# Duplicate checking function
findduplicates <- function(df) {
  if ("peak" %in% colnames(df)) {
    duplicate_peaks <- df$peak[duplicated(df$peak)]
    if (length(duplicate_peaks) > 0) {
      return(paste(
        "The following peaks are duplicated:",
        duplicate_peaks,
        ". "
      ))
    } else {
      return(paste("Hooray, no duplicate peaks :-) "))
    }
  } else {
    return("The peak column missing! ")
  }
}

message(sapply(datalist, findduplicates))

# SCRIPT =======================================================================
## Modifying data ----

#Function to filter out rows of unassigned peaks or side-chain peaks
filter_data <- function(df) {
  df %>% filter(grepl("^[A-Z]\\d+N-H$", peak))
}


datalist <- lapply(datalist, filter_data)

#Function to extract the position from the peak name
extract_ogpos <- function(df) {
  df$ogpos <- as.numeric(str_extract(df$peak, "\\d+(?=N-H)"))
  return(df)
}

datalist <- lapply(datalist, extract_ogpos)


## Extracting PRE_ox AND PRE_red dataframes ----
# Separating the dataframes in datalist into separate dataframes
for (i in seq_along(datalist)) {
  dataframe_name <- names(datalist[i])
  assign(dataframe_name, datalist[[i]])
}

# Generalise name, so it works for any PRE-DNA
PRE_ox <- get(paste0(PREDNA, "_ox"))
PRE_red <- get(paste0(PREDNA, "_red"))

# Calculate the real positions
PRE_ox$pos <- PRE_ox$ogpos + positionshift
PRE_red$pos <- PRE_red$ogpos + positionshift


## Combining to yield preratio dataframe ----
# Make a new dataset which merges PRE_ox and PRE_red based on pos.
# Sort by pos and remove negatives.
preratio <- full_join(PRE_ox[, c("pos", "Data.Height")], PRE_red[, c("pos", "Data.Height")], by = "pos")
preratio <- subset(arrange(preratio, pos), pos > 0)

# Calculate ratio; if peak in PRE completely disappeared, set ratio to 0
preratio$ratio <- with(preratio, ifelse(is.na(Data.Height.x), 0, Data.Height.x / Data.Height.y))
preratio$isitassigned <- "assigned"

maximum <- max(preratio$ratio)

message(paste("The maximum ratio is", maximum))

##Scale so max ratio = 1???
if (scale_ratio == TRUE) {
  preratio$ratio <- preratio$ratio / maximum
  
  # Redefine maximum to its new value of 1, for downstream PyMOL script.
  maximum <- max(preratio$ratio)
  message(paste("The maximum ratio after scaling is", maximum))
}


## Assigned and affected residues ----
# Make a searchable list of amino acids
a <- unlist(strsplit(aseq, ""))
aa <- paste0(a, 1:nchar(aseq))
pos <- 1:nchar(aseq)

# Cutoff
cutoff <- mean(preratio$ratio)
cutoffstrict <- mean(preratio$ratio) - sd(preratio$ratio)
standarddeviation <- sd(preratio$ratio)

#Assigned residues
assignedaa <- aa[pos %in% preratio$pos]
unassignedaa <- aa[!(pos %in% preratio$pos)]

# Filter amino acids based on Ratio (assuming numeric Ratio)
broadaa <- assignedaa[preratio$ratio < cutoff]
strictaa <- assignedaa[preratio$ratio < cutoffstrict]
zeroaa <- assignedaa[preratio$ratio == 0]

# Affected levels===============================================================
assignedpos <- pos[pos %in% preratio$pos]
unassignedpos <- pos[!(pos %in% preratio$pos)]

# DISTANCE DATA ================================================================
if (include_distance_curves) {
  # Load data based on name of file
  distance <- read.csv(here('Data', 'Distance_Data', paste0(protein, '_distance.csv')))
  
  # Identify amino acid number
  distance$pos <- as.integer(str_extract(distance$Residue, '[0-9]+'))
  
  # Join to preratio dataframe, by the amino acid number
  preratio <- preratio %>%
    left_join(distance, by = "pos")
}


# GGPLOT2 BARPLOT ==============================================================
# Set graph title
title <- paste0(protein, " ", PREDNA)

# Barplot graph using ggplot2
ratio_graph <- ggplot(preratio, aes(x = pos, y = ratio))

# Optionally, add cutoff lines
if (plot_cutoffs) {
  ratio_graph <- ratio_graph +
    geom_hline(yintercept = cutoff, colour = "#4178D2") +
    geom_hline(yintercept = cutoffstrict, colour = "#9EDA4B")
}

# Make main bars, coloured by ratio value
ratio_graph <- ratio_graph +
  geom_bar(
    stat = "identity",
    aes(fill = ifelse(
      ratio < cutoffstrict,
      "< 1 S.D. below mean",
      ifelse(ratio < cutoff, "< mean", "≥ mean")
    )),
    color = "black",
    show.legend = FALSE,
    width = 1
  ) +
  scale_fill_manual(
    values = c(
      "< 1 S.D. below mean" = "#9EDA4B",
      "< mean" = "#4178D2",
      "≥ mean" = "grey50"
    ),
    name = "Signal intensity ratio"
  )

# Graph themeing and labels
ratio_graph <- ratio_graph +
  theme_minimal() +
  
  # Adding x and y axis
  theme(axis.line.x = element_line(), axis.line.y = element_line()) +
  
  labs(title = title, x = "Residue", y = "Signal intensity ratio") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.ticks = element_line(linewidth = 0.5),
    plot.title = element_text(hjust = 0.5)
  )


# Turning off labels for every second tick
ratio_graph <- ratio_graph +
  theme(axis.text.x = element_text(colour = c(NA, "black")),
        axis.text.y = element_text(colour = c("black", NA)))

# Axes
ratio_graph <- ratio_graph +
  scale_x_continuous(breaks = seq(5, ceiling(nchar(aseq) / 5) * 5, by = 5),
                     expand = c(0, 0),
  ) +
  coord_cartesian(clip = 'off', expand = F)

if(standardise_y_axis){
ratio_graph <- ratio_graph +
  scale_y_continuous(
    breaks = seq(0, 1.35, by = 0.1),
    expand = c(0, 0),
    limits = c(0, 1),
    oob = scales::squish
  )
}

# Add circles where peaks completely disappeared
ratio_graph <- ratio_graph +
  geom_point(
    data = preratio[preratio$ratio == 0, ],
    aes(x = pos, y = ratio),
    shape = 21,
    colour = "black",
    fill = "#9EDA4B",
    stroke = 1,
    size = 2
  )

# DISTANCE ON SECONDARY AXIS ===================================================
if (include_distance_curves) {

sec_axis_scaling_factor <- max(preratio$ratio) / max(preratio$mean_angstrom_5,
                                                     preratio$mean_angstrom_3,
                                                     na.rm = TRUE)

# Select which distance lines to graph
graph_5 <- graph_both_distances || PREDNA == "5-PRE"
graph_3 <- graph_both_distances || PREDNA == "3-PRE"


# Add ribbon graph
if (graph_ribbons) {
  if (graph_5) {
    ratio_graph <- ratio_graph +
      geom_ribbon(
        data = preratio,
        aes(
          x = pos,
          ymin = (mean_angstrom_5 - sd_5) * sec_axis_scaling_factor,
          ymax = (mean_angstrom_5 + sd_5) * sec_axis_scaling_factor
        ),
        fill = '#f8902b',
        alpha = 0.5
      )
  }
  
  if (graph_3) {
    ratio_graph <- ratio_graph +
      geom_ribbon(
        data = preratio,
        aes(
          x = pos,
          ymin = (mean_angstrom_3 - sd_3) * sec_axis_scaling_factor,
          ymax = (mean_angstrom_3 + sd_3) * sec_axis_scaling_factor
        ),
        fill = '#fb3b58',
        alpha = 0.5
      )
  }
}

# Add line graph

if (graph_5) {
ratio_graph <- ratio_graph +
  geom_line(aes(x = pos, y = mean_angstrom_5 * sec_axis_scaling_factor),
            colour = '#f8902b',
            linewidth = 1) }

if (graph_3) {
ratio_graph <- ratio_graph +
  geom_line(aes(x = pos, y = mean_angstrom_3 * sec_axis_scaling_factor),
            colour = '#fb3b58',
            linewidth = 1) }

ratio_graph <- ratio_graph +
  scale_y_continuous(
    name = "Signal intensity ratio",
    sec.axis = sec_axis(
      ~./sec_axis_scaling_factor,
      name = "Distance to DNA end (Å)"
    )
  )

}
ratio_graph

# GRAPH SAVING =================================================================

print(ratio_graph)



if (save_ratio_graph == TRUE) {
  graph_filepath <- here('Output',
                         'Graphs',
                         paste0(protein, "_", PREDNA, if (scale_ratio == TRUE) {
                           "_scaled"
                         }, if (scale_ratio == FALSE) {
                           "_unscaled"
                         }, ".pdf"))
  
  ggsave(
    filename = graph_filepath,
    device = "pdf",
    create.dir = TRUE,
    width = 75,
    height = 45,
    units = "mm",
    scale = 2
  )
  
  message("Graph written to ", graph_filepath)
}

#CREATING CSV FOR PYMOL SCRIPT COLOURING =======================================
if (export_pymol_txt) {
  unassignedvalueforpymol <- -1
  
  unassignedpos <- setdiff((1:nchar(aseq)), preratio$pos)
  unassigneddf <- data.frame(unassignedpos,
                             rep(unassignedvalueforpymol, length(unassignedpos)),
                             rep("unassigned", length(unassignedpos)))
  names(unassigneddf) <- c("pos", "ratio", "isitassigned")
  
  pymolscriptdfcheck <- preratio[, -c(2, 3)]
  pymolscriptdfcheck <- arrange(rbind(pymolscriptdfcheck, unassigneddf), pos)
  pymolscriptdf <- subset(pymolscriptdfcheck, select = -c(isitassigned))
  
  pymol_output <- here('Output',
                       'PyMOL',
                       'Ratio_Files',
                       paste0(protein, "_", PREDNA, if (scale_ratio == TRUE) {
                         "_scaled"
                       }, if (scale_ratio == FALSE) {
                         "_unscaled"
                       }, ".txt"))
  
  write.table(
    pymolscriptdf,
    file = pymol_output,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t"
  )
  
  message(".txt file written to ", pymol_output)
}

# WRITE INPUT FOR PYMOL ========================================================
if (export_pymol_txt) {
source(here('Scripts', 'PyMOL_Scripts', 'pymol_writer.R'))

write_pymol_colouring_script()
message(paste0('Entry for PyMOL command line saved to: PyMOL_Structures/',
               protein, '/', protein, '_', PREDNA, '_colouring.txt'))
}

# FORMULA FOR MANUAL PYMOL INPUT ===============================================
## Enter as... "select SSB and (pymolpos)"
##             "color salmon, sele"

broadposlist <- preratio$pos[preratio$ratio < cutoff]
pymolbroad <- paste0("resi ", broadposlist, collapse = " or ")

strictposlist <- preratio$pos[preratio$ratio < cutoffstrict]
pymolstrict <- paste0("resi ", strictposlist, collapse = " or ")

pymolunassigned <- paste0("resi ", unassignedpos, collapse = " or ")

# Return how affected it is
level.aa <- function(aa) {
  if (aa %in% unassignedaa)
    return("unassigned")
  if (aa %in% zeroaa)
    return("strictly affected, peak completely disappeared!")
  if (aa %in% strictaa)
    return("strictly affected")
  if (aa %in% broadaa)
    return("broadly affected")
  if (aa %in% assignedaa)
    return("unaffected")
  else
    (return("not part of protein!"))
}

aa.pos <- function(pos) {
  aa[pos]
}
level.pos <- function(pos) {
  level.aa(aa.pos(pos))
}

ogpos <- function(pos) {
  paste0(a[pos], pos - positionshift)
}
realpos <- function(pos) {
  paste0(a[pos + positionshift], pos + positionshift)
}

# SAVE PLAIN CSV OF preratio ===================================================
if (export_preratio_csv) {
  csv_output <- here('Output',
                     'CSV',
                     paste0(protein, "_", PREDNA, if (scale_ratio == TRUE) {
                       "_scaled"
                     }, if (scale_ratio == FALSE) {
                       "_unscaled"
                     }, ".csv"))
  
  write.csv(preratio, file = csv_output)
  
  message("Ratio file written to ", csv_output)
}