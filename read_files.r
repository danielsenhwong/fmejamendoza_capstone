# read_files.r
#
# The goal of this file will be to read in GDC data files to a single dataframe

# Set up environment
## Define packages to use, then download and load them to the library
### Set repository for packages
r_repo <- "https://cloud.r-project.org"

### We need these packages
packages <- c(
  "data.table"
)

### Install packages
install.packages(packages, repo = r_repo)

### Load all packages into library
lapply(packages,
  library,
  character.only = TRUE)

## Working Directory
cwd <- getwd()

## Define location of data files
data_root <- "gdc_data_files/"

## Build a list of folders in the data folder root to support navigation later
data_folder_list <- list.dirs(
  path = data_root,
  full.names = TRUE,
  recursive = FALSE
)
data_folder_list

# Peek inside one subfolder
data_file_list <- list.files(data_folder_list[1])
data_file_list

## Try to read in the compressed file
# Example: gdc_data_files/00326327-8da5-4226-a894-2e4e1d5acd45/2562547c-7d88-45e1-8ff4-d23faff4669d.wxs.aliquot_ensemble_masked.maf.gz
gdc_df <- fread(
  paste(
    data_folder_list[1],
    data_file_list[1],
    sep = "/"
  )

  # Things to do:
  # 1. Add a new column with the folder (gdc_sample sheet File ID column) OR file name (gdc_sample_sheet File Name column)
  # 2. Iterate (i.e., repeat, loop through) through every folder in the data_root (gdc_data_files folder)
  # 2a. When you iterate through, add the data from each new *.maf file as new rows to the existing dataframe gdc_df
  # 2b. When you pull in new data from each new file as rows in the dataframe, also add in the folder (or file) name column
)

# ===== HOMEWORK STOPS HERE =====

# Connect *.maf file file and folder names to Case ID
bridge_df <- fread(
  "gdc_data_files/gdc_sample_sheet.2024-01-12.tsv"
)

# Look up case clinical annotations by Case ID
clinical_annot_df <- fread(
  "gdc_data_files/clinical.cart.2024-01-12/clinical.tsv"
)

# Next, to support analysis:
# 1. Write code to look up the clinical annotations (e.g,. sex, age, etc.) for each Case ID based on the *.maf file folder (or file) name
# 1a. This will allow us to see who was male or female and actually support that analysis
# 2. Build a list of genes of interest (i.e., which are important in T-cell exhaustion, stem cell pluripotency, etc.)


# Example future list of genes of interest, categorized by position in relevant pathway(s)
gene_df <- data.frame(
  gene_name =       c("NOTCH1",
                      "X",
                      "Y",
                      "b",
                      "c"),
  pathway_step =    c("step 1",
                      "step 1",
                      "step 1",
                      "step 2",
                      "step 3")
)
