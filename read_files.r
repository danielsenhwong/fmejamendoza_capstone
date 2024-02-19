# read_files.r
#
# The goal of this file will be to read in GDC data files to a single dataframe

# Set up environment
## Define packages to use, then download and load them to the library
### Set repository for packages
r_repo <- "https://cloud.r-project.org"

### We need these packages
packages <- c(
  "purrr",
  "dplyr",
  "stringr",
  "data.table",
  "MASS",
  "BiocManager",
  "curl"
)

### Install packages
install.packages(packages, repo = r_repo)

### Load all packages into library
lapply(packages,
  library,
  character.only = TRUE)

### Also install a specific package for analysis of MAF files through Biocunductor
BiocManager::install("maftools")
library(maftools)

## Working Directory
cwd <- getwd()

## Define location of data files
maf_data_root <- "gdc_maf_data_files/"

# Start loading the mutation data from GDC MAF files
# Format details: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/
## Build a list of folders in the data folder root to support navigation later
maf_data_folder_list <- list.dirs(
  path = maf_data_root,
  full.names = TRUE,
  recursive = FALSE
)
maf_data_folder_list

## Use the list of folders to build a list of data files, use full names so we do not have to keep appending the folder name
maf_data_file_list <- list.files(
  maf_data_folder_list,
  pattern = "\\.maf.gz",
  full.names = TRUE
)
maf_data_file_list

## Use the path to the data files as the names so we can find the sources later if needed
## Unnecessary, MAF file and clinical annotations are joined using 'case_id' column
# names(data_file_list) <- str_remove(data_file_list, "\\.maf.gz")

## A list of all the column names, but there seems to be some variability between files
# gdc_maf_colnames <- c(
#   "Hugo_Symbol",
#   "Entrez_Gene_Id",
#   "Center",
#   "NCBI_Build",
#   "Chromosome",
#   "Start_Position",
#   "End_Position",
#   "Strand",
#   "Variant_Classification",
#   "Variant_Type",
#   "Reference_Allele",
#   "Tumor_Seq_Allele1",
#   "Tumor_Seq_Allele2",
#   "dbSNP_RS",
#   "dbSNP_Val_Status",
#   "Tumor_Sample_Barcode",
#   "Matched_Norm_Sample_Barcode",
#   "Match_Norm_Seq_Allele1",
#   "Match_Norm_Seq_Allele2",
#   "Tumor_Validation_Allele1",
#   "Tumor_Validation_Allele2",
#   "Match_Norm_Validation_Allele1",
#   "Match_Norm_Validation_Allele2",
#   "Verification_Status",
#   "Validation_Status",
#   "Mutation_Status",
#   "Sequencing_Phase",
#   "Sequence_Source",
#   "Validation_Method",
#   "Score",
#   "BAM_File",
#   "Sequencer",
#   "Tumor_Sample_UUID",
#   "Matched_Norm_Sample_UUID",
#   "HGVSc",
#   "HGVSp",
#   "HGVSp_Short",
#   "Transcript_ID",
#   "Exon_Number",
#   "t_depth",
#   "t_ref_count",
#   "t_alt_count",
#   "n_depth",
#   "n_ref_count",
#   "n_alt_count",
#   "all_effects",
#   "Allele",
#   "Gene",
#   "Feature",
#   "Feature_type",
#   "One_Consequence",
#   "Consequence",
#   "cDNA_position",
#   "CDS_position",
#   "Protein_position",
#   "Amino_acids",
#   "Codons",
#   "Existing_variation",
#   "DISTANCE",
#   "TRANSCRIPT_STRAND",
#   "SYMBOL",
#   "SYMBOL_SOURCE",
#   "HGNC_ID",
#   "BIOTYPE",
#   "CANONICAL",
#   "CCDS",
#   "ENSP",
#   "SWISSPROT",
#   "TREMBL",
#   "UNIPARC",
#   "UNIPROT_ISOFORM",
#   "RefSeq",
#   "MANE",
#   "APPRIS",
#   "FLAGS",
#   "SIFT",
#   "PolyPhen",
#   "EXON",
#   "INTRON",
#   "DOMAINS",
#   "1000G_AF",
#   "1000G_AFR_AF",
#   "1000G_AMR_AF",
#   "1000G_EAS_AF",
#   "1000G_EUR_AF",
#   "1000G_SAS_AF",
#   "ESP_AA_AF",
#   "ESP_EA_AF",
#   "gnomAD_AF",
#   "gnomAD_AFR_AF",
#   "gnomAD_AMR_AF",
#   "gnomAD_ASJ_AF",
#   "gnomAD_EAS_AF",
#   "gnomAD_FIN_AF",
#   "gnomAD_NFE_AF",
#   "gnomAD_OTH_AF",
#   "gnomAD_SAS_AF",
#   "MAX_AF",
#   "MAX_AF_POPS",
#   "gnomAD_non_cancer_AF",
#   "gnomAD_non_cancer_AFR_AF",
#   "gnomAD_non_cancer_AMI_AF",
#   "gnomAD_non_cancer_AMR_AF",
#   "gnomAD_non_cancer_ASJ_AF",
#   "gnomAD_non_cancer_EAS_AF",
#   "gnomAD_non_cancer_FIN_AF",
#   "gnomAD_non_cancer_MID_AF",
#   "gnomAD_non_cancer_NFE_AF",
#   "gnomAD_non_cancer_OTH_AF",
#   "gnomAD_non_cancer_SAS_AF",
#   "gnomAD_non_cancer_MAX_AF_adj",
#   "gnomAD_non_cancer_MAX_AF_POPS_adj",
#   "CLIN_SIG",
#   "SOMATIC",
#   "PUBMED",
#   "TRANSCRIPTION_FACTORS",
#   "MOTIF_NAME",
#   "MOTIF_POS",
#   "HIGH_INF_POS",
#   "MOTIF_SCORE_CHANGE",
#   "miRNA",
#   "IMPACT",
#   "PICK",
#   "VARIANT_CLASS",
#   "TSL",
#   "HGVS_OFFSET",
#   "PHENO",
#   "GENE_PHENO",
#   "CONTEXT",
#   "tumor_bam_uuid",
#   "normal_bam_uuid",
#   "case_id",
#   "GDC_FILTER",
#   "COSMIC",
#   "hotspot",
#   "RNA_Support",
#   "RNA_depth",
#   "RNA_ref_count",
#   "RNA_alt_count",
#   "callers"
# )

## Try to read in the compressed file
## Example: gdc_data_files/00326327-8da5-4226-a894-2e4e1d5acd45/2562547c-7d88-45e1-8ff4-d23faff4669d.wxs.aliquot_ensemble_masked.maf.gz
## The following works for one file, the goal is to also add the file name and folder name this data came from as a column
#
# gdc_df <- as.data.frame(fread(
#   paste(
#     data_folder_list[1],
#     data_file_list[1],
#     sep = "/"
#   )
# ))

## Read all of the *.maf files to a dataframe
maf_data <- map_dfr(
  .x = maf_data_file_list,
  .f = ~fread(.x, fill=TRUE, skip=7) %>% 
    mutate(
      .x,
      across(
        everything(),
        as.character)
    )
)

# Start loading the clinical data
# Clinical annotations (e.g., age, sex) for each of the cases a MAF file was retrieved
# can be downloaded as a "clinical cart" bundle of tab-separated value files (*.tsv)
# https://docs.gdc.cancer.gov/Data/Data_Model/GDC_Data_Model/
clinical_annot_file = "gdc_clinical_cart_files/clinical.tsv"
clinical_annot <- as.data.frame(
  fread(
    clinical_annot_file,
    colClasses = c(
      "days_to_last_follow_up" = "numeric",
      "vital_status" = "factor"
    ),
    na.strings = c(
      "'--"
    )
  )
)

analysis_df <- maf_data %>% inner_join(clinical_annot,
  by = c("case_id" = "case_id"),
  relationship = "many-to-many"
)

# Next, to support analysis:
# 1. Write code to look up the clinical annotations (e.g,. sex, age, etc.) for each Case ID based on the *.maf file folder (or file) name
# 1a. This will allow us to see who was male or female and actually support that analysis
# 2. Build a list of genes of interest (i.e., which are important in T-cell exhaustion, stem cell pluripotency, etc.)


# # Example future list of genes of interest, categorized by position in relevant pathway(s)
# gene_df <- data.frame(
#   gene_name =       c("NOTCH1",
#                       "X",
#                       "Y",
#                       "b",
#                       "c"),
#   pathway_step =    c("step 1",
#                       "step 1",
#                       "step 1",
#                       "step 2",
#                       "step 3")
# )

# Discriminant Analysis
lda_result <- lda(Hugo_Symbol ~ factor(gender) + as.numeric(age_at_index), data = analysis_df)
lda_result

###########
# Try this with maftools package
# https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.htm
# What TCGA projects are available
tcga_avail <- tcgaAvailable()
tcga_avail

# Grab the GBM project, it is only 400 patients from the MC3 source, fewer from the Broad
# https://github.com/PoisonAlien/TCGAmutations
tcga_gbm <- tcgaLoad(study = "GBM")
tcga_gbm

# Try making some plots
plotmafSummary(
  maf = tcga_gbm,
  rmOutlier = TRUE,
  addStat = 'median',
  dashboard = TRUE,
  titvRaw = FALSE
)

# Plot mutations with a gender annotation
oncoplot(
  maf = tcga_gbm,
  clinicalFeatures = "gender",
  sortByAnnotation = TRUE
)

# Plot with pathways
pathways = "smgbp"
oncoplot(
  maf = tcga_gbm,
  pathways = pathways,
  clinicalFeatures = "gender",
  gene_mar = 10
)
oncoplot(
  maf = tcga_gbm,
  pathways = pathways,
  clinicalFeatures = "gender",
  gene_mar = 10,
  sortByAnnotation = TRUE
)
oncoplot(
  maf = tcga_gbm,
  pathways = pathways,
  clinicalFeatures = "gender",
  gene_mar = 10,
  sortByAnnotation = TRUE,
  collapsePathway = TRUE
)

# Focus on a few genes
oncoplot(
  maf = tcga_gbm,
  genes = c("EZH2", "KDM6A", "S1PR1"),
  clinicalFeatures = "gender",
  sortByAnnotation = TRUE
)

# Try some analysis
# not sure why this doesn't work
mafSurvival(
  maf = tcga_gbm,
  genes = c("KDM6A", "EZH2", "PTEN"),
  time = "days_to_last_followup",
  Status = "vital_status",
  isTCGA = TRUE
)

# Let's use the data we pulled
# Bring the Tumor Sample Barcode column into the clinical data, and set the days_to_last_follow_up as numeric
clinical_annot_merge <- left_join(
  clinical_annot,
  maf_data[, c("case_id", "Tumor_Sample_Barcode")],
  by = "case_id"
)

# Merge all of the MAF files we downloaded ourselves, including the clinical data
gdc_gbm <- merge_mafs(
  maf_data_file_list,
  clinicalData = clinical_annot_merge
)

plotmafSummary(
  maf = gdc_gbm,
  rmOutlier = TRUE,
  addStat = 'median',
  dashboard = TRUE,
  titvRaw = FALSE
)

oncoplot(
  maf = gdc_gbm,
  clinicalFeatures = "gender",
  sortByAnnotation = TRUE
)

oncoplot(
  maf = gdc_gbm,
  pathways = pathways,
  clinicalFeatures = "gender",
  gene_mar = 10,
  sortByAnnotation = TRUE,
  collapsePathway = TRUE
)

oncoplot(
  maf = gdc_gbm,
  pathways = pathways,
  clinicalFeatures = "gender",
  gene_mar = 10,
  sortByAnnotation = TRUE
)

# Focus on a few genes
oncoplot(
  maf = gdc_gbm,
  genes = c("EZH2", "KDM6A", "S1PR1"),
  clinicalFeatures = "gender",
  sortByAnnotation = TRUE
)

# Not sure why this doesn't work
mafSurvival(
  maf = gdc_gbm,
  genes = c("KDM6A", "EZH2", "PTEN"),
  time = "days_to_last_follow_up",
  Status = "vital_status"
)

gdc_gbm.titv = titv(maf = gdc_gbm, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = gdc_gbm.titv)

lollipopPlot(
  maf = gdc_gbm,
  gene = "KDM6A",
  AACol = "HGVSp_Short",
  showMutationRate = TRUE,
  labelPos = "all"
)

lollipopPlot(
  maf = tcga_gbm,
  gene = "KDM6A",
  AACol = "HGVSp_Short",
  showMutationRate = TRUE,
  labelPos = "all"
)

rainfallPlot(
  maf = gdc_gbm,
  detectChangePoints = TRUE,
  pointSize = 0.4
)

#exclusive/co-occurence event analysis on top 10 mutated genes. 
somaticInteractions(maf = gdc_gbm, top = 25, pvalue = c(0.05, 0.1))

# Detecting cancer driver genes based on positional clustering
gdc_gbm.sig = oncodrive(maf = gdc_gbm, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = gdc_gbm.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)

#Using top 20 mutated genes to identify a set of genes (of size 2) to predict poor prognostic groups
# this doesn't work either
prog_geneset = survGroup(
  maf = gdc_gbm,
  top = 20,
  geneSetSize = 2,
  time = "days_to_last_follow_up",
  Status = "vital_status",
  verbose = FALSE
)
