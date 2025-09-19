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

count(maf_data[maf_data$Chromosome == "chrX"])
unique(maf_data$Chromosome)
unique(maf_data$Hugo_Symbol[maf_data$Chromosome == "chrX"])

# Notes for Flor's analysis
# Fisher's test? Chi sq n must be >5
# What is the frequency of mutation for these genes? The concern is the dataset is not large enough and insufficiently powered
# Can Flor just do analysis on the X chromosome genes? How do we make a heatmap with just these
# Is there a difference in the frequency of disease between men vs. women
# What is the formula to determine sample size required for adequate power?
## ? https://csg.sph.umich.edu/abecasis/cats/gas_power_calculator/
## ? https://www.ejeph.com/download/the-association-between-snps-and-a-quantitative-trait-power-calculation-3925.pdf
## ? https://www.semanticscholar.org/paper/Sample-Size-and-Statistical-Power-Calculation-in-Hong-Park/e476b89ad6acc97b06307549a58b9b81034b4c2b

analysis_kdm6a <- (analysis_df[analysis_df$Hugo_Symbol == "KDM6A"])
analysis_kdm6a_male <- analysis_kdm6a[analysis_kdm6a$gender == "male"]
analysis_kdm6a_female <- analysis_kdm6a[analysis_kdm6a$gender == "female"]

# Next, to support analysis:
# 1. Write code to look up the clinical annotations (e.g,. sex, age, etc.) for each Case ID based on the *.maf file folder (or file) name
# 1a. This will allow us to see who was male or female and actually support that analysis
# 2. Build a list of genes of interest (i.e., which are important in T-cell exhaustion, stem cell pluripotency, etc.)

# all mutated genes, count of unqiue patients
# 14,465 genes total
analysis_summary_all <- analysis_df %>%
  group_by(Hugo_Symbol) %>%
  summarize(
    n_pts = n_distinct(case_id),
    n_tumor_samples = n_distinct(Tumor_Sample_Barcode)
  )

# mutated genes on X chromosome, unique patients
# 628 genes total
analysis_summary_chrX <- analysis_df[analysis_df$Chromosome == "chrX"] %>%
  group_by(Hugo_Symbol) %>%
  summarize(
    n_pts = n_distinct(case_id),
    n_tumor_samples = n_distinct(Tumor_Sample_Barcode)
  )

# mutated genes in males, unique patients
# 8,045 total
analysis_summary_male <- analysis_df[analysis_df$gender == "male"] %>%
  group_by(Hugo_Symbol) %>%
  summarize(
    n_pts = n_distinct(case_id),
    n_tumor_samples = n_distinct(Tumor_Sample_Barcode)
  )

# mutated genes in females, unique patients
# 12,899 total
analysis_summary_female <- analysis_df[analysis_df$gender == "female"] %>%
  group_by(Hugo_Symbol) %>%
  summarize(
    n_pts = n_distinct(case_id),
    n_tumor_samples = n_distinct(Tumor_Sample_Barcode)
  )

# mutated genes on X chromosome, by gender, unique patients
# 859 genes total
analysis_summary_chrX_gender_group <- analysis_df[analysis_df$Chromosome == "chrX"] %>%
  group_by(Hugo_Symbol, gender) %>%
  summarize(
    n_pts = n_distinct(case_id),
    n_tumor_samples = n_distinct(Tumor_Sample_Barcode)
  )

# mutated genes on X chromosome in males, unique patients
# 261 total
analysis_summary_chrX_male <- analysis_df[
  (analysis_df$Chromosome == "chrX") &
    (analysis_df$gender == "male")
] %>%
  group_by(Hugo_Symbol) %>%
  summarize(
    n_pts = n_distinct(case_id),
    n_tumor_samples = n_distinct(Tumor_Sample_Barcode)
  )

# mutated genes on X chromosome in females, unique patients
# 598 total
analysis_summary_chrX_female <- analysis_df[
  (analysis_df$Chromosome == "chrX") &
    (analysis_df$gender == "female")
] %>%
  group_by(Hugo_Symbol) %>%
  summarize(
    n_pts = n_distinct(case_id),
    n_tumor_samples = n_distinct(Tumor_Sample_Barcode)
  )

# all mutated genes, count of unqiue patients
# only 4 patients? this does not agree with above, seems like women are missing
# -- no, this is right.
# there are only four unique case ids, there are 6 unique tumor ids
analysis_summary_kdm6a <- analysis_df[(analysis_df$Hugo_Symbol == "KDM6A")] %>%
  group_by(Hugo_Symbol, gender) %>%
  summarize(
    n_pts = n_distinct(case_id),
    n_tumor_samples = n_distinct(Tumor_Sample_Barcode)
  )

write.csv(analysis_summary_all,
          "export_csv\\summary_all_mutated_tcga_gbm_genes.csv")
write.csv(analysis_summary_chrX,
          "export_csv\\summary_chrX_mutated_tcga_gbm_genes.csv")
write.csv(analysis_summary_male,
          "export_csv\\summary_male_mutated_tcga_gbm_genes.csv")
write.csv(analysis_summary_female,
          "export_csv\\summary_female_mutated_tcga_gbm_genes.csv")
write.csv(analysis_summary_chrX_gender_group,
          "export_csv\\summary_chrX_mutated_tcga_gbm_genes_by_gender.csv")
write.csv(analysis_summary_chrX_male,
          "export_csv\\summary_chrX_mutated_tcga_gbm_genes_male.csv")
write.csv(analysis_summary_chrX_female,
          "export_csv\\summary_chrX_mutated_tcga_gbm_genes_female.csv")
write.csv(analysis_summary_kdm6a,
          "export_csv\\summary_mutated_tcga_gbm_kdm6a.csv")

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

# # Discriminant Analysis
# lda_result <- lda(Hugo_Symbol ~ factor(gender) + as.numeric(age_at_index), data = analysis_df)
# lda_result

###########
# Try this with maftools package
# https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.htm
# # What TCGA projects are available
# tcga_avail <- tcgaAvailable()
# tcga_avail

# # Grab the GBM project, it is only 400 patients from the MC3 source, fewer from the Broad
# # https://github.com/PoisonAlien/TCGAmutations
# tcga_gbm <- tcgaLoad(study = "GBM")
# tcga_gbm

# # Try making some plots
# plotmafSummary(
#   maf = tcga_gbm,
#   rmOutlier = TRUE,
#   addStat = 'median',
#   dashboard = TRUE,
#   titvRaw = FALSE
# )

# oncoplot(
#   maf = tcga_gbm
# )

# # Plot mutations with a gender annotation
# oncoplot(
#   maf = tcga_gbm,
#   clinicalFeatures = "gender",
#   sortByAnnotation = TRUE
# )

# # Plot with pathways
# pathways = "smgbp"
# oncoplot(
#   maf = tcga_gbm,
#   pathways = pathways,
#   clinicalFeatures = "gender",
#   gene_mar = 10
# )
# oncoplot(
#   maf = tcga_gbm,
#   pathways = pathways,
#   clinicalFeatures = "gender",
#   gene_mar = 10,
#   sortByAnnotation = TRUE
# )
# oncoplot(
#   maf = tcga_gbm,
#   pathways = pathways,
#   clinicalFeatures = "gender",
#   gene_mar = 10,
#   sortByAnnotation = TRUE,
#   collapsePathway = TRUE
# )

# # Focus on a few genes
# oncoplot(
#   maf = tcga_gbm,
#   genes = c("EZH2", "KDM6A", "S1PR1"),
#   clinicalFeatures = "gender",
#   sortByAnnotation = TRUE
# )

# # Try some analysis
# # not sure why this doesn't work
# mafSurvival(
#   maf = tcga_gbm,
#   genes = "KDM6A",
#   time = "days_to_last_followup",
#   Status = "vital_status",
#   isTCGA = TRUE
# )

# Let's use the data we pulled
# Bring the Tumor Sample Barcode column into the clinical data, and set the days_to_last_follow_up as numeric
clinical_annot_merge <- merge(
  clinical_annot,
  maf_data[, c("case_id", "Tumor_Sample_Barcode")],
  by = "case_id"
)
colnames(clinical_annot_merge)[colnames(clinical_annot_merge)=="vital_status"] <- "Overall_Survival_Status"
clinical_annot_merge$days_to_last_follow_up <- as.numeric(as.character(clinical_annot_merge$days_to_last_follow_up))
clinical_annot_merge <- subset(clinical_annot_merge, Overall_Survival_Status==c("Alive", "Dead"))
levels(clinical_annot_merge$Overall_Survival_Status)
clinical_annot_merge$Overall_Survival_Status <- droplevels(clinical_annot_merge$Overall_Survival_Status)
levels(clinical_annot_merge$Overall_Survival_Status)
clinical_annot_merge$survival_cens <- as.logical(clinical_annot_merge$Overall_Survival_Status == "Alive")

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
  maf = gdc_gbm
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
  genes = c("PTEN"),
  time = "days_to_last_follow_up",
  Status = "survival_cens"
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

oncoplot(
  maf = gdc_gbm,
  pathways = "smgbp",
  selectedPathways = "Immune signaling",
  clinicalFeatures = "gender",
  gene_mar = 10,
  sortByAnnotation = TRUE
)

somaticInteractions(
  maf = gdc_gbm,
  top = 25,
  pvalue = c(0.05, 0.1),
  fontSize = 0.5
)

mafSurvival(maf = gdc_gbm, gene="PTEN", time = 'days_to_last_follow_up', Status = 'Overall_Survival_Status', isTCGA = TRUE)

#########
# Subset GDC GBM data by MAF column
gdc_gbm_Xlinked <- subsetMaf(maf = gdc_gbm, query = "Chromosome == 'chrX'")
gdc_gbm_kdm6a <- subsetMaf(maf = gdc_gbm, query = "Hugo_Symbol == 'KDM6A'")

# Subset GDC GBM data by clinical annotations
# gender is male
gdc_gbm_male <- subsetMaf(maf = gdc_gbm, clinQuery = "gender %in% 'male'")
gdc_gbm_Xlinked_male <- subsetMaf(
  maf = gdc_gbm,
  query = "Chromosome == 'chrX'",
  clinQuery = "gender %in% 'male'"
)
# gender is female
gdc_gbm_female <- subsetMaf(maf = gdc_gbm, clinQuery = "gender %in% 'female'")
gdc_gbm_Xlinked_female <- subsetMaf(
  maf = gdc_gbm,
  query = "Chromosome == 'chrX'",
  clinQuery = "gender %in% 'female'"
)

oncoplot(
  maf = gdc_gbm,
  top = 25
)

oncoplot(
  maf = gdc_gbm_Xlinked,
  top = 25
)

plotmafSummary(
  maf = gdc_gbm_kdm6a
)
plotmafSummary(
  maf = gdc_gbm_male
)
plotmafSummary(
  maf = gdc_gbm_female
)

# Plot mutations with a gender annotation for amle subset
oncoplot(
  maf = gdc_gbm_male,
  top = 25,
  clinicalFeatures = "gender",
  sortByAnnotation = TRUE
)
oncoplot(
  maf = gdc_gbm_Xlinked_male,
  top = 25,
  clinicalFeatures = "gender",
  sortByAnnotation = TRUE
)

# Plot mutations with a gender annotation for female subset
oncoplot(
  maf = gdc_gbm_female,
  top = 25,
  clinicalFeatures = "gender",
  sortByAnnotation = TRUE
)
oncoplot(
  maf = gdc_gbm_Xlinked_female,
  top = 25,
  clinicalFeatures = "gender",
  sortByAnnotation = TRUE
)

# Not useful
# laml.titv_gdc_gbm_kdm6a <- titv(
#   maf = gdc_gbm_kdm6a,
#   plot = FALSE,
#   useSyn = TRUE
# )
# plotTiTv(
#   res = laml.titv_gdc_gbm_kdm6a
# )

lollipopPlot(
  maf = gdc_gbm_kdm6a,
  gene = "KDM6A",
  showMutationRate = TRUE
)

rainfallPlot(
  maf = gdc_gbm
)

rainfallPlot(
  maf = gdc_gbm_male
)

# does not work
# laml.mutload_gdc_gbm = tcgaCompare(
#   maf = gdc_gbm,
#   cohortName = 'GDC-GBM',
#   logscale = TRUE,
#   capture_size = 50
# )

# plotVaf(
#   maf = gdc_gbm
# )

# maftools survival functions look for Tumor_Sample_Barcode in clinical data
# But GBM project data only has case_id, so use our merged dataframe to provide
prog_geneset = survGroup(
  maf = gdc_gbm,
  top = 20,
  geneSetSize = 2,
  time = "days_to_last_follow_up",
  Status = "vital_status",
  verbose = FALSE,
  clinicalData = analysis_df
)
print(prog_geneset)

tsb_kdm6a <- unique(analysis_df$Tumor_Sample_Barcode[analysis_df$Hugo_Symbol == "KDM6A"])
tsb_not_kdm6a <- unique(analysis_df$Tumor_Sample_Barcode[!(analysis_df$Tumor_Sample_Barcode %in% tsb_kdm6a)])

mafSurvival(
  maf = gdc_gbm,
  genes = "KDM6A",
  clinicalData = analysis_df,
  time = "days_to_last_follow_up",
  Status = "vital_status",
  isTCGA = TRUE
)

# Plot with pathways
pathways = "smgbp"
# Top pathways
oncoplot(
  maf = gdc_gbm,
  pathways = pathways,
  gene_mar = 10
)
# Summarizing known driver genes [Bailey et al., https://doi.org/10.1016/j.cell.2018.02.060]
# Drawing upto top 3 mutated pathways
#                                Pathway  N n_affected_genes fraction_affected Mutated_samples Fraction_mutated_samples
#  1:                     PI3K signaling  9                9         1.0000000             229              0.500000000
#  2:                   Genome integrity 14               14         1.0000000             190              0.414847162
#  3:                      RTK signaling 16               16         1.0000000             143              0.312227074
#  4:                              Other 22               19         0.8636364             121              0.264192140
#  5:               Transcription factor 39               29         0.7435897              87              0.189956332
#  6:                     MAPK signaling  9                7         0.7777778              66              0.144104803
#  7:        Chromatin histone modifiers 15               14         0.9333333              51              0.111353712
#  8:          Chromatin SWI/SNF complex  8                7         0.8750000              50              0.109170306
#  9:                         Cell cycle  8                5         0.6250000              49              0.106986900
# 10:                    Chromatin other 14               11         0.7857143              41              0.089519651
# 11:                    Other signaling 28               22         0.7857143              41              0.089519651
# 12: Protein homeostasis/ubiquitination 15               11         0.7333333              40              0.087336245
# 13:                         Metabolism  2                2         1.0000000              27              0.058951965
# 14:                      RNA abundance 15               11         0.7333333              22              0.048034934
# 15:            Wnt/B-catenin signaling  8                7         0.8750000              18              0.039301310
# 16:               Histone modification  3                2         0.6666667              13              0.028384279
# 17:                   Immune signaling 10                5         0.5000000              12              0.026200873
# 18:                     TGFB signaling  7                6         0.8571429              10              0.021834061
# 19:                      TOR signaling  3                2         0.6666667               8              0.017467249
# 20:                     NFKB signaling  2                2         1.0000000               8              0.017467249
# 21:                          Apoptosis  3                2         0.6666667               4              0.008733624
# 22:                           Splicing  6                4         0.6666667               4              0.008733624
# 23:          Epigenetics DNA modifiers  1                1         1.0000000               3              0.006550218
# 24:                    NOTCH signaling  1                0         0.0000000               0              0.000000000
#                                Pathway  N n_affected_genes fraction_affected Mutated_samples Fraction_mutated_samples
# Specific pathways
oncoplot(
  maf = gdc_gbm,
  pathways = pathways,
  gene_mar = 10,
  selectedPathways = "Transcription factor"
)
oncoplot(
  maf = gdc_gbm,
  pathways = pathways,
  gene_mar = 10,
  selectedPathways = "Chromatin histone modifiers"
)
oncoplot(
  maf = gdc_gbm,
  pathways = pathways,
  gene_mar = 10,
  selectedPathways = "Histone modification"
)
oncoplot(
  maf = gdc_gbm,
  pathways = pathways,
  gene_mar = 10,
  selectedPathways = "Immune signaling"
)

# Oncogenic signaling pathways
pws = pathways(
  maf = gdc_gbm,
  pathdb = "smgbp"
)
write.csv(pws,
          "export_csv\\pathway_summary.csv")

pws_xlinked = pathways(
  maf = gdc_gbm_Xlinked,
  pathdb = "smgbp"
)
write.csv(pws_xlinked,
          "export_csv\\pathway_summary_xlinked.csv")

pws_male = pathways(
  maf = gdc_gbm_male,
  pathdb = "smgbp"
)
write.csv(pws_male,
          "export_csv\\pathway_summary_male.csv")

pws_female = pathways(
  maf = gdc_gbm_female,
  pathdb = "smgbp"
)
write.csv(pws_female,
          "export_csv\\pathway_summary_female.csv")

pws_xlinked_male = pathways(
  maf = gdc_gbm_Xlinked_male,
  pathdb = "smgbp"
)
write.csv(pws_xlinked_male,
          "export_csv\\pathway_summary_xlinked_male.csv")

pws_xlinked_female = pathways(
  maf = gdc_gbm_Xlinked_female,
  pathdb = "smgbp"
)
write.csv(pws_xlinked_female,
          "export_csv\\pathway_summary_xlinked_female.csv")