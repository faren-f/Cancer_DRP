rm(list = ls())
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")


# Read Data ---------------------------------------------------------------
# PRISM tissue type
cellline_info = read.csv("Raw_data/PRISM/Secondary/secondary-screen-cell-line-info.csv")
cellline_info = cellline_info[complete.cases(cellline_info), ]
tissue_types_PRISM = unique(cellline_info$primary_tissue)
tissue_types_PRISM = data.frame(tissue_types_PRISM)
Sample_Tissue = readRDS("Processed_data/S19/sample_tissue_types.rds")
tissue_types_PRISM$NO_Tissues = apply(Sample_Tissue,2,sum)


# TCGA tissue types
NO_tissues_TCGA = readRDS("Processed_data/S27/NO_tissues.rds")   

TCGA_Patients = readRDS("Processed_data/S22/TCGA_Patients.rds")
TCGA_Patients = data.frame(TCGA_Patients)
colnames(TCGA_Patients) = c("Cancer_type", "Patient", "Drug",
                            "Response","Cancer_type_perfect_name")

tissue_types_TCGA = data.frame(unique(TCGA_Patients[,5]))



