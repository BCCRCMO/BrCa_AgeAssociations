Lung cancer age associated genes - TCGA lung adenocarcinoma data

Example is in
main/code/Analysis_TCGA_BrCa_RNASeq_Expression_Trends
TCGAStudy.R

with report in
TCGA_BrCa_AgeAssociated_ERBinding_Genes.docx

Set up a similar script in

main/code/Analysis_TCGA_LungCa_RNASeq_Expression_Trends
called
TCGA_Lung_Study.R

using data files as follows:

### Clinical data for lung:

main/data/Data_TCGA_LungCa_RNASeq_Expression_Trends
data_bcr_clinical_data_patient.txt


### Expression data for lung:

main/data/Data_TCGA_LungCa_RNASeq_Expression_Trends
data_RNA_Seq_v2_expression_median.txt


### Annotation data:

main/data/Data_METABRIC_Expression_Trends
Annotation_Illumina_Human-WG-V3_hg18_V1.0.0_Aug09.txt


### Age associated gene list from Tomo Osako (in BrCa data dir already):

main/data/Data_TomoOsako_BrCa_AgeRelated_EZH2Related
AgeRelated_Accession_Numbersv02.csv
EZH2relatedGenesv08.csv


### ER binding genes from ChIP-Seq:

main/data/Data_JasonCarroll
ER_binding_gene_hg19_threshold_1e-5.txt


Goals:

Obtain regression fits for all 20000 genes for cases <=60, >60, all cases

Plots of data for all genes in gene name union of
AgeRelated_Accession_Numbersv02.csv
EZH2relatedGenesv08.csv

(All plots in one PDF, and also each plot in its own PDF)

Enrichment assessment via Fisher's exact test






