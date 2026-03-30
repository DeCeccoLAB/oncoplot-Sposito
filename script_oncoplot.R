library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(maftools)

# read the mutation file with allele frequency (af)
file_af <- read_xlsx('Mutations_af.xlsx')

file_mut <- subset(file_af, select = -c(SAMPLE, `SUBJECT ID`, VAF))
dim(file_mut) #119

file_mut <- file_mut %>%
  pivot_longer(
    cols = ends_with("af"),
    names_to = "Hugo_Symbol",
    values_to = "af"
  ) %>%
  filter(af > 0) %>%
  mutate(
    Hugo_Symbol = str_remove(Hugo_Symbol, " af")
  )
dim(file_mut)

file_mut <- file_mut %>%
  mutate(
    cDNA = str_extract(ID, "c\\..*(?=_p\\.)"),
    Protein_Change = str_extract(ID, "p\\..+")
  )

file_mut <- file_mut %>%
  mutate(
    Variant_Classification = case_when(
      str_detect(Protein_Change, "fs") & str_detect(cDNA, "del") ~ "Frame_Shift_Del",
      str_detect(Protein_Change, "fs") & str_detect(cDNA, "ins") ~ "Frame_Shift_Ins",
      str_detect(cDNA, "del") & str_detect(Protein_Change,"del") ~ "In_Frame_Del",
      str_detect(cDNA, "ins") & str_detect(Protein_Change,"ins") ~ "In_Frame_Ins",
      str_detect(Protein_Change, "\\*") ~ "Nonsense_Mutation",
      str_detect(Protein_Change, "^p\\.[A-Z][0-9]+[A-Z]$") ~ "Missense_Mutation"
    )
  )

file_mut <- file_mut%>%
  mutate(
    Variant_Type = case_when(
      str_detect(cDNA, "del") ~ "DEL",
      str_detect(cDNA, "ins") ~ "INS",
      str_detect(cDNA, ">") ~ "SNP",
      TRUE ~ "SNP"
    )
  )

file_mut <- file_mut %>%
  mutate(
    Reference_Allele = case_when(
      str_detect(Variant_Type,"INS") ~ "-",
      str_detect(Variant_Type,"DEL") ~ str_extract(cDNA, "(?<=del)[A-Z]+"),
      str_detect(Variant_Type,"SNP") ~ str_extract(cDNA, "(?<=\\d)[A-Z](?=>)")
    ),
    Tumor_Seq_Allele2 = case_when(
      str_detect(Variant_Type,"INS") ~ str_extract(cDNA, "(?<=ins)[A-Z]+"),
      str_detect(Variant_Type,"DEL") ~ "-",
      str_detect(Variant_Type,"SNP") ~ str_extract(cDNA, "(?<=>)[A-Z]")
    ),
    Chromosome = NA,   # or real data if available
    Start_Position = NA, # ideally from annotation
    End_Position = NA,
    Tumor_Sample_Barcode = `SAMPLE ID`,
    Tumor_Seq_Allele1 = Reference_Allele
  )
### get the chromosome
library(biomaRt) 
# Connect to Ensembl
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
# Get chromosome name, start, and end position
gene_locs <- getBM(attributes = c("hgnc_symbol", "chromosome_name", 
                                  "start_position", "end_position"),
                   filters = "hgnc_symbol",
                   values = unique(file_mut$Hugo_Symbol),
                   mart = mart)
file_mut$Chromosome <- gene_locs$chromosome_name[match(file_mut$Hugo_Symbol,gene_locs$hgnc_symbol)]

# not biologically correct, but useful for the visualization
file_mut <- file_mut %>%
  mutate(
    Start_Position = as.numeric(str_extract(cDNA, "\\d+")),
    End_Position = Start_Position
  )

maf <- read.maf(maf = file_mut)

pdf('oncoplot.pdf',height = 10,width = 10)
oncoplot(maf = maf,top = 50,draw_titv = T,fontSize = 0.7)
dev.off()
