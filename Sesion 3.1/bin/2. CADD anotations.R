library(tidyverse)

vep <-read.delim("../data/S13_PASS_VEP.tsv")

rank_impact <- c(HIGH = 1, MODERATE = 2, LOW = 3, MODIFIER = 4)

##Colapsamos a una unica fila por variante 

vep_unica <- vep %>%
  mutate(impact_rank = rank_impact[IMPACT]) %>%
  distinct(Location, Ref, Allele, .keep_all = TRUE) %>%
  select(-impact_rank)

## Filtro de las variantes

filtradas <- vep_unica %>%
  filter(
    (!is.na(ClinVar_CLNSIG) & ClinVar_CLNSIG != "" & !grepl("benign", ClinVar_CLNSIG, ignore.case = TRUE)) | 
      (!is.na(CADD_PHRED) & CADD_PHRED > 20))

##Tabla final  
tabla_final <- filtradas %>%
  select(Location, Ref, Allele, SYMBOL, Consequence, IMPACT,HGVSc, HGVSp, CADD_PHRED, CADD_RAW, ClinVar_CLNSIG, PHENOTYPES)

tabla_final

