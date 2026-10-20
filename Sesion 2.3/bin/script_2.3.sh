#!/bin/bash

# Variantes eliminadas por tasa de datos perdidos > 0.2
grep "variants removed" chilean_all48_hg19_2.log

# Individuos con tasa de datos perdidos > 0.02
grep "people removed" chilean_all48_hg19_5.log

# Individuos eliminados por discrepancia de sexo
wc -l sex_discrepancy.txt
wc -l chilean_all48_hg19_5.fam
wc -l chilean_all48_hg19_6.fam

# SNPs en cromosomas sexuales
wc -l chilean_all48_hg19_6.bim
wc -l snp_1_22.txt

# Número de cromosomas que porta cada alelo por SNP
plink --bfile chilean_all48_hg19_7 --freq counts --out allele_counts
awk '{print $2,$3,$4,$5,$6}' allele_counts.frq.counts | head

# SNPs en aparente equilibrio de ligamiento
wc -l indepSNP.prune.in

# SNPs eliminados por estar en regiones de inversiones conocidas
grep "variants removed" indepSNP.log

# Individuos restantes tras el filtro de parentesco
wc -l chilean_all48_hg19_10.fam

# SNPs restantes tras el filtro de r2 < 0.2
wc -l MDS_merge_r2_lt_0.2.bim

# Valores de error de validación cruzada de ADMIXTURE
grep "CV error" MDS_merge_r2_lt_0.2.K*.log