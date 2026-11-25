#!/bin/bash

cd ~/pvalenzuela-carvajal   # Ir a la carpeta de trabajo donde estan los datos filtrados

module load bwa samtools qualimap   # Cargar programas

# Alineamiento
bwa mem -t 4 -M /datos/reference/genomes/hg19_reference/hg19.fasta S13_R1_filter.fastq.gz S13_R2_filter.fastq.gz > S13.sam   # Alinear contra hg19
java -jar /opt/picard/picard-2.25.2/picard.jar SamFormatConverter -I S13.sam -O S13.bam   # SAM a BAM
java -jar /opt/picard/picard-2.25.2/picard.jar SortSam -I S13.bam -O S13_sorted.bam -SO coordinate   # Ordenar por posición
java -jar /opt/picard/picard-2.25.2/picard.jar AddOrReplaceReadGroups -I S13_sorted.bam -O S13_sorted_RG.bam -ID S13 -LB Paired-end -PL Illumina -PU Unknown -SM S13   # Agregar Read Group
samtools index S13_sorted_RG.bam   # Indexar BAM
ls S13_sorted_RG.bam*   # Verificar archivos

# Primera lectura con soft-clipping 
grep -v '^@' S13.sam | awk '$6 ~ /S/ {print; exit}' > lectura_softclip.txt   # Primera lectura con "S" en el CIGAR
cat lectura_softclip.txt   # Mostrar registro
tr '\t' '\n' < lectura_softclip.txt | head -11 | nl   # Mostrar campos numerados

# Reporte de calidad 
java -Xmx4g -cp "/opt/qualimap/qualimap_v2.2.1/qualimap.jar:/opt/qualimap/qualimap_v2.2.1/lib/*" \
org.bioinfo.ngs.qc.qualimap.main.NgsSmartMain bamqc \
-bam S13_sorted_RG.bam -gff regiones_blanco.bed -outdir ./S13_sorted_RG   # Qualimap sobre regiones blanco

#  Variables 
REF=/datos/reference/genomes/hg19_reference/hg19.fasta
DBSNP=/datos/reference/genomes/hg19_reference/dbSNP_hg19.vcf

# Recalibración de calidades
java -jar /opt/GenomeAnalysisTK-3.7-0/GenomeAnalysisTK.jar -T BaseRecalibrator -R $REF -I S13_sorted_RG.bam -knownSites $DBSNP -o S13_recal_data.table   # Modelo de errores

java -jar /opt/GenomeAnalysisTK-3.7-0/GenomeAnalysisTK.jar -T PrintReads -R $REF -I S13_sorted_RG.bam --BQSR S13_recal_data.table -o S13_recal_reads.bam   # Aplicar recalibración

# Llamado de variantes sobre las regiones blanco
java -jar /opt/GenomeAnalysisTK-3.7-0/GenomeAnalysisTK.jar -T HaplotypeCaller -R $REF -I S13_recal_reads.bam --dbsnp $DBSNP -stand_call_conf 30 -L regiones_blanco.bed -o S13_raw_variants.vcf   # Llamar variantes

# Separar SNPs e INDELs
java -jar /opt/GenomeAnalysisTK-3.7-0/GenomeAnalysisTK.jar -T SelectVariants -R $REF -V S13_raw_variants.vcf -selectType SNP -o S13_RAW_SNP.vcf   # Solo SNPs

java -jar /opt/GenomeAnalysisTK-3.7-0/GenomeAnalysisTK.jar -T SelectVariants -R $REF -V S13_raw_variants.vcf -selectType INDEL -o S13_RAW_INDEL.vcf   # Solo INDELs

# Filtros: DP < 10 y QD < 2.0
java -jar /opt/GenomeAnalysisTK-3.7-0/GenomeAnalysisTK.jar -T VariantFiltration -R $REF -V S13_RAW_SNP.vcf --filterExpression "DP < 10" --filterName "LowDP" --filterExpression "QD < 2.0" --filterName "LowQD" -o S13_FILTERED_SNP.vcf   # Filtrar SNPs

java -jar /opt/GenomeAnalysisTK-3.7-0/GenomeAnalysisTK.jar -T VariantFiltration -R $REF -V S13_RAW_INDEL.vcf --filterExpression "DP < 10" --filterName "LowDP" --filterExpression "QD < 2.0" --filterName "LowQD" -o S13_FILTERED_INDEL.vcf   # Filtrar INDELs

# Combinar SNPs e INDELs filtrados en un VCF
java -jar /opt/GenomeAnalysisTK-3.7-0/GenomeAnalysisTK.jar -T CombineVariants -R $REF --variant:foo S13_FILTERED_SNP.vcf --variant:bar S13_FILTERED_INDEL.vcf -o S13_FILTER_VARIANTS.vcf -genotypeMergeOptions PRIORITIZE -priority foo,bar

# Conteo por etiqueta de filtro (PASS, LowDP, LowQD, ambos)
for f in S13_FILTER_VARIANTS S13_FILTERED_SNP S13_FILTERED_INDEL; do
  echo "== $f"
  grep -v '^#' $f.vcf | cut -f7 | sort | uniq -c
done

# VCF solo con variantes que pasan ambos filtros
java -jar /opt/GenomeAnalysisTK-3.7-0/GenomeAnalysisTK.jar -T SelectVariants -R /datos/reference/genomes/hg19_reference/hg19.fasta -V S13_FILTER_VARIANTS.vcf --excludeFiltered -o S13_PASS.vcf