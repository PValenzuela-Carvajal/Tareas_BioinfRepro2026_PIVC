VCF="/datos/compartido/ChileGenomico/GATK_ChGdb_recalibrated.autosomes.12262013.snps.known.vcf"

## Count variants

N_VAR=$(grep -v '^#' "$VCF" | wc -l)
N_IND=$(awk '{if ($1 == "#CHROM"){print NF-9; exit}}' "$VCF")

echo "Number of individuals: $N_IND"
echo "Number of variants:    $N_VAR"

## Count sites without missing data

NO_MISSING=$(vcftools --vcf $VCF --max-missing 1  2>&1 | awk '/After filtering, kept/ && /Sites/ {print $4}')
echo "Number of sites without missing data: $NO_MISSING"

## Create 2 Mb window

CHR="1"
START_MB=1
END_MB=3

START_BP=$((START_MB * 1000000))
END_BP=$((END_MB * 1000000))

WINDOW_NAME="CLG_Chr${CHR}_${START_MB}-${END_MB}Mb"

vcftools \
    --vcf "$VCF" \
    --chr "$CHR" \
    --from-bp "$START_BP" \
    --to-bp "$END_BP" \
    --recode \
    --recode-INFO-all \
    --out "./$WINDOW_NAME"

mv "./${WINDOW_NAME}.recode.vcf" \
   "./${WINDOW_NAME}.vcf"

WINDOW_VCF="./${WINDOW_NAME}.vcf"
echo "The output file is: $WINDOW_VCF"


# Number of variants in window

N_WINDOW=$(grep -v '^#' "$WINDOW_VCF" | wc -l)
echo "Variants in window: $N_WINDOW"


# Calculate the mean depth

vcftools \
    --vcf "$VCF" \
    --depth \
    --out depth

echo -e "Mean depth per individual:\n$(head depth.idepth)"


# Calculate allele frequencies

vcftools \
    --vcf "$VCF" \
    --freq \
    --out all_allele_frequencies

echo -e "Allele frequencies per site:\n$(head all_allele_frequencies.frq)"


## Calculate frequencies from biallelic sites only

vcftools \
    --vcf "$VCF" \
    --min-alleles 2 \
    --max-alleles 2 \
    --freq \
    --out biallelic_freq_output

echo -e "Allele frequencies for biallelic sites:\n$(head biallelic_freq_output.frq)"


## Generate MAF graphic

Rscript script_maf_graph.R


## Count sites with MAF < 0.05

vcftools \
    --vcf "$VCF" \
    --max-maf 0.05 \
    2>&1 | awk '/After filtering, kept/ && /Sites/ {print $4}'


## Calculate heterozygosity

vcftools \
    --vcf "$VCF" \
    --het \
    --out heterozygosity

echo -e "Heterozygosity per individual:\n$(head heterozygosity.het)"


## Pi calculations

vcftools \
    --vcf "$VCF" \
    --site-pi \
    --out nucleotide_diversity

echo -e "Nucleotide diversity per site:\n$(head nucleotide_diversity.sites.pi)"


## Filter the file

vcftools \
    --vcf "$VCF" \
    --maf 0.05 \
    --recode \
    --recode-INFO-all \
    --out GATK_MAF05


## Convert to PLINK

vcftools \
    --vcf "$VCF" \
    --plink \
    --out GATK_ChGdb_recalibrated
