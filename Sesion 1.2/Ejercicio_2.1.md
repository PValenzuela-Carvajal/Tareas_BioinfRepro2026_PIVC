# 4. Running the pipeline

## 4.1 Clean the data

In a typical analysis, data will be received from an Illumina sequencer, or some other type of sequencer as FASTQ files. The first requirement is to demultiplex, or sort, the raw data to recover the individual samples in the Illumina library. While doing this, we will use the Phred scores provided in the FASTQ files to discard sequencing reads of low quality. These tasks are accomplished using <font style="color: green;"> process_radtags </font> program. 

Some things to consider when running this program: 

* <font style="color: green;"> process_radtags </font> can handle both single-end or paired-end Illumina sequencing.

* The raw data can be compressed, or gzipped (files end with a `.gz` suffix).

* You can supply a list of barcodes, or indexes, to <font style="color: green;"> process_radtags </font> in order for it to demultiplex your samples. These barcodes can be single-end barcodes or combinatorial barcodes (pairs of barcodes, one on each of the paired reads). Barcodes are specified, one per line (or in tab separated pairs per line), in a text file.
  
  * If, in addition to your barcodes, you also supply a sample name in an extra column within the barcodes file, <span style="color: green;"> process_radtags </span> will name your output files according to sample name instead of barcode.
  
  * If you supply the same sample name for multiple barcodes, reads containing those barcodes will be consolidated into a single output file, merging them

* If you believe your reads may contain adapter contamination, <font style="color: green;"> process_radtags </font> can filter it out.

* You can supply the restriction enzyme used to construct the library. In the case of double-digest RAD, you can supply both restriction enzymes.

* If instructed, (`-r` command line option), <font style="color: green;"> process_radtags </font> will correct barcodes and restriction enzyme sites that are within a certain distance from the true barcode or restriction enzyme cutsite.

* By default, <font style="color: green;"> process_radtags </font> will identify and filter reads containing runs of poly-Gs. These runs are often indicative of synethesis termination in two-channel sequencing platforms &mdash they represent the absence of sequence once the DNA molecule "runs out" during sequencing. The presence of these poly-G runs is platform dependent, thus <font style="color: green;"> process_radtags </font> will use the FASTQ headers to identify the correct platorm whenever possible and adjust filters accordingly. This behavior can be modified in the advanced options. 

### 4.1.1 Understanding barcodes/indexes and specifying the barcode type

Genotype by sequencing libraries sample the genome by selecting DNA adjacent to one or more restriction enzyme cutsites. By reducing the amount of total DNA sampled, most researchers will multiplex many samples into one molecular library. Individual samples are demarcated in the library by ligating an oligo barcode onto the restriction enzyme-associated DNA for each sample. Alternatively, an index barcode is used, where the barcode is located upstream of the sample DNA within the sequencing adaptor. Regardless of the type of barcode used, after sequencing, the data must be demultiplexed so the samples are again separated. The process_radtags program will perform this task, but we much specify the type of barcodes used, and where to find them in the sequencing data.

There are a number of different configurations possible, each of them is detailed below. 

1. If your data are **single-end** or **paired-end**, with an inline barcode present only on the single-end (marked in red): 
```
@HWI-ST0747:188:C09HWACXX:1:1101:2968:2083 1:N:0:
TTATGATGCAGGACCAGGATGACGTCAGCACAGTGCGGGTCCTCCATGGATGCTCCTCGGTCGTGGTTGGGGGAGGAGGCA
+
@@@DDDDDBHHFBF@CCAGEHHHBFGIIFGIIGIEDBBGFHCGIIGAEEEDCC;A?;;5,:@A?=B5559999B@BBBBBA
@HWI-ST0747:188:C09HWACXX:1:1101:2863:2096 1:N:0:
TTATGATGCAGGCAAATAGAGTTGGATTTTGTGTCAGTAGGCGGTTAATCCCATACAATTTTACACTTTATTCAAGGTGGA
+
CCCFFFFFHHHHHJJGHIGGAHHIIGGIIJDHIGCEGHIFIJIH7DGIIIAHIJGEDHIDEHJJHFEEECEFEFFDECDDD
@HWI-ST0747:188:C09HWACXX:1:1101:2837:2098 1:N:0:
GTGCCTTGCAGGCAATTAAGTTAGCCGAGATTAAGCGAAGGTTGAAAATGTCGGATGGAGTCCGGCAGCAGCGAATGTAAA
```
Then you can specify the `--inline_null` flag to <font style="color: green;"> process_radtags </font>. This is also the default behavior and the flag can be ommitted in this case. 
2. If your data are **single-end** or **paired-end**, with a single index barcode (in blue):
```
@9432NS1:54:C1K8JACXX:8:1101:6912:1869 1:N:0:ATGACT
TCAGGCATGCTTTCGACTATTATTGCATCAATGTTCTTTGCGTAATCAGCTACAATATCAGGTAATATCAGGCGCA
+
CCCFFFFFHHHHHJJJJJJJJIJJJJJJJJJJJHIIJJJJJJIJJJJJJJJJJJJJJJJJJJGIJJJJJJJHHHFF
@9432NS1:54:C1K8JACXX:8:1101:6822:1873 1:N:0:ATGACT
CAGCGCATGAGCTAATGTATGTTTTACATTCCAGAAAGAGAGCTACTGCTGCAGGTTGTGATAAAATAAAGTAAGA
+
B@@FFFFFHFFHHJJJJFHIJHGGGHIJIIJIJCHJIIGGIIIGGIJEHIJJHII?FFHICHFFGGHIIGG@DEHH
@9432NS1:54:C1K8JACXX:8:1101:6793:1916 1:N:0:ATGACT
TTTCGCATGCCCTATCCTTTTATCACTCTGTCATTCAGTGTGGCAGCGGCCATAGTGTATGGCGTACTAAGCGAAA
+
@C@DFFFFHGHHHGIGHHJJJJJJJGIJIJJIGIJJJJHIGGGHGII@GEHIGGHDHEHIHD6?493;AAA?;=;=
```
Then yoy can specify the `--index_null` flag to process_radtags.

3. If your data are **single-end** with both an inline barcode (in red) and an index barcode (in blue):
```
@9432NS1:54:C1K8JACXX:8:1101:6912:1869 1:N:0:ATCACG
TCACGCATGCTTTCGACTATTATTGCATCAATGTTCTTTGCGTAATCAGCTACAATATCAGGTAATATCAGGCGCA
+
CCCFFFFFHHHHHJJJJJJJJIJJJJJJJJJJJHIIJJJJJJIJJJJJJJJJJJJJJJJJJJGIJJJJJJJHHHFF
@9432NS1:54:C1K8JACXX:8:1101:6822:1873 1:N:0:ATCACG
GTCCGCATGAGCTAATGTATGTTTTACATTCCAGAAAGAGAGCTACTGCTGCAGGTTGTGATAAAATAAAGTAAGA
+
B@@FFFFFHFFHHJJJJFHIJHGGGHIJIIJIJCHJIIGGIIIGGIJEHIJJHII?FFHICHFFGGHIIGG@DEHH
@9432NS1:54:C1K8JACXX:8:1101:6793:1916 1:N:0:ATCACG
GTCCGCATGCCCTATCCTTTTATCACTCTGTCATTCAGTGTGGCAGCGGCCATAGTGTATGGCGTACTAAGCGAAA
+
@C@DFFFFHGHHHGIGHHJJJJJJJGIJIJJIGIJJJJHIGGGHGII@GEHIGGHDHEHIHD6?493;AAA?;=;=
```
Then yoy can specify the `--inline_null` flag to process_radtags.
4. If  your data are **paired-end** with an inline barcode on the single-end (in red) and an index barcode (in blue):
```
@9432NS1:54:C1K8JACXX:7:1101:5584:1725 1:N:0:CGATGT
ACTGGCATGATGATCATAGTATAACGTGGGATACATATGCCTAAGGCTAAAGATGCCTTGAAGCTTGGCTTATGTT
+
#1=DDDFFHFHFHIFGIEHIEHGIIHFFHICGGGIIIIIIIIAEIGIGHAHIEGHHIHIIGFFFGGIIIGIIIEE7
@9432NS1:54:C1K8JACXX:7:1101:5708:1737 1:N:0:CGATGT
TTCGACATGTGTTTACAACGCGAACGGACAAAGCATTGAAAATCCTTGTTTTGGTTTCGTTACTCTCTCCTAGCAT
+
#1=DFFFFHHHHHJJJJJJJJJJJJJJJJJIIJIJJJJJJJJJJIIJJHHHHHFEFEEDDDDDDDDDDDDDDDDD@
```
```
@9432NS1:54:C1K8JACXX:7:1101:5584:1725 2:N:0:CGATGT
AATTTACTTTGATAGAAGAACAACATAAGCCAAGCTTCAAGGCATCTTTAGCCTTAGGCATATGTATCCCACGTTA
+
@@@DFFFFHGHDHIIJJJGGIIIEJJJCHIIIGIJGGEGGIIGGGIJIJIHIIJJJJIJJJIIIGGIIJJJIICEH
@9432NS1:54:C1K8JACXX:7:1101:5708:1737 2:N:0:CGATGT
AGTCTTGTGAAAAACGAAATCTTCCAAAATGCTAGGAGAGAGTAACGAAACCAAAACAAGGATTTTCAATGCTTTG
+
C@CFFFFFHHHHHJJJJJJIJJJJJJJJJJJJJJIJJJHIJJFHIIJJJJIIJJJJJJJJJHGHHHHFFFFFFFED
```
Then yoy can specify the `--inline_index` flag to process_radtags.
5. If your data are **paired-end** with indexed barcodes on the single and paired-ends (in blue):
```
@9432NS1:54:C1K8JACXX:7:1101:5584:1725 1:N:0:ATCACG+CGATGT
ACTGGCATGATGATCATAGTATAACGTGGGATACATATGCCTAAGGCTAAAGATGCCTTGAAGCTTGGCTTATGTT
+
#1=DDDFFHFHFHIFGIEHIEHGIIHFFHICGGGIIIIIIIIAEIGIGHAHIEGHHIHIIGFFFGGIIIGIIIEE7
@9432NS1:54:C1K8JACXX:7:1101:5708:1737 1:N:0:ATCACG+CGATGT
TTCGACATGTGTTTACAACGCGAACGGACAAAGCATTGAAAATCCTTGTTTTGGTTTCGTTACTCTCTCCTAGCAT
+
#1=DFFFFHHHHHJJJJJJJJJJJJJJJJJIIJIJJJJJJJJJJIIJJHHHHHFEFEEDDDDDDDDDDDDDDDDD@
```
```
@9432NS1:54:C1K8JACXX:7:1101:5584:1725 2:N:0:ATCACG+CGATGT
AATTTACTTTGATAGAAGAACAACATAAGCCAAGCTTCAAGGCATCTTTAGCCTTAGGCATATGTATCCCACGTTA
+
@@@DFFFFHGHDHIIJJJGGIIIEJJJCHIIIGIJGGEGGIIGGGIJIJIHIIJJJJIJJJIIIGGIIJJJIICEH
@9432NS1:54:C1K8JACXX:7:1101:5708:1737 2:N:0:ATCACG+CGATGT
AGTCTTGTGAAAAACGAAATCTTCCAAAATGCTAGGAGAGAGTAACGAAACCAAAACAAGGATTTTCAATGCTTTG
+
C@CFFFFFHHHHHJJJJJJIJJJJJJJJJJJJJJIJJJHIJJFHIIJJJJIIJJJJJJJJJHGHHHHFFFFFFFED
```
Then you can specify the `--index_index` flag to process_radtags. 
6. If your data are **paired-end** with inline barcodes on the single and paired-ends (in red): 
```
@9432NS1:54:C1K8JACXX:7:1101:5584:1725 1:N:0:
ACTGGCATGATGATCATAGTATAACGTGGGATACATATGCCTAAGGCTAAAGATGCCTTGAAGCTTGGCTTATGTT
+
#1=DDDFFHFHFHIFGIEHIEHGIIHFFHICGGGIIIIIIIIAEIGIGHAHIEGHHIHIIGFFFGGIIIGIIIEE7
@9432NS1:54:C1K8JACXX:7:1101:5708:1737 1:N:0:
TTCGACATGTGTTTACAACGCGAACGGACAAAGCATTGAAAATCCTTGTTTTGGTTTCGTTACTCTCTCCTAGCAT
+
#1=DFFFFHHHHHJJJJJJJJJJJJJJJJJIIJIJJJJJJJJJJIIJJHHHHHFEFEEDDDDDDDDDDDDDDDDD@
```
```
@9432NS1:54:C1K8JACXX:7:1101:5584:1725 2:N:0:
AATTTACTTTGATAGAAGAACAACATAAGCCAAGCTTCAAGGCATCTTTAGCCTTAGGCATATGTATCCCACGTTA
+
@@@DFFFFHGHDHIIJJJGGIIIEJJJCHIIIGIJGGEGGIIGGGIJIJIHIIJJJJIJJJIIIGGIIJJJIICEH
@9432NS1:54:C1K8JACXX:7:1101:5708:1737 2:N:0:
AGTCTTGTGAAAAACGAAATCTTCCAAAATGCTAGGAGAGAGTAACGAAACCAAAACAAGGATTTTCAATGCTTTG
+
C@CFFFFFHHHHHJJJJJJIJJJJJJJJJJJJJJIJJJHIJJFHIIJJJJIIJJJJJJJJJHGHHHHFFFFFFFED
```
Then you can specify the `--inline_inline` flag to process_radtags. 
### 4.1.2 Specifying the barcodes

The process_radtags program will demultiplex data if it is told which barcodes/indexes to expect in the data. It will also properly name the output files if the user specifies how to translate a particular barcode to a specific output file neam. This is done with the barcodes file, which we provide to process_radtags program. The barcode file is a very simple format — one barcode per line; if you want to rename the output files, the sample name prefix is provided in the second column. 

```
cat barcodes_lane3
CGATA<tab>sample_01
CGGCG     sample_02
GAAGC     sample_03
GAGAT     sample_04
TAATG     sample_05
TAGCA     sample_06
AAGGG     sample_07
ACACG     sample_08
ACGTA     sample_09
```

The sample names can be whatever is meaningful for your project: 

```
more barcodes_run01_lane01
CGATA<tab>spruce_site_12-01
CGGCG     spruce_site_12-02
GAAGC     spruce_site_12-03
GAGAT     spruce_site_12-04
TAATG     spruce_site_06-01
TAGCA     spruce_site_06-02
AAGGG     spruce_site_06-03
ACACG     spruce_site_06-04
```

Combinatorial barcodes are specified, one per column, separated by a tab: 

```
cat barcodes_lane07
CGATA<tab>ACGTA<tab>sample_01
CGGCG     ACGTA     sample_02
GAAGC     ACGTA     sample_03
GAGAT     ACGTA     sample_04
CGATA     TAGCA     sample_05
CGGCG     TAGCA     sample_06
GAAGC     TAGCA     sample_07
GAGAT     TAGCA     sample_08
```

 Merging the first three sets of combinatorial barcodes into a single output file: 

```
cat barcodes_lane07
CGATA<tab>ACGTA<tab>sample_01
CGGCG     ACGTA     sample_01
GAAGC     ACGTA     sample_01
GAGAT     ACGTA     sample_02
CGATA     TAGCA     sample_03
```

1. If you don’t want process_radtags to rename your samples, simply do not specify the last column in the barcodes file, and the output files will instead be named after the barcode.
2. Often, sequencing centers will return data from indexed libraries already demultiplexed. In this case, omit the barcodes file and process_radtags will not attempt to demulitplex the data, but can still be used to clean the data.

### 4.1.3 Running process_radtags

 Here is how single-end data received from an Illumina sequencer might look: 

```
ls ./raw
lane3_NoIndex_L003_R1_001.fastq.gz  lane3_NoIndex_L003_R1_006.fastq.gz  lane3_NoIndex_L003_R1_011.fastq.gz
lane3_NoIndex_L003_R1_002.fastq.gz  lane3_NoIndex_L003_R1_007.fastq.gz  lane3_NoIndex_L003_R1_012.fastq.gz
lane3_NoIndex_L003_R1_003.fastq.gz  lane3_NoIndex_L003_R1_008.fastq.gz  lane3_NoIndex_L003_R1_013.fastq.gz
lane3_NoIndex_L003_R1_004.fastq.gz  lane3_NoIndex_L003_R1_009.fastq.gz
lane3_NoIndex_L003_R1_005.fastq.gz  lane3_NoIndex_L003_R1_010.fastq.gz
```

I specify the directory containing the input files, `./raw`, the directory I want process_radtags to enter the output files, `./samples`, and a file containing my barcodes, `./barcodes/barcodes_lane3`, along with the restrction enzyme I used and instructions to clean the data and correct barcodes and restriction enzyme cutsites (`-r`, `-c`, `-q`).

Here is a more complex example, using paired-end double-digested data (two restriction enzymes) with combinatorial barcodes, and gzipped input files. Here is what the raw Illumina files may look like: 

```
ls ./raw
GfddRAD1_005_ATCACG_L007_R1_001.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_001.fastq.gz
GfddRAD1_005_ATCACG_L007_R1_002.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_002.fastq.gz
GfddRAD1_005_ATCACG_L007_R1_003.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_003.fastq.gz
GfddRAD1_005_ATCACG_L007_R1_004.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_004.fastq.gz
GfddRAD1_005_ATCACG_L007_R1_005.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_005.fastq.gz
GfddRAD1_005_ATCACG_L007_R1_006.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_006.fastq.gz
GfddRAD1_005_ATCACG_L007_R1_007.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_007.fastq.gz
GfddRAD1_005_ATCACG_L007_R1_008.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_008.fastq.gz
GfddRAD1_005_ATCACG_L007_R1_009.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_009.fastq.gz
```

Now we specify both restriction enzymes using the --renz_1 and --renz_2 flags along with the type combinatorial barcoding used. Here is the command:

```
process_radtags -P -p ./raw  -b ./barcodes/barcodes -o ./samples/ \
                  -c -q -r --inline_index --renz_1 nlaIII --renz_2 mluCI
```

##### The output of process_radtags

The output of the process_radtags differs depending if you are processing single-end or paired-end data. In the case of single-end reads, the program will output one file per barcode into the output directory you specify. If the data do not have barcodes, then the file will retain its original name.

If you are processing paired-end reads, then you will get four files per barcode, two for the single-end read and two for the paired-end read. For example, given barcode ACTCG, you would see the following four files: 

```
sample_ACTCG.1.fq
sample_ACTCG.rem.1.fq
sample_ACTCG.2.fq
sample_ACTCG.rem.2.fq
```

The process_radtags program wants to keep the reads in phase, so that the first read in the `sample_XXX.1.fq` file is the mate of the first read in the `sample_XXX.2.fq` file. Likewise for the second pair of reads being the second record in each of the two files and so on. When one read in a pair is discarded due to low quality or a missing restriction enzyme cut site, the remaining read can’t simply be output to the `sample_XXX.1.fq` or `sample_XXX.2.fq` files as it would cause the remaining reads to fall out of phase. Instead, this read is considered a *remainder* read and is output into the `sample_XXX.rem.1.fq` file if the paired-end was discarded, or the `sample_XXX.rem.2.fq` file if the single-end was discarded. 

##### Modifying how process_radtags executes

 The process_radtags program can be modified in several ways. If your data do not have barcodes, omit the barcodes file and the program will not try to demultiplex the data. You can also disable the checking of the restriction enzyme cut site, or modify what types of quality are checked for. So, the program can be modified to only demultiplex and not clean, clean but not demultiplex, or some combination.

There is additional information available in process_radtags [manual page](https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php). 

##### Example: Processing data from NCBI’s Short Read Archive

process_radtags can be used to process public data available in SRA and other short-read sequence databases; however, this often requires providing specific parameters when running the program. 

- SRA data is often available separately for each sample, therefore de-multiplexing may not be required. In such cases, process_radtags is run separately for each sample, omitting the barcodes `(-b/--barcodes)` option. 

- FASTQ files downloaded from SRA often do not follow Illumina’s file naming convention and should therefore be specified as generic FASTQ file inputs to process_radtags. This can be done using the `-f` and `-1/-2` options, depending on if the downloaded data is single- or paired-end, respectively. In such cases, the `--basename` option can be used to provide a specific prefix name to the output files and logs. Output files and logs are otherwise named according to the specified inputs. 

- By default, process_radtags will use the information available on the input file’s FASTQ headers to detect the specific sequencing platform. This information is used to identify poly-G artefacts when filtering reads. These intact FASTQ headers might not be available in FASTQ files downloaded from SRA, which will lead to poly-G runs being ignored during filtering. The user can specify the `--force-poly-g-check` option, which will result in the filtering of reads with poly-G runs, independent of the sequencing platform information. 

- Additional parameters might be needed depending on how the data was processed prior to submission to SRA. For example, the data might contain variable read lengths (which can be made uniform using `-t/--truncate along` with `--len-limit`) or the restriction enzyme cut site was previously removed, which may require disabling the cut site check (`--disable-rad-check`). 

Here, we provide an example of how to process RADseq data available on SRA (accession number [SRR20082702](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR20082702&display=metadata)). The data is for an single mackarel icefish individual, paired-end sequenced on a single-digest RADseq library using the enzyme SbfI ([Rivera-Colón et al. 2023](doi.org/10.1093/molbev/msad029)). 

First, we download the data from NCBI using the fastq-dump program, available as part of NCBI’s [*SRA Toolkit*](https://github.com/ncbi/sra-tools).

```
fastq-dump --accession SRR20082702 --gzip --split-files \
    --defline-seq '@$ac.$si/$ri' --defline-qual '+$ac.$si/$ri' --outdir ./raw
```

After downloading the paired files from this accession, named `SRR20082702_1.fastq.gz` and `SRR20082702_2.fastq.gz`, we can proceed with process_radtags. We want to apply quality filters to this data (`--clean` and `--quality`), specifying and rescueing SbfI cut site (`--renz-1 sbfI` and `--rescue`), and manually enabling the filtering of poly-G runs (`--force-poly-g-check`). Also, we want to apply a more informative name to the resulting files, which we can specify using the `--basename` option. 

```
process_radtags -1 ./raw/SRR20082702_1.fastq.gz -2 ./raw/SRR20082702_2.fastq.gz \
    --out-path ./processed --clean --quality --rescue --renz-1 sbfI \
    --force-poly-g-check --basename cgunnari_49
```

Upon completion, process_radtags will generate the usual outputs, which will be named according to the value specified by the `--basename option`. In this case, `cgunnari_49`: 

```
cgunnari_49.1.fq
    cgunnari_49.rem.1.fq
    cgunnari_49.2.fq
    cgunnari_49.rem.2.fq
    cgunnari_49.log
```
