Human ATACSEQ
================
RO
4/28/2025

# Objective: Determine if doxycycline exposure changes chromatin accessibility in humans

# Approach:

## 1.1 Loading in ATACseq peak files with custom function import_human_peaks (list of GRanges output)

``` r
# establishing peak path to the dir with MACS2 output peak files from NF_CORE ATACseq pipeline

human_peak_path <- "/scratch/Shares/rinnclass/MASTER_CLASS/DATA/human_atacseq"

# creating a file list also needed for import_peaks function to get sample name associated with file
human_fl <- list.files(human_peak_path, full.names = TRUE, pattern = ".broadPeak")

# running import_human_peaks
my_human_peaks <- import_human_peaks(human_consensus_file_path = human_peak_path)

print("here are the number of peaks for each human sample")
```

    ## [1] "here are the number of peaks for each human sample"

``` r
print(num_peaks <- sapply(my_human_peaks, length) %>% as.data.frame)
```

    ##                            .
    ## novageneGFP_0.5h_REP1 148391
    ## novageneGFP_0.5h_REP2 136522
    ## novageneGFP_0.5h_REP3 154628
    ## novageneGFP_0h_REP1   148083
    ## novageneGFP_0h_REP2   153262
    ## novageneGFP_0h_REP3   140370
    ## novageneGFP_1.5h_REP1 115685
    ## novageneGFP_1.5h_REP2 119006
    ## novageneGFP_1.5h_REP3 104847
    ## novageneGFP_1h_REP1    89923
    ## novageneGFP_1h_REP2   130800
    ## novageneGFP_1h_REP3   131043
    ## novageneGFP_2.5h_REP1 144102
    ## novageneGFP_2.5h_REP2 154542
    ## novageneGFP_2.5h_REP3 148161
    ## novageneGFP_2h_REP1   136505
    ## novageneGFP_2h_REP2   143190
    ## novageneGFP_2h_REP3   150659

### Result: The number of peaks in these samples range from 90k-155k

## 1.2 Finding number of peaks common in all samples using find_common_peaks custom function

``` r
# run find_common_peaks function (no alterations needed)
common_peaks <-  suppressWarnings(find_common_peaks(my_human_peaks))

length(common_peaks)
```

    ## [1] 58371

### Result: There are 58,371 common peaks across all human samples. Visualization in IGV seems to confirm (see below).

<figure>
<img
src="/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/reanna/MASTER_CLASS/finalproject2/images/igv_snapshot.png"
alt="Human peaks visualized in IGV." />
<figcaption aria-hidden="true">Human peaks visualized in
IGV.</figcaption>
</figure>

## 1.3 Finding peaks unique to non-doxycyline and doxycycline conditions

Non-doxycycline conditions (time point 0h):

``` r
non_dox_human_samples <- my_human_peaks[grepl("novageneGFP_0h_REP", names(my_human_peaks))]

non_dox_common_human_peaks <- suppressWarnings(find_common_peaks(non_dox_human_samples))

print(c("This is how many peaks are common in non-doxycline conditions:", length(non_dox_common_human_peaks)))
```

    ## [1] "This is how many peaks are common in non-doxycline conditions:"
    ## [2] "101526"

Doxycycline conditions (all other time points):

``` r
dox_human_samples <- names(my_human_peaks)[!grepl("_0$", names(my_human_peaks))]

dox_human_peaks <- my_human_peaks[dox_human_samples]

dox_common_human_peaks <- suppressWarnings(find_common_peaks(dox_human_peaks))
print(c("This is how many peaks are common in doxycyline conditions",length(dox_common_human_peaks)))
```

    ## [1] "This is how many peaks are common in doxycyline conditions"
    ## [2] "58371"

Overlap of non-doxycyline and doxycycline conditions (to find peaks
unique to each condition)

``` r
human_dox_compare_list <- list(non_dox = non_dox_common_human_peaks, dox = dox_common_human_peaks)

human_dox_non_dox_ov <- suppressWarnings(find_common_peaks(human_dox_compare_list))

print(c("This is how many peaks are common in both non-doxycyline and doxycyline conditions", length(human_dox_non_dox_ov)))
```

    ## [1] "This is how many peaks are common in both non-doxycyline and doxycyline conditions"
    ## [2] "58654"

``` r
# extracting peaks unique to each condition (non-dox vs. dox)

# Peaks unique to non_dox:
unique_to_non_dox_human <-suppressWarnings(find_my_peaks(human_dox_non_dox_ov, non_dox_common_human_peaks))

print(c("This is how many peaks are unique to non-dox condition",length(unique_to_non_dox_human)))
```

    ## [1] "This is how many peaks are unique to non-dox condition"
    ## [2] "42872"

``` r
# Peaks unique to dox:
unique_to_dox_human <- suppressWarnings(find_my_peaks(human_dox_non_dox_ov, dox_common_human_peaks))

print(c("This is how many peaks are unique to dox condition", length(unique_to_dox_human)))
```

    ## [1] "This is how many peaks are unique to dox condition"
    ## [2] "13"

## 2 Creating mouse gene, lincrna, mRNA annotation GRange objects

Here I am going to create GRange objects of genome annotations from
Gencode hg38. Specifically we will create gene annotation GRanges and
their corresponding promoter region. These objects will be used for
overlaps with ATAC peaks.

``` r
# Loading gencode genome annotation as GRanges 
gencode_gr_human <- rtracklayer::import("/scratch/Shares/rinnclass/MASTER_CLASS/DATA/genomes/Homo_sapiens/Gencode/v38/gencode.v38.annotation.gtf")

# all gene promoters
gencode_genes_human <- gencode_gr_human[gencode_gr_human$type == "gene"] 
human_gene_promoters <- promoters(gencode_genes_human, upstream = 2000, downstream = 2000)

# mrna_genes
human_mrna_genes <- gencode_genes_human[gencode_genes_human$gene_type %in% "protein_coding"]
human_mrna_promoters <- promoters(human_mrna_genes, upstream = 2000, downstream = 2000)

# lincrna_genes
human_lincrna_genes <- gencode_genes_human[gencode_genes_human$gene_type %in% "lincRNA"]
human_lincrna_gene_promoters <- promoters(human_lincrna_genes, upstream = 2000, downstream = 2000)

# This is more reliable than looking at all peaks because these regions in the genome are more stable
```

## 2.1 Compare overlaps of dox and non-dox peaks with gene annotations

Now we will overlap our dox and non-dox unique peaks with genome
annotations (gene promoters). First we will find number of overlaps with
gene promoters and then genes that had changed in RNAseq

``` r
# gr_list of promoters and peaks unique to non_dox condition
human_gr_list_gene_promoter_non_dox_ov <- list(human_gene_promoters = human_gene_promoters, non_dox_peaks = unique_to_non_dox_human)

human_non_dox_gene_promoter_ov <- suppressWarnings(find_common_peaks(human_gr_list_gene_promoter_non_dox_ov))

print("This is how many non-dox_unique peaks overlapped gene promoters")
```

    ## [1] "This is how many non-dox_unique peaks overlapped gene promoters"

``` r
length(human_non_dox_gene_promoter_ov)
```

    ## [1] 4086

``` r
# 4086 genes had an overlap with an ATAC peak at the zero time point but not in dox treated samples

# peaks unique to dox condition overlapped with gene promoters
human_gr_list_gene_promoter_dox_ov <- list(human_gene_promoters = human_gene_promoters, dox_human_peaks = unique_to_dox_human)

human_dox_gene_promoter_ov <- suppressWarnings(find_common_peaks(human_gr_list_gene_promoter_dox_ov))

print(c("This is how many dox_unique peaks overlapped gene promoters", length(human_dox_gene_promoter_ov)))
```

    ## [1] "This is how many dox_unique peaks overlapped gene promoters"
    ## [2] "2"

``` r
# 2 genes are unique to only dox

# Now find same gene_id in RNAseq resuts
# Loading RNAseq results from 06_Differential_expression_analyses
load("/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/RESULTS/HUMAN/human_DESEQ_results.RData")

# filter significant genes from RNAseq "filtered_res_df" to non-dox unique promoter overlaps
# res df and filtered res df are the results we need to use (filtered by padj)
human_sig_rnaseq_atac_non_dox <- human_non_dox_gene_promoter_ov[human_non_dox_gene_promoter_ov$gene_id %in% human_filtered_res_df$gene_id]

print(c("this is how many genes overlap between RNAseq and ATACseq non-dox peaks", length(human_sig_rnaseq_atac_non_dox)))
```

    ## [1] "this is how many genes overlap between RNAseq and ATACseq non-dox peaks"
    ## [2] "208"

``` r
# 208 genes have no peak at zero time point but peaks appeared after dox treatment

print(c("Here are the gene names that overlap", human_sig_rnaseq_atac_non_dox$gene_name))
```

    ##   [1] "Here are the gene names that overlap"
    ##   [2] "SAMD11"                              
    ##   [3] "HES4"                                
    ##   [4] "MXRA8"                               
    ##   [5] "AJAP1"                               
    ##   [6] "HES3"                                
    ##   [7] "MIR34AHG"                            
    ##   [8] "AADACL3"                             
    ##   [9] "JUN"                                 
    ##  [10] "IL23R"                               
    ##  [11] "LPAR3"                               
    ##  [12] "MCOLN3"                              
    ##  [13] "LINC02795"                           
    ##  [14] "LINC02884"                           
    ##  [15] "RP4-590F24.2"                        
    ##  [16] "SRGAP2-AS1"                          
    ##  [17] "RP11-404F10.2"                       
    ##  [18] "RP11-384L19.1"                       
    ##  [19] "RP1-45C12.1"                         
    ##  [20] "CRB1"                                
    ##  [21] "OPTC"                                
    ##  [22] "KISS1"                               
    ##  [23] "RD3"                                 
    ##  [24] "GJC2"                                
    ##  [25] "OTOF"                                
    ##  [26] "CAPN13"                              
    ##  [27] "CEBPZOS"                             
    ##  [28] "RMDN2-AS1"                           
    ##  [29] "CYP1B1-AS1"                          
    ##  [30] "ATP6V1B1"                            
    ##  [31] "LRRTM4-AS1"                          
    ##  [32] "ATOH8"                               
    ##  [33] "IL1R1"                               
    ##  [34] "PTCHD3P2"                            
    ##  [35] "ITGA6-AS1"                           
    ##  [36] "RP11-171I2.1"                        
    ##  [37] "IMPDH1P10"                           
    ##  [38] "NMUR1"                               
    ##  [39] "AC096574.5"                          
    ##  [40] "HRH1"                                
    ##  [41] "DAZL"                                
    ##  [42] "DPPA2"                               
    ##  [43] "RP11-375I20.6"                       
    ##  [44] "RP11-271K21.12"                      
    ##  [45] "TM4SF18"                             
    ##  [46] "CHRD"                                
    ##  [47] "RTP1"                                
    ##  [48] "TM4SF19"                             
    ##  [49] "DOK7"                                
    ##  [50] "RPL17P19"                            
    ##  [51] "SPARCL1"                             
    ##  [52] "RP11-115D19.1"                       
    ##  [53] "RP11-562F9.3"                        
    ##  [54] "CICP16"                              
    ##  [55] "FAT4"                                
    ##  [56] "SORBS2"                              
    ##  [57] "HTR1A"                               
    ##  [58] "LINC01338"                           
    ##  [59] "PITX1"                               
    ##  [60] "WNT8A"                               
    ##  [61] "CTB-180C19.1"                        
    ##  [62] "CTB-114C7.4"                         
    ##  [63] "LINC00518"                           
    ##  [64] "EDN1"                                
    ##  [65] "RIPOR2"                              
    ##  [66] "IP6K3"                               
    ##  [67] "RP3-462C17.1"                        
    ##  [68] "FILIP1"                              
    ##  [69] "RP11-474L11.5"                       
    ##  [70] "TRDN"                                
    ##  [71] "HEY2"                                
    ##  [72] "RP11-69I8.2"                         
    ##  [73] "CCN2"                                
    ##  [74] "ALDH8A1"                             
    ##  [75] "OLIG3"                               
    ##  [76] "TXLNB"                               
    ##  [77] "RP11-30B1.1"                         
    ##  [78] "RP11-32P3.1"                         
    ##  [79] "PPP1R17"                             
    ##  [80] "GCK"                                 
    ##  [81] "CALN1"                               
    ##  [82] "VGF"                                 
    ##  [83] "LRRC17"                              
    ##  [84] "EFCAB10"                             
    ##  [85] "RP11-274B21.12"                      
    ##  [86] "CICP14"                              
    ##  [87] "DOK2"                                
    ##  [88] "EGR3"                                
    ##  [89] "IDO1"                                
    ##  [90] "RP11-697N18.3"                       
    ##  [91] "PRDM14"                              
    ##  [92] "RP11-531A24.5"                       
    ##  [93] "CTD-2501M5.1"                        
    ##  [94] "RP11-240B13.2"                       
    ##  [95] "LY6L"                                
    ##  [96] "ALDH1A1"                             
    ##  [97] "ANXA1"                               
    ##  [98] "LMX1B"                               
    ##  [99] "DBH"                                 
    ## [100] "DKFZP434A062"                        
    ## [101] "CLIC3"                               
    ## [102] "AKR1C1"                              
    ## [103] "AKR1C2"                              
    ## [104] "LINP1"                               
    ## [105] "NEUROG3"                             
    ## [106] "LRMDA"                               
    ## [107] "FAM24B"                              
    ## [108] "FAM53B-AS1"                          
    ## [109] "C10orf90"                            
    ## [110] "RP11-288A5.2"                        
    ## [111] "CDKN1C"                              
    ## [112] "TSSC2"                               
    ## [113] "TPH1"                                
    ## [114] "RP11-1110F20.1"                      
    ## [115] "KCNK7"                               
    ## [116] "EFEMP2"                              
    ## [117] "RP11-826F13.1"                       
    ## [118] "GRIA4"                               
    ## [119] "IL18"                                
    ## [120] "C11orf45"                            
    ## [121] "RP11-861E23.2"                       
    ## [122] "RP11-234B24.4"                       
    ## [123] "LAG3"                                
    ## [124] "DPPA3"                               
    ## [125] "PHC1"                                
    ## [126] "LMO3"                                
    ## [127] "YAF2"                                
    ## [128] "RP11-793H13.11"                      
    ## [129] "LINC01490"                           
    ## [130] "TMEM119"                             
    ## [131] "TBX5-AS1"                            
    ## [132] "TMEM132B"                            
    ## [133] "TNFRSF19"                            
    ## [134] "RLIMP1"                              
    ## [135] "RP11-457D13.4"                       
    ## [136] "SOX21-AS1"                           
    ## [137] "SOX21"                               
    ## [138] "KARS1P2"                             
    ## [139] "LMLN2"                               
    ## [140] "RPL10L"                              
    ## [141] "RP11-1046B16.4"                      
    ## [142] "PROX2"                               
    ## [143] "ESRRB"                               
    ## [144] "RP11-7F17.8"                         
    ## [145] "RP11-497E19.1"                       
    ## [146] "SNURF"                               
    ## [147] "RP11-133K1.8"                        
    ## [148] "CYP1A1"                              
    ## [149] "GOLGA6L4"                            
    ## [150] "PRR35"                               
    ## [151] "RP11-20I23.10"                       
    ## [152] "SULT1A2"                             
    ## [153] "PYDC1"                               
    ## [154] "ADGRG1"                              
    ## [155] "CDH5"                                
    ## [156] "CMTM1"                               
    ## [157] "NECAB2"                              
    ## [158] "LINC02188"                           
    ## [159] "TRARG1"                              
    ## [160] "RP1-59D14.1"                         
    ## [161] "RP11-81A22.4"                        
    ## [162] "RP11-47L3.1"                         
    ## [163] "HNF1B"                               
    ## [164] "PNMT"                                
    ## [165] "KRT222"                              
    ## [166] "RP11-304F15.3"                       
    ## [167] "RP11-304F15.4"                       
    ## [168] "CSH2"                                
    ## [169] "ICAM2"                               
    ## [170] "AC144831.1"                          
    ## [171] "RP11-621L6.2"                        
    ## [172] "LINC01901"                           
    ## [173] "SYT4"                                
    ## [174] "RP11-9H20.2"                         
    ## [175] "ALPK2"                               
    ## [176] "TSHZ1"                               
    ## [177] "LLNLR-304A6.2"                       
    ## [178] "GNA15"                               
    ## [179] "EBI3"                                
    ## [180] "F2RL3"                               
    ## [181] "CTC-439O9.1"                         
    ## [182] "MAG"                                 
    ## [183] "UPK1A"                               
    ## [184] "KLK6"                                
    ## [185] "KLK7"                                
    ## [186] "NLRP12"                              
    ## [187] "PDYN"                                
    ## [188] "COX4I2"                              
    ## [189] "SLC32A1"                             
    ## [190] "CTCFL"                               
    ## [191] "DSCAM"                               
    ## [192] "KB-68A7.2"                           
    ## [193] "S100B"                               
    ## [194] "CRYBB1"                              
    ## [195] "CABP7"                               
    ## [196] "RFPL3"                               
    ## [197] "MB"                                  
    ## [198] "LGALS1"                              
    ## [199] "RP3-508I15.21"                       
    ## [200] "DGKK"                                
    ## [201] "IQSEC2"                              
    ## [202] "RP3-326L13.2"                        
    ## [203] "SERTM2"                              
    ## [204] "LONRF3"                              
    ## [205] "XPNPEP2"                             
    ## [206] "BCORL1"                              
    ## [207] "PNCK"                                
    ## [208] "SRY"                                 
    ## [209] "RP11-400O10.1"

``` r
# filter significant genes from RNAseq "filtered_res_df" to non-dox unique promoter overlaps
human_sig_rnaseq_atac_dox <- human_non_dox_gene_promoter_ov[human_dox_gene_promoter_ov$gene_id %in% human_filtered_res_df$gene_id]

print(c("this is how many genes overlap between RNAseq and ATACseq non-dox peaks", length(human_sig_rnaseq_atac_dox)))
```

    ## [1] "this is how many genes overlap between RNAseq and ATACseq non-dox peaks"
    ## [2] "2043"

``` r
# 2043 genes had a peak at zero time point but peaks disappeared after dox treatment
```

print(c(“Number of genes overlap between RNAseq and ATACseq non-dox
peaks:”,length(human_sig_rnaseq_atac_dox)))

### WARNING: Visual inspection of peaks does not find this approach convincing

Importing bigWig files and peak calls showed very little convincing data
that peaks were changing in dox and non-dox conditions. Does show good
overlaps with gene annotations so approach is working well.

## 3.1 Use DESEQ2 to find out if peaks are changing in dox and non-dox conditions

The overlap analysis did not contain statistical analysis and is a
logical approach that ended up not being that compelling based on raw
data. So now we will use a statistical approach DESEQ2 to compare peak
read counts across samples of dox and non-dox conditions. To do so we
will use consensus peaks (any peak called in any condition) and feature
counts of each consensus peak to be used as input into DESEQ2 for
differential expression of peak counts.

``` r
# consensus peaks
human_broad_consensus_peaks <- read.table("/scratch/Shares/rinnclass/MASTER_CLASS/DATA/human_atacseq/consensus_peaks.mLb.clN.annotatePeaks.txt",
                             sep = "\t", header = TRUE)

# consensus peak counts
human_broad_consensus_counts <- read.table("/scratch/Shares/rinnclass/MASTER_CLASS/DATA/human_atacseq/consensus_peaks.mLb.clN.featureCounts.txt",
                             sep = "\t", header = TRUE)

# creating sample sheet "atac_samples" from file names (col names of consensus):
rownames(human_broad_consensus_counts) <- human_broad_consensus_counts$Geneid

human_broad_consensus_counts <- human_broad_consensus_counts %>%
  dplyr::select(-c(Geneid, Chr, Start, End, Strand, Length))

# cleaning up column names
colnames(human_broad_consensus_counts) <- gsub(
  "\\.mLb\\.clN\\.sorted\\.bam", 
  "", 
  colnames(human_broad_consensus_counts)
)

human_count_columns <- colnames(human_broad_consensus_counts)

# creating the sample sheet for DESEQ by definition non-dox and dox conditions, and switching hours to minutes
human_atac_samples <- data.frame(
  sample = human_count_columns,
  condition = ifelse(
    grepl("_0h_", human_count_columns),
    "non-dox",
    ifelse(grepl("non-dox", human_count_columns), "non-dox", "dox")
  ),
  timepoint_minutes = as.numeric(sub(".*_(\\d+(\\.\\d+)?)h_.*", "\\1", human_count_columns)) * 60
)

# human_atac_samples is a sample sheet for DESEQ now!

# Factor condition for DESEQ2 !!
human_atac_samples <- human_atac_samples %>%
  mutate(condition = factor(condition, levels = c("non-dox", "dox")))

# matrix for Deseq2
human_atac_dds_condition <- suppressWarnings(DESeqDataSetFromMatrix(countData = human_broad_consensus_counts, 
                                   colData = human_atac_samples, 
                                   design = ~ condition))
```

    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not a warning or an error]

``` r
# Run DESeq2 condition model
human_atac_dds_condition <- suppressWarnings(DESeq(human_atac_dds_condition))
```

    ## estimating size factors
    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not a warning or an error]

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not a warning or an error]

    ## final dispersion estimates

    ## fitting model and testing

    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not a warning or an error]

## 3.2 Extracting and analyzing DESEQ2 results

Creating a dataframe with p-vlaues and lfc etc

``` r
# Extract DESeq2 results
# print the results from DESEQ which is the intercept file in the comparison and make it a data frame
human_atac_lfc_condition <- results(human_atac_dds_condition) %>%
  as.data.frame() %>%
  rownames_to_column("interval_id")

# Merge with broad_consensus_peaks info

colnames(human_broad_consensus_peaks)[1] <- "interval_id"

human_atac_lfc_condition <- merge(human_atac_lfc_condition, 
                  human_broad_consensus_peaks %>%
                    dplyr::select(interval_id, Gene.Name, Nearest.PromoterID, 
                                  Distance.to.TSS, Chr, Start, End),
                  by = "interval_id")
```

## 3.3 Analyzing results of DESEQ on ATACseq peak counts

Now we will see how many peaks had a padj \< 0.05!

``` r
# padj histogram
ggplot(human_atac_lfc_condition, aes(x = padj)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "white") +
  labs(title = "Distribution of Human DESEQ2/ATAC padj Values",
       x = "Adjusted p-value",
       y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
```

![](humanATACSEQ_files/figure-gfm/analysis%20of%20DESEQ2%20on%20ATAC%20peak%20counts-1.png)<!-- -->

``` r
# p value histogram
ggplot(human_atac_lfc_condition, aes(x = pvalue)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "white") +
  labs(title = "Distribution of Human DESEQ2/ATAC p Values",
       x = "p-value",
       y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
```

![](humanATACSEQ_files/figure-gfm/analysis%20of%20DESEQ2%20on%20ATAC%20peak%20counts-2.png)<!-- -->

![histogram of padj
values](/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/reanna/MASTER_CLASS/finalproject2/images/humandeseq2atacpadj.png)
\# no padj values at or under 0.05 :(

![histogram of p
values](/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/reanna/MASTER_CLASS/finalproject2/images/humandeseq2atacpvalues.png)
\# some p values are significant but no padj values are…

### Result: No peaks are significant by padj value.

### Result: Doxycyline does not significantly affect chromatin accessibility.

#### Conclusion: There are ~60k peaks common across all samples, and those peaks did not signficantly change from doxycyline exposure. However, RNA gene expression changed for 27 genes, so doxycyline could be influencing post transcriptional or transcriptional changes but not epigenetic.

## Bonus question: Are the RNA regulatory events (e.g. differential expression) biased toward genes with open chromatin? Will be answered soon…
