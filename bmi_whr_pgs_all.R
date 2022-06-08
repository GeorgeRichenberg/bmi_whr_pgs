### Setting working directory and loading required packages:
setwd("/Users/cs21140/Downloads/main_thesis")
library(vroom)
library(tidyverse)
library(TwoSampleMR)

### Uploading bmi data:
bmi <- vroom(file = "Meta-analysis_Locke_et_al+UKBiobank_2018_top_941_from_COJO_analysis_UPDATED.txt")

### Formatting and clumping bmi data:
dim(bmi)
colnames(bmi)
bmi_selected <- bmi %>%
  select(SNP, CHR, POS, Tested_Allele, Other_Allele, Freq_Tested_Allele_in_HRS, BETA, SE, P, N)
bmi_arranged <- bmi_selected %>% 
  arrange(CHR, POS)
bmi_renamed <- bmi_arranged %>%
  rename(chr_name = CHR,
         chrom_start = POS,
         Effect_allele = Tested_Allele,
         Other_allele = Other_Allele,
         Freq_effect_allele = Freq_Tested_Allele_in_HRS,
         pval.exposure = P)
bmi_clumped <- clump_data(bmi_renamed,
                          clump_kb = 10000,
                          clump_r2 = 0.001,
                          clump_p1 = 1,
                          clump_p2 = 1)
dim(bmi_clumped)
write.csv(bmi_renamed, file = "bmi_all_unclumped_pgs(n=941).csv", row.names = FALSE, quote = FALSE)
write.csv(bmi_clumped, file = "bmi_all_clumped_pgs(n=333).csv", row.names = FALSE, quote = FALSE)

### Uploading whr data:
whr_primary <- vroom(file = "whradjbmi.giant-ukbb.meta.1.merged.indexSnps.combined.parsed.txt")
whr_secondary <- vroom(file = "whradjbmi.giant-ukbb.meta.1.merged.secondarySnps.combined.parsed.txt")

### Combining whr datasets:
dim(whr_primary)
dim(whr_secondary)
whr <- rbind(whr_primary, whr_secondary)
sum(duplicated(whr))
dim(whr)

### Formatting and clumping whr data:
whr_selected <- whr %>%
  select(SNP, Chr.ref.males, Pos.ref.males, A1.combined, A2.combined, frqA1.combined, beta.combined, se.combined, pval.combined, nmeta.combined)
whr_arranged <- whr_selected %>% 
  arrange(Chr.ref.males, Pos.ref.males)
whr_arranged$SNP <- gsub(":.*", "", whr_arranged$SNP)
whr_renamed <- whr_arranged %>%
  rename(chr_name = Chr.ref.males,
         chrom_start = Pos.ref.males,
         Effect_allele = A1.combined,
         Other_allele = A2.combined,
         Freq_effect_allele = frqA1.combined,
         BETA = beta.combined,
         SE = se.combined,
         pval.exposure = pval.combined,
         N = nmeta.combined)
dim(whr_renamed)
whr_clumped <- clump_data(whr_renamed,
                          clump_kb = 10000,
                          clump_r2 = 0.001,
                          clump_p1 = 1,
                          clump_p2 = 1)
dim(whr_clumped)
write.csv(whr_renamed, file = "whr_all_unclumped_pgs(n=463).csv", row.names = FALSE, quote = FALSE)
write.csv(whr_clumped, file = "whr_all_clumped_pgs(n=217).csv", row.names = FALSE, quote = FALSE)