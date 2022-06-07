### Set working directory and load required packages:
setwd("/Users/cs21140/Downloads/bmi_cancer")
library(vroom)
library(tidyverse)
library(TwoSampleMR)

### Uploading bmi data:
bmi_females_primary <- vroom(file = "bmi.giant-ukbb.meta.1.merged.indexSnps.females.parsed.txt")
bmi_females_secondary <- vroom(file = "bmi.giant-ukbb.meta.1.merged.secondarySnps.females.parsed.txt")
bmi_males_primary <- vroom(file = "bmi.giant-ukbb.meta.1.merged.indexSnps.males.parsed.txt")
bmi_males_secondary <- vroom(file = "bmi.giant-ukbb.meta.1.merged.secondarySnps.males.parsed.txt")

### Combining primary and secondary bmi datasets per sex:
bmi_females <- rbind(bmi_females_primary, bmi_females_secondary)
sum(duplicated(bmi_females))
dim(bmi_females)
bmi_males <- rbind(bmi_males_primary, bmi_males_secondary)
sum(duplicated(bmi_males))
dim(bmi_males)

### Formatting and clumping bmi_females data:
bmi_females_selected <- bmi_females %>%
  select(SNP, Chr.ref.males, Pos.ref.males, A1.combined, A2.combined, frqA1.females, beta.females, se.females, pval.females, nmeta.females)
bmi_females_arranged <- bmi_females_selected %>% 
  arrange(Chr.ref.males, Pos.ref.males)
bmi_females_arranged$SNP <- gsub(":.*", "", bmi_females_arranged$SNP)
bmi_females_renamed <- bmi_females_arranged %>%
  rename(chr_name = Chr.ref.males,
         chrom_start = Pos.ref.males,
         Effect_allele = A1.combined,
         Other_allele = A2.combined,
         Freq_effect_allele = frqA1.females,
         BETA = beta.females,
         SE = se.females,
         pval.exposure = pval.females,
         N = nmeta.females)
dim(bmi_females_renamed)
bmi_females_clumped <- clump_data(bmi_females_renamed,
                                  clump_kb = 10000,
                                  clump_r2 = 0.001,
                                  clump_p1 = 1,
                                  clump_p2 = 1)
dim(bmi_females_clumped)
write.csv(bmi_females_renamed, file = "bmi_females_unclumped_pgs.csv", row.names = FALSE, quote = FALSE)
write.csv(whr_females_clumped, file = "bmi_females_clumped_pgs.csv", row.names = FALSE, quote = FALSE)

### Formatting and clumping bmi_males data:
bmi_males_selected <- bmi_males %>%
  select(SNP, Chr.ref.males, Pos.ref.males, A1.combined, A2.combined, frqA1.males, beta.males, se.males, pval.males, nmeta.males)
bmi_males_arranged <- bmi_males_selected %>% 
  arrange(Chr.ref.males, Pos.ref.males)
bmi_males_arranged$SNP <- gsub(":.*", "", bmi_males_arranged$SNP)
bmi_males_renamed <- bmi_males_arranged %>%
  rename(chr_name = Chr.ref.males,
         chrom_start = Pos.ref.males,
         Effect_allele = A1.combined,
         Other_allele = A2.combined,
         Freq_effect_allele = frqA1.males,
         BETA = beta.males,
         SE = se.males,
         pval.exposure = pval.males,
         N = nmeta.males)
dim(bmi_males_renamed)
bmi_males_clumped <- clump_data(bmi_males_renamed,
                                clump_kb = 10000,
                                clump_r2 = 0.001,
                                clump_p1 = 1,
                                clump_p2 = 1)
dim(bmi_males_clumped)
write.csv(bmi_males_renamed, file = "bmi_males_unclumped_pgs.csv", row.names = FALSE, quote = FALSE)
write.csv(bmi_males_clumped, file = "bmi_males_clumped_pgs.csv", row.names = FALSE, quote = FALSE)

### Uploading whr data:
whr_females_primary <- vroom(file = "whradjbmi.giant-ukbb.meta.1.merged.indexSnps.females.parsed.txt")
whr_females_secondary <- vroom(file = "whradjbmi.giant-ukbb.meta.1.merged.secondarySnps.females.parsed.txt")
whr_males_primary <- vroom(file = "whradjbmi.giant-ukbb.meta.1.merged.indexSnps.males.parsed.txt")
whr_males_secondary <- vroom(file = "whradjbmi.giant-ukbb.meta.1.merged.secondarySnps.males.parsed.txt")

### Combining primary and secondary whr datasets per sex:
whr_females <- rbind(whr_females_primary, whr_females_secondary)
sum(duplicated(whr_females))
dim(whr_females)
whr_males <- rbind(whr_males_primary, whr_males_secondary)
sum(duplicated(whr_males))
dim(whr_males)

### Formatting and clumping whr_females data:
whr_females_selected <- whr_females %>%
  select(SNP, Chr.ref.males, Pos.ref.males, A1.combined, A2.combined, frqA1.females, beta.females, se.females, pval.females, nmeta.females)
whr_females_arranged <- whr_females_selected %>% 
  arrange(Chr.ref.males, Pos.ref.males)
whr_females_arranged$SNP <- gsub(":.*", "", whr_females_arranged$SNP)
whr_females_renamed <- whr_females_arranged %>%
  rename(chr_name = Chr.ref.males,
         chrom_start = Pos.ref.males,
         Effect_allele = A1.combined,
         Other_allele = A2.combined,
         Freq_effect_allele = frqA1.females,
         BETA = beta.females,
         SE = se.females,
         pval.exposure = pval.females,
         N = nmeta.females)
dim(whr_females_renamed)
whr_females_clumped <- clump_data(whr_females_renamed,
                          clump_kb = 10000,
                          clump_r2 = 0.001,
                          clump_p1 = 1,
                          clump_p2 = 1)
dim(whr_females_clumped)
write.csv(whr_females_renamed, file = "whr_females_unclumped_pgs.csv", row.names = FALSE, quote = FALSE)
write.csv(whr_females_clumped, file = "whr_females_clumped_pgs.csv", row.names = FALSE, quote = FALSE)

### Formatting and clumping whr_males data:
whr_males_selected <- whr_males %>%
  select(SNP, Chr.ref.males, Pos.ref.males, A1.combined, A2.combined, frqA1.males, beta.males, se.males, pval.males, nmeta.males)
whr_males_arranged <- whr_males_selected %>% 
  arrange(Chr.ref.males, Pos.ref.males)
whr_males_arranged$SNP <- gsub(":.*", "", whr_males_arranged$SNP)
whr_males_renamed <- whr_males_arranged %>%
  rename(chr_name = Chr.ref.males,
         chrom_start = Pos.ref.males,
         Effect_allele = A1.combined,
         Other_allele = A2.combined,
         Freq_effect_allele = frqA1.males,
         BETA = beta.males,
         SE = se.males,
         pval.exposure = pval.males,
         N = nmeta.males)
dim(whr_males_renamed)
whr_males_clumped <- clump_data(whr_males_renamed,
                                  clump_kb = 10000,
                                  clump_r2 = 0.001,
                                  clump_p1 = 1,
                                  clump_p2 = 1)
dim(whr_males_clumped)
write.csv(whr_males_renamed, file = "whr_males_unclumped_pgs.csv", row.names = FALSE, quote = FALSE)
write.csv(whr_males_clumped, file = "whr_males_clumped_pgs.csv", row.names = FALSE, quote = FALSE)