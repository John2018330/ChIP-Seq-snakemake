## Rscript for 'Integration with RNA-Seq' Section of project 2

library(tidyverse)
library(scales)

deseq_results <- read_tsv('results/GSE75070_MCF7_shRUNX1_shNS_RNAseq_log2_foldchange.txt')

filtered_deseq <- deseq_results %>%
  filter(abs(log2FoldChange) > 1 & padj < 0.01)

annotated_peaks <- read_tsv('results/RUNX1_noBlacklist_annotated.txt')

fig2f_data <- filtered_deseq %>%
  left_join(y=annotated_peaks, by = c('genename' = 'Gene Name')) %>%
  select(genename, log2FoldChange, padj, 'Distance to TSS') %>%
  filter(!duplicated(genename))

fig2f_5k <- fig2f_data %>%
  mutate(bound_status = case_when(
                          (abs(`Distance to TSS`) >= 5000 | is.na(`Distance to TSS`)) ~ 'Not bound',
                          .default = 'RUNX1 bound'), 
         regulation = case_when(
                          log2FoldChange > 0 ~ '5k-Up-regulated',
                          log2FoldChange < 0 ~ '5k-Down-regulated'))

fig2f_20k <- fig2f_data %>%
  mutate(bound_status = case_when(
                            (abs(`Distance to TSS`) >= 20000 | is.na(`Distance to TSS`)) ~ 'Not bound',
                            .default = 'RUNX1 bound'), 
         regulation = case_when(
                            log2FoldChange > 0 ~ '20k-Up-regulated',
                            log2FoldChange < 0 ~ '20k-Down-regulated'))

fig2f_100k <- fig2f_data %>%
  filter(abs(`Distance to TSS`) <= 100000 | is.na(`Distance to TSS`)) %>%
  mutate(bound_status = case_when(
                            (abs(`Distance to TSS`) >= 100000 | is.na(`Distance to TSS`)) ~ 'Not bound',
                            .default = 'RUNX1 bound'), 
         regulation = case_when(
                            log2FoldChange > 0 ~ '100k-Up-regulated',
                            log2FoldChange < 0 ~ '100k-Down-regulated'))


xlab <- c('5k-Up-regulated', '5k-Down-regulated', '20k-Up-regulated', '20k-Down-regulated',
          '100k-Up-regulated', '100k-Down-regulated')

counts_5k <- fig2f_5k %>%
  group_by(regulation, bound_status) %>%
  count()

counts_20k <- fig2f_20k %>%
  group_by(regulation, bound_status) %>%
  count()

counts_100k <- fig2f_100k %>%
  group_by(regulation, bound_status) %>%
  count()


fig2f <- ggplot(fig2f_5k) +
  geom_bar(aes(fill=bound_status, x=regulation), position='fill', width=0.5) +
  geom_bar(data=fig2f_20k, aes(fill=bound_status, x=regulation), position='fill', width=0.5) +
  geom_bar(data=fig2f_100k, aes(fill=bound_status, x=regulation), position='fill', width=0.5) +
  scale_x_discrete(limits = xlab, guide = guide_axis(n.dodge=2)) +
  theme_classic()

ggsave('2f.png', plot=fig2f, dpi=300)
