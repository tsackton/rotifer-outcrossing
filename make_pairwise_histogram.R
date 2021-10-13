library(tidyverse)
setwd("~/Projects/genomics/rotifer/rotifer-outcrossing/")

get_dist_mat_MACR <- function(pairs, file, elem) {
  read_tsv(str_c("data/raw_distance_matrices/", pairs, "/", file), show_col_types = FALSE) %>% 
  pivot_longer(CR1:MA2) %>% 
  rename_at(vars(names(.)), ~c("seq1", "seq2", "diff")) %>% 
  mutate(element = elem, source = pairs) %>%
  filter(seq1 != seq2) %>%
  rowwise() %>% 
  mutate(comp = paste(sort(unlist(strsplit(paste(seq1, seq2, sep = " "), " "))), collapse = "-")) %>%
  distinct(comp, .keep_all = TRUE)
}

get_dist_mat_MAMM <- function(pairs, file, elem) {
  read_tsv(str_c("data/raw_distance_matrices/", pairs, "/", file), show_col_types = FALSE) %>% 
    pivot_longer(MA1:MM2) %>% 
    rename_at(vars(names(.)), ~c("seq1", "seq2", "diff")) %>% 
    mutate(element = elem, source = pairs) %>%
    filter(seq1 != seq2) %>%
    rowwise() %>% 
    mutate(comp = paste(sort(unlist(strsplit(paste(seq1, seq2, sep = " "), " "))), collapse = "-")) %>%
    distinct(comp, .keep_all = TRUE)
}

get_dist_mat <- function(pairs, file, elem) {
  read_tsv(str_c("data/raw_distance_matrices/", pairs, "/", file), show_col_types = FALSE) %>% 
    pivot_longer(CR1:MM2) %>% 
    rename_at(vars(names(.)), ~c("seq1", "seq2", "diff")) %>% 
    mutate(element = elem, source = pairs) %>%
    filter(seq1 != seq2) %>%
    rowwise() %>% 
    mutate(comp = paste(sort(unlist(strsplit(paste(seq1, seq2, sep = " "), " "))), collapse = "-")) %>%
    distinct(comp, .keep_all = TRUE)
}


pairs<-"MACR"
macr <- tibble(file = list.files(path = str_c("data/raw_distance_matrices/", pairs), pattern = "*.tab")) %>%
  separate(file, into = c("elem"), sep = "N", remove=FALSE, extra="drop") %>%
  rowwise() %>%
  mutate(distmat = list(get_dist_mat_MACR(pairs, file, elem)))

pairs<-"MAMM"
mamm <- tibble(file = list.files(path = str_c("data/raw_distance_matrices/", pairs), pattern = "*.tab")) %>%
  separate(file, into = c("elem"), sep = "N", remove=FALSE, extra="drop") %>%
  rowwise() %>%
  mutate(distmat = list(get_dist_mat_MAMM(pairs, file, elem)))


pairs<-"MMCR"
mmcr<- tibble(file = list.files(path = str_c("data/raw_distance_matrices/", pairs), pattern = "*.tab")) %>%
  separate(file, into = c("elem"), sep = "N", remove=FALSE, extra="drop") %>%
  rowwise() %>%
  mutate(distmat = list(get_dist_mat(pairs, file, elem)))


pairs<-"MAMMCR"
mammcr <- tibble(file = list.files(path = str_c("data/raw_distance_matrices/", pairs), pattern = "*.tab")) %>%
  separate(file, into = c("elem"), sep = "N", remove=FALSE, extra="drop") %>%
  rowwise() %>%
  mutate(distmat = list(get_dist_mat(pairs, file, elem)))

allpairs<-rbind(macr, mamm, mmcr, mammcr) %>% unnest(cols=c("distmat")) %>% mutate(file = str_remove(file, "_dist.tab"))

### get final sets ##

rtb<-read_csv("data/combining_3_and_2_ways-final-between.csv", col_names=c("nanopore", "region_3", "length_3", "MAMM_3", "MACR_3", "MMCR_3", "region_MAMM", "length_MAMM", "MAMM_2", "region_MACR", "length_MACR", "MACR_2", "region_MMCR", "length_MMCR", "MMCR_2"), skip=1)

#make nanopore <-> region filter
region_key <- rtb %>% select(nanopore, region_3, region_MAMM, region_MMCR, region_MACR)


#clean up to get final set

rtb_final <- rtb %>% mutate(num_na = as.numeric(paste0(as.numeric(!is.na(length_3)),
                                                       as.numeric(!is.na(length_MAMM)),
                                                       as.numeric(!is.na(length_MACR)),
                                                       as.numeric(!is.na(length_MMCR))))) %>% 
  mutate(length = case_when(num_na >= 1000 ~ length_3,
                            num_na == 1 ~ length_MMCR,
                            num_na == 10 ~ length_MACR,
                            num_na == 100 ~ length_MAMM)) %>%
  filter(!is.na(length)) %>% 
  mutate(MAMM = case_when(num_na >= 1000 ~ MAMM_3,
                          num_na == 100 ~ MAMM_2)) %>%
  mutate(MACR = case_when(num_na >= 1000 ~ MACR_3,
                          num_na == 10 ~ MACR_2)) %>%
  mutate(MMCR = case_when(num_na >= 1000 ~ MMCR_3,
                          num_na == 1 ~ MMCR_2)) %>%
  mutate(comparison = case_when(num_na >= 1000 ~ "MAMMCR",
                                num_na == 100 ~ "MAMM",
                                num_na == 10 ~ "MACR",
                                num_na == 1 ~ "MMCR")) %>%
  mutate(region = as.character(case_when(num_na >= 1000 ~ region_3,
                            num_na == 100 ~ region_MAMM,
                            num_na == 10 ~ region_MACR,
                            num_na == 1 ~ region_MMCR))) %>%
  select(nanopore, length, MAMM, MACR, MMCR, comparison, region)

#now filter/process allpairs

pairdiff <- left_join(rtb_final, allpairs, by=c("region" = "elem", "comparison" = "source")) %>% 
  select(comparison, region, length, seq1, seq2, comp, diff) %>%
  mutate(perdiff = diff/length)

pairdiff %>% ggplot(aes(x=perdiff)) + geom_histogram(bins=50)

#unfiltered histogram

#first get lengths

macr_len <- read_tsv("data/aligns/MA_CR/MACR_seq_lengths_gb_trim.txt", col_names = c("file", "allele", "length")) %>% 
  select(file, length) %>% distinct() %>% mutate(source="MACR", file = str_remove(file, ".fasta"))
mamm_len <- read_tsv("data/aligns/MA_MM/MAMM_all_seq_lengths_gb_trim.txt", col_names = c("file", "allele", "length")) %>% 
  select(file, length) %>% distinct() %>% mutate(source="MAMM", file = str_remove(file, ".fasta"))
mmcr_len <- read_tsv("data/aligns/MM_CR/MMCR_all_seq_lengths_gb_trim.txt", col_names = c("file", "allele", "length")) %>% 
  select(file, length) %>% distinct() %>% mutate(source="MMCR", file = str_remove(file, ".fasta"))
mammcr_len <- read_tsv("data/aligns/three-way/all_seq_lengths_gb.txt", col_names = c("file", "allele", "length")) %>% 
  select(file, length) %>% distinct() %>% mutate(source="MAMMCR", file = str_remove(file, ".fasta"))

lengths<-rbind(macr_len, mamm_len, mmcr_len, mammcr_len)

pairdiff_uf <- left_join(allpairs,lengths, by=c("file" = "file", "source" = "source")) %>%
  select(source, elem, length, seq1, seq2, comp, diff) %>%
  mutate(perdiff = diff/length)

pairdiff_uf %>% ggplot(aes(x=perdiff)) + geom_histogram(bins=100)

