library(tidyverse)

cr<-read_tsv("~/Projects/genomics/rotifer/PNAS/vcfs/CR.ad.out", col_types = "cnccc") %>% 
  separate(CHROM, into=c("SET", "CHROM"), sep="\\.", extra="merge") %>%
  separate(CR, into=c("CR.1", "CR.2")) %>% 
  separate(MA, into=c("MA.1", "MA.2")) %>%
  separate(MM, into=c("MM.1", "MM.2")) %>%
  mutate(across(CR.1:MM.2, as.numeric)) %>% 
  mutate(CRfreq = CR.2/(CR.1+CR.2),
         MAfreq = MA.2/(MA.1+MA.2),
         MMfreq = MM.2/(MM.1+MM.2))

ma<-read_tsv("~/Projects/genomics/rotifer/PNAS/vcfs/MA.ad.out", col_types = "cnccc") %>% 
  separate(CHROM, into=c("SET", "CHROM"), sep="\\.", extra="merge") %>%
  separate(CR, into=c("CR.1", "CR.2")) %>% 
  separate(MA, into=c("MA.1", "MA.2")) %>%
  separate(MM, into=c("MM.1", "MM.2")) %>%
  mutate(across(CR.1:MM.2, as.numeric)) %>% 
  mutate(CRfreq = CR.2/(CR.1+CR.2),
         MAfreq = MA.2/(MA.1+MA.2),
         MMfreq = MM.2/(MM.1+MM.2))


mm<-read_tsv("~/Projects/genomics/rotifer/PNAS/vcfs/MM.ad.out", col_types = "cnccc") %>% 
  separate(CHROM, into=c("SET", "CHROM"), sep="\\.", extra="merge") %>%
  separate(CR, into=c("CR.1", "CR.2")) %>% 
  separate(MA, into=c("MA.1", "MA.2")) %>%
  separate(MM, into=c("MM.1", "MM.2")) %>%
  mutate(across(CR.1:MM.2, as.numeric)) %>% 
  mutate(CRfreq = CR.2/(CR.1+CR.2),
         MAfreq = MA.2/(MA.1+MA.2),
         MMfreq = MM.2/(MM.1+MM.2))

ma %>% 
  rename(CR=CRfreq, MM=MMfreq, MA=MAfreq) %>%
  pivot_longer(CR:MM, names_to ="isolate", values_to = "alt_freq") %>% 
  ggplot(aes(alt_freq, fill=isolate)) + geom_histogram(bins=50)

ma %>% 
  rename(CR=CRfreq, MM=MMfreq, MA=MAfreq) %>% select(-CR) %>%
  pivot_longer(MA:MM, names_to ="isolate", values_to = "alt_freq") %>% 
  ggplot(aes(alt_freq, fill=isolate)) + geom_histogram(bins=50) 
                                                                                            
cr %>% ggplot(aes(MAfreq)) + geom_histogram(bins=50)
cr %>% ggplot(aes(MMfreq)) + geom_histogram(bins=50)
