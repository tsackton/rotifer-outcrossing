library(tidyverse)
setwd("~/Projects/genomics/rotifer/manuscript/")

rtb<-read_csv("~/Projects/genomics/rotifer/manuscript/analysis/combining_3_and_2_ways-final.csv", col_names=c("nanopore", "region_3", "length_3", "MAMM_3", "MACR_3", "MMCR_3", "region_MAMM", "length_MAMM", "MAMM_2", "region_MACR", "length_MACR", "MACR_2", "region_MMCR", "length_MMCR", "MMCR_2"), skip=1)

rtw<-read_csv("~/Projects/genomics/rotifer/manuscript/analysis/combining_3_and_2_ways-final-within.csv", col_names=c("nanopore", "region_3", "length_3", "MA_3way", "CR_3way", "MM_3way", "region_MAMM", "length_MAMM", "MA_2", "MM_2", "region_MACR", "length_MACR", "MA_3", "CR_3", "region_MMCR", "length_MMCR", "MM_4", "CR_4"), skip=1)

#make nanopore <-> region filter
region_key <- rtb %>% select(nanopore, region_3, region_MAMM, region_MMCR, region_MACR)


#### BETWEEN ###

#some data checks
rtb %>% mutate(num_na = paste0(as.numeric(!is.na(length_3)),
                         as.numeric(!is.na(length_MAMM)),
                         as.numeric(!is.na(length_MACR)),
                         as.numeric(!is.na(length_MMCR)))) %>%
  with(., table(num_na))


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
  mutate(comparison = case_when(num_na >= 1000 ~ "MA/MM/CR",
                                num_na == 100 ~ "MA/MM",
                                num_na == 10 ~ "MA/CR",
                                num_na == 1 ~ "MM/CR")) %>%
  select(nanopore, length, MAMM, MACR, MMCR, comparison)

#check

table(rtb_final$comparison)

#### NEW PLOTS 6/16/2020 ##

rtb_final %>% filter(MACR <= 50) %>% count(MACR) %>% 
  full_join(tibble(MACR=seq(0,50,1), count=0)) %>% 
  arrange(MACR) %>% mutate(n = ifelse(is.na(n), 0, n)) %>%
  ggplot(aes(y=n, x=MACR)) + geom_point() + scale_y_continuous(name="Number of regions", breaks=c(-1,0,1,2,3,4,5,6,7,8,9,10,11,12)) + scale_x_continuous(name="Number of differences (MA vs CR)", breaks=seq(0,50,2))

rtb_final %>% filter(MMCR <= 50) %>% count(MMCR) %>% 
  full_join(tibble(MMCR=seq(0,50,1), count=0)) %>% 
  arrange(MMCR) %>% mutate(n = ifelse(is.na(n), 0, n)) %>%
  ggplot(aes(y=n, x=MMCR)) + geom_point() + scale_y_continuous(name="Number of regions", breaks=seq(0, 18, 1)) + scale_x_continuous(name="Number of differences (MM vs CR)", breaks=seq(0,50,2))

rtb_final %>% filter(MAMM <= 50, MAMM != 0) %>% count(MAMM) %>% 
  full_join(tibble(MAMM=seq(0,50,1), count=0)) %>% 
  arrange(MAMM) %>% mutate(n = ifelse(is.na(n), 0, n)) %>%
  ggplot(aes(y=n, x=MAMM)) + geom_point() + scale_y_continuous(name="Number of regions", breaks=seq(0, 4, 1)) + scale_x_continuous(name="Number of differences (MA vs MM)", breaks=seq(0,50,2))

#percent, rounded

rtb_final %>% filter(MACR <= 50) %>% mutate(norm_diff = round((MACR/length)*10000)) %>% filter(norm_diff <= 50) %>% count(norm_diff) %>% 
  full_join(tibble(norm_diff=seq(0,50,1), count=0)) %>% 
  arrange(norm_diff) %>% mutate(n = ifelse(is.na(n), 0, n)) %>%
  ggplot(aes(y=n, x=norm_diff)) + geom_point() + scale_y_continuous(name="Number of regions", breaks=c(-1,0,1,2,3,4,5,6,7,8,9,10,11,12)) + scale_x_continuous(name="Number of differences (MA vs CR)", breaks=seq(0,50,2))


rtb_final %>% filter(MMCR <= 50) %>% mutate(norm_diff = round((MMCR/length)*10000)) %>% filter(norm_diff <= 50) %>% count(norm_diff) %>% 
  full_join(tibble(norm_diff=seq(0,50,1), count=0)) %>% 
  arrange(norm_diff) %>% mutate(n = ifelse(is.na(n), 0, n)) %>%
  ggplot(aes(y=n, x=norm_diff)) + geom_point() + scale_y_continuous(name="Number of regions", breaks=seq(0, 18, 1)) + scale_x_continuous(name="Number of differences (MM vs CR)", breaks=seq(0,50,2))


rtb_final %>% filter(MAMM <= 50) %>% mutate(norm_diff = round((MAMM/length)*10000)) %>% filter(norm_diff <= 50) %>% count(norm_diff) %>% 
  full_join(tibble(norm_diff=seq(0,50,1), count=0)) %>% 
  arrange(norm_diff) %>% mutate(n = ifelse(is.na(n), 0, n)) %>%
  ggplot(aes(y=n, x=norm_diff)) + geom_point() + scale_y_continuous(name="Number of regions", breaks=seq(0, 5, 1)) + scale_x_continuous(name="Number of differences (MA vs CR)", breaks=seq(0,50,2))

#calculating age

rtb_final %>% filter(MAMM <= 50) %>% mutate(norm_diff = round((MAMM/length)*10000)) %>% filter(norm_diff <= 1) %>% summarize(length = sum(length), diffs = sum(MAMM)) %>%
  mutate(perbp = diffs/length, age = diffs/(3.9e-9*length))

rtb_final %>% filter(MAMM <= 50) %>% mutate(norm_diff = round((MAMM/length)*10000)) %>% filter(norm_diff <= 2) %>% summarize(length = sum(length), diffs = sum(MAMM)) %>%
  mutate(perbp = diffs/length, age = diffs/(3.9e-9*length))

rtb_final %>% filter(MAMM <= 50) %>% mutate(norm_diff = round((MAMM/length)*10000)) %>% filter(norm_diff <= 3) %>% summarize(length = sum(length), diffs = sum(MAMM)) %>%
  mutate(perbp = diffs/length, age = diffs/(3.9e-9*length))

rtb_final %>% filter(MAMM <= 50) %>% mutate(norm_diff = round((MAMM/length)*10000)) %>% filter(norm_diff <= 4) %>% summarize(length = sum(length), diffs = sum(MAMM)) %>%
  mutate(perbp = diffs/length, age = diffs/(3.9e-9*length))

rtb_final %>% filter(MAMM <= 50) %>% mutate(norm_diff = round((MAMM/length)*10000)) %>% filter(norm_diff <= 5) %>% summarize(length = sum(length), diffs = sum(MAMM)) %>%
  mutate(perbp = diffs/length, age = diffs/(3.9e-9*length))

rtb_final %>% filter(MAMM <= 50) %>% mutate(norm_diff = round((MAMM/length)*10000)) %>% filter(norm_diff <= 6) %>% summarize(length = sum(length), diffs = sum(MAMM)) %>%
  mutate(perbp = diffs/length, age = diffs/(3.9e-9*length))

#alternate approach
3.9e-9*10000

sim_res_1<-numeric(100)
sim_res_2<-numeric(100)
sim_res_3<-numeric(100)

for (i in 1:100) {
  gens=200
  lambda=(gens*3.9e-5)
  n <- rpois(311, lambda) %>% as_tibble() %>% filter(value > 0) %>% count(value)
  sim_res_1[i] <- ifelse(length(n$n[n$value==1]) > 0, n$n[n$value==1], 0)
  sim_res_2[i] <- ifelse(length(n$n[n$value==2]) > 0, n$n[n$value==2], 0)
  sim_res_3[i] <- ifelse(length(n$n[n$value==3]) > 0, n$n[n$value==3], 0)
}

mean(sim_res_1)
mean(sim_res_2)
mean(sim_res_3)


rtb_final %>% filter(MAMM <= 50) %>% mutate(norm_diff = round((MAMM/length)*10000)) %>% filter(norm_diff <= 7) %>% summarize(length = sum(length), diffs = sum(MAMM)) %>%
  mutate(perbp = diffs/length, age = diffs/(3.9e-9*length))


#### WITHIN ###

#some data checks
rtw %>% mutate(num_na = paste0(as.numeric(!is.na(length_3)),
                               as.numeric(!is.na(length_MAMM)),
                               as.numeric(!is.na(length_MACR)),
                               as.numeric(!is.na(length_MMCR)))) %>%
  with(., table(num_na))


#clean up to get final set

rtw_final <- rtw %>% mutate(num_na = as.numeric(paste0(as.numeric(!is.na(length_3)),
                                                       as.numeric(!is.na(length_MAMM)),
                                                       as.numeric(!is.na(length_MACR)),
                                                       as.numeric(!is.na(length_MMCR))))) %>% 
  mutate(length = case_when(num_na >= 1000 ~ length_3,
                            num_na == 1 ~ length_MMCR,
                            num_na == 10 ~ length_MACR,
                            num_na == 100 ~ length_MAMM)) %>%
  filter(!is.na(length)) %>% 
  mutate(MM = case_when(num_na >= 1000 ~ MM_3way,
                        num_na == 100 ~ MM_2,
                        num_na == 1 ~ MM_4)) %>%
  mutate(CR = case_when(num_na >= 1000 ~ CR_3way,
                        num_na == 10 ~ CR_3,
                        num_na == 1 ~ CR_4)) %>%
  mutate(MA = case_when(num_na >= 1000 ~ MA_3way,
                        num_na == 10 ~ MA_3,
                        num_na == 100 ~ MA_2)) %>%
  mutate(comparison = case_when(num_na >= 1000 ~ "MA/MM/CR",
                                num_na == 100 ~ "MA/MM",
                                num_na == 10 ~ "MA/CR",
                                num_na == 1 ~ "MM/CR")) %>%
  select(nanopore, length, MM, MA, CR, comparison)

#check

table(rtw_final$comparison)

#### END DATA PREP ###

##ADD CODE TO GET NEAR-PERFECT SHARERS

rtb_final %>% inner_join(region_key, by=c("nanopore" = "nanopore")) %>% 
  filter(MMCR/length <= 0.0014969) %>% write_csv(path="mmcr_near_sharers.csv")

rtb_final %>% inner_join(region_key, by=c("nanopore" = "nanopore")) %>% 
  filter(MAMM/length <= 0.0014969) %>% write_csv(path="mamm_near_sharers.csv")

rtb_final %>% inner_join(region_key, by=c("nanopore" = "nanopore")) %>% 
  filter(MACR/length <= 0.0014969) %>% write_csv(path="macr_near_sharers.csv")

rtw_final %>% inner_join(region_key, by=c("nanopore" = "nanopore")) %>% 
  filter(MM/length <= 0.0014969) %>% write_csv(path="mm_near_sharers.csv")

rtw_final %>% inner_join(region_key, by=c("nanopore" = "nanopore")) %>% 
  filter(MA/length <= 0.0014969) %>% write_csv(path="ma_near_sharers.csv")

rtw_final %>% inner_join(region_key, by=c("nanopore" = "nanopore")) %>% 
  filter(CR/length <= 0.0014969) %>% write_csv(path="cr_near_sharers.csv")


#number of each type
sum(!is.na(rtb_final$MAMM))
sum(!is.na(rtb_final$MACR))
sum(!is.na(rtb_final$MMCR))

#length
min(rtb_final$length)
max(rtb_final$length)
mean(rtb_final$length)

#looking at allele sharing
table(rtb_final$MAMM == 0)
table(rtb_final$MAMM == 0, rtb_final$comparison)

#conversion
table(rtw_final$MM == 0)
table(rtw_final$MA == 0)
table(rtw_final$CR == 0)

rtw_final %>% filter(MM == 0) %>% write_csv("MM_sharing.csv")
rtw_final %>% filter(MA == 0) %>% write_csv("MA_sharing.csv")
rtw_final %>% filter(CR == 0) %>% write_csv("CR_sharing.csv")

rtw_final %>% filter(MM <= 10) %>% write_csv("MM_near_sharing.csv")
rtw_final %>% filter(MA <= 10) %>% write_csv("MA_near_sharing.csv")
rtw_final %>% filter(CR <= 10) %>% write_csv("CR_near_sharing.csv")

rtw_final %>% filter(MM != 0 | MA != 0) %>% left_join(rtb_final, by=c("nanopore" = "nanopore")) %>% group_by(MAMM == 0) %>% tally()

conv_rate = 3.3e-5
conv_rate2 = 4.8e-5
bp_total = rtw_final %>% filter(!is.na(MM)) %>% summarize(total = sum(length))
bp_conv = rtw_final %>% filter(MM == 0) %>% summarize(total = sum(length)) 
bp_conv / (conv_rate*bp_total)
bp_conv / (conv_rate2*bp_total)


bp_total = rtw_final %>% filter(!is.na(MA)) %>% summarize(total = sum(length))
bp_conv = rtw_final %>% filter(MA == 0) %>% summarize(total = sum(length))
bp_conv / (conv_rate*bp_total)
bp_conv / (conv_rate2*bp_total)

#using non-perfect sharers to get age

cutoff = 2
perfect_len = sum(rtb_final %>% filter(MAMM == 0) %>% pull(length))
non_perf_len = sum(rtb_final %>% filter(MAMM > 0 & MAMM <= cutoff) %>% pull(length))
total_len = perfect_len + non_perf_len
subs = sum(rtb_final %>% filter(MAMM > 0 & MAMM <= cutoff) %>% pull(MAMM))
mut_rate_gen = 3e-9
gens=subs/(mut_rate_gen*total_len)

#total elements 
tot_elem <- 622
sharers <- rtb_final %>% filter(MAMM <= cutoff) %>% select(length) %>% mutate(diffs = length*mut_rate_gen*gens)
summarize(sharers, (311-sum(diffs))/tot_elem)

#frac
sharers/tot_elem

#expected muts per 12kb
(subs/(mut_rate_gen*total_len)) * 10000 * mut_rate_gen * 622
length(rtb_final %>% filter(MAMM > 0 & MAMM <= cutoff) %>% pull(length))

#mtdna
2.5e-5/1.5e-7


#mean MA/MM outside perfect regions
rtb_final %>% filter(MAMM > 0 & !is.na(MAMM)) %>%
  mutate(diff = MAMM/length) %>%
  summarize(mean = mean(diff), sd = sqrt(var(diff)))

rtb_final %>% filter(MAMM > 3 & !is.na(MAMM)) %>%
  mutate(diff = MAMM/length) %>%
  summarize(mean = mean(diff), sd = sqrt(var(diff)))

#CR sharing

table(rtb_final$MACR == 0) 
table(rtb_final$MMCR == 0)

rtb_final %>% filter(MACR > 0 & !is.na(MACR)) %>%
  mutate(diff = MACR/length) %>%
  summarize(mean = mean(diff), sd = sqrt(var(diff)))

rtb_final %>% filter(MMCR > 0 & !is.na(MMCR)) %>%
  mutate(diff = MMCR/length) %>%
  summarize(mean = mean(diff), sd = sqrt(var(diff)))

#double sharing

table(rtb_final$MACR == 0 & rtb_final$MAMM == 0 & rtb_final$MMCR == 0, rtb_final$comparison) 

#table 2

rtb_final %>% filter(comparison == "MA/MM/CR") %>% summarize(mlen = min(length), maxlen = max(length), sharing = sum(MACR == 0), total = n())
rtb_final %>% filter(comparison == "MA/MM/CR") %>% summarize(mlen = min(length), maxlen = max(length), sharing = sum(MAMM == 0), total = n())
rtb_final %>% filter(comparison == "MA/MM/CR") %>% summarize(mlen = min(length), maxlen = max(length), sharing = sum(MMCR == 0), total = n())

rtb_final %>% filter(comparison == "MA/CR") %>% summarize(mlen = min(length), maxlen = max(length), sharing = sum(MACR == 0), total = n())
rtb_final %>% filter(comparison == "MA/MM") %>% summarize(mlen = min(length), maxlen = max(length), sharing = sum(MAMM == 0), total = n())
rtb_final %>% filter(comparison == "MM/CR") %>% summarize(mlen = min(length), maxlen = max(length), sharing = sum(MMCR == 0), total = n())

rtw_final %>% filter(comparison == "MA/MM/CR") %>% summarize(mlen = min(length), maxlen = max(length), sharing = sum(MA == 0), total = n())
rtw_final %>% filter(comparison == "MA/MM/CR") %>% summarize(mlen = min(length), maxlen = max(length), sharing = sum(MM == 0), total = n())
rtw_final %>% filter(comparison == "MA/MM/CR") %>% summarize(mlen = min(length), maxlen = max(length), sharing = sum(CR == 0), total = n())

rtw_final %>% filter(comparison == "MA/CR") %>% summarize(mlen = min(length), maxlen = max(length), sharing = sum(MA == 0), total = n())
rtw_final %>% filter(comparison == "MA/CR") %>% summarize(mlen = min(length), maxlen = max(length), sharing = sum(CR == 0), total = n())
rtw_final %>% filter(comparison == "MA/MM") %>% summarize(mlen = min(length), maxlen = max(length), sharing = sum(MA == 0), total = n())
rtw_final %>% filter(comparison == "MA/MM") %>% summarize(mlen = min(length), maxlen = max(length), sharing = sum(MM == 0), total = n())
rtw_final %>% filter(comparison == "MM/CR") %>% summarize(mlen = min(length), maxlen = max(length), sharing = sum(MM == 0), total = n())
rtw_final %>% filter(comparison == "MM/CR") %>% summarize(mlen = min(length), maxlen = max(length), sharing = sum(CR == 0), total = n())

#het
rtw_final %>% filter(comparison == "MA/MM/CR", MA > 0) %>% mutate(diff = MA/length) %>% summarize(avghet = mean(diff))
rtw_final %>% filter(comparison == "MA/MM/CR", MM > 0) %>% mutate(diff = MM/length) %>% summarize(avghet = mean(diff))
rtw_final %>% filter(comparison == "MA/MM/CR", CR > 0) %>% mutate(diff = CR/length) %>% summarize(avghet = mean(diff))
rtw_final %>% filter(comparison == "MA/CR", MA >0) %>% mutate(diff = MA/length) %>% summarize(avghet = mean(diff))
rtw_final %>% filter(comparison == "MA/CR", CR>0) %>% mutate(diff = CR/length) %>% summarize(avghet = mean(diff))
rtw_final %>% filter(comparison == "MA/MM", MA>0) %>% mutate(diff = MA/length) %>% summarize(avghet = mean(diff))
rtw_final %>% filter(comparison == "MA/MM", MM>0) %>% mutate(diff = MM/length) %>% summarize(avghet = mean(diff))
rtw_final %>% filter(comparison == "MM/CR", MM>0) %>% mutate(diff = MM/length) %>% summarize(avghet = mean(diff))
rtw_final %>% filter(comparison == "MM/CR", CR>0) %>% mutate(diff = CR/length) %>% summarize(avghet = mean(diff))

##plots

##FIG 4
interval_num=8
rtb_final %>% filter(!is.na(MAMM)) %>% 
    mutate(length_kb = length/1000, diffs = MAMM) %>% 
    mutate(length_bin = cut(length_kb, 
                            breaks=quantile(length_kb, probs=seq(0,1,length.out = interval_num+1), na.rm=T),
                            include.lowest=TRUE, dig.lab=2)) %>%
    group_by(length_bin) %>% 
    summarize(prop_zero = sum(diffs == 0)/length(diffs), num_obs = n(), se=sqrt(((prop_zero*(1-prop_zero))/n()))) %>% 
    ggplot(aes(x=length_bin, y=prop_zero)) + geom_point(size=3) + 
    geom_errorbar(aes(ymin=prop_zero-se, ymax=prop_zero+se), width=.2, position=position_dodge(0.05)) +
#    geom_smooth(mapping = aes(y=prop_zero, x=as.numeric(length_bin)), method=lm, formula = y~x) +
    theme_classic() + 
    coord_cartesian(ylim=c(0,1)) + 
    geom_hline(yintercept=0.5, color="black", linetype="dashed") + 
    ylab ("Fraction Identical in Bin") + 
    xlab("Length Interval (kb)") + 
    theme(axis.text.x = element_text(face="bold", color="#993333", size=12, angle=45, vjust=1, hjust=1), 
          axis.text.y = element_text(face="bold", color="#993333", size=12, angle=0))



rtb_final %>% filter(!is.na(MAMM)) %>% 
  mutate(length_kb = length/1000, diffs = MAMM) %>% 
  mutate(length_bin = cut(length_kb, 
                          breaks=quantile(length_kb, probs=seq(0,1,length.out = interval_num+1), na.rm=T),
                          include.lowest=TRUE, dig.lab=2)) %>%
  group_by(length_bin) %>% summarize(n())

#Fig 2
fig2 <- rtb_final %>% mutate(MAMM = MAMM/length, MMCR  = MMCR/length, MACR = MACR/length) %>%
  select(MACR, MMCR, MAMM) %>% gather(key = "comp", value = "diff") %>%
  ggplot(aes(x=diff, y=..count..)) + 
  geom_histogram(color="black", fill="white", bins=50) +
  facet_grid(comp ~ .)

#zoomed, with 0.01% differences
library(ggforce)
rtb_final %>% mutate(MAMM = MAMM/length, MMCR  = MMCR/length, MACR = MACR/length) %>%
  select(MACR) %>% gather(key = "comp", value = "diff") %>% 
  mutate(diff = diff*100) %>%
  filter(!is.na(diff), diff <= 0.5) %>%
  ggplot(aes(x=diff, y=..count..)) + 
  geom_histogram(color="black", fill="gray", binwidth=0.01) + labs(x="Percent difference (MA vs CR)")

rtb_final %>% mutate(MAMM = MAMM/length, MMCR  = MMCR/length, MACR = MACR/length) %>%
  select(MMCR) %>% gather(key = "comp", value = "diff") %>% 
  mutate(diff = diff*100) %>%
  filter(!is.na(diff), diff <= 0.5) %>%
  ggplot(aes(x=diff, y=..count..)) + 
  geom_histogram(color="black", fill="gray", binwidth=0.01) + labs(x="Percent difference (MM vs CR)")

rtb_final %>% mutate(MAMM = MAMM/length, MMCR  = MMCR/length, MACR = MACR/length) %>%
  select(MAMM) %>% gather(key = "comp", value = "diff") %>% 
  mutate(diff = diff*100) %>%
  filter(!is.na(diff), diff <= 0.5) %>%
  ggplot(aes(x=diff, y=..count..)) + 
  geom_histogram(color="black", fill="gray", binwidth=0.01) + labs(x="Percent difference (MA vs MM)") + facet_zoom(ylim=c(0,5))

rtw_final %>% mutate(MA = MA/length, MM = MM/length, CR = CR/length) %>%
  select(MA, MM, CR) %>% gather(key = "comp", value = "diff") %>%
  mutate(diff = diff*100) %>%
  filter(!is.na(diff), diff <= 0.5) %>%
  ggplot(aes(x=diff, y=..count..)) + 
  geom_histogram(color="black", fill="white", binwidth=0.01) +
  labs(x="Percent difference") +
  facet_grid(comp ~ .)


#Fig 3
##to do##
fig3 <- rtw_final %>% mutate(MA = MA/length, MM = MM/length, CR = CR/length) %>%
  select(MA, MM, CR) %>% gather(key = "comp", value = "diff") %>%
  ggplot(aes(x=diff, y=..count..)) + 
  geom_histogram(color="black", fill="white", bins=50) +
  facet_grid(comp ~ .)

#####
##extract counts

fig2_counts <- ggplot_build(fig3)$data[[1]] %>% as_tibble() %>% 
  mutate(comparison = case_when(PANEL == 1 ~ "MA/CR",
                                PANEL == 2 ~ "MA/MM",
                                PANEL == 3 ~ "MM/CR")) %>% 
  select(comparison, xmin, xmax, count)

fig3_counts <- ggplot_build(fig4)$data[[1]] %>% as_tibble() %>% 
  mutate(comparison = case_when(PANEL == 1 ~ "CR",
                                PANEL == 2 ~ "MA",
                                PANEL == 3 ~ "MM")) %>% 
  select(comparison, xmin, xmax, count)

write_csv(fig2_counts, "fig2_counts.csv")
write_csv(fig3_counts, "fig3_counts.csv")

###########

mamm<-  rt %>% filter(!is.na(MAMM_2) | !is.na(MAMM_3)) %>% select(nanopore, length_3, MAMM_3, length_MAMM, MAMM_2) %>%
  mutate(which_align = ifelse(is.na(length_3), "2", "3"),
         length = ifelse(which_align == "3", length_3, length_MAMM),
         diffs = ifelse(which_align == "3", MAMM_3, MAMM_2)) %>%  select(nanopore, length, diffs)

mamm_3<-rt %>% filter(!is.na(MAMM_3)) %>% select(nanopore, length_3, MAMM_3)

table(mamm$diffs)

#1
mamm %>% filter(mamm$diffs <= 2) %>% summarize(ls = mean(length), mut = sum(diffs))

table(mamm_3$MAMM_3 == 0)

for (interval_num in 5:25) {
  rt %>% filter(!is.na(MAMM_2) | !is.na(MAMM_3)) %>% select(nanopore, length_3, MAMM_3, length_MAMM, MAMM_2) %>%
    mutate(length_3 = ifelse(is.na(length_3), 0, length_3),
           length_MAMM = ifelse(is.na(length_MAMM), 0, length_MAMM)) %>% 
    mutate(which_align = ifelse(length_3 > length_MAMM, "3", "2"),
           length = ifelse(which_align == "3", length_3, length_MAMM),
           diffs = ifelse(which_align == "3", MAMM_3, MAMM_2)) %>%
    select(nanopore, length, diffs) %>% filter(diffs == 0 | diffs > 5) %>%
    mutate(length_kb = length/1000) %>% 
    mutate(length_bin = cut(length_kb, 
                          breaks=quantile(length_kb, probs=seq(0,1,length.out = interval_num+1), na.rm=T),
                          include.lowest=TRUE, dig.lab=2)) %>%
    group_by(length_bin) %>% 
    summarize(prop_zero = sum(diffs == 0)/length(diffs), num_obs = n(), se=sqrt(((prop_zero*(1-prop_zero))/n()))) %>% 
    ggplot(aes(x=length_bin, y=prop_zero)) + geom_point(size=3) + 
    geom_errorbar(aes(ymin=prop_zero-se, ymax=prop_zero+se), width=.2, position=position_dodge(0.05)) +
    theme_classic() + 
    coord_cartesian(ylim=c(0,1)) + 
    geom_hline(yintercept=0.5, color="black", linetype="dashed") + 
    ylab ("Fraction Identical in Bin") + 
    xlab("Length Interval (kb)") + 
    theme(axis.text.x = element_text(face="bold", color="#993333", size=12, angle=45, vjust=1, hjust=1), 
          axis.text.y = element_text(face="bold", color="#993333", size=12, angle=0))
  ggsave(filename = paste0("Length_plot_", interval_num, "filtered_intervals.pdf"))
}


for (interval_num in 5:25) {
  rt %>% filter(!is.na(MAMM_2) | !is.na(MAMM_3)) %>% select(nanopore, length_3, MAMM_3, length_MAMM, MAMM_2) %>%
    mutate(length_3 = ifelse(is.na(length_3), 0, length_3),
           length_MAMM = ifelse(is.na(length_MAMM), 0, length_MAMM)) %>% 
    mutate(which_align = ifelse(length_3 > length_MAMM, "3", "2"),
           length = ifelse(which_align == "3", length_3, length_MAMM),
           diffs = ifelse(which_align == "3", MAMM_3, MAMM_2)) %>%
    select(nanopore, length, diffs) %>% filter(diffs == 0 | diffs > 5) %>%
    mutate(length_kb = length/1000) %>% 
    mutate(length_bin = cut(length_kb, 
                            breaks=quantile(length_kb, probs=seq(0,1,length.out = interval_num+1), na.rm=T),
                            include.lowest=TRUE, dig.lab=2)) %>%
    group_by(length_bin) %>% 
    summarize(prop_zero = sum(diffs == 0)/length(diffs), num_obs = n(), se=sqrt(((prop_zero*(1-prop_zero))/n()))) %>% 
    ggplot(aes(x=length_bin, y=prop_zero)) + geom_point(size=3) + 
    geom_smooth(mapping = aes(y=prop_zero, x=as.numeric(length_bin)), method=lm, formula = y~x) +
    geom_errorbar(aes(ymin=prop_zero-se, ymax=prop_zero+se), width=.2, position=position_dodge(0.05)) +
    theme_classic() + 
    coord_cartesian(ylim=c(0,1)) + 
    geom_hline(yintercept=0.5, color="black", linetype="dashed") + 
    ylab ("Fraction Identical in Bin") + 
    xlab("Length Interval (kb)") + 
    theme(axis.text.x = element_text(face="bold", color="#993333", size=12, angle=45, vjust=1, hjust=1), 
          axis.text.y = element_text(face="bold", color="#993333", size=12, angle=0))
  ggsave(filename = paste0("Length_plot_", interval_num, "_filtered_intervals_linefit.pdf"))
}


num_interval = 15
model<-rt %>% filter(!is.na(MAMM_2) | !is.na(MAMM_3)) %>% select(nanopore, length_3, MAMM_3, length_MAMM, MAMM_2) %>%
  mutate(length_3 = ifelse(is.na(length_3), 0, length_3),
         length_MAMM = ifelse(is.na(length_MAMM), 0, length_MAMM)) %>% 
  mutate(which_align = ifelse(length_3 > length_MAMM, "3", "2"),
         length = ifelse(which_align == "3", length_3, length_MAMM),
         diffs = ifelse(which_align == "3", MAMM_3, MAMM_2)) %>%
  select(nanopore, length, diffs) %>%
  mutate(length_kb = length/1000) %>% 
  mutate(length_bin = cut(length_kb, 
                          breaks=interval_num,
                          include.lowest=TRUE, labels=FALSE)) %>%
  group_by(length_bin) %>% 
  summarize(prop_zero = sum(diffs == 0)/length(diffs), num_obs = n(), se=sqrt(((prop_zero*(1-prop_zero))/n()))) %>% 
  filter(num_obs >= 5) %>%
  lm(prop_zero ~ length_bin, data=.)
  
length_interval<-rt %>% filter(!is.na(MAMM_2) | !is.na(MAMM_3)) %>% select(nanopore, length_3, MAMM_3, length_MAMM, MAMM_2) %>%
  mutate(length_3 = ifelse(is.na(length_3), 0, length_3),
         length_MAMM = ifelse(is.na(length_MAMM), 0, length_MAMM)) %>% 
  mutate(which_align = ifelse(length_3 > length_MAMM, "3", "2"),
         length = ifelse(which_align == "3", length_3, length_MAMM),
         diffs = ifelse(which_align == "3", MAMM_3, MAMM_2)) %>%
  select(nanopore, length, diffs) %>%
  mutate(length_kb = length/1000) %>% summarize(length_interval = max(length_kb) - min(length_kb))

length_per_interval = length_interval / num_interval
length_per_interval * coef(model)[2]


rt %>% filter(!is.na(MAMM_2) | !is.na(MAMM_3)) %>% select(nanopore, length_3, MAMM_3, length_MAMM, MAMM_2) %>%
  mutate(length_3 = ifelse(is.na(length_3), 0, length_3),
         length_MAMM = ifelse(is.na(length_MAMM), 0, length_MAMM)) %>% 
  mutate(which_align = ifelse(length_3 > length_MAMM, "3", "2"),
         length = ifelse(which_align == "3", length_3, length_MAMM),
         diffs = ifelse(which_align == "3", MAMM_3, MAMM_2)) %>%
  select(nanopore, length, diffs) %>%
  mutate(length_kb = length/1000) %>% 
  mutate(length_bin = cut(length_kb, 
                          breaks=quantile(length_kb, probs=seq(0,1,length.out = interval_num+1), na.rm=T),
                          include.lowest=TRUE, labels=FALSE),
         diff_bin = cut(diffs, breaks=c(0,1,5,1000), include.lowest = TRUE, labels=c(0,1,2))) %>%
  with(., table(diff_bin, length_bin))



##### plotting trees ####
library(ape)
r372<-read.tree("~/Projects/genomics/rotifer/manuscript/trees/RAxML_bipartitionsBranchLabels.372Nsremoved_rename_npspace_gb_ML_boots_noremoval.phy")  
r382<-read.tree("~/Projects/genomics/rotifer/manuscript/trees/RAxML_bipartitionsBranchLabels.382Nsremoved_rename_npspace_gb_ML_boots_noremoval.phy")  


## making filters for plotting ###
rtb_filter <- inner_join(rtb_final, region_key, by = c("nanopore" = "nanopore")) %>% select(comparison, region_3, region_MAMM, region_MACR, region_MMCR)

rtb_filter %>% filter(comparison == "MA/MM/CR") %>% select(region_3) %>% write_tsv("../tick_plots/three_way_filter.txt", col_names = FALSE)
rtb_filter %>% filter(comparison == "MA/MM") %>% select(region_MAMM) %>% write_tsv("../tick_plots/MAMM_filter.txt", col_names = FALSE)
rtb_filter %>% filter(comparison == "MA/CR") %>% select(region_MACR) %>% write_tsv("../tick_plots/MACR_filter.txt", col_names = FALSE)
rtb_filter %>% filter(comparison == "MM/CR") %>% select(region_MMCR) %>% write_tsv("../tick_plots/MMCR_filter.txt", col_names = FALSE)


### make NJ trees for 2-way comparisons ###

### From Saitou and Nei (1987, Table 1):
library(ape)
#region 557
x <- c(296,256,246,413,394,206)
M <- matrix(0, 4, 4)
M[lower.tri(M)] <- x
M <- t(M)
M[lower.tri(M)] <- x
dimnames(M) <- list(c("MA1", "MA2", "MM1", "MM2"), c("MA1", "MA2", "MM1", "MM2"))
tr <- nj(M)
plot(tr, "u")

#region 650
x <- c(241,0,241,134,271,134)
M <- matrix(0, 4, 4)
M[lower.tri(M)] <- x
M <- t(M)
M[lower.tri(M)] <- x
dimnames(M) <- list(c("MA1", "MA2", "MM1", "MM2"), c("MA1", "MA2", "MM1", "MM2"))
tr <- nj(M)
tr$edge.length[4] = 0
tr$edge.length[5] = 0
plot(tr, "u")

#region 747
x <- c(216,115,217,115,217,0)
M <- matrix(0, 4, 4)
M[lower.tri(M)] <- x
M <- t(M)
M[lower.tri(M)] <- x
dimnames(M) <- list(c("MA1", "MA2", "MM1", "MM2"), c("MA1", "MA2", "MM1", "MM2"))
tr <- nj(M)
tr$edge.length[1] = 0
tr$edge.length[2] = 0
plot(tr, "u")

##additional analyses -- getting counts
rtb_final %>% filter(MAMM <= 10) %>% rename(diffs=MAMM) %>% group_by(diffs) %>% summarize(number=n()) %>%
  write_excel_csv("MAMM_near_id.csv")

rtb_final %>% filter(MACR <= 10) %>% rename(diffs=MACR) %>% group_by(diffs) %>% summarize(number=n()) %>%
  write_excel_csv("MACR_near_id.csv")

rtb_final %>% filter(MMCR <= 10) %>% rename(diffs=MMCR) %>% group_by(diffs) %>% summarize(number=n()) %>%
  write_excel_csv("MMCR_near_id.csv")
