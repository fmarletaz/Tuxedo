library(tidyverse)
library(zoo)

setwd("~/Dropbox/Projects/Tuxedo/Tuxedo_urchin/GenomeStats")

univ_col=c("#A6CEE3","#1F78B4","#1B9E77","#33A02C","#B2DF8A","#E31A1C","#FB9A99","#FFFF99","#B15928","#D95F02",
           "#7570B3","#E7298A","#66A61E","#CAB2D6","#6A3D9A","#E6AB02","#A6761D","#8DD3C7","#666666","#FF7F00",
           "#FDBF6F","#BEBADA","#ffd92f","#e9521e","lightgrey")


# Blue male genome 
male=read.table("Mg_BlMal_BlMalRd_10kov.regions.bed")
female=read.table("Mg_BlMal_BlFemRd_10kov.regions.bed")
names(male)=c('chr','start','end','m_cov')
names(female)=c('chr','start','end','f_cov')
levels(as.factor(male$chr))
#male %>% group_by(chr) %>% summarise(max=max(end)) %>% arrange(-max) %>% head(n=25) -> chrsize#,decreasing=T)
#chroms=c("scaffold_1","scaffold_2","scaffold_3","scaffold_4","scaffold_5","scaffold_6","scaffold_7","scaffold_8","scaffold_9","scaffold_10","scaffold_11","scaffold_12","scaffold_13","scaffold_14","scaffold_15","scaffold_16","scaffold_17","scaffold_18","scaffold_19","scaffold_20","scaffold_21")

inner_join(male,female)-> sex_depth
sex_depth %>% #filter(chr %in% chroms) %>% 
  filter(m_cov<200) %>%
  filter(f_cov<200) %>%
  mutate(nm_cov=m_cov/median(m_cov)) %>%
  mutate(nf_cov=f_cov/median(f_cov)) %>% 
  mutate(ratio=nm_cov/nf_cov) %>%
  mutate(log2ratio=log2(ratio))-> sex_depth
sex_depth %>% arrange(chr, start) %>%
  group_by(chr) %>%
  mutate(smth_ratio = rollmedian(ratio, k = 21, fill = NA))%>%
  #mutate(ksm_ratio = ksmooth(start, ratio, kernel = "normal", bandwidth = 10000))%>%
  ungroup() -> sex_depth_blmale
write_tsv(sex_depth_blmale,file="Mg_BlueMale_SexCov.tsv")
sex_depth_blmale$chr <- factor(sex_depth_blmale$chr,levels = str_sort(unique(sex_depth_blmale$chr), numeric = TRUE))
ggplot(filter(sex_depth_blmale,grepl('Chr',chr)), aes(start, log2ratio)) +
  geom_point(size = 0.2, alpha = 0.3) +
  geom_smooth(span = 0.05, se = FALSE) +
  geom_hline(yintercept = 0, linetype = 1, colour = "grey50") +
  geom_hline(yintercept = c(-1, 1), linetype = 2, colour = "grey50") +
  facet_wrap(~ chr, scales = "free_x") + ylim(-2,2)+
  labs(y = expression(log[2](male/female)), x = "position")+theme_bw()


ggplot(filt.male,aes(x=m_cov,col=chr))+ geom_density(alpha=.2,lwd=0.75)  +theme_bw()
ggplot(filt.female,aes(x=f_cov,col=chr))+ geom_density(alpha=.2,lwd=0.75)  +theme_bw()
ggplot(filter(sex_depth_blmale,smth_ratio<2),aes(x=smth_ratio,col=chr))+ geom_density(alpha=.2,lwd=0.75)  +theme_bw()
ggplot(filter(sex_depth_blmale,ratio<2),aes(x=ratio,col=chr))+ geom_density(alpha=.2,lwd=0.75)  +theme_bw()
ggplot(filter(sex_depth_blmale,smth_ratio<2),aes(x=smth_ratio,col=chr))+ geom_density(alpha=.2,lwd=0.75)  +theme_bw()

sex_depth_blmale %>% filter(chr=='scaffold_4',smth_ratio<1.5) %>% ggplot(aes(start,smth_ratio))+geom_line()+
  geom_hline(yintercept = c(0.75,1,1.25), color='red') + theme_classic()+ggtitle("M/F ration in scaf_4 (BlueMale)")




# Red Female genome UPDATED!!
male=read.table("Mg_RdFem_BlMalRd_10kov.regions.bed")
female=read.table("Mg_RdFem_BlFemRd_10kov.regions.bed")
names(male)=c('chr','start','end','m_cov')
names(female)=c('chr','start','end','f_cov')
levels(as.factor(male$chr))

inner_join(male,female)-> sex_depth

sex_depth %>% filter(grepl('Chr',chr)) %>% 
  filter(m_cov<200) %>%
  filter(f_cov<200) %>%
  mutate(nm_cov=m_cov/median(m_cov)) %>%
  mutate(nf_cov=f_cov/median(f_cov)) %>% 
  mutate(ratio=nm_cov/nf_cov)  -> sex_depth_redfem

sex_depth_redfem %>% arrange(chr, start) %>%
  group_by(chr) %>%
  mutate(smth_ratio = rollmedian(ratio, k = 21, fill = NA))%>%
  #mutate(ksm_ratio = ksmooth(start, ratio, kernel = "normal", bandwidth = 10000))%>%
  ungroup() -> sex_depth_redfem
write_tsv(sex_depth_redfem,file="Mg_RedFem_SexCov.tsv")

# testing 
require(rstatix)   # effect sizes + tidy post-hoc
require(FSA)       # dunnTest
require(purrr)
sex_depth_redfem <- sex_depth_redfem %>% mutate(log2ratio=log2(ratio)) 
kruskal.test(log2(ratio) ~ chr, data=sex_depth_redfem)
sex_depth_redfem %>% kruskal_effsize(log2(ratio) ~ chr)
sex_depth_redfem %>% dunn_test(log2ratio ~ chr, p.adjust.method = "fdr")

res <- map_dfr(unique(sex_depth_blmale$chr), function(c) {
  d <- sex_depth_blmale %>% mutate(grp = if_else(chr == c, "this", "rest"))
  w   <- wilcox.test(log2ratio ~ grp, data = d)
  eff <- d %>% wilcox_effsize(log2ratio ~ grp)
  tibble(
    chr         = c,
    median_log2 = median(d$log2ratio[d$grp == "this"]),
    p           = w$p.value,
    effsize_r   = eff$effsize
  )
}) %>%
  mutate(p_adj = p.adjust(p, method = "fdr")) %>%
  arrange(desc(abs(effsize_r)))
sex_depth_redfem$chr <- factor(
  sex_depth_redfem$chr,
  levels = str_sort(unique(sex_depth_redfem$chr), numeric = TRUE)
)
write_tsv(≈,file='res_wilcox_redfem.txt')
write_tsv(res,file='res_wilcox_blumale.txt')

ggplot(sex_depth_redfem, aes(start, log2ratio)) +
  geom_point(size = 0.2, alpha = 0.3) +
  geom_smooth(span = 0.05, se = FALSE) +
  geom_hline(yintercept = 0, linetype = 1, colour = "grey50") +
  geom_hline(yintercept = c(-1, 1), linetype = 2, colour = "grey50") +
  facet_wrap(~ chr, scales = "free_x") + ylim(-2,2)+
  labs(y = expression(log[2](male/female)), x = "position")+theme_bw()


ggplot(filter(sex_depth_redfem,ratio<2),aes(x=log2(ratio),col=chr))+ 
  geom_density(alpha=.2,lwd=0.75) +
  scale_color_manual(values=univ_col)+theme_bw()
ggplot(filter(sex_depth_redfem,ratio<2),aes(x=ratio,col=chr))+ 
  stat_ecdf() +  scale_color_manual(values=univ_col)+theme_bw()
sex_depth_redfem %>% filter(chr=='Chr07' & smth_ratio<2) %>% ggplot(aes(start,smth_ratio))+geom_line()+
  geom_smooth(aes(start,ratio))+  geom_hline(yintercept = c(1), color='red') + theme_classic()+ggtitle("M/F ration chr07 (Red female)")

######  Red Male genome #####

male=read.table("BluMale-RedMale_cov.regions.bed")
female=read.table("BluFemale-RedMale.regions.bed")
names(male)=c('chr','start','end','m_cov')
names(female)=c('chr','start','end','f_cov')
#male %>% mutate(nm_cov=m_cov/median(filter(male,grepl('Chr',chr))$m_cov)) -> male
#female %>% mutate(nf_cov=f_cov/median(filter(female,grepl('Chr',chr))$f_cov)) -> female

inner_join(male,female)-> sex_depth_redmale

sex_depth_redmale %>% filter(grepl('Chr',chr)) %>% 
  filter(m_cov<200) %>%
  filter(f_cov<200) %>%
  mutate(nm_cov=m_cov/median(m_cov)) %>%
  mutate(nf_cov=f_cov/median(f_cov)) %>% 
  mutate(ratio=nm_cov/nf_cov)  -> sex_depth_redmale

sex_depth_redmale %>% arrange(chr, start) %>%
  group_by(chr) %>%
  mutate(smth_ratio = rollmedian(ratio, k = 21, fill = NA))%>%
  ungroup() -> sex_depth_redmale
write_tsv(sex_depth_redmale,file="RedMale_SexCov.tsv")

ggplot(filter(sex_depth_redmale,ratio<2),aes(x=ratio,col=chr))+ 
  geom_density(alpha=.5,lwd=0.75) + 
  scale_color_manual(values=univ_col)+theme_bw()
ggplot(filter(sex_depth_redmale,ratio<2),aes(x=ratio,col=chr))+ 
  stat_ecdf() +  scale_color_manual(values=univ_col)+theme_bw()
sex_depth_redmale %>% filter(chr=='Chr15' & smth_ratio<2) %>% 
  ggplot(aes(start,smth_ratio))+geom_line()+geom_smooth(aes(start,ratio))+ 
  geom_hline(yintercept = c(1), color='red')+ theme_classic()
sex_depth_redmale %>% filter(chr=='Chr06' & smth_ratio<2) %>% 
  ggplot(aes(start,smth_ratio))+geom_line()+geom_smooth(aes(start,ratio))+ 
  geom_hline(yintercept = c(1), color='red')+ theme_classic()

length(univ_col)




#regcov %>% filter(scaf %in% cs$V1 & cov<200) %>% ggplot(aes(x=cov,color=scaf)) + geom_density(alpha=.2,lwd=0.75) +
#  geom_vline(xintercept = c(50,80), color='darkgrey') +theme_bw()
# male: 348827063
# female: 370426615
