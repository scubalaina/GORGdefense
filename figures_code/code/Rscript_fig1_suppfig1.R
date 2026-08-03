setwd("/Users/aweinheimer/Documents/GORG_defense/gorgdef_man/v3_files/code_v3/fig1_suppfig1")
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(vegan)
library(reshape2)
library(maps)

dataset_col = c("#0C8AF3","#00437A")
gorg_sys = read.csv("../gorgB_sys_50c_wide_meta_nobs.tsv",sep="\t",header=T,check.names = F)


# Figure 1a

gorg_sys$totaldefenses_5plus = ifelse(gorg_sys$totaldefenses >= 5, 5,
                  gorg_sys$totaldefenses)

gorg_dataset_syscount_sagcount = gorg_sys %>% group_by(dataset,totaldefenses_5plus) %>% 
  summarise(sags = n()) %>% mutate(prop_sags = sags/sum(sags)) %>% as.data.frame()
gorg_dataset_syscount_sagcount$dataset = ifelse(gorg_dataset_syscount_sagcount$dataset == "Dark","zDark","Tropics")

gorg_dataset_syscount_bar = gorg_dataset_syscount_sagcount %>% 
  ggplot(aes(x=as.factor(totaldefenses_5plus),y=prop_sags)) + 
  geom_bar(stat="identity",aes(fill=dataset),position="dodge",width=0.75) + theme_minimal() +
  scale_fill_manual(values=dataset_col,labels=c("Tropics","Dark")) + 
  theme(legend.position="none",axis.text = element_text(color="black",size=10),
        axis.title.y = element_text(size=10),axis.title.x = element_blank()) +
  xlab("Total defenses") + scale_x_discrete(labels=c("0","1","2","3","4","5+")) + ylab("Proportion") + ggtitle("SAGs")
gorg_dataset_syscount_bar


# MAG info
mag_sys_meta = read.csv("magB_50c_sys_meta.tsv",sep="\t",header=T,check.names = F)

mag_defcount = mag_sys_meta %>% group_by(dataset,totaldefenses_5plus) %>% summarise(mags= n()) %>% 
  mutate(prop_mags = mags/sum(mags)) %>% as.data.frame()
mag_defcount$dataset = ifelse(mag_defcount$dataset == "MAG_dark","zMAG_dark",mag_defcount$dataset)
mag_defcount %>% ggplot(aes(x=as.character(totaldefenses_5plus),y=prop_mags,fill=dataset)) +
  geom_bar(stat="identity") + theme_minimal() + scale_fill_manual(values=dataset_col)

mag_dataset_syscount_bar = mag_defcount %>% 
  ggplot(aes(x=as.factor(totaldefenses_5plus),y=prop_mags)) + 
  geom_bar(stat="identity",aes(fill=dataset),position="dodge",width=0.75) + theme_minimal() +
  scale_fill_manual(values=dataset_col,labels=c("Tropics","Dark")) + 
  theme(legend.position="none",axis.text = element_text(color="black",size=10),
        axis.title.y = element_text(size=10),axis.title.x = element_blank()) +
  xlab("Total defenses") + scale_x_discrete(labels=c("0","1","2","3","4","5+")) + 
  ylab("Proportion") + ggtitle("MAGs")
mag_dataset_syscount_bar

tesson_meltALL =read.csv("tesson2022_full_sys_info.csv",sep=",",header = T,check.names = F)
tesson_meltALL$type = ifelse(tesson_meltALL$type == "Lamassu-Fam","Lamassu",
                             tesson_meltALL$type)
tesson_meltALL$type = ifelse(tesson_meltALL$type == "Gao_Mza","Mza",
                             tesson_meltALL$type)
tesson_meta = tesson_meltALL %>% select(Assembly,Superkingdom,phylum,class,order,family,genus,species) %>%
  unique() %>% as.data.frame()
tesson_meltALL$presence = 1
nodefs_tess = tesson_meltALL %>% filter(type == "No system found") %>% as.data.frame()
defs_tess = tesson_meltALL %>% filter(type != "No system found") %>% as.data.frame()
nodefs_tess_gen = as.vector(unique(nodefs_tess$Assembly))
defs_tess_gen = as.vector(unique(defs_tess$Assembly))

tesson_melt = tesson_meltALL %>% filter(type!= "No system found") %>% as.data.frame()
tesson_melt2 = tesson_melt %>% group_by(Assembly,type) %>% summarise(instances = n()) %>% as.data.frame()
tesson_sys = unique(tesson_melt2$type)
tesson_wide = tesson_melt2 %>% pivot_wider(names_from = type,values_from = instances) %>% as.data.frame()
tesson_wide = tesson_wide %>% replace(is.na(.),0) %>% as.data.frame()
tesson_wide$totalsystems = rowSums(tesson_wide[,tesson_sys])
tesson_wide$totaldefenses_5plus = ifelse(tesson_wide$totalsystems >=5, 5, tesson_wide$totalsystems)
tesson_defcount_tbl = tesson_wide %>% group_by(totaldefenses_5plus) %>% summarise(gens=n()) %>% 
  as.data.frame()
tesson_defcount_tbl = tesson_defcount_tbl %>% add_row(totaldefenses_5plus = 0, gens = 1450) %>% as.data.frame()
tesson_defcount_tbl$dataset = "RefSeq"
tesson_defcount_tbl$prop_gens = tesson_defcount_tbl$gens/sum(tesson_defcount_tbl$gens)

refseq_dataset_syscount_bar = tesson_defcount_tbl %>%
  ggplot(aes(x=as.factor(totaldefenses_5plus),y=prop_gens)) +
  geom_bar(stat="identity",width=0.75) + theme_minimal() +
  theme(legend.position="none",axis.text = element_text(color="black",size=10),
        axis.title = element_text(size=10)) +
  xlab("Total defenses") + scale_x_discrete(labels=c("0","1","2","3","4","5+")) + 
  ylab("Proportion") + ggtitle("RefSeq assemblies")
refseq_dataset_syscount_bar

ggarrange(gorg_dataset_syscount_bar,mag_dataset_syscount_bar,refseq_dataset_syscount_bar,
          nrow=3,ncol=1)


# Figure 1b, c prep

gorg_datset_order_defstats = gorg_sys %>% group_by(dataset,Order) %>% 
  summarise(mean_def = mean(totaldefenses),sd_def = sd(totaldefenses),sags = n()) %>%
  replace(is.na(.),0) %>% as.data.frame()
gorg_datset_order_defstats$Order_count = paste(gorg_datset_order_defstats$Order,gorg_datset_order_defstats$sags,sep="\n(")
gorg_datset_order_defstats$Order_count = paste(gorg_datset_order_defstats$Order_count,"",sep=")")
gorg_datset_order_defstats$se_def = gorg_datset_order_defstats$sd_def / sqrt(gorg_datset_order_defstats$sags)
gorg_datset_order_defstats$minse_df = gorg_datset_order_defstats$mean_def - gorg_datset_order_defstats$se_def - gorg_datset_order_defstats$se_def
gorg_datset_order_defstats$maxse_df = gorg_datset_order_defstats$mean_def + gorg_datset_order_defstats$se_def + gorg_datset_order_defstats$se_def

# Figure 1b
gorgT_top10_ord = gorg_datset_order_defstats %>% filter(dataset == "Tropics") %>% arrange(desc(sags)) %>% head(n=10) %>% select(Order_count)
gorgT_top10_ord = as.vector(gorgT_top10_ord$Order)
gorg_sys %>% filter(dataset == "Tropics") %>% summarise(meandef = mean(totaldefenses))
gorgT_order_def_bar = gorg_datset_order_defstats %>% 
  filter(Order_count %in% gorgT_top10_ord) %>% filter(dataset=="Tropics") %>% 
  ggplot(aes(y=Order_count,x=mean_def)) + geom_bar(fill="#0C8AF3",stat="identity",width=0.5) + 
  theme_minimal()  +
  geom_errorbar(aes(xmin=minse_df, xmax=maxse_df), width=.2, color="gray20") +
  xlab("Total defenses") + ggtitle("Top 10 Orders in GORG Tropics") +
  theme(axis.title.x = element_text(size=12), 
        axis.text = element_text(size=10,color="black"),axis.title.y = element_blank()) +
  scale_y_discrete(limits=rev(gorgT_top10_ord)) +
  geom_vline(xintercept = 1.089,color="red") + xlim(0,2.5)
gorgT_order_def_bar

# Figure 1c
gorgD_top10_ord = gorg_datset_order_defstats %>% filter(dataset == "Dark") %>% arrange(desc(sags)) %>% head(n=10) %>% select(Order_count)
gorgD_top10_ord = as.vector(gorgD_top10_ord$Order)
gorg_sys %>% filter(dataset == "Dark") %>% summarise(meandef = mean(totaldefenses))
gorgD_order_def_bar = gorg_datset_order_defstats %>% filter(Order_count %in% gorgD_top10_ord) %>% filter(dataset=="Dark") %>% 
  ggplot(aes(y=Order_count,x=mean_def)) + geom_bar(fill="#00437A",stat="identity",width=0.5) + 
  theme_minimal()  + ggtitle("Top 10 Orders in GORG Dark") +
  geom_errorbar(aes(xmin=minse_df, xmax=maxse_df), width=.2, color="gray50") +
  xlab("Total defenses") + ylab("Top 10 Orders in GORG Dark") + 
  theme(axis.title.x = element_text(size=12), axis.title.y = element_blank(),
        axis.text = element_text(size=10,color="black")) +
  scale_y_discrete(limits=rev(gorgD_top10_ord)) +
  geom_vline(xintercept = 1.213,color="red") + xlim(0,2.5)
gorgD_order_def_bar

ggarrange(gorgT_order_def_bar,gorgD_order_def_bar,nrow=1,ncol=2)


# Figure S1

gorgB_stat_sum = gorg_sys %>% group_by(dataset,defpres) %>% summarise(sags = n()) %>% as.data.frame()
gorgB_compl_sum = gorg_sys %>% group_by(dataset) %>% summarise(mean_compl = mean(completeness)) %>% as.data.frame()
gorgB_stat_sum$dataset = ifelse(gorgB_stat_sum$dataset == "Tropics","gorgT","gorgD")
gorgB_compl_sum$dataset = ifelse(gorgB_compl_sum$dataset == "Tropics","gorgT","gorgD")
gorgB_stat_sum$prop_sags = ifelse(gorgB_stat_sum$dataset == "gorgD",gorgB_stat_sum$sags / 2156, 
                                  gorgB_stat_sum$sags / 3880)

gorgB_stat_sum = merge(gorgB_stat_sum,gorgB_compl_sum,all.x=T)

mag_sys_meta$defpres = ifelse(mag_sys_meta$totaldefenses > 0, "yes","no")
magB_compl = mag_sys_meta %>% group_by(dataset) %>% summarise(mean_compl = mean(CheckM1_completeness)) %>% as.data.frame()
magB_stat_sum = mag_sys_meta %>% group_by(dataset,defpres) %>% summarise(sags = n()) %>% as.data.frame()
magB_stat_sum$prop_sags = ifelse(magB_stat_sum$dataset == "MAG_dark",
                                 magB_stat_sum$sags/5468, magB_stat_sum$sags/4588)
magB_stat_sum = merge(magB_stat_sum,magB_compl,all.x=T)
sagmag_stats = rbind(gorgB_stat_sum,magB_stat_sum)
sagmag_stats$combo_label  = paste(sagmag_stats$dataset,sagmag_stats$defpres,sep="_")
sagmag_stats
sagmag_defcomp_bar = sagmag_stats %>% ggplot(aes(x=dataset,y=prop_sags*100)) + 
  geom_bar(stat="identity",aes(fill=combo_label,alpha=as.factor(defpres)),width=0.5) + 
  theme_minimal() +
  ylab("Percent of assemblies\nwith a defense") + xlab("Dataset") + 
  geom_point(data=sagmag_stats, aes(dataset, mean_compl), shape=95, size=20,color="red") +
  theme(legend.position = "none") + 
  scale_fill_manual(values=c("white","#00437A","white","#0C8AF3","white","#556B2F","white","#A2CD5A"))+
  scale_x_discrete(limits=c("gorgT","gorgD","MAG_sunlit","MAG_dark"),
                   labels = c("GORG\nTropics","GORG\nDark","MAG\nSunlit","MAG\nDark"))
sagmag_defcomp_bar

