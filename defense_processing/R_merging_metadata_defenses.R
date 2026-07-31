setwd("/Users/aweinheimer/Documents/GORG_defense/gorgdef_man/v3_files/code_v3/")

library(ggplot2)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(vegan)
library(reshape2)
library(maps)

# set GORG Dark and Tropics colors
dataset_col = c("#0C8AF3","#00437A")

# load in table where DefenseFinder and Padloc system terminology has been manually reconciled for consistency
sysinfo = read.csv("sys_df_ploc_common_name_utf_v3.tsv",
                   sep="\t",header=T,check.names = F)
# subset  columns needed for downstream analyses: system name and subtype name
sysinfosub = sysinfo %>% select(common_system,common_subtype) %>% as.data.frame()

# load in system count and substype count tables generated from step XX 
# each genome is a line and the columns are system or subtype and the number of instances in that genome
gorgD_sys = read.csv("gorgD_system_count.tsv",sep="\t",header=T,check.names = F)
gorgT_sys = read.csv("gorgT_system_count.tsv",sep="\t",header=T,check.names = F)
gorgD_sub = read.csv("gorgD_subtype_count.tsv",sep="\t",header=T,check.names = F)
gorgT_sub = read.csv("gorgT_subtype_count.tsv",sep="\t",header=T,check.names = F)

# load in metadata for GORG Tropics and Dark
saginfo = read.csv("GORG_SAG_metadata_combined_v3.tsv",sep="\t",header=T,check.names = F)
# filter table for genomes with a minimum completeness of 50%
saginfo_50c = saginfo %>% filter(completeness >= 50) %>% as.data.frame()

# add dataset information to system table
gorgD_sys$dataset="Dark"
gorgT_sys$dataset = "Tropics"
# join Dark and Tropics system tables
gorgB_sys = rbind(gorgD_sys,gorgT_sys)
# add dataset information to subtype table
gorgD_sub$dataset="Dark"
gorgT_sub$dataset = "Tropics"
# join Dark and Tropics subtype tables
gorgB_sub = rbind(gorgD_sub,gorgT_sub)


# Add SAG genomic information to system dataframe
gorgB_sys_meta = merge(gorgB_sys,saginfo,by=c("genome","dataset"),all.x = T)
# filter for genomes with a minimum of 50% completeness
gorgB_sys_50c_meta = gorgB_sys_meta %>% filter(completeness >= 50) %>% as.data.frame()
# remove SAGs from special environments that skew comparisons - identified by Chang et al 
gorgB_sys_50c_meta_nobs = gorgB_sys_50c_meta %>% 
  filter(! niche_latest %in% c("Black","Baltic","Ross", "Arctic")) %>%
  as.data.frame()
# get list of SAGs with minimum 50% completeness for other filtering
saginfo_50c = saginfo %>% filter(completeness >= 50) %>% as.data.frame()
sag50c_list = as.vector(saginfo_50c$genome)
gorgB_sys_50c = gorgB_sys %>% filter(genome %in% sag50c_list) %>% as.data.frame()
gorgB_sub_50c = gorgB_sub %>% filter(genome %in% sag50c_list) %>% as.data.frame()

# pivot dataframe so each row is a genome and the columns are defenses with values as the number of instances of that defense
gorgB_sys_50c_wide = gorgB_sys_50c %>%  pivot_wider(names_from = system,values_from = instances) %>% replace(is.na(.), 0)  %>% as.data.frame()
gorgB_sub_50c_wide = gorgB_sub_50c %>% pivot_wider(names_from = subtype,values_from = instances) %>% replace(is.na(.), 0)  %>% as.data.frame()

# add in SAG metadata and fill empty fields with 0s in case
gorgB_sys_50c_wide_meta = merge(saginfo_50c,gorgB_sys_50c_wide,by=c("genome","dataset"),all.x = T)
gorgB_sys_50c_wide_meta = gorgB_sys_50c_wide_meta %>% replace(is.na(.),0) %>% as.data.frame()
gorgB_sub_50c_wide_meta = merge(saginfo_50c,gorgB_sub_50c_wide,by=c("genome","dataset"),all.x = T)
gorgB_sub_50c_wide_meta = gorgB_sub_50c_wide_meta %>% replace(is.na(.),0) %>% as.data.frame()

# calculate number of defenses in a genome (totaldefenses) by summing the columns of a defense
gorgB_sys_50c_list = as.vector(unique(gorgB_sys_50c$system))
gorgB_sys_50c_wide_meta$totaldefenses = rowSums(gorgB_sys_50c_wide_meta[,gorgB_sys_50c_list])
gorgB_sub_50c_list = as.vector(unique(gorgB_sub_50c$subtype))
gorgB_sub_50c_wide_meta$totaldefenses = rowSums(gorgB_sub_50c_wide_meta[,gorgB_sub_50c_list])
# create column (defpres) thet specifies whether a SAG has a defense 
gorgB_sys_50c_wide_meta$defpres = ifelse(gorgB_sys_50c_wide_meta$totaldefenses > 0, "yes","no")
gorgB_sub_50c_wide_meta$defpres = ifelse(gorgB_sub_50c_wide_meta$totaldefenses > 0, "yes","no")

# get list of SAGs that have a defense to add defpres columns to the melted tables
gorgB_50c_defpres_sags = gorgB_sys_50c_wide_meta %>% filter(defpres == "yes") %>% as.data.frame()
gorgB_50c_defpres_sags = as.vector(gorgB_50c_defpres_sags$genome)

# add SAG metadata to the "melted" version of the dataframes (gorgB_sys_50c, gorgB_sub_50c) for other types of analyses
gorgB_sys_50c_meta = merge(saginfo_50c,gorgB_sys_50c,by=c("genome","dataset"),all.x=T)
gorgB_sys_50c_meta$defpres = ifelse(gorgB_sys_50c_meta$genome %in% gorgB_50c_defpres_sags, "yes","no")
# for SAGs that lack a defense, put "0" in the defense column - to retain all genomes
gorgB_sys_50c_meta = gorgB_sys_50c_meta %>% replace(is.na(.),0)  %>% as.data.frame()

# repeats the above lines on sub dataframe (gorgB_sub_50c)
gorgB_sub_50c_meta = merge(saginfo_50c,gorgB_sub_50c,by=c("genome","dataset"),all.x=T)
gorgB_sub_50c_meta = gorgB_sub_50c_meta %>% replace(is.na(.),0)  %>% as.data.frame()
gorgB_sub_50c_meta$defpres = ifelse(gorgB_sub_50c_meta$genome %in% gorgB_50c_defpres_sags, "yes","no")

# remove SAGs from "unique" environments (designated by Chang et al) filtered label 'nobs'
gorgB_sys_50c_wide_meta_nobs = gorgB_sys_50c_wide_meta %>% 
  filter(! niche_latest %in% c("Black","Baltic","Ross", "Arctic")) %>%  as.data.frame()
gorgB_sub_50c_wide_meta_nobs = gorgB_sub_50c_wide_meta %>% 
  filter(! niche_latest %in% c("Black","Baltic","Ross", "Arctic")) %>%   as.data.frame()
gorgB_sys_50c_meta_nobs = gorgB_sys_50c_meta %>% filter(! niche_latest %in% c("Black","Baltic","Ross", "Arctic")) %>%
  as.data.frame()
gorgB_sub_50c_meta_nobs = gorgB_sub_50c_meta %>% 
  filter(! niche_latest %in% c("Black","Baltic","Ross", "Arctic")) %>%   as.data.frame()

# calculate number of defense systems or subtypes in a SAG to add to gorgB_sys_50c_wide_meta_nobs / gorgB_sub_50c_wide_meta_nobs tables
gorgB_sys_50c_systype_count = gorgB_sys_50c_meta_nobs %>% filter(system != 0) %>% group_by(genome) %>% summarise(sys_types = n()) %>% as.data.frame()
gorgB_sys_50c_subtype_count = gorgB_sub_50c_meta_nobs %>% filter(subtype != 0) %>% group_by(genome) %>% summarise(sub_types = n()) %>% as.data.frame()
gorgB_sys_50c_sysdiv_df = merge(gorgB_sys_50c_systype_count,gorgB_sys_50c_subtype_count,on=c("genome"))
gorgB_sys_50c_wide_meta_nobs = merge(gorgB_sys_50c_wide_meta_nobs,gorgB_sys_50c_sysdiv_df,on=c("genome"),all.x = T)
gorgB_sub_50c_wide_meta_nobs = merge(gorgB_sub_50c_wide_meta_nobs,gorgB_sys_50c_sysdiv_df,on=c("genome"),all.x = T)

# Load GORG Tropics SAG information seperately - includes fields absent in GORG Dark information
gorgt_init = read.csv("gorg-tropics_sag_metadata.csv",sep=",",header=T,check.names = F)
# get list of SAGs sorted with Syto9 prokaryote stain and gate
gorgt_keep = gorgt_init %>% filter(`FACS mode` == "Syto9-Prok") %>% select(genome) %>% as.data.frame()
# make list of these SAGs to filter downstream dataframes
gorgt_keep = as.vector(gorgt_keep$genome)


# compile list of SAGs to retain with ("_keep")
gorgD_keepgen =gorgB_sub_50c_wide_meta_nobs %>% filter(dataset == "Dark") %>% select(genome) %>%
  as.data.frame()
gorgD_keepgen = as.vector(gorgD_keepgen$genome)
gorgT_keepgen = gorgB_sys_50c_wide_meta_nobs %>% filter(genome %in% gorgt_keep) %>% select(genome) %>% as.data.frame()
gorgT_keepgen = as.vector(gorgT_keepgen$genome)
gorgB_keep = c(gorgD_keepgen,gorgT_keepgen)

# filter dataframes to retain the keep list
gorgB_sys_50c_wide_meta_nobs = gorgB_sys_50c_wide_meta_nobs %>% filter(genome %in% gorgB_keep) %>% as.data.frame()
gorgB_sub_50c_wide_meta_nobs = gorgB_sub_50c_wide_meta_nobs %>% filter(genome %in% gorgB_keep) %>% as.data.frame()
gorgB_sys_50c_meta_nobs = gorgB_sys_50c_meta_nobs %>% filter(genome %in% gorgB_keep) %>% as.data.frame()
gorgB_sub_50c_meta_nobs = gorgB_sub_50c_meta_nobs %>% filter(genome %in% gorgB_keep) %>% as.data.frame()

# create dataframes that seperate GORG Tropics and Dark for downstream analyses
gorgD_sys_50c_wide_meta_nobs = gorgB_sys_50c_wide_meta_nobs %>% filter(dataset == "Dark") %>% as.data.frame()
gorgT_sys_50c_wide_meta_nobs = gorgB_sys_50c_wide_meta_nobs %>% filter(dataset == "Tropics") %>% as.data.frame()
gorgD_sub_50c_wide_meta_nobs = gorgB_sub_50c_wide_meta_nobs %>% filter(dataset == "Dark") %>% as.data.frame()
gorgT_sub_50c_wide_meta_nobs = gorgB_sub_50c_wide_meta_nobs %>% filter(dataset == "Tropics") %>% as.data.frame()
gorgD_sub_50c_meta_nobs = gorgB_sub_50c_meta_nobs %>% filter(dataset == "Dark") %>% as.data.frame()
gorgT_sub_50c_meta_nobs = gorgB_sub_50c_meta_nobs %>% filter(dataset == "Tropics") %>% as.data.frame()
gorgD_sys_50c_meta_nobs = gorgB_sys_50c_meta_nobs %>% filter(dataset == "Dark") %>% as.data.frame()
gorgT_sys_50c_meta_nobs = gorgB_sys_50c_meta_nobs %>% filter(dataset == "Tropics") %>% as.data.frame()

# write files for input in other analyses 
# Supplemental Table S1
write.table(gorgB_sys_50c_wide_meta_nobs,"gorgB_sys_50c_wide_meta_nobs.tsv",sep="\t",quote=F,row.names=F)
# Supplemental Table S2
write.table(gorgB_sub_50c_wide_meta_nobs,"gorgB_sub_50c_wide_meta_nobs.tsv",sep="\t",quote=F,row.names=F)
