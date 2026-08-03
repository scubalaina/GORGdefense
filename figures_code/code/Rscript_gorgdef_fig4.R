setwd("/Users/aweinheimer/Documents/GORG_defense/gorgdef_man/v3_files/code_v3/fig4/")
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(vegan)
library(reshape2)
library(maps)

# specify the dataset colors
dataset_col = c("#0C8AF3","#00437A")
# load in the wide and melted forms of the GORG defenses system info
gorg_wide = read.csv("../gorgB_sys_50c_wide_meta_nobs.tsv",sep="\t",header=T,check.names = F)
gorg_melt = read.csv("../gorgB_sys_50c_meta_nobs.tsv",sep="\t",header=T,check.names = F)

# separate Tropics and Dark
gorgT_wide = gorg_wide %>% filter(dataset == "Tropics") %>% as.data.frame()
gorgD_wide = gorg_wide %>% filter(dataset == "Dark") %>% as.data.frame()
gorgT_melt = gorg_melt %>% filter(dataset == "Tropics") %>% as.data.frame()
gorgD_melt = gorg_melt %>% filter(dataset == "Dark") %>% as.data.frame()

# count the number of SAGs in each order
gorgT_ord_sagcount =  gorgT_wide %>% group_by(Order) %>% summarise(sags = n()) %>% 
  arrange(desc(sags)) %>% as.data.frame()
gorgD_ord_sagcount =  gorgD_wide %>% group_by(Order) %>% summarise(sags = n()) %>% 
  arrange(desc(sags)) %>% as.data.frame()

# count the number of SAGs with each system 
gorgT_sys_sagcount = gorgT_melt %>% group_by(system) %>% summarise(sys_sags = n()) %>% as.data.frame()
gorgD_sys_sagcount = gorgD_melt %>% group_by(system) %>% summarise(sys_sags = n()) %>% as.data.frame()

# count the number of SAGs with in each Order with each system
gorgT_ord_sys_sagcount = gorgT_melt %>% group_by(Order,system) %>% 
  summarise(ordsys_sags = n()) %>% as.data.frame()
gorgD_ord_sys_sagcount = gorgD_melt %>% group_by(Order,system) %>% 
  summarise(ordsys_sags = n()) %>% as.data.frame()

# combine order system count and order count tables to get proportion of SAGs in an order with a given system
gorgT_ord_sys_stats = merge(gorgT_ord_sys_sagcount,gorgT_ord_sagcount,on=c("Order"),all.x=T)
gorgT_ord_sys_stats$prop_ordsys_sags = gorgT_ord_sys_stats$ordsys_sags / gorgT_ord_sys_stats$sags
gorgT_ord_sys_stats = merge(gorgT_ord_sys_stats,gorgT_sys_sagcount,on=c("system"),all.x=T)
gorgD_ord_sys_stats = merge(gorgD_ord_sys_sagcount,gorgD_ord_sagcount,on=c("Order"),all.x=T)
gorgD_ord_sys_stats$prop_ordsys_sags = gorgD_ord_sys_stats$ordsys_sags / gorgD_ord_sys_stats$sags
gorgD_ord_sys_stats = merge(gorgD_ord_sys_stats,gorgD_sys_sagcount,on=c("system"),all.x=T)

# make bubble plot of proportion of SAGs in an order with each system for orders that contain at least
# four SAGs and systems that are present in more than 1 SAG
gorgT_ord_sys_stats2 = gorgT_ord_sys_stats %>% add_row(system = "RM",Order="AXxx", 
                                                       ordsys_sags = 30,sags = 5,prop_ordsys_sags=1,
                                sys_sags = 50) %>% as.data.frame()
gorgT_ord_sys_bubble = gorgT_ord_sys_stats %>% filter(sags > 4)  %>% filter(sys_sags > 1)  %>% filter(system != 0) %>% 
  ggplot(aes(y=reorder(Order,sags),x=reorder(system,-sys_sags))) + 
  geom_point(aes(alpha=prop_ordsys_sags),color="dodgerblue",size=3) +
  theme_minimal() + geom_point(shape=1,color="black",size=3) + 
  theme(axis.text.x = element_text(angle=90, size = 11),axis.text.y=element_text(size = 11),axis.title = element_blank(),legend.position = "top") +
  scale_alpha(range = c(0, 1),limits = c(0, 1)) + guides(alpha = guide_legend(title = "Proportion of SAGs in GORG Tropics Order"))
gorgT_ord_sys_bubble
ggsave("/Users/aweinheimer/Documents/GORG_defense/gorgdef_man/v3_files/v3_figs/v3_svgs/gorgT_ord_sys_bubble.svg",gorgT_ord_sys_bubble,width = 6.8,height = 6.4)

gorgD_ord_sys_bubble = gorgD_ord_sys_stats %>% filter(sags > 4) %>% filter(sys_sags > 1) %>% filter(system != 0) %>% 
  ggplot(aes(y=reorder(Order,sags),x=reorder(system,-sys_sags))) + 
  geom_point(aes(alpha=prop_ordsys_sags),color="navy",size=3) +
  theme_minimal() + geom_point(shape=1,color="black",size=3) + 
  theme(axis.text.x = element_text(angle=90, size = 11),axis.text.y=element_text(size = 11),axis.title = element_blank(),legend.position = "top") +
  scale_alpha(range = c(0, 1)) + guides(alpha = guide_legend(title = "Proportion of SAGs in Order"))
gorgD_ord_sys_bubble
ggsave("/Users/aweinheimer/Documents/GORG_defense/gorgdef_man/v3_files/v3_figs/v3_svgs/gorgD_ord_sys_bubble.svg",gorgD_ord_sys_bubble,width = 6.8,height = 6.4)

# subset for SAGs in the Tropics bubble plot
gorgT_ord_sagcount_inplot = gorgT_ord_sys_stats %>% filter(sags > 4)  %>% filter(sys_sags > 1)  %>% filter(system != 0)  %>% 
  select(Order,sags) %>% unique() %>% as.data.frame()
gorgT_ord_sagcount_bar = gorgT_ord_sagcount_inplot %>% ggplot(aes(y=reorder(Order,log10(sags)),x=log10(sags))) + 
  geom_bar(stat="identity",fill="gray40") + theme_minimal()  + scale_x_reverse()
gorgT_ord_sagcount_bar
ggsave("/Users/aweinheimer/Documents/GORG_defense/gorgdef_man/v3_files/v3_figs/v3_svgs/gorgT_ord_sagcount.svg",gorgT_ord_sagcount_bar,width = 6.8,height = 6.4)

# subset for systems in the Tropics bubble plot
gorgT_syscount_inplot = gorgT_ord_sys_stats %>% filter(sags > 4)  %>% filter(sys_sags > 1)  %>% filter(system != 0)  %>% 
  select(system,sys_sags) %>% unique() %>% as.data.frame()
gorgT_syscount_bar = gorgT_syscount_inplot %>% filter(sys_sags > 1) %>% filter(system != 0) %>% ggplot(aes(x=reorder(system,-sys_sags),y=log10(sys_sags))) + 
  geom_bar(stat="identity",fill="gray40") + theme_minimal()  + scale_y_reverse() + theme(axis.text.x = element_text(angle = 90))
gorgT_syscount_inplot
ggsave("/Users/aweinheimer/Documents/GORG_defense/gorgdef_man/v3_files/v3_figs/v3_svgs/gorgT_sys_sagcount.svg",gorgT_syscount_bar,width = 6.8,height = 6.4)


# repeat above steps with Dark
gorgD_ord_sagcount_inplot = gorgD_ord_sys_stats %>% filter(sags > 4)  %>% filter(sys_sags > 1)  %>% filter(system != 0)  %>% 
  select(Order,sags) %>% unique() %>% as.data.frame()
gorgD_ord_sagcount_bar = gorgD_ord_sagcount_inplot %>% ggplot(aes(y=reorder(Order,log10(sags)),x=log10(sags))) + 
  geom_bar(stat="identity",fill="gray40") + theme_minimal()  + scale_x_reverse()
gorgD_ord_sagcount_bar
ggsave("/Users/aweinheimer/Documents/GORG_defense/gorgdef_man/v3_files/v3_figs/v3_svgs/gorgD_ord_sagcount_bar.svg",gorgD_ord_sagcount_bar,width=2.5,height=5.17)

gorgD_syscount_inplot = gorgD_ord_sys_stats %>% filter(sags > 4)  %>% filter(sys_sags > 1)  %>% filter(system != 0)  %>% 
  select(system,sys_sags) %>% unique() %>% as.data.frame()


gorgD_syscount_bar = gorgD_syscount_inplot %>% filter(sys_sags > 1) %>% filter(system != 0) %>% ggplot(aes(x=reorder(system,-sys_sags),y=log10(sys_sags))) + 
  geom_bar(stat="identity",fill="gray40") + theme_minimal()  + scale_y_reverse() + theme(axis.text.x = element_text(angle = 90))
gorgD_syscount_bar

ggsave("/Users/aweinheimer/Documents/GORG_defense/gorgdef_man/v3_files/v3_figs/v3_svgs/gorgD_sys_sagcount_bar.svg",gorgD_syscount_bar,width = 5.6, height  = 1.6, units="in" )



