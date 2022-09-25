setwd("/data/thomas-covid/testing_canada/")
library(dplyr)
library(zoo)
library(data.table)
library(colorspace)
library(rcartocolor)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(ggpubr)
library(data.table)

population = 14826276 
pal_col = carto_pal(n = 8, name = "Bold")[c(1, 3, 7, 2,8)]
# Functions and parameters ---------------------------------------------------------------

source("functions_paper.R")
folder = "fmild_0.0_fwork_0.5/"

# Create plot time and box ------------------------------------------------
type_d = "lat"


create_plots("lat",80,folder,2)
create_plots("lat",20,folder,2)
create_plots("lat",0,folder,2)
# 
# create_plots("lat",80,folder,3)
# create_plots("lat",20,folder,3)
# create_plots("lat",0,folder,3)
# create_plots("lat",80,folder,4)
# create_plots("lat",20,folder,4)
# create_plots("lat",0,folder,4)

# Tables outcomes -------------------------------------------------------------------
folderss = c("fmild_0.0_fwork_0.5/","fmild_0.5_fwork_0.5/", "fmild_1.0_fwork_0.5/","fmild_1.0_fwork_1.0/", "fmild_0.5_fwork_1.0/", "fmild_0.0_fwork_1.0/")

for(folder in folderss){
  
  save_tables(80,folder,2)
  save_tables(20,folder,2)
  save_tables(0,folder,2)
}
# save_tables(80,folder,3)
# save_tables(20,folder,3)
# save_tables(0,folder,3)
# save_tables(80,folder,4)
# save_tables(20,folder,4)
# save_tables(0,folder,4)



# Costs -------------------------------------------------------------------
#folderss = c("fmild_0.0_fwork_0.5/","fmild_0.5_fwork_0.5/", "fmild_1.0_fwork_0.5/","fmild_1.0_fwork_1.0/", "fmild_0.5_fwork_1.0/", "fmild_0.0_fwork_1.0/")

folderss = c("fmild_0.0_fwork_0.5/", "fmild_0.0_fwork_1.0/")

for(folder in folderss){
  print(folder)
  #booster 0
  #scenario A
  save_table_cost(1,0,2,folder)
  save_table_cost(1,0,3,folder)
  save_table_cost(1,0,4,folder)
  #scenario B
  save_table_cost(2,0,2,folder)
  save_table_cost(2,0,3,folder)
  save_table_cost(2,0,4,folder)
  #scenario C
  save_table_cost(3,0,2,folder)
  save_table_cost(3,0,3,folder)
  save_table_cost(3,0,4,folder)
  #scenario D
  save_table_cost(4,0,2,folder)
  save_table_cost(4,0,3,folder)
  save_table_cost(4,0,4,folder)
  #scenario E
  save_table_cost(5,0,2,folder)
  save_table_cost(5,0,3,folder)
  save_table_cost(5,0,4,folder)
  
  #booster 20
  #scenario A
  save_table_cost(1,20,2,folder)
  save_table_cost(1,20,3,folder)
  save_table_cost(1,20,4,folder)
  #scenario B
  save_table_cost(2,20,2,folder)
  save_table_cost(2,20,3,folder)
  save_table_cost(2,20,4,folder)
  #scenario C
  save_table_cost(3,20,2,folder)
  save_table_cost(3,20,3,folder)
  save_table_cost(3,20,4,folder)
  #scenario D
  save_table_cost(4,20,2,folder)
  save_table_cost(4,20,3,folder)
  save_table_cost(4,20,4,folder)
  #scenario E
  save_table_cost(5,20,2,folder)
  save_table_cost(5,20,3,folder)
  save_table_cost(5,20,4,folder)
  
  #booster 80
  #scenario A
  save_table_cost(1,80,2,folder)
  save_table_cost(1,80,3,folder)
  save_table_cost(1,80,4,folder)
  #scenario B
  save_table_cost(2,80,2,folder)
  save_table_cost(2,80,3,folder)
  save_table_cost(2,80,4,folder)
  #scenario C
  save_table_cost(3,80,2,folder)
  save_table_cost(3,80,3,folder)
  save_table_cost(3,80,4,folder)
  #scenario D
  save_table_cost(4,80,2,folder)
  save_table_cost(4,80,3,folder)
  save_table_cost(4,80,4,folder)
  #scenario E
  save_table_cost(5,80,2,folder)
  save_table_cost(5,80,3,folder)
  save_table_cost(5,80,4,folder)

}



# Plot ICER ---------------------------------------------------------------


folder = "fmild_0.0_fwork_0.5/"

folder2 = "paper_sims/"
test_idx = 2
booster = 80

df = create_df_4_ICER(1,0,test_idx,folder,folder2)

for(idx in c(2,3,4,5)){
  df = df %>%
    bind_rows(create_df_4_ICER(idx,0,test_idx,folder,folder2))
}

df = df %>% bind_rows(create_df_4_ICER(1,20,test_idx,folder,folder2))

for(idx in c(2,3,4,5)){
  df = df %>%
    bind_rows(create_df_4_ICER(idx,20,test_idx,folder,folder2))
}


df = df %>% bind_rows(create_df_4_ICER(1,80,test_idx,folder,folder2))

for(idx in c(2,3,4,5)){
  df = df %>%
    bind_rows(create_df_4_ICER(idx,80,test_idx,folder,folder2))
}

#df  = df %>% filter(scen != 0)
df$ICER = -df$cost/df$qaly

df$idx = factor(df$idx,levels = c("A","B","C","D","E"))
df$booster = factor(df$booster,levels = c(0,20,80))
df$scen = factor(df$scen,levels = c(0,3,4))
label_scen = c("SC1","SC2","SC3","SC4")

C0=rgb(0.0,0.0,0.0)
C1=rgb(0,0.45,0.74)
C2=rgb(0.9,0.65,0.13)
C3=rgb(0.49,0.18,0.56)
C4=rgb(0.47,0.67,0.19)
# 
# ggplot()+
#   geom_boxplot(data = df, aes(x = scen, y = ICER, color = scen))+
#   scale_color_manual(values = c(C0,C1,C2,C3,C4))+
#   theme_bw()+facet_wrap(booster~idx,scales = "free",ncol = 5)
# 
# ggplot(df,aes(x = scen, y = ICER, color = scen))+
#   geom_boxplot(data = df, aes(x = scen, y = ICER, color = scen))+
#   stat_summary(fun=mean,geom = "point")+
#   stat_summary(fun.data=mean_cl_normal,geom = "errorbar")+
#   scale_color_manual(values = c(C0,C1,C2,C3,C4))+
#   theme_bw()+facet_wrap(booster~idx,scales = "free",ncol = 5)


df %>% group_by(booster,idx,scen) %>%
  summarise(qq1 = quantile(ICER, 0.025,names=F,type=5,na.rm=T),
            qq2 = quantile(ICER, 0.975,names=F,type=5,na.rm=T),
            mm = mean(ICER,na.rm=T),
            mm2 = quantile(ICER, 0.5,names=F,type=5,na.rm=T),
  ) %>% filter(scen != 0) %>%  view
  
df %>% filter(booster == 80) %>% group_by(booster,idx,scen) %>% summarize(m_q = mean(qaly), mc = mean(cost)) %>%
  view

# Costs2 -------------------------------------------------------------------
folderss = c("fmild_0.0_fwork_0.5/","fmild_0.5_fwork_0.5/", "fmild_1.0_fwork_0.5/","fmild_1.0_fwork_1.0/", "fmild_0.5_fwork_1.0/", "fmild_0.0_fwork_1.0/")

for(folder in folderss){
  print(folder)
  #booster 0
  #scenario A
  save_table_cost2(1,0,2,folder)
  #save_table_cost2(1,0,3,folder)
  #save_table_cost2(1,0,4,folder)
  #scenario B
  save_table_cost2(2,0,2,folder)
  #save_table_cost2(2,0,3,folder)
  #save_table_cost2(2,0,4,folder)
  #scenario C
  save_table_cost2(3,0,2,folder)
  #save_table_cost2(3,0,3,folder)
  #save_table_cost2(3,0,4,folder)
  #scenario D
  save_table_cost2(4,0,2,folder)
  #save_table_cost2(4,0,3,folder)
  #save_table_cost2(4,0,4,folder)
  #scenario E
  save_table_cost2(5,0,2,folder)
  #save_table_cost2(5,0,3,folder)
  #save_table_cost2(5,0,4,folder)
  
  #booster 20
  #scenario A
  save_table_cost2(1,20,2,folder)
  #save_table_cost2(1,20,3,folder)
  #save_table_cost2(1,20,4,folder)
  #scenario B
  save_table_cost2(2,20,2,folder)
  #save_table_cost2(2,20,3,folder)
  #save_table_cost2(2,20,4,folder)
  #scenario C
  save_table_cost2(3,20,2,folder)
  #save_table_cost2(3,20,3,folder)
  #save_table_cost2(3,20,4,folder)
  #scenario D
  save_table_cost2(4,20,2,folder)
  #save_table_cost2(4,20,3,folder)
  #save_table_cost2(4,20,4,folder)
  #scenario E
  save_table_cost2(5,20,2,folder)
  #save_table_cost2(5,20,3,folder)
  #save_table_cost2(5,20,4,folder)
  
  #booster 80
  #scenario A
  save_table_cost2(1,80,2,folder)
  #save_table_cost2(1,80,3,folder)
  #save_table_cost2(1,80,4,folder)
  #scenario B
  save_table_cost2(2,80,2,folder)
  #save_table_cost2(2,80,3,folder)
  #save_table_cost2(2,80,4,folder)
  #scenario C
  save_table_cost2(3,80,2,folder)
  #save_table_cost2(3,80,3,folder)
  #save_table_cost2(3,80,4,folder)
  #scenario D
  save_table_cost2(4,80,2,folder)
  #save_table_cost2(4,80,3,folder)
  #save_table_cost2(4,80,4,folder)
  #scenario E
  save_table_cost2(5,80,2,folder)
  #save_table_cost2(5,80,3,folder)
  #save_table_cost2(5,80,4,folder)
  
}


# NMB -------------------------------------------------------------------

 folderss = c("fmild_0.0_fwork_0.5/","fmild_0.5_fwork_0.5/", "fmild_1.0_fwork_0.5/","fmild_1.0_fwork_1.0/", "fmild_0.5_fwork_1.0/", "fmild_0.0_fwork_1.0/")
#folderss = c("fmild_0.0_fwork_0.5/","fmild_0.0_fwork_1.0/")
for(folder in folderss){
  print(folder)
  #booster 0
  #scenario A
  save_table_NMB(1,0,2,folder)
  #save_table_NMB(1,0,3,folder)
  #save_table_NMB(1,0,4,folder)
  #scenario B
  save_table_NMB(2,0,2,folder)
  #save_table_NMB(2,0,3,folder)
  #save_table_NMB(2,0,4,folder)
  #scenario C
  save_table_NMB(3,0,2,folder)
  #save_table_NMB(3,0,3,folder)
  #save_table_NMB(3,0,4,folder)
  #scenario D
  save_table_NMB(4,0,2,folder)
  #save_table_NMB(4,0,3,folder)
  #save_table_NMB(4,0,4,folder)
  #scenario E
  save_table_NMB(5,0,2,folder)
  #save_table_NMB(5,0,3,folder)
  #save_table_NMB(5,0,4,folder)
  
  #booster 20
  #scenario A
  save_table_NMB(1,20,2,folder)
  #save_table_NMB(1,20,3,folder)
  #save_table_NMB(1,20,4,folder)
  #scenario B
  save_table_NMB(2,20,2,folder)
  #save_table_NMB(2,20,3,folder)
  #save_table_NMB(2,20,4,folder)
  #scenario C
  save_table_NMB(3,20,2,folder)
  #save_table_NMB(3,20,3,folder)
  #save_table_NMB(3,20,4,folder)
  #scenario D
  save_table_NMB(4,20,2,folder)
  #save_table_NMB(4,20,3,folder)
  #save_table_NMB(4,20,4,folder)
  #scenario E
  save_table_NMB(5,20,2,folder)
  #save_table_NMB(5,20,3,folder)
  #save_table_NMB(5,20,4,folder)
  
  #booster 80
  #scenario A
  save_table_NMB(1,80,2,folder)
  #save_table_NMB(1,80,3,folder)
  #save_table_NMB(1,80,4,folder)
  #scenario B
  save_table_NMB(2,80,2,folder)
  #save_table_NMB(2,80,3,folder)
  #save_table_NMB(2,80,4,folder)
  #scenario C
  save_table_NMB(3,80,2,folder)
  #save_table_NMB(3,80,3,folder)
  #save_table_NMB(3,80,4,folder)
  #scenario D
  save_table_NMB(4,80,2,folder)
  #save_table_NMB(4,80,3,folder)
  #save_table_NMB(4,80,4,folder)
  #scenario E
  save_table_NMB(5,80,2,folder)
  #save_table_NMB(5,80,3,folder)
  #save_table_NMB(5,80,4,folder)
  
}


# NMB2 -------------------------------------------------------------------

#folderss = c("fmild_0.0_fwork_0.5/","fmild_0.5_fwork_0.5/", "fmild_1.0_fwork_0.5/","fmild_1.0_fwork_1.0/", "fmild_0.5_fwork_1.0/", "fmild_0.0_fwork_1.0/")



folderss = c("fmild_0.0_fwork_0.5/", "fmild_0.0_fwork_1.0/")

for(folder in folderss){
  print(folder)
  #booster 0
  #scenario A
  save_table_NMB2(1,0,2,folder)
  save_table_NMB2(1,0,3,folder)
  save_table_NMB2(1,0,4,folder)
  #scenario B
  save_table_NMB2(2,0,2,folder)
  save_table_NMB2(2,0,3,folder)
  save_table_NMB2(2,0,4,folder)
  #scenario C
  save_table_NMB2(3,0,2,folder)
  save_table_NMB2(3,0,3,folder)
  save_table_NMB2(3,0,4,folder)
  #scenario D
  save_table_NMB2(4,0,2,folder)
  save_table_NMB2(4,0,3,folder)
  save_table_NMB2(4,0,4,folder)
  #scenario E
  save_table_NMB2(5,0,2,folder)
  save_table_NMB2(5,0,3,folder)
  save_table_NMB2(5,0,4,folder)
  
  #booster 20
  #scenario A
  save_table_NMB2(1,20,2,folder)
  save_table_NMB2(1,20,3,folder)
  save_table_NMB2(1,20,4,folder)
  #scenario B
  save_table_NMB2(2,20,2,folder)
  save_table_NMB2(2,20,3,folder)
  save_table_NMB2(2,20,4,folder)
  #scenario C
  save_table_NMB2(3,20,2,folder)
  save_table_NMB2(3,20,3,folder)
  save_table_NMB2(3,20,4,folder)
  #scenario D
  save_table_NMB2(4,20,2,folder)
  save_table_NMB2(4,20,3,folder)
  save_table_NMB2(4,20,4,folder)
  #scenario E
  save_table_NMB2(5,20,2,folder)
  save_table_NMB2(5,20,3,folder)
  save_table_NMB2(5,20,4,folder)
  
  #booster 80
  #scenario A
  save_table_NMB2(1,80,2,folder)
  save_table_NMB2(1,80,3,folder)
  save_table_NMB2(1,80,4,folder)
  #scenario B
  save_table_NMB2(2,80,2,folder)
  save_table_NMB2(2,80,3,folder)
  save_table_NMB2(2,80,4,folder)
  #scenario C
  save_table_NMB2(3,80,2,folder)
  save_table_NMB2(3,80,3,folder)
  save_table_NMB2(3,80,4,folder)
  #scenario D
  save_table_NMB2(4,80,2,folder)
  save_table_NMB2(4,80,3,folder)
  save_table_NMB2(4,80,4,folder)
  #scenario E
  save_table_NMB2(5,80,2,folder)
  save_table_NMB2(5,80,3,folder)
  save_table_NMB2(5,80,4,folder)
  
}




# Plot NMB ---------------------------------------------------------------


folder = "fmild_0.0_fwork_0.5/"

folder2 = "paper_sims/"
test_idx = 2
booster = 0

df = create_df_4_NMB(1,0,test_idx,folder,folder2)

for(idx in c(2,3,4,5)){
  df = df %>%
    bind_rows(create_df_4_NMB(idx,0,test_idx,folder,folder2))
}

df = df %>% bind_rows(create_df_4_NMB(1,20,test_idx,folder,folder2))

for(idx in c(2,3,4,5)){
  df = df %>%
    bind_rows(create_df_4_NMB(idx,20,test_idx,folder,folder2))
}


df = df %>% bind_rows(create_df_4_NMB(1,80,test_idx,folder,folder2))

for(idx in c(2,3,4,5)){
  df = df %>%
    bind_rows(create_df_4_NMB(idx,80,test_idx,folder,folder2))
}

#df  = df %>% filter(scen != 0)
df$ICER = df$cost+30000*df$qaly

df$idx = factor(df$idx,levels = c("A","B","C","D","E"))
df$booster = factor(df$booster,levels = c(0,20,80))
df$scen = factor(df$scen,levels = c(0,1,2,3,4))
label_scen = c("SC1","SC2","SC3","SC4")

C0=rgb(0.0,0.0,0.0)
C1=rgb(0,0.45,0.74)
C2=rgb(0.9,0.65,0.13)
C3=rgb(0.49,0.18,0.56)
C4=rgb(0.47,0.67,0.19)
# 
# ggplot()+
#   geom_boxplot(data = df, aes(x = scen, y = ICER, color = scen))+
#   scale_color_manual(values = c(C0,C1,C2,C3,C4))+
#   theme_bw()+facet_wrap(booster~idx,scales = "free",ncol = 5)
# 
# ggplot(df,aes(x = scen, y = ICER, color = scen))+
#   geom_boxplot(data = df, aes(x = scen, y = ICER, color = scen))+
#   stat_summary(fun=mean,geom = "point")+
#   stat_summary(fun.data=mean_cl_normal,geom = "errorbar")+
#   scale_color_manual(values = c(C0,C1,C2,C3,C4))+
#   theme_bw()+facet_wrap(booster~idx,scales = "free",ncol = 5)

df0 = df %>% filter(scen == 0, idx == "A",booster == 0)

df %>% group_by(booster,idx,scen) %>%
  mutate(incremental = (df0$ICER-ICER)/1e6) %>%
  summarise(qq1 = quantile(incremental, 0.025,names=F,type=5,na.rm=T),
            qq2 = quantile(incremental, 0.975,names=F,type=5,na.rm=T),
            mm = mean(incremental,na.rm=T),
            mm2 = quantile(incremental, 0.5,names=F,type=5,na.rm=T),
  ) %>% view




dd = df %>% group_by(booster,idx,scen) %>%
  mutate(incremental = (df0$ICER-ICER)/1e6)

dd %>% filter(booster == 80, idx == "C",scen == 3) %>% select(incremental) %>%
  ggplot(data = .,aes(x = incremental))+geom_histogram()+
  geom_vline(aes(xintercept = mean(incremental)))

dd %>% filter(booster == 0, idx == "B",scen == 3) %>% pull(incremental) %>%
  wilcox.test(alternative = "two.sided")



x = dd %>% filter(booster == 80, idx == "C",scen == 3) %>% pull(incremental)
y = dd %>% filter(booster == 80, idx == "C",scen == 4) %>% pull(incremental)
wilcox.test(x,y,alternative = "two.sided")


x = dd %>% filter(booster == 80, idx == "C",scen == 3) %>% pull(incremental)
y = dd %>% filter(booster == 80, idx == "C",scen == 4) %>% pull(incremental)
t.test(x,y)

dd %>% filter(booster == 80, idx == "C") %>%
  ggplot(aes(x = incremental, color = scen, fill = scen))+
  geom_density(alpha = 0.3)+theme_bw()+
  labs(x = "iNMB", y = "Density",color = "Scenario", fill = "Scenario",
       title = "Panel C3")



folder = "fmild_0.0_fwork_0.5/"

folder2 = "paper_data/"
test_idx = 2
booster = 0

df = create_df_4_NMB(1,3,test_idx,folder,folder2)


#df  = df %>% filter(scen != 0)
df$ICER = df$cost+30000*df$qaly

df$idx = factor(df$idx,levels = c("A","B","C","D","E"))
df$booster = factor(df$booster,levels = c(0,20,80))
df$scen = factor(df$scen,levels = c(0,1,2,3,4))
label_scen = c("SC1","SC2","SC3","SC4")
