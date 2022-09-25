setwd("/data/thomas-covid/testing_canada/")
library(dplyr)
library(zoo)
library(data.table)
library(colorspace)
library(rcartocolor)
library(ggplot2)
library(tidyverse)
library(ggpubr)
population = 14826276 
pal_col = carto_pal(n = 8, name = "Bold")[c(1, 3, 7, 2,8)]
# Functions ---------------------------------------------------------------


read_file_incidence <- function(index,type,strain = 1,scen = 0,test = 0,eb = 0,size = 50,beta = "0515",folder = ".",st2 = "ontario",hi="10",ag="all"){
  
  data.cases1 = read.table(paste0(folder,"/results_prob_0_",beta,"_herd_immu_",hi,"_idx_",index,"_",st2,"_strain_",strain,"_scen_",scen,"_test_",test,"_eb_",eb,"_size_",size,"/simlevel_",type,"_inc_",ag,".dat"),',',h = T)
  data.cases1 = data.cases1[,-1]

  data.cases2 = read.table(paste0(folder,"/results_prob_0_",beta,"_herd_immu_",hi,"_idx_",index,"_",st2,"_strain_",strain,"_scen_",scen,"_test_",test,"_eb_",eb,"_size_",size,"/simlevel_",type,"2_inc_",ag,".dat"),',',h = T)
  data.cases2 = data.cases2[,-1]

  data.cases3 = read.table(paste0(folder,"/results_prob_0_",beta,"_herd_immu_",hi,"_idx_",index,"_",st2,"_strain_",strain,"_scen_",scen,"_test_",test,"_eb_",eb,"_size_",size,"/simlevel_",type,"3_inc_",ag,".dat"),',',h = T)
  data.cases3 = data.cases3[,-1]


  l = list(data.cases1,data.cases2,data.cases3)

  return(l[[strain]])
  #return(xx)
}


fc <- function(d, i){
  return(mean(d[i],na.rm=T))
}
# R0 ----------------------------------------------------------------------
R0 = read.table("results_prob_0_053_herd_immu_10_idx_98_ontario_strain_3_scen_0_test_0_eb_0_size_50/R01.dat")
mean(R0$V1)
fc <- function(d, i){
  return(mean(d[i],na.rm=T))
}

bb = boot::boot(R0$V1,fc,R=500)
mean(bb$t[,1])
boot::boot.ci(bb,0.95)


# Isolation ---------------------------------------------------------------
Iso = read.table("results_prob_0_042_herd_immu_10_idx_1_ontario_strain_3_scen_4_test_2_eb_0_size_50/niso_f_w.dat")
sum(Iso)
iso=colSums(Iso)
mean(iso)
fc <- function(d, i){
  return(mean(d[i],na.rm=T))
}

bb = boot::boot(iso,fc,R=1000)
mean(bb$t[,1])
boot::boot.ci(bb,0.95)

# Simple plot -------------------------------------------------------------

type_d = "lat"
xx = read_file_incidence(98,type_d,3,0,0,booster,50,"053",folder)
nr = nrow(xx)
sum(xx)/500

df = data.frame(days = seq(1,nr),mm = rowMeans(xx))
df$roll7 = frollmean(df$mm,7)



df %>% filter(days <= 500) %>% ggplot()+
  geom_line(aes(x=days,y=mm),color = "black", size = 1.2)+
  geom_line(aes(x=days,y=roll7),color = "red", size = 1.2)+
  theme_bw()

sum(df$mm)



# PCR ----------------------------------------------------------------------
R0 = read.table("results_prob_0_042_herd_immu_10_idx_1_ontario_strain_3_scen_1_test_1_eb_0_size_100/npcr.dat")
mean(colSums(R0))
fc <- function(d, i){
  return(mean(d[i],na.rm=T))
}

bb = boot::boot(R0$V1,fc,R=500)
mean(bb$t[,1])
boot::boot.ci(bb,0.95)



# plot incidence ----------------------------------------------------------

type_d = "lat"
beta = "042"
x4.2 = read_file_incidence(1,type_d,3,4,2,50,beta,"weird_run")

x2.2 = read_file_incidence(1,type_d,3,2,2,50,beta,"weird_run")
x3.2 = read_file_incidence(1,type_d,3,3,2,50,beta,"weird_run")

x0 = read_file_incidence(1,type_d,3,0,0,50,beta,"weird_run")
x1 = read_file_incidence(2,type_d,3,0,0,50,beta,"weird_run")

nr = nrow(x4.2)


df4.2 = data.frame(mm = frollmean(rowMeans(x4.2),7),days = seq(1,nr),scen = rep("4",nr),test = rep(tests[2],nr))
df2.2 = data.frame(mm = frollmean(rowMeans(x2.2),7),days = seq(1,nr),scen = rep("2",nr),test = rep(tests[2],nr))
df3.2 = data.frame(mm = frollmean(rowMeans(x3.2),7),days = seq(1,nr),scen = rep("3",nr),test = rep(tests[2],nr))

df = rbind(df4.2,df2.2,df3.2)


df1 = data.frame(mm = frollmean(rowMeans(x1),7),days = seq(1,nr),scen = rep("1",nr),test = rep(tests[1],nr))
df0 = data.frame(mm = frollmean(rowMeans(x0),7),days = seq(1,nr),scen = rep("3",nr),test = rep(tests[1],nr))


ggplot()+
  geom_line(data = df,aes(x = days,y=mm,color = scen),size = 1.2)+
  geom_line(data = df0,aes(x = days,y=mm,color = "0"),size = 1.2)+
  geom_line(data = df1,aes(x = days,y=mm,color = "1"),size = 1.2)+
  scale_color_manual(values = pal_col)+theme_bw()

  

# Comparing ---------------------------------------------------------------
tests = c("PCR","Abbott_PanBio","BTNX_Rapid_Response","Artron")

type_d = "lat"
booster = 0
beta = "053"
folder = "fmild_1.0_fwork_1.0/"
test_idx = 2
idx = 1



df_f$scen = factor(df_f$scen,levels = c("0","1","2","3","4"))

# ggplot()+
#   geom_line(data = df_f,aes(x = days,y=mm,color = scen),size = 1.2)+
#   scale_x_continuous(limits = c(0,250))+
#   scale_y_continuous(limits = c(0,5))+
#   facet_grid(.~idx)+
#   geom_vline(xintercept=c(112,224),linetype="dashed")+
#   scale_color_manual(values = pal_col)+theme_bw()


p1 = ggplot()+
  geom_line(data = df_f,aes(x = days,y=mm,color = scen),size = 1.2)+
  facet_grid(.~idx)+
  geom_vline(xintercept=c(112,224),linetype="dashed")+
  scale_y_continuous(name = "Incidence")+
  scale_x_continuous(name = "Days")+
  scale_color_manual(values = pal_col,name = NULL)+theme_bw()
#p1


# Boxplot -----------------------------------------------------------------


######## 
x4.2 = read_file_incidence(2,type_d,3,4,2,booster,50,beta,folder)
x2.2 = read_file_incidence(2,type_d,3,2,2,booster,50,beta,folder)
x3.2 = read_file_incidence(2,type_d,3,3,2,booster,50,beta,folder)
x0 = read_file_incidence(2,type_d,3,0,0,booster,50,beta,folder)
x1 = read_file_incidence(2,type_d,3,1,0,booster,50,beta,folder)

nr = 500

df4.2 = data.frame(mm = boot::boot(colSums(x4.2),fc,R=nr)$t[,1],days = seq(1,nr),scen = rep("4",nr),test = rep(tests[2],nr))
df2.2 = data.frame(mm = boot::boot(colSums(x2.2),fc,R=nr)$t[,1],days = seq(1,nr),scen = rep("2",nr),test = rep(tests[2],nr))
df3.2 = data.frame(mm = boot::boot(colSums(x3.2),fc,R=nr)$t[,1],days = seq(1,nr),scen = rep("3",nr),test = rep(tests[2],nr))

df1 = data.frame(mm = boot::boot(colSums(x1),fc,R=nr)$t[,1],days = seq(1,nr),scen = rep("1",nr),test = rep(tests[1],nr))
df0 = data.frame(mm = boot::boot(colSums(x0),fc,R=nr)$t[,1],days = seq(1,nr),scen = rep("0",nr),test = rep(tests[1],nr))

df = rbind(df4.2,df2.2,df3.2,df0,df1)

df$idx = rep(2,nrow(df))

df_f = rbind(df_f,df)

#########

x4.2 = read_file_incidence(3,type_d,3,4,2,booster,50,beta,folder)
x2.2 = read_file_incidence(3,type_d,3,2,2,booster,50,beta,folder)
x3.2 = read_file_incidence(3,type_d,3,3,2,booster,50,beta,folder)
x0 = read_file_incidence(3,type_d,3,0,0,booster,50,beta,folder)
x1 = read_file_incidence(3,type_d,3,1,0,booster,50,beta,folder)

nr = 500

df4.2 = data.frame(mm = boot::boot(colSums(x4.2),fc,R=nr)$t[,1],days = seq(1,nr),scen = rep("4",nr),test = rep(tests[2],nr))
df2.2 = data.frame(mm = boot::boot(colSums(x2.2),fc,R=nr)$t[,1],days = seq(1,nr),scen = rep("2",nr),test = rep(tests[2],nr))
df3.2 = data.frame(mm = boot::boot(colSums(x3.2),fc,R=nr)$t[,1],days = seq(1,nr),scen = rep("3",nr),test = rep(tests[2],nr))

df1 = data.frame(mm = boot::boot(colSums(x1),fc,R=nr)$t[,1],days = seq(1,nr),scen = rep("1",nr),test = rep(tests[1],nr))
df0 = data.frame(mm = boot::boot(colSums(x0),fc,R=nr)$t[,1],days = seq(1,nr),scen = rep("0",nr),test = rep(tests[1],nr))

df = rbind(df4.2,df2.2,df3.2,df0,df1)

df$idx = rep(3,nrow(df))

df_f = rbind(df_f,df)


#########

x4.2 = read_file_incidence(4,type_d,3,4,2,booster,50,beta,folder)
x2.2 = read_file_incidence(4,type_d,3,2,2,booster,50,beta,folder)
x3.2 = read_file_incidence(4,type_d,3,3,2,booster,50,beta,folder)
x0 = read_file_incidence(4,type_d,3,0,0,booster,50,beta,folder)
x1 = read_file_incidence(4,type_d,3,1,0,booster,50,beta,folder)

nr = 500

df4.2 = data.frame(mm = boot::boot(colSums(x4.2),fc,R=nr)$t[,1],days = seq(1,nr),scen = rep("4",nr),test = rep(tests[2],nr))
df2.2 = data.frame(mm = boot::boot(colSums(x2.2),fc,R=nr)$t[,1],days = seq(1,nr),scen = rep("2",nr),test = rep(tests[2],nr))
df3.2 = data.frame(mm = boot::boot(colSums(x3.2),fc,R=nr)$t[,1],days = seq(1,nr),scen = rep("3",nr),test = rep(tests[2],nr))

df1 = data.frame(mm = boot::boot(colSums(x1),fc,R=nr)$t[,1],days = seq(1,nr),scen = rep("1",nr),test = rep(tests[1],nr))
df0 = data.frame(mm = boot::boot(colSums(x0),fc,R=nr)$t[,1],days = seq(1,nr),scen = rep("0",nr),test = rep(tests[1],nr))

df = rbind(df4.2,df2.2,df3.2,df0,df1)

df$idx = rep(4,nrow(df))

df_f = rbind(df_f,df)


#########

x4.2 = read_file_incidence(5,type_d,3,4,2,booster,50,beta,folder)
x2.2 = read_file_incidence(5,type_d,3,2,2,booster,50,beta,folder)
x3.2 = read_file_incidence(5,type_d,3,3,2,booster,50,beta,folder)
x0 = read_file_incidence(5,type_d,3,0,0,booster,50,beta,folder)
x1 = read_file_incidence(5,type_d,3,1,0,booster,50,beta,folder)

nr = 500

df4.2 = data.frame(mm = boot::boot(colSums(x4.2),fc,R=nr)$t[,1],days = seq(1,nr),scen = rep("4",nr),test = rep(tests[2],nr))
df2.2 = data.frame(mm = boot::boot(colSums(x2.2),fc,R=nr)$t[,1],days = seq(1,nr),scen = rep("2",nr),test = rep(tests[2],nr))
df3.2 = data.frame(mm = boot::boot(colSums(x3.2),fc,R=nr)$t[,1],days = seq(1,nr),scen = rep("3",nr),test = rep(tests[2],nr))

df1 = data.frame(mm = boot::boot(colSums(x1),fc,R=nr)$t[,1],days = seq(1,nr),scen = rep("1",nr),test = rep(tests[1],nr))
df0 = data.frame(mm = boot::boot(colSums(x0),fc,R=nr)$t[,1],days = seq(1,nr),scen = rep("0",nr),test = rep(tests[1],nr))

df = rbind(df4.2,df2.2,df3.2,df0,df1)

df$idx = rep(5,nrow(df))

df_f = rbind(df_f,df)

df_f$scen = factor(df_f$scen,levels = c("0","1","2","3","4"))


p2 = ggplot()+
  geom_boxplot(data = df_f,aes(x = scen,y=mm,color = scen,fill=after_scale(colorspace::lighten(color, .5))),notch = T,size = 1.2)+
  facet_grid(.~idx)+
  scale_y_continuous(name = "Total infections")+
  scale_x_discrete(name = "Scenario")+
  scale_color_manual(values = pal_col,name = NULL)+theme_bw()+
  scale_fill_manual(values = pal_col,name = NULL)+theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.y = element_text(angle = 90,hjust=0.5))
#p2


# Arrange Plots -----------------------------------------------------------

ggarrange(p1,p2,nrow=2)

# plot incidence compare two----------------------------------------------------------

type_d = "lat"

x0 = read_file_incidence(0,type_d,3,999,0,80,50,"049",".")
x1 = read_file_incidence(1,type_d,3,999,0,80,50,"049",".")

nr = nrow(x1)


df1 = data.frame(mm = frollmean(rowMeans(x1),7),days = seq(1,nr),scen = rep("1",nr),test = rep(tests[1],nr))
df0 = data.frame(mm = frollmean(rowMeans(x0),7),days = seq(1,nr),scen = rep("0",nr),test = rep(tests[1],nr))


ggplot()+
  geom_line(data = df0,aes(x = days,y=mm,color = scen),size = 1.2)+
  geom_line(data = df1,aes(x = days,y=mm,color = scen),size = 1.2)+
  scale_color_manual(values = pal_col)+theme_bw()

sum(df1$mm,na.rm=T)
sum(df0$mm,na.rm=T)

# Some tests --------------------------------------------------------------
library(dplyr)
library(ggplot2)

type_d = "lat"
booster = 80
beta = "049"
folder = "."

R0 = read.table("results_prob_0_049_herd_immu_10_idx_4_ontario_strain_3_scen_2_test_2_eb_80_size_50/npcr.dat")
mean(colSums(R0))

R0 = read.table("results_prob_0_049_herd_immu_10_idx_4_ontario_strain_3_scen_2_test_2_eb_80_size_50/nra.dat")
mean(colSums(R0[-nrow(R0),]))



type_d = "inf"
#index,type,strain = 1,scen = 0,test = 0,eb = 0,size = 50,beta = "045",st2 = "ontario",hi="10",ag="all")
x0 = read_file_incidence(4,type_d,3,2,2,booster,50,beta,folder)
a=mean(colSums(x0))
x1 =read_file_incidence(4,type_d,3,2,2,booster,50,beta,folder,"ontario","10","ag1")
b = mean(colSums(x1))
a-b


aux = seq(112,224)
type_d = "mild"
#index,type,strain = 1,scen = 0,test = 0,eb = 0,size = 50,beta = "045",st2 = "ontario",hi="10",ag="all")
x0 = read_file_incidence(4,type_d,3,1,0,booster,50,beta,folder)
a=mean(colSums(x0[aux,]))
x1 =read_file_incidence(4,type_d,3,1,0,booster,50,beta,folder,"ontario","10","ag1")
b = mean(colSums(x1[aux,]))
a-b


xx = x0-x1
aa = seq(1,500)[colSums(xx == R0) != 365]

nn = xx[,aa]
rr = R0[,aa]

colnames(nn) = c("sim1","sim2","sim3")
colnames(rr) = c("sim1","sim2","sim3")

nnn = data.frame(ii = rep(seq(1,nrow(nn)),3),stack(nn))
rrr = data.frame(ii = rep(seq(1,nrow(rr)),3),stack(rr))

df = inner_join(nnn,rrr,by = c("ind","ii"))

df %>% filter(ind == "sim1") %>% ggplot()+geom_line(aes(x=ii,y=values.y,color = ind),linetype = "solid")+
  geom_line(aes(x=ii,y=values.x,color = ind),linetype = "dashed")+
  geom_point(data = df%>%filter(ind == "sim1",values.x != values.y),
             aes(x=ii,y = values.y),shape=1,size=3)+
  scale_y_continuous(limits=c(75,100))
  

type_d = "mild"
#index,type,strain = 1,scen = 0,test = 0,eb = 0,size = 50,beta = "045",st2 = "ontario",hi="10",ag="all")
x0 = read_file_incidence(1,type_d,3,1,0,0,50,"045","ontario","10","all")
a=mean(colSums(x0[-nrow(x0),]))
x1 = read_file_incidence(1,type_d,3,1,0,0,50,"045","ontario","10","ag1")
b = mean(colSums(x1[-nrow(x1),]))
a-b




# Plot 1 scen -------------------------------------------------------------


type_d = "lat"
booster = 80
beta = "049"
idx = 4
folder = "."
x4.2 = read_file_incidence(idx,type_d,3,4,2,booster,50,beta,folder)
x2.2 = read_file_incidence(idx,type_d,3,2,2,booster,50,beta,folder)
x3.2 = read_file_incidence(idx,type_d,3,3,2,booster,50,beta,folder)
x0 = read_file_incidence(idx,type_d,3,0,0,booster,50,beta,folder)
x1 = read_file_incidence(idx,type_d,3,1,0,booster,50,beta,folder)
nr = nrow(x4.2)


df4.2 = data.frame(mm = frollmean(rowMeans(x4.2),7),days = seq(1,nr),scen = rep("4",nr),test = rep(tests[2],nr))
df2.2 = data.frame(mm = frollmean(rowMeans(x2.2),7),days = seq(1,nr),scen = rep("2",nr),test = rep(tests[2],nr))
df3.2 = data.frame(mm = frollmean(rowMeans(x3.2),7),days = seq(1,nr),scen = rep("3",nr),test = rep(tests[2],nr))

df1 = data.frame(mm = frollmean(rowMeans(x1),7),days = seq(1,nr),scen = rep("1",nr),test = rep(tests[1],nr))
df0 = data.frame(mm = frollmean(rowMeans(x0),7),days = seq(1,nr),scen = rep("0",nr),test = rep(tests[1],nr))

df = rbind(df4.2,df2.2,df3.2,df0,df1)

df$idx = rep(idx,nrow(df))

ggplot()+
  geom_line(data = df,aes(x = days,y=mm,color = scen),size = 1.2)+
  scale_x_continuous(limits = c(0,250))+
  scale_y_continuous(limits = c(0,5))+
  facet_grid(.~idx)+
  geom_vline(xintercept=c(112,224),linetype="dashed")+
  scale_color_manual(values = pal_col)+theme_bw()


ggplot()+
  geom_line(data = df,aes(x = days,y=mm,color = scen),size = 1.2)+
  facet_grid(.~idx)+
  geom_vline(xintercept=c(112,224),linetype="dashed")+
  scale_color_manual(values = pal_col)+theme_bw()

df %>% group_by(scen) %>% summarize(total = sum(mm,na.rm=T))

aux = seq(112,224)
df %>% filter(days %in% aux) %>% group_by(scen) %>% summarize(total = sum(mm,na.rm=T))

