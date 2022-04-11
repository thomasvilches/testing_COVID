setwd("~/PosDoc/Coronavirus/Cost_eff/")
library(ggplot2)
library(dplyr)
library(rcartocolor)
library(scales)
library(colorspace) 

# parameters --------------------------------------------------------------

pal_col = carto_pal(n = 8, name = "Bold")[c(1, 3, 7, 2)]
tests = c("RT-PCR","Abbott_PanBio","BTNX_Rapid_Response",	"Artron")
scenario_labs = c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4")
population = 14826276 #Ontario https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=1710000501 on March 24, 2022


# Functions ---------------------------------------------------------------
theme_flip <-
  theme(
    axis.text.x = element_text(face = "plain", family = "Helvetica", size = 16),
    #axis.text.y = element_text(face = "bold", family = "Arial", size = 26),
    axis.title.x = element_text(face = "bold", family = "Helvetica", size = 16),
    axis.text.y = element_text(face = "plain", family = "Helvetica", size = 16),
    #axis.text.y = element_text(face = "bold", family = "Arial", size = 26),
    axis.title.y = element_text(face = "bold", family = "Helvetica", size = 16),
    axis.ticks.y = element_line(color = "black", size = 1.0),
    legend.position = "none", 
    legend.text = element_text(family = "Helvetica", size = 18),
    legend.title = element_text(face = "bold", size = 18)
  )

read_file_incidence <- function(index,type,strain = 1,scen = 0,test = 0,eb = 0,size = 100,st2 = "ontario",beta = "121",hi="25",ag="all"){
  
  data.cases1 = read.table(paste0("data_cluster/results_prob_0_",beta,"_herd_immu_",hi,"_idx_",index,"_",st2,"_strain_",strain,"_scen_",scen,"_test_",test,"_eb_",eb,"_size_",size,"/simlevel_",type,"_inc_",ag,".dat"),',',h = T) 
  data.cases1 = data.cases1[,-1]
  
  data.cases2 = read.table(paste0("data_cluster/results_prob_0_",beta,"_herd_immu_",hi,"_idx_",index,"_",st2,"_strain_",strain,"_scen_",scen,"_test_",test,"_eb_",eb,"_size_",size,"/simlevel_",type,"2_inc_",ag,".dat"),',',h = T) 
  data.cases2 = data.cases2[,-1]
  
  data.cases3 = read.table(paste0("data_cluster/results_prob_0_",beta,"_herd_immu_",hi,"_idx_",index,"_",st2,"_strain_",strain,"_scen_",scen,"_test_",test,"_eb_",eb,"_size_",size,"/simlevel_",type,"3_inc_",ag,".dat"),',',h = T) 
  data.cases3 = data.cases3[,-1]
  
  data.cases4 = read.table(paste0("data_cluster/results_prob_0_",beta,"_herd_immu_",hi,"_idx_",index,"_",st2,"_strain_",strain,"_scen_",scen,"_test_",test,"_eb_",eb,"_size_",size,"/simlevel_",type,"4_inc_",ag,".dat"),',',h = T) 
  data.cases4 = data.cases4[,-1]
  
  data.cases5 = read.table(paste0("data_cluster/results_prob_0_",beta,"_herd_immu_",hi,"_idx_",index,"_",st2,"_strain_",strain,"_scen_",scen,"_test_",test,"_eb_",eb,"_size_",size,"/simlevel_",type,"5_inc_",ag,".dat"),',',h = T) 
  data.cases5 = data.cases5[,-1]
  
  data.cases6 = read.table(paste0("data_cluster/results_prob_0_",beta,"_herd_immu_",hi,"_idx_",index,"_",st2,"_strain_",strain,"_scen_",scen,"_test_",test,"_eb_",eb,"_size_",size,"/simlevel_",type,"6_inc_",ag,".dat"),',',h = T) 
  data.cases6 = data.cases6[,-1]
  
  l = list(data.cases1,data.cases2,data.cases3,data.cases4,data.cases5,data.cases6)
  
  return(l[[strain]])
}

# And a function to bootstrap 

fc <- function(d, i){
  return(mean(d[i],na.rm=T))
}

bb_ci <- function(x){
  bx = boot::boot(as.vector(t(x)),fc,R = 500)
  bci=boot::boot.ci(bx,conf=0.95,type = "norm")
  return(c(mean(bx$t[,1]),bci$normal[,c(2,3)]))
}


# Simple plot of incidence ------------------------------------------------
type_d = "lat"
x4.2 = read_file_incidence(1,type_d,1,4,2)
x4.3 = read_file_incidence(1,type_d,1,4,3)
x4.4 = read_file_incidence(1,type_d,1,4,4)

x2.2 = read_file_incidence(1,type_d,1,2,2)
x2.3 = read_file_incidence(1,type_d,1,2,3)
x2.4 = read_file_incidence(1,type_d,1,2,4)

x1 = read_file_incidence(1,type_d,1,1,1)
x3 = read_file_incidence(1,type_d,1,3,1)

nr = nrow(x4.2)


df4.2 = data.frame(mm = rowMeans(x4.2),days = seq(1,nr),scen = rep("4",nr),test = rep(tests[2],nr))
df4.3 = data.frame(mm = rowMeans(x4.3),days = seq(1,nr),scen = rep("4",nr),test = rep(tests[3],nr))
df4.4 = data.frame(mm = rowMeans(x4.4),days = seq(1,nr),scen = rep("4",nr),test = rep(tests[4],nr))


df2.2 = data.frame(mm = rowMeans(x2.2),days = seq(1,nr),scen = rep("2",nr),test = rep(tests[2],nr))
df2.3 = data.frame(mm = rowMeans(x2.3),days = seq(1,nr),scen = rep("2",nr),test = rep(tests[3],nr))
df2.4 = data.frame(mm = rowMeans(x2.4),days = seq(1,nr),scen = rep("2",nr),test = rep(tests[4],nr))

df = rbind(df4.2,df4.3,df4.4,df2.2,df2.3,df2.4)


df1 = data.frame(mm = rowMeans(x1),days = seq(1,nr),scen = rep("1",nr),test = rep(tests[1],nr))
df3 = data.frame(mm = rowMeans(x3),days = seq(1,nr),scen = rep("3",nr),test = rep(tests[1],nr))


ggplot()+
  geom_line(data = df,aes(x = days,y=mm,color = scen),size = 1.2)+facet_grid(.~test)+
  geom_line(data = df1 %>% select(mm,days),aes(x = days,y=mm,color = "1"),size = 1.2)+
  geom_line(data = df3 %>% select(mm,days),aes(x = days,y=mm,color = "3"),size = 1.2)+
  scale_color_manual(values = pal_col)+theme_bw()
  

# Plots with CI -----------------------------------------------------------

create_df <- function(M,scen,tt){
  nr = ncol(M)
  df = data.frame(time = seq(1,nr),mm = M[1,],ci1 = M[2,],ci2 = M[3,],scen = rep(scen,nr),test = rep(tt,nr))
  return(df)
}

type_d = "lat"
index = 4

x4.2 = read_file_incidence(index,type_d,1,4,2)
x4.3 = read_file_incidence(index,type_d,1,4,3)
x4.4 = read_file_incidence(index,type_d,1,4,4)

x2.2 = read_file_incidence(index,type_d,1,2,2)
x2.3 = read_file_incidence(index,type_d,1,2,3)
x2.4 = read_file_incidence(index,type_d,1,2,4)

x1 = read_file_incidence(index,type_d,1,1,1)
x3 = read_file_incidence(index,type_d,1,3,1)


bbx4.2 = create_df(apply(x4.2,1,bb_ci),"4",tests[2])
bbx4.3 = create_df(apply(x4.3,1,bb_ci),"4",tests[3])
bbx4.4 = create_df(apply(x4.4,1,bb_ci),"4",tests[4])

bbx2.2 = create_df(apply(x2.2,1,bb_ci),"2",tests[2])
bbx2.3 = create_df(apply(x2.3,1,bb_ci),"2",tests[3])
bbx2.4 = create_df(apply(x2.4,1,bb_ci),"2",tests[4])

bbx1 = create_df(apply(x1,1,bb_ci),"1",tests[1])
bbx3 = create_df(apply(x3,1,bb_ci),"3",tests[1])


bbx = rbind(bbx4.2,bbx4.3,bbx4.4,bbx2.2,bbx2.3,bbx2.4)


ggplot()+
  geom_line(data = bbx4.2,aes(x = time,y=mm),color = pal_col[1])+
  geom_ribbon(data = bbx4.2,aes(x = time,ymin=ci1,ymax = ci2),fill = pal_col[1],alpha = 0.5)+
  theme_bw()


ggplot()+
  geom_line(data = bbx,aes(x = time,y=mm,color = scen),size = 1.2)+
  geom_ribbon(data = bbx,aes(x = time,ymin=ci1,ymax = ci2,fill = scen),alpha = 0.5)+facet_grid(.~test)+
  geom_line(data = bbx1 %>% select(mm,time,ci1,ci2),aes(x = time,y=mm,color = "1"),size = 1.2)+
  geom_ribbon(data = bbx1 %>% select(mm,time,ci1,ci2),aes(x = time,ymin=ci1,ymax = ci2,fill = "1"),alpha = 0.5)+
  geom_line(data = bbx3 %>% select(mm,time,ci1,ci2),aes(x = time,y=mm,color = "3"),size = 1.2)+
  geom_ribbon(data = bbx3 %>% select(mm,time,ci1,ci2),aes(x = time,ymin=ci1,ymax = ci2,fill = "3"),alpha = 0.5)+
  scale_color_manual(values = pal_col,name = NULL,breaks = c("1","2","3","4"),label = scenario_labs)+scale_fill_manual(values = pal_col,name = NULL,breaks = c("1","2","3","4"),label = scenario_labs)+
  scale_x_continuous(limits = c(0,220),expand=expansion(mult=c(0,0)))+
  scale_y_continuous(expand=expansion(mult=c(0,0)))+
  labs(y = "Incidence of infections",x = "Days")+theme_bw()+theme_flip+
  theme(strip.text.x = element_text(size = 13,face = "bold", hjust = 0),
        strip.background = element_blank(),
        legend.position = "bottom")


ggsave(
  paste0("figures/time-series-",index,"-",type_d,".pdf"),
  device = "pdf",
  width = 9,
  height = 4,
  dpi = 300
)



# Boxplot -----------------------------------------------------------------


create_df <- function(M,scen,tt){
  nr = length(M)
  df = data.frame(mm=M,scen = rep(scen,nr),test = rep(tt,nr))
  return(df)
}

type_d = "lat"
index = 4
x4.2 = read_file_incidence(index,type_d,1,4,2)
x4.3 = read_file_incidence(index,type_d,1,4,3)
x4.4 = read_file_incidence(index,type_d,1,4,4)

x2.2 = read_file_incidence(index,type_d,1,2,2)
x2.3 = read_file_incidence(index,type_d,1,2,3)
x2.4 = read_file_incidence(index,type_d,1,2,4)

x1 = read_file_incidence(index,type_d,1,1,1)
x3 = read_file_incidence(index,type_d,1,3,1)

xx = boot::boot(colSums(x4.2),fc,500)$t[,1]
bbx4.2 = create_df(xx,"4",tests[2])
xx = boot::boot(colSums(x4.3),fc,500)$t[,1]
bbx4.3 = create_df(xx,"4",tests[3])
xx = boot::boot(colSums(x4.4),fc,500)$t[,1]
bbx4.4 = create_df(xx,"4",tests[4])

xx = boot::boot(colSums(x2.2),fc,500)$t[,1]
bbx2.2 = create_df(xx,"2",tests[2])
xx = boot::boot(colSums(x2.3),fc,500)$t[,1]
bbx2.3 = create_df(xx,"2",tests[3])
xx = boot::boot(colSums(x2.4),fc,500)$t[,1]
bbx2.4 = create_df(xx,"2",tests[4])

xx = boot::boot(colSums(x1),fc,500)$t[,1]
bbx1 = create_df(xx,"1",tests[1])
xx = boot::boot(colSums(x3),fc,500)$t[,1]
bbx3 = create_df(xx,"3",tests[1])



bbx = rbind(bbx4.2,bbx4.3,bbx4.4,bbx2.2,bbx2.3,bbx2.4)



ggplot()+
  geom_boxplot(data = bbx,aes(x = scen,y=mm,color = scen,fill=scen,fill=after_scale(colorspace::lighten(fill, .5))),size = 1.2)+facet_grid(.~test)+
  geom_boxplot(data = bbx1 %>% select(mm,scen),aes(x = "1",y=mm,color = "1",fill="1",fill=after_scale(colorspace::lighten(fill, .5))),size = 1.2)+
  geom_boxplot(data = bbx3 %>% select(mm,scen),aes(x = "3",y=mm,color = "3",fill="3",fill=after_scale(colorspace::lighten(fill, .5))),size = 1.2)+
  scale_color_manual(values = pal_col,name = NULL,breaks = c("1","2","3","4"),label = scenario_labs)+scale_fill_manual(values = pal_col,name = NULL,breaks = c("1","2","3","4"),label = scenario_labs)+
  scale_y_continuous(expand=expansion(mult=c(0,0)))+
  labs(y = "Number of infections",x = "Scenarios")+theme_bw()+theme_flip+
  theme(strip.text.x = element_text(size = 13,face = "bold", hjust = 0),
        strip.background = element_blank(),
        legend.position = "none")


ggsave(
  paste0("figures/boxplot-",index,"-",type_d,".pdf"),
  device = "pdf",
  width = 9,
  height = 4,
  dpi = 300
)


# Figure ------------------------------------------------------------------


# Functions ---------------------------------------------------------------
theme_flip <-
  theme(
    axis.text.x = element_text(face = "plain", family = "Helvetica", size = 14,angle=45,hjust=1.0),
    #axis.text.y = element_text(face = "bold", family = "Arial", size = 26),
    axis.title.x = element_text(face = "bold", family = "Helvetica", size = 16),
    axis.text.y = element_text(face = "plain", family = "Helvetica", size = 16),
    #axis.text.y = element_text(face = "bold", family = "Arial", size = 26),
    axis.title.y = element_text(face = "bold", family = "Helvetica", size = 16),
    axis.ticks.y = element_line(color = "black", size = 1.0),
    legend.position = "none", 
    legend.text = element_text(family = "Helvetica", size = 18),
    legend.title = element_text(face = "bold", size = 18)
  )

bb = c("1-4","5-9","10-19","20-49","50-99","100-199","200-499","500-1000")
breaks = bb %>% factor() %>% factor(levels = bb)

aux = c(0.5795052593106524,0.18068968828208243,0.11820672374061501,0.07722637956240554,0.025467236167207672,0.01122172113883589,0.00551251951334341,0.0021704722848576454)

ages = c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75-84","85-99")
ages = ages %>% factor() %>% factor(levels = ages)

agep=c(0.04807822, 0.10498712, 0.12470340, 0.14498051, 0.13137129, 0.12679091, 0.13804896, 0.10292032, 0.05484776, 0.02327152)



df1 = data.frame(breaks,aux)
df2 = data.frame(ages,agep)


ggplot()+
  geom_col(data=df1,aes(x=breaks,y=aux),color = pal_col[1],fill=pal_col[1])+
  scale_x_discrete(name="Size")+
  scale_y_continuous(name="Probability",expand = c(0,0.025,0.0,0.05))+
  theme_bw()+theme_flip

ggsave(
  paste0("workplace_size.pdf"),
  device = "pdf",
  width = 4,
  height = 5,
  dpi = 300
)


ggplot()+
  geom_col(data=df2,aes(x=ages,y=agep),color = pal_col[4],fill=pal_col[4])+
  scale_x_discrete(name="Age group")+
  scale_y_continuous(name="Probability",expand = c(0,0.005,0.0,0.01))+
  theme_bw()+theme_flip

ggsave(
  paste0("agegroup.pdf"),
  device = "pdf",
  width = 4,
  height = 5,
  dpi = 300
)

dev.new(width=4, height=5, unit="in")
