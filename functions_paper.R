
source("cost_parameters.R")
tests = c("PCR","Abbott_PanBio","BTNX_Rapid_Response","Artron")

idx_letter = c("C","A","B","D","E")

df_time = data.frame(x=c(1,365,1,112,1,224,112,224,112,365),
                     idx=c(1,1,2,2,3,3,4,4,5,5),
                     type=c("Beginning","End","Beginning","End","Beginning","End",
                            "Beginning","End","Beginning","End"))

df_time$idx = factor(df_time$idx,levels=c(1,2,3,4,5))
df_time$idx_letter = idx_letter[df_time$idx]

figure_type = "pdf"


if(!file.exists("figures")){
  dir.create("figures")
}


if(!file.exists("ICER_tables")){
  dir.create("ICER_tables")
}


if(!file.exists("tables")){
  dir.create("tables")
}

# function that reads the data
read_file_incidence <- function(index,type,strain = 1,scen = 0,test = 0,eb = 0,size = 50,beta = "053",
                                folder = ".",ag="all",st2 = "ontario",hi="10"){
  
  data.cases1 = read.table(paste0(folder,"/results_prob_0_",beta,"_herd_immu_",hi,"_idx_",index,
                                  "_",st2,"_strain_",strain,"_scen_",scen,"_test_",test,"_eb_",eb,
                                  "_size_",size,"/simlevel_",type,"_inc_",ag,".dat"),',',h = T)
  data.cases1 = data.cases1[,-1]
  
  data.cases2 = read.table(paste0(folder,"/results_prob_0_",beta,"_herd_immu_",hi,"_idx_",index,"_",
                                  st2,"_strain_",strain,"_scen_",scen,"_test_",test,"_eb_",eb,"_size_",
                                  size,"/simlevel_",type,"2_inc_",ag,".dat"),',',h = T)
  data.cases2 = data.cases2[,-1]
  
  data.cases3 = read.table(paste0(folder,"/results_prob_0_",beta,"_herd_immu_",hi,"_idx_",index,"_",st2,
                                  "_strain_",strain,"_scen_",scen,"_test_",test,"_eb_",eb,"_size_",size,
                                  "/simlevel_",type,"3_inc_",ag,".dat"),',',h = T)
  data.cases3 = data.cases3[,-1]
  
  
  l = list(data.cases1,data.cases2,data.cases3)
  
  return(l[[strain]])
  #return(xx)
}

# mean for bootstrap
fc <- function(d, i){
  return(mean(d[i],na.rm=T))
}

#used for temporal series plot
time_average_return <- function(idx,type_d,test_idx,booster,folder,ag="all",size_w=50,beta="053"){
  
  x4.2 = read_file_incidence(idx,type_d,3,4,test_idx,booster,size_w,beta,folder)
  x2.2 = read_file_incidence(idx,type_d,3,2,test_idx,booster,size_w,beta,folder)
  x3.2 = read_file_incidence(idx,type_d,3,3,test_idx,booster,size_w,beta,folder)
  x0 = read_file_incidence(1,type_d,3,0,0,booster,size_w,beta,folder)
  x1 = read_file_incidence(idx,type_d,3,1,0,booster,size_w,beta,folder)
  nr = nrow(x4.2)
  
  
  df4.2 = data.frame(mm = frollmean(rowMeans(x4.2),7),days = seq(1,nr),scen = rep("4",nr),test = rep(tests[2],nr))
  df2.2 = data.frame(mm = frollmean(rowMeans(x2.2),7),days = seq(1,nr),scen = rep("2",nr),test = rep(tests[2],nr))
  df3.2 = data.frame(mm = frollmean(rowMeans(x3.2),7),days = seq(1,nr),scen = rep("3",nr),test = rep(tests[2],nr))
  
  df1 = data.frame(mm = frollmean(rowMeans(x1),7),days = seq(1,nr),scen = rep("1",nr),test = rep(tests[1],nr))
  df0 = data.frame(mm = frollmean(rowMeans(x0),7),days = seq(1,nr),scen = rep("0",nr),test = rep(tests[1],nr))
  
  df = rbind(df4.2,df2.2,df3.2,df0,df1)
  #df = rbind(df4.2,df3.2,df0)
  
  df$idx = rep(idx,nrow(df))
  df$scen = factor(df$scen,levels = c("0","1","2","3","4"))
  
  return(df)
}

# read and return total outcome in each simulation with type, scen and idx

total_outcome <- function(idx,type_d,strain,scen,test_idx,booster,folder,ag="all",size_w=50,beta="053"){
  
  x = colSums(read_file_incidence(idx,type_d,strain,scen,test_idx,booster,size_w,beta,folder))
 
  nr = length(x)
  
  df = data.frame(total = x,sim = seq(1,nr),scen = rep(scen,nr),type = rep(type_d,nr))
  
  df$idx = rep(idx_letter[idx],nrow(df))
  
  return(df)
}

# Function with boostrap, can be used later

boot_scenarios <- function(idx,type_d,test_idx,booster,folder,ag="all",size_w=50,beta="053"){
  
  x4.2 = read_file_incidence(idx,type_d,3,4,test_idx,booster,size_w,beta,folder,ag)
  x2.2 = read_file_incidence(idx,type_d,3,2,test_idx,booster,size_w,beta,folder,ag)
  x3.2 = read_file_incidence(idx,type_d,3,3,test_idx,booster,size_w,beta,folder,ag)
  x0 = read_file_incidence(1,type_d,3,0,0,booster,size_w,beta,folder,ag)
  x1 = read_file_incidence(idx,type_d,3,1,0,booster,size_w,beta,folder,ag)
  
  nr = 500
  
  df4.2 = data.frame(mm = boot::boot(colSums(x4.2),fc,R=nr)$t[,1],days = seq(1,nr),scen = rep("4",nr),test = rep(tests[2],nr))
  df2.2 = data.frame(mm = boot::boot(colSums(x2.2),fc,R=nr)$t[,1],days = seq(1,nr),scen = rep("2",nr),test = rep(tests[2],nr))
  df3.2 = data.frame(mm = boot::boot(colSums(x3.2),fc,R=nr)$t[,1],days = seq(1,nr),scen = rep("3",nr),test = rep(tests[2],nr))
  
  df1 = data.frame(mm = boot::boot(colSums(x1),fc,R=nr)$t[,1],days = seq(1,nr),scen = rep("1",nr),test = rep(tests[1],nr))
  df0 = data.frame(mm = boot::boot(colSums(x0),fc,R=nr)$t[,1],days = seq(1,nr),scen = rep("0",nr),test = rep(tests[1],nr))
  
  df = rbind(df4.2,df2.2,df3.2,df0,df1)
  #df = rbind(df4.2,df3.2,df0)
  df$idx = rep(idx,nrow(df))
  df$scen = factor(df$scen,levels = c("0","1","2","3","4"))
  return(df)
}

# Boostrap any vector and return the MEAN and CI

bootstrapci <- function(x){
  
  bx = boot::boot(x,fc,R=500)
  
  med = mean(bx$t[,1])
  ci = boot::boot.ci(bx,0.95,"bca")$bca[c(4,5)]
  
  return(c(med,ci))
}

create_text <- function(x){
  texto=paste0(format(x[1],scientific = T,digits=4)," (CI95%: ",format(x[2],scientific = T,digits=4)," to ",format(x[3],scientific = T,digits=4),")")
  return(texto)
}

# Function that creates the two plots (time-series and boxplot)

create_plots <- function(type_d,booster,folder,test_idx,beta="053"){
  
  set.seed(12534)
  df1 = time_average_return(1,type_d,test_idx,booster,folder) %>% 
    bind_rows(time_average_return(2,type_d,test_idx,booster,folder)) %>%
    bind_rows(time_average_return(3,type_d,test_idx,booster,folder)) %>%
    bind_rows(time_average_return(4,type_d,test_idx,booster,folder)) %>%  
    bind_rows(time_average_return(5,type_d,test_idx,booster,folder)) 
  
  
  df1$idx = factor(df1$idx,levels = c(1,2,3,4,5))
  df1$idx_letter = idx_letter[df1$idx]
  
  p1 = ggplot()+
    geom_line(data = df1,aes(x = days,y=mm,color = scen),size = 1.2)+
    geom_vline(data=df_time,aes(xintercept = x,linetype=type))+
    facet_grid(~idx_letter)+
    scale_y_continuous(name = "Incidence")+
    scale_x_continuous(name = "Days")+
    scale_linetype_discrete(name = "Testing")+
    scale_color_manual(values = pal_col,name ="Scenario")+theme_bw()+
    theme(strip.background = element_blank(),
          strip.text = element_text(size=14,face="bold",hjust=0.0),
          axis.text.y = element_text(angle = 90,hjust=0.5,size=12),
          axis.text.x = element_text(size=12),
          axis.title.y = element_text(size=14,face="bold"),
          axis.title.x = element_text(size=14,face="bold"))
  
  
  
  df2 = boot_scenarios(1,type_d,test_idx,booster,folder) %>%
    bind_rows(boot_scenarios(2,type_d,test_idx,booster,folder)) %>%
    bind_rows(boot_scenarios(3,type_d,test_idx,booster,folder)) %>%
    bind_rows(boot_scenarios(4,type_d,test_idx,booster,folder)) %>%
    bind_rows(boot_scenarios(5,type_d,test_idx,booster,folder))
  
  df2$idx_letter = idx_letter[df2$idx]
  
  p2 = ggplot()+
    geom_boxplot(data = df2,aes(x = scen,y=mm,color = scen,fill=after_scale(colorspace::lighten(color, .5))),notch = T,size = 1.2)+
    facet_grid(.~idx_letter)+
    scale_y_continuous(name = "Total infections")+
    scale_x_discrete()+
    scale_color_manual(values = pal_col,name = "Scenario")+theme_bw()+
    scale_fill_manual(values = pal_col,name = "Scenario")+theme_bw()+
    theme(strip.background = element_blank(),
          strip.text = element_text(size=14,face="bold",hjust=0.0),
          axis.text.y = element_text(angle = 90,hjust=0.5,size=10),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=14,face="bold"),
          axis.ticks.length.x.bottom = unit(0.0,"mm"),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),)
  
  
  
  if(!file.exists(paste0("figures/",folder))){
    dir.create(paste0("figures/",folder))
  }
  
  path1 = paste0("figures/",folder,"/time-series_",type_d,"_",booster,"_",test_idx,".",figure_type)
  
  ggsave(path1,p1,device = figure_type,dpi = 300,width = 23,height = 9,units = "cm")
  
  
  path2 = paste0("figures/",folder,"/boxplot_",type_d,"_",booster,"_",test_idx,".",figure_type)
  
  ggsave(path2,p2,device = figure_type,dpi = 300,width = 23,height = 9,units = "cm")
}


# Function that calculates the cost

costs <- function(idx,scen,test_idx,booster,folder,size_w=50,strain = 3,beta="053",st2="ontario",hi="10"){
  set.seed(15234)
  #isolation cost
  isodays = read.table(paste0(folder,"results_prob_0_",beta,"_herd_immu_",hi,"_idx_",idx,"_",st2,"_strain_",strain,"_scen_",scen,"_test_",test_idx,"_eb_",booster,"_size_",size_w,"/totalisog.dat"),h = F)
  cost_iso = (as.matrix(isodays) %*% income_per_day)[,1]
  #cost RA
  cost_ra = read.table(paste0(folder,"results_prob_0_",beta,"_herd_immu_",hi,"_idx_",idx,"_",st2,"_strain_",strain,"_scen_",scen,"_test_",test_idx,"_eb_",booster,"_size_",size_w,"/nra.dat"),
                       h = F) %>% colSums()*ra_cost
                        
  #cost PCR
  cost_pcr = read.table(paste0(folder,"results_prob_0_",beta,"_herd_immu_",hi,"_idx_",idx,"_",st2,"_strain_",strain,"_scen_",scen,"_test_",test_idx,"_eb_",booster,"_size_",size_w,"/npcr.dat"),
                       h = F) %>% colSums()*pcr_cost
    
  #read outcomes files
  sev = colSums(read_file_incidence(idx,"inf",strain,scen,test_idx,booster,size_w,beta,folder))
  mild = colSums(read_file_incidence(idx,"mild",strain,scen,test_idx,booster,size_w,beta,folder))
  
  hos = colSums(read_file_incidence(idx,"hos",strain,scen,test_idx,booster,size_w,beta,folder))
  icu = colSums(read_file_incidence(idx,"icu",strain,scen,test_idx,booster,size_w,beta,folder))
  
  
  emerg_cost = sev*emergency_per_severe_cost # cost emergency room per simulation
  outp_cost = (mild+sev)*outpatient_cost # cost outpatient visit per simulation
  hos_cost = hos*cost_hosp #cost hospitalization per simulation
  icu_cost = icu*cost_icu # cost icu per simulation
  long_cost = (hos+icu)*cost_long_covid #cost long covid per simulation
  
  #loss of productivity due to hospitalization is per age group
  
  ag_groups = c("ag1","ag2","ag3","ag4","ag5","ag6","ag7")
  xx = seq(1,length(ag_groups))
  
  dum = lapply(xx, function(x) colSums(read_file_incidence(idx,"hos",strain,scen,test_idx,booster,size_w,beta,folder,ag_groups[x]))*prod_loss_hosp[x])
  hos_prod_cost = Reduce("+",dum)
  
  dum = lapply(xx, function(x) colSums(read_file_incidence(idx,"icu",strain,scen,test_idx,booster,size_w,beta,folder,ag_groups[x]))*prod_loss_hosp[x])
  icu_prod_cost = Reduce("+",dum)
  
  # Cost loss of productivity for death
  dum = lapply(xx, function(x) colSums(read_file_incidence(idx,"ded",strain,scen,test_idx,booster,size_w,beta,folder,ag_groups[x]))*prod_loss_death[x])
  death_prod_cost = Reduce("+",dum)
  
  total_cost = cost_iso+cost_ra+cost_pcr+emerg_cost+outp_cost+hos_cost+icu_cost+long_cost+
    hos_prod_cost+icu_prod_cost+death_prod_cost
  
  #calculate QALY Loss
  
  qaly_symp = (mild+sev)*qaly_symp_case #symp case
  qaly_hos = (hos+icu)*qaly_hosp_case #hospitalized case
  #QALY due to death
  dum = lapply(xx, function(x) colSums(read_file_incidence(idx,"ded",strain,scen,test_idx,booster,size_w,beta,folder,ag_groups[x]))*qaly_death_case[x])
  qaly_death = Reduce("+",dum)
  
  #Total QALY loss
  total_qaly = qaly_symp+qaly_hos+qaly_death
  
  nr = length(total_cost)
  df = data.frame(cost = total_cost,qaly = total_qaly, sim = seq(1,nr), scen = rep(scen,nr),idx = rep(idx,nr))
  return(df)
}


# create the table for a given scenario (A,B,C,D,E), booster, strain, folder, and test

create_table <- function(idx,test_idx,booster,folder,size_w=50,strain = 3,beta="053",st2="ontario",hi="10"){
  
  set.seed(12534)
  
  df = total_outcome(1,"lat",strain,0,0,booster,folder) %>%
    bind_rows(total_outcome(idx,"lat",strain,1,0,booster,folder)) %>%
    bind_rows(total_outcome(idx,"lat",strain,2,test_idx,booster,folder)) %>%
    bind_rows(total_outcome(idx,"lat",strain,3,test_idx,booster,folder)) %>%
    bind_rows(total_outcome(idx,"lat",strain,4,test_idx,booster,folder)) %>%
    ##
    bind_rows(total_outcome(1,"hos",strain,0,0,booster,folder)) %>%
    bind_rows(total_outcome(idx,"hos",strain,1,0,booster,folder)) %>%
    bind_rows(total_outcome(idx,"hos",strain,2,test_idx,booster,folder)) %>%
    bind_rows(total_outcome(idx,"hos",strain,3,test_idx,booster,folder)) %>%
    bind_rows(total_outcome(idx,"hos",strain,4,test_idx,booster,folder)) %>%
    ##
    bind_rows(total_outcome(1,"icu",strain,0,0,booster,folder)) %>%
    bind_rows(total_outcome(idx,"icu",strain,1,0,booster,folder)) %>%
    bind_rows(total_outcome(idx,"icu",strain,2,test_idx,booster,folder)) %>%
    bind_rows(total_outcome(idx,"icu",strain,3,test_idx,booster,folder)) %>%
    bind_rows(total_outcome(idx,"icu",strain,4,test_idx,booster,folder)) %>%
    ##
    bind_rows(total_outcome(1,"ded",strain,0,0,booster,folder)) %>%
    bind_rows(total_outcome(idx,"ded",strain,1,0,booster,folder)) %>%
    bind_rows(total_outcome(idx,"ded",strain,2,test_idx,booster,folder)) %>%
    bind_rows(total_outcome(idx,"ded",strain,3,test_idx,booster,folder)) %>%
    bind_rows(total_outcome(idx,"ded",strain,4,test_idx,booster,folder)) 
  
  
  
  df  = df %>% group_by(idx,scen,type) %>% summarise(result = create_text(bootstrapci(total))) %>%
    pivot_wider(names_from = type, values_from = result)
  
  return(df)
}


save_tables <- function(booster,folder,test_idx){
  
  dff = create_table(1,test_idx,booster,folder) %>%
    bind_rows(create_table(2,test_idx,booster,folder)) %>%
    bind_rows(create_table(3,test_idx,booster,folder)) %>%
    bind_rows(create_table(4,test_idx,booster,folder)) %>%
    bind_rows(create_table(5,test_idx,booster,folder)) 
  
  names.c = c("Testing timeline","Testing scenario","Infections","Hospitalizations","ICU","Deaths")
  colnames(dff) = names.c
  
  
  if(!file.exists(paste0("tables"))){
    dir.create(paste0("tables/"))
  }
  if(!file.exists(paste0("tables/",folder))){
    dir.create(paste0("tables/",folder))
  }
  
  write.csv(dff,paste0("tables/",folder,"table_outcomes_",test_idx,"_",booster,".csv"),row.names = F)
}



save_table_cost <- function(idx,booster,test_idx,folder){
  
  folder2 = paste0("",folder)
  df_base = costs(1,0,0,booster,folder2)
  
  df0 = costs(1,0,0,booster,folder2)
  #df1 = costs(idx,1,0,booster,folder)
  #df2 = costs(idx,2,test_idx,booster,folder)
  df3 = costs(idx,3,test_idx,booster,folder2)
  df4 = costs(idx,4,test_idx,booster,folder2)
  df5 = costs(idx,5,test_idx,booster,folder2)
  
  cost0 = boot::boot(df0$cost-df_base$cost,fc,R=500)$t[,1]
  qaly0 = boot::boot(df0$qaly-df_base$qaly,fc,R=500)$t[,1]
  
  #cost1 = boot::boot(df1$cost-df_base$cost,fc,R=500)$t[,1]
  #qaly1 = boot::boot(df1$qaly-df_base$qaly,fc,R=500)$t[,1]
  
  #cost2 = boot::boot(df2$cost-df_base$cost,fc,R=500)$t[,1]
  #qaly2 = boot::boot(df2$qaly-df_base$qaly,fc,R=500)$t[,1]
  
  cost3 = boot::boot(df3$cost-df_base$cost,fc,R=500)$t[,1]
  qaly3 = boot::boot(df3$qaly-df_base$qaly,fc,R=500)$t[,1]
  
  cost4 = boot::boot(df4$cost-df_base$cost,fc,R=500)$t[,1]
  qaly4 = boot::boot(df4$qaly-df_base$qaly,fc,R=500)$t[,1]
  
  cost5 = boot::boot(df5$cost-df_base$cost,fc,R=500)$t[,1]
  qaly5 = boot::boot(df5$qaly-df_base$qaly,fc,R=500)$t[,1]
  
  if(!file.exists(paste0("ICER_tables"))){
    dir.create(paste0("ICER_tables"))
  }
  
  if(!file.exists(paste0("ICER_tables/",folder))){
    dir.create(paste0("ICER_tables/",folder))
  }
  
  
  if(!file.exists(paste0("ICER_tables/",folder,"tables_",booster,"_",test_idx))){
    dir.create(paste0("ICER_tables/",folder,"tables_",booster,"_",test_idx))
  }
  
  #M = cbind(cost0,qaly0,cost1,qaly1,cost2,qaly2,cost3,qaly3,cost4,qaly4)
  M = cbind(cost0,qaly0,cost3,qaly3,cost4,qaly4, cost5, qaly5)
  
  write.table(M,paste0("ICER_tables/",folder,"tables_",booster,"_",test_idx,"/table_",idx_letter[idx],".csv"),row.names = F,col.names = F)
}

# function for plot

create_df_4_ICER <- function(idx,booster,test_idx,folder,folder2="./"){
  
  M = read.table(paste0(folder2,"ICER_tables/",folder,"tables_",booster,"_",test_idx,"/table_",idx_letter[idx],".csv"),h=F)
  
  out <- data.frame(cost = unlist(M[c(TRUE, FALSE)]),
                    qaly = unlist(M[c(FALSE, TRUE)]),
                    scen = c(rep(0,nrow(M)),rep(3,nrow(M)),
                             rep(4,nrow(M))))
  out$booster = rep(booster,nrow(out))
  out$idx = rep(idx_letter[idx],nrow(out))
  return(out)
  
}


# function for plot

create_df_4_NMB <- function(idx,booster,test_idx,folder,folder2="./"){
  
  M = read.table(paste0(folder2,"NMB_tables2/",folder,"tables_",booster,"_",test_idx,"/table_",idx_letter[idx],".csv"),h=F)
  
  out <- data.frame(cost = unlist(M[c(TRUE, FALSE)]),
                    qaly = unlist(M[c(FALSE, TRUE)]),
                    scen = c(rep(0,nrow(M)),rep(3,nrow(M)),
                             rep(4,nrow(M))))
  out$booster = rep(booster,nrow(out))
  out$idx = rep(idx_letter[idx],nrow(out))
  return(out)
  
}



save_table_cost2 <- function(idx,booster,test_idx,folder){
  
  folder2 = paste0(folder)
  df_base = costs(1,3,test_idx,0,folder2)
  
  df0 = costs(1,0,0,booster,folder2)
  df1 = costs(idx,1,0,booster,folder2)
  df2 = costs(idx,2,test_idx,booster,folder2)
  df3 = costs(idx,3,test_idx,booster,folder2)
  df4 = costs(idx,4,test_idx,booster,folder2)
  
  cost0 = boot::boot(df0$cost-df_base$cost,fc,R=500)$t[,1]
  qaly0 = boot::boot(df0$qaly-df_base$qaly,fc,R=500)$t[,1]
  
  cost1 = boot::boot(df1$cost-df_base$cost,fc,R=500)$t[,1]
  qaly1 = boot::boot(df1$qaly-df_base$qaly,fc,R=500)$t[,1]
  
  cost2 = boot::boot(df2$cost-df_base$cost,fc,R=500)$t[,1]
  qaly2 = boot::boot(df2$qaly-df_base$qaly,fc,R=500)$t[,1]
  
  cost3 = boot::boot(df3$cost-df_base$cost,fc,R=500)$t[,1]
  qaly3 = boot::boot(df3$qaly-df_base$qaly,fc,R=500)$t[,1]
  
  cost4 = boot::boot(df4$cost-df_base$cost,fc,R=500)$t[,1]
  qaly4 = boot::boot(df4$qaly-df_base$qaly,fc,R=500)$t[,1]
  
  
  if(!file.exists(paste0("ICER_tables_2"))){
    dir.create(paste0("ICER_tables_2"))
  }
  if(!file.exists(paste0("ICER_tables_2/",folder))){
    dir.create(paste0("ICER_tables_2/",folder))
  }
  
  
  if(!file.exists(paste0("ICER_tables_2/",folder,"tables_",booster,"_",test_idx))){
    dir.create(paste0("ICER_tables_2/",folder,"tables_",booster,"_",test_idx))
  }
  
  M = cbind(cost0,qaly0,cost1,qaly1,cost2,qaly2,cost3,qaly3,cost4,qaly4)
  #M = cbind(cost0,qaly0,cost3,qaly3,cost4,qaly4)
  
  write.table(M,paste0("ICER_tables_2/",folder,"tables_",booster,"_",test_idx,"/table_",idx_letter[idx],".csv"),row.names = F,col.names = F)
}




fv_qaly <- function(exp){
  rate = 0.015
  weigth_r = 0.87
  xx = seq(0,round(exp))
  sum(unlist(lapply(xx, function(x){ weigth_r*(1/(1+rate)^x)})))
}


total_lifeexp = sum(unlist(lapply(life_exp,fv_qaly))*n_age_group)


save_table_NMB <- function(idx,booster,test_idx,folder){
  
  folder2 = paste0(folder)
  df0 = costs(1,0,0,booster,folder)
  df1 = costs(idx,1,0,booster,folder)
  df2 = costs(idx,2,test_idx,booster,folder)
  df3 = costs(idx,3,test_idx,booster,folder)
  df4 = costs(idx,4,test_idx,booster,folder)
  
  cost0 = boot::boot(df0$cost,fc,R=500)$t[,1]
  qaly0 = boot::boot(total_lifeexp-df0$qaly,fc,R=500)$t[,1]
  
  cost1 = boot::boot(df1$cost,fc,R=500)$t[,1]
  qaly1 = boot::boot(total_lifeexp-df1$qaly,fc,R=500)$t[,1]

  cost2 = boot::boot(df2$cost,fc,R=500)$t[,1]
  qaly2 = boot::boot(total_lifeexp-df2$qaly,fc,R=500)$t[,1]
  
  cost3 = boot::boot(df3$cost,fc,R=500)$t[,1]
  qaly3 = boot::boot(total_lifeexp-df3$qaly,fc,R=500)$t[,1]
  
  cost4 = boot::boot(df4$cost,fc,R=500)$t[,1]
  qaly4 = boot::boot(total_lifeexp-df4$qaly,fc,R=500)$t[,1]
  
  
  if(!file.exists(paste0("NMB_tables"))){
    dir.create(paste0("NMB_tables/"))
  }
  if(!file.exists(paste0("NMB_tables/",folder))){
    dir.create(paste0("NMB_tables/",folder))
  }
  
  
  if(!file.exists(paste0("NMB_tables/",folder,"tables_",booster,"_",test_idx))){
    dir.create(paste0("NMB_tables/",folder,"tables_",booster,"_",test_idx))
  }
  
  M = cbind(cost0,qaly0,cost1,qaly1,cost2,qaly2,cost3,qaly3,cost4,qaly4)
  #M = cbind(cost0,qaly0,cost3,qaly3,cost4,qaly4)
  
  write.table(M,paste0("NMB_tables/",folder,"tables_",booster,"_",test_idx,"/table_",idx_letter[idx],".csv"),row.names = F,col.names = F)
}


save_table_NMB2 <- function(idx,booster,test_idx,folder){
  
  folder2 = paste0("",folder)
  df0 = costs(1,0,0,booster,folder2)
  #df1 = costs(idx,1,0,booster,folder2)
  #df2 = costs(idx,2,test_idx,booster,folder2)
  df3 = costs(idx,3,test_idx,booster,folder2)
  df4 = costs(idx,4,test_idx,booster,folder2)
  df5 = costs(idx,5,test_idx,booster,folder2)
  
  set.seed(12342)
  cost0 = boot::boot(df0$cost,fc,R=500)$t[,1]
  qaly0 = boot::boot(df0$qaly,fc,R=500)$t[,1]
  # set.seed(12342)
  # cost1 = boot::boot(df1$cost,fc,R=500)$t[,1]
  # qaly1 = boot::boot(df1$qaly,fc,R=500)$t[,1]
  # set.seed(12342)
  # cost2 = boot::boot(df2$cost,fc,R=500)$t[,1]
  # qaly2 = boot::boot(df2$qaly,fc,R=500)$t[,1]
  set.seed(12342)
  cost3 = boot::boot(df3$cost,fc,R=500)$t[,1]
  qaly3 = boot::boot(df3$qaly,fc,R=500)$t[,1]
  set.seed(12342)
  cost4 = boot::boot(df4$cost,fc,R=500)$t[,1]
  qaly4 = boot::boot(df4$qaly,fc,R=500)$t[,1]
  set.seed(12342)
  cost5 = boot::boot(df5$cost,fc,R=500)$t[,1]
  qaly5 = boot::boot(df5$qaly,fc,R=500)$t[,1]
  
  
  if(!file.exists(paste0("NMB_tables2"))){
    dir.create(paste0("NMB_tables2/"))
  }
  if(!file.exists(paste0("NMB_tables2/",folder))){
    dir.create(paste0("NMB_tables2/",folder))
  }
  
  
  if(!file.exists(paste0("NMB_tables2/",folder,"tables_",booster,"_",test_idx))){
    dir.create(paste0("NMB_tables2/",folder,"tables_",booster,"_",test_idx))
  }
  
  #M = cbind(cost0,qaly0,cost1,qaly1,cost2,qaly2,cost3,qaly3,cost4,qaly4)
  M = cbind(cost0,qaly0,cost3,qaly3,cost4,qaly4,cost5,qaly5)
  
  write.table(M,paste0("NMB_tables2/",folder,"tables_",booster,"_",test_idx,"/table_",idx_letter[idx],".csv"),row.names = F,col.names = F)
}


