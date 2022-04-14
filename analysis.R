setwd("D:/PosDoc/Coronavirus/Testing/")
library(readxl)
library(dplyr)
library(ggplot2)
library(zoo)
library(lubridate)
library(data.table)
library(tidyverse)

population = 14826276 #Ontario https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=1710000501 on March 24, 2022
enddate = as.Date("2022-03-31")

# Cases -------------------------------------------------------------------
#https://www.canada.ca/en/public-health/services/diseases/2019-novel-coronavirus-infection.html#a1
data = read.csv("covid19-download.csv")
head(data)
data = data[data$prname == "Ontario",]
data$date = as.Date(data$date)
head(data)
data = data[,c(4,6,8)]
head(data)

data$inccases = diff(c(0,data$numconf))
data$incdeaths = diff(c(0,data$numdeaths))


data$roll7cases = frollmean(data$inccases,7)
data$roll7deaths = frollmean(data$incdeaths,7)
  
ggplot()+geom_col(data=data,aes(x=date,y=inccases),fill="grey")+
  geom_line(data=data,aes(x = date, y = roll7cases),color = "navy")+
  #geom_vline(aes(xintercept = d1), linetype = "dashed")+
  theme_bw()


ggplot()+geom_col(data=data,aes(x=date,y=incdeaths),fill="grey")+
  geom_line(data=data,aes(x = date, y = roll7deaths),color = "red")+
  #geom_vline(aes(xintercept = d2), linetype = "dashed")+
  theme_bw()


aa = max(data$roll7cases[data$date <= enddate & data$date >= as.Date("2020-11-01") & !is.na(data$roll7cases)])
aa

d1 = data$date[data$roll7cases == aa & !is.na(data$roll7cases) ]


aa = max(data$roll7deaths[data$date <= as.Date("2021-04-01") & data$date >= as.Date("2020-11-01")  & !is.na(data$roll7deaths)])
aa

d2 = data$date[data$roll7deaths == aa & !is.na(data$roll7deaths) ]
d2

d2-d1

data$dated = data$date-16




ggplot()+#geom_col(data=data,aes(x=date,y=inccases),fill="grey")+
  geom_line(data=data,aes(x = date, y = roll7cases),size=1.5,color = "navy")+
  geom_line(data=data,aes(x = date, y = roll7deaths),size=1.5,color = "red")+
  scale_y_continuous(trans="log")+
  #geom_vline(aes(xintercept = data$dated[length(data$dated)]),linetype = "dashed")+
  #geom_vline(aes(xintercept = d1),linetype = "dashed")+
  theme_bw()




calc_cor <- function(k){
  x1 = x[!is.na(x)]
  y1 = y[!is.na(y)]
  
  x1 = x1[1:(length(x1)-k)]
  y1 = y1[(k+1):length(y1)]
  
  cc = cor(x1,y1,method="spearman")
  
  par(mfrow=c(1,2))
  plot(x1,y1, main = paste("k=",k," cor=", cc))
  plot(rank(x1),rank(y1),main = paste(cc))
  
  return(cc)
}

x = data$roll7cases
y = data$roll7deaths

kk = seq(0,30)

ccc = sapply(kk, calc_cor)
par(mfrow=c(1,1))
plot(ccc)

l = length(data$date)
aaa = data$roll7deaths[16:l]/data$roll7cases[1:(l-15)]
plot(aaa,ylim = c(0,0.04))

ll = length(aaa)
m = mean(aaa[(ll-10):ll])

m*data$roll7cases[nrow(data)]


# Draw first plot using axis y1
par(mar = c(7, 3, 5, 4) + 0.3)              
plot(data$date, data$roll7cases, type = "l", col = "navy")  
abline(v = d1)
# set parameter new=True for a new axis
par(new = TRUE)         

# Draw second plot using axis y2
plot(data$dated, data$roll7deaths, type = "l", col = "red", axes = FALSE)

axis(side = 4, at = pretty(range(data$roll7deaths[!is.na(data$roll7deaths)])))      
mtext("y2", side = 4, line = 3)
abline(v = d1)




##############################################

setwd("~/PosDoc/Coronavirus/Cost_eff/")

data.file = read.csv("data/age-distribution.csv") #https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=1710000501

y = data.file[,6]
y
age = as.numeric(gsub(",", "", y))

data = data.frame(count = age)
rownames(data) = data.file[,1]

g1 = data["0 to 4 years","count"]
g2 = sum(data[c("5 to 9 years","10 to 14 years"),"count"])
g3 = sum(data[c("15 to 19 years","20 to 24 years"),"count"])
g4 = sum(data[c("25 to 29 years","30 to 34 years"),"count"])
g5 = sum(data[c("35 to 39 years","40 to 44 years"),"count"])
g6 = sum(data[c("45 to 49 years","50 to 54 years"),"count"])
g7 = sum(data[c("55 to 59 years","60 to 64 years"),"count"])
g8 = sum(data[c("65 to 69 years","70 to 74 years"),"count"])
g9 = sum(data[c("75 to 79 years","80 to 84 years"),"count"])
g10 = sum(data[c("85 to 89 years","90 to 94 years","95 to 99 years","100 years and over"),"count"])

vector = c(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10)
sum(vector) == data["All ages","count"]

vector/sum(vector)

population = sum(vector)
population

##############################################

setwd("~/PosDoc/Coronavirus/Cost_eff/")

data = read.csv("covid19-download.csv")
data = data %>% filter(prname=="Ontario") %>% mutate(date = as.Date(date)) %>% select(date,numconf)

head(data)
library(dplyr)
library(tidyr)

df = data %>%
  complete(date = seq.Date(min(date), max(date), by="day")) %>%
  #group_by(Product) %>%
  fill(`numconf`)
df$inccases = diff(c(0,df$numconf))
df = df %>% filter(date <= enddate)


ggplot()+geom_col(data=df,aes(x=date,y=inccases),fill="grey")+
  #geom_vline(aes(xintercept = d1), linetype = "dashed")+
  theme_bw()


vector_aux = df$inccases/sum(df$inccases)
va = rev(vector_aux)

write.table(va,"vector_prob_ontario.dat",row.names= F,col.names = F)



# Vaccination Coverage ----------------------------------------------------

library(ggplot2)
library(reshape2)
# https://health-infobase.canada.ca/covid-19/vaccination-coverage/

data = read.csv("vaccination-coverage-byAgeAndSex-overTimeDownload.csv", encoding = "UTF-8")
head(data)
data %>% pull(prfname) %>% unique()
data = data %>% filter(prfname == "Ontario") %>% rename(date = week_end) %>% mutate(date = ymd(date)) %>% filter(date <= enddate)

data %>% pull(sex) %>% unique()
data %>% pull(age) %>% unique()


df = data %>%
  filter(!(age %in% c("Not reported","Unknown", "All ages")),sex == "All sexes") %>%
  mutate(numtotal_fully = as.numeric(numtotal_fully),numtotal_atleast1dose = as.numeric(numtotal_atleast1dose),
         numtotal_additional = as.numeric(numtotal_additional)) %>% 
  select(date,age,numtotal_fully,numtotal_atleast1dose,numtotal_additional)%>%
  rename(groups=age,fully=numtotal_fully,partial=numtotal_atleast1dose,booster=numtotal_additional) %>%
  mutate(booster=replace_na(booster,0))
head(df)


df_f = df %>% select(date,groups,fully) %>%   
  pivot_wider(names_from = c(groups), values_from = c(fully))
df_f[is.na(df_f)]=0
df_f
order = c("date","0–4","05–11","12–17","18–29","30–39","40–49","50–59","60–69","70–79","80+")
df_f = as.data.frame(df_f)
df_f = df_f[,order]
head(df_f)


aa = df_f %>%
  complete(date = seq.Date(min(date), enddate, by="day")) %>%
  #group_by(Product) %>%
  fill(colnames(df_f)[-1])


M2 = sapply(aa[,3:11], diff) ###starting at 2020-12-20
#M2 = abs(M2)
M2[M2<0] = 0

ggplot(melt(abs(M2)), aes(x=Var1, y=value, col=Var2))+
  geom_line()


sum(M2)/population

M22 = round(abs(M2)/population*100000)
M22

write.table(M22,"dose2.dat",row.names = F,col.names = F)

###dose 1


df_f = df %>% select(date,groups,partial) %>% as_tibble() %>%
  pivot_wider(names_from = c(groups), values_from = c(partial))
df_f[is.na(df_f)]=0
df_f

order = c("date","0–4","05–11","12–17","18–29","30–39","40–49","50–59","60–69","70–79","80+")
df_f = as.data.frame(df_f)
df_f = df_f[,order]
head(df_f)


aa = df_f %>%
  complete(date = seq.Date(min(date),enddate, by="day")) %>%
  #group_by(Product) %>%
  fill(colnames(df_f)[-1])


M1 = sapply(aa[,3:11], diff) ###starting at 2020-12-20
#M1 = abs(M1)

ggplot(melt(M1), aes(x=Var1, y=value, col=Var2))+
  geom_line()

sum(abs(M1))/population

M1[M1<0] = 0

M12 = round(abs(M1)/population*100000)
sum(M12)/100000

write.table(M12,"dose1.dat",row.names = F,col.names = F)



###booster


df_f = df %>% select(date,groups,booster)  %>% group_by(date) %>% summarise(booster = sum(booster))
head(df_f)

df_f = df %>% filter(date <= enddate) %>% group_by(date) %>% summarise(booster = sum(booster))


ggplot(df_f)+geom_point(aes(x = date,y=booster/population))+scale_x_date(date_breaks = "1 week")+
  theme(axis.text.x = element_text(angle = 45,hjust=1.0))

head(df_f)
  
  
  
  aa = df_f %>%
  complete(date = seq.Date(min(date), enddate, by="day")) %>%
  #group_by(Product) %>%
  fill(colnames(df_f)[-1])
  
aa$boosterd =diff(c(0,aa$booster)) ###starting at 2020-12-20
B = aa$boosterd  
B[B<0]
  
  
  
  ggplot(aa)+geom_col(aes(x = date,y=boosterd/population))+scale_x_date(date_breaks = "1 week")+
    theme(axis.text.x = element_text(angle = 45,hjust=1.0))
  
  sum(B)/population
  
  B2 = round(B/population*100000)
  sum(B2)/100000
  
  write.table(B2,"booster_dose.dat",row.names = F,col.names = F)
  
  
  ## vaccine type https://health-infobase.canada.ca/covid-19/vaccination-coverage/
  pfizer = 6299358
  moderna = 1760955
  
  v = c(pfizer,moderna)
  
  v/sum(v)

  

# Contact Distribution ----------------------------------------------------
  
  data.file = read.csv("data/age-distribution.csv") #https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=1710000501
  
  y = data.file[,6]
  y
  age = as.numeric(gsub(",", "", y))
  
  data = data.frame(group=data.file[,1],count = age)
  
  
  together_group = c("80 to 84 years","85 to 89 years","90 to 94 years","95 to 99 years","100 years and over")
  g10 = data %>% filter(group %in% together_group) %>% pull(count) %>% sum()
  
  df = data %>% filter(!(group %in% together_group)) %>% rbind(data.frame(group=c("80 and over"),count = c(g10))) %>%
    filter(group != "All ages")
  
  
  
  data.contact = readxl::read_excel("data/CONNECT1_Mixing matrices_20210512 copy.xlsx",sheet = "Total",skip = 1)
  data.contact = data.contact[,-1]
  #[0:4, 5:19, 20:49, 50:64, 65:99]
  group_idx = list(1,c(2,3,4),c(5,6,7,8,9,10),c(11,12,13),c(14,15,16,17))
  
  #Let's create a matrix with the weighted average of the columns based on age groups
  
  fwm <- function(x){
    pp = df$count[x]/sum(df$count[x])
    cc = as.matrix(data.contact[,x]) %*% t(t(pp))
    ct = lapply(group_idx, function(y) sum(cc[y,])) %>% do.call("rbind",.)
    return(ct)
  }
  
  lcontacts = lapply(group_idx, fwm) %>% do.call("cbind",.)
  
  ncontacts = colSums(lcontacts)
  pmatrix = lapply(1:5, function(x) (lcontacts[,x]/ncontacts[x])) %>% do.call("cbind",.)
  t(pmatrix)
  colSums(pmatrix)  
  
  mean0 = c(10.21,
            16.79,
            13.79,
            11.26,
            8.00)
  sd_0 = c(7.65,
           11.72,
           10.50,
           9.59,
           6.96)
  
  
  (sd = ncontacts/mean0*sd_0)
  shelter = c(2.86, 4.7, 3.86, 3.15, 2.24)
  shelter/mean0  
  
  sd_shelter = c(2.14, 3.28, 2.94, 2.66, 1.95)
  
  sd_shelter/sd_0
  
  sd*0.28
  ncontacts*0.28
  
  