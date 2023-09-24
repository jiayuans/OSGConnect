#!/usr/bin/env Rscript

#library(dplyr)
#df <- data.frame(
#  id = c(10,11,12,13,14,15,16,17),
#  name = c('sai','ram','deepika','sahithi','kumar','scott','Don','Lin'),
#  gender = c('M','M','F','F','M','M','M','F'),
#  dob = as.Date(c('1990-10-02','1981-3-24','1987-6-14','1985-8-16',
#                  '1995-03-02','1991-6-21','1986-3-24','1990-8-26')),
#  state = c('CA','NY',NA,NA,'DC','DW','AZ','PH'),
#  row.names=c('r1','r2','r3','r4','r5','r6','r7','r8')
#)
#df
#df %>% filter(rownames(df) == 'r3')
#df %>% select(c(1,2))


library(tidyverse)
library(nlme)

cortdata <- read.csv(file="cortdata.csv")
cortdatalong <-cortdata %>% pivot_longer(cols = 'Time1':'Time6', names_to = "Time", values_to = "Cortisollevel") %>%
  mutate(casecontrol = factor(casecontrol, labels = c("Control","Patient"))) ## Reshape the data to a long format using tidyverse

head(cortdatalong)
cortdatalong$Time <- as.factor(cortdatalong$Time)
cortdatalong <- within(cortdatalong, Time <- relevel(Time, ref = "Time6")) ## change the reference level
cortdatalong$casecontrol <- as.factor(cortdatalong$casecontrol)

lmm.un.fit<-gls(Cortisollevel~Time + casecontrol, data=cortdatalong, correlation = corSymm(form = ~ 1 | SubjectId))
lmm.un.summary <- summary(lmm.un.fit)
lmm.un.summary