library(tidyverse)
library(car)
library(lme4)
theme_set(theme_classic())
### create mock dataset with t technical replicates

### parameters
t=8        ### techical reps (i.e. number of leaves measured)
l=t*1024   ### total number of leaves measured during experiment

### Mean treatment effect (units are same as response variable)
T1_effect=0
T2_effect=0.5
D1_effect=-1
IR_effect=0
C_effect=0
B_effect=0
Low_effect=0
High_effect=2
Sep_1_effect=-0.25
Sep_7_effect=0.25

### create treatment structure and a random response variable
dataset=data.frame(Treatment=rep(c(rep('T1',64),rep('T2',64)),l/128),
               Replicate=rep(c(rep('R1',16),rep('R2',16),rep('R3',16),rep('R4',16)),l/64),
               Footprint=rep(c(rep('D1',8),rep('IR',8)),l/16),
               Subplot=rep(c(rep('C',4),rep('B',4)),l/8),
               CO2level=rep(c(rep('Low',2),rep('High',2)),l/4),
               Day=rep(c(rep('Sep_1',1),rep('Sep_7',1)),l/2),
               A=rnorm(l,5,3))

### add treatment effects to response variable
dataset$A[dataset$Treatment=='T1']=T1_effect+dataset$A[dataset$Treatment=='T1']
dataset$A[dataset$Treatment=='T2']=T2_effect+dataset$A[dataset$Treatment=='T2']
dataset$A[dataset$Footprint=='D1']=D1_effect+dataset$A[dataset$Footprint=='D1']
dataset$A[dataset$Footprint=='IR']=IR_effect+dataset$A[dataset$Footprint=='IR']
dataset$A[dataset$Subplot=='C']=C_effect+dataset$A[dataset$Subplot=='C']
dataset$A[dataset$Subplot=='B']=B_effect+dataset$A[dataset$Subplot=='B']
dataset$A[dataset$CO2level=='Low']=Low_effect+dataset$A[dataset$CO2level=='Low']
dataset$A[dataset$CO2level=='High']=High_effect+dataset$A[dataset$CO2level=='High']
dataset$A[dataset$Day=='Sep_1']=Sep_1_effect+dataset$A[dataset$Day=='Sep_1']
dataset$A[dataset$Day=='Sep_7']=Sep_7_effect+dataset$A[dataset$Day=='Sep_7']

head(dataset)
str(dataset)

### model construction
datlm <- lmer((A) ~ Treatment+Footprint+Subplot+CO2level+Day+
                (1:Replicate)+                             # block
                (1|Treatment:Replicate)+                   # whole plot
                (1|Footprint:Treatment:Replicate)+         # footprint
                (1|Subplot:Footprint:Treatment:Replicate)+ # subplot 
                (1|CO2level:Subplot:Footprint:Treatment:Replicate), #CO2 level
              data = dataset%>%
                group_by(Treatment,Replicate,Footprint,Subplot,CO2level,Day)%>%
                summarise(A=mean(A,na.rm=))) 

### model evaluation
joint_tests(datlm)
hist(resid(datlm),20,col='Blue',xlab="Model Residuals",main="")
shapiro.test((resid(datlm)))
boxplot(resid(datlm)~Treatment+Footprint+Subplot+CO2level+Day,
        data = dataset%>%
          group_by(Treatment,Replicate,Footprint,Subplot,CO2level,Day)%>%
          summarise(A=mean(A,na.rm=)),
        col='Blue',xlab="Treatment",ylab="Model Residuals")
plot(datlm,pch=20,cex=5,col='Blue',xlab="Fitted Values",ylab="Model Residuals")
qqnorm(resid(datlm))
qqline(resid(datlm))
leveneTest((A)~Treatment*Footprint*Subplot*CO2level*Day,
           data = dataset%>%
             group_by(Treatment,Replicate,Footprint,Subplot,CO2level,Day)%>%
             summarise(A=mean(A,na.rm=)))

### Multple Comparisons
datlm.lsd=emmeans(datlm, specs=c( 'Treatment', 'Footprint','Subplot','Day'),
                  by=c('CO2level'),
                  adjust='Tukey',cov.reduce = F)
datlm.lsd_output=multcomp::cld(datlm.lsd,alpha=0.05,Letters=letters,
                               adjust='none',type = "response")
datlm.lsd_output$.group=trimws(datlm.lsd_output$.group)
datlm.lsd_output

### Graphing
ggplot(datlm.lsd_output,
       aes(x=Treatment:Footprint:Subplot,y=response,fill = Day,
           label = .group))+
  facet_grid(vars(CO2level))+
  geom_bar(stat="identity", position=position_dodge(), show.legend = TRUE)+
  geom_errorbar(aes(ymin  =  response-SE, ymax  =  response+SE),
                width =  0.3, size  =  1, position = position_dodge(0.9))+
  geom_text(position = position_dodge(0.9), 
            aes(y=response+SE),vjust=-0.8, size = 5)+
  ylab('A mol co2 / m^2 / s')+
  coord_cartesian(ylim=c(0,1.15*max(datlm.lsd_output$response+datlm.lsd_output$SE)))+
  ggtitle('Leaf Level Instantaneous CO2 Assimilation')+
  xlab('')+
  theme(
    legend.position=c(0.1, 0.1),
    legend.title = element_blank(),
    legend.box = 1,
    plot.title = element_text(face= "bold", size=25,hjust = 0.5),
    axis.title = element_text(face= "bold", size=25,color = 'Black'),
    text=element_text(face= "bold",size=25,color = 'Black'),
    plot.margin = margin(1,1,0,1,"cm"),
    axis.text = element_text(face= "bold",size=25,color = 'Black')
  )



