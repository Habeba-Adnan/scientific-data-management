

install.packages("ggplot2")
library(ggplot2)
library(reshape2)
install.packages('corrplot')
library(corrplot)
library(dplyr)
library(tidyr)
##################################
df=read.csv("METABRIC-RNA-Mutation.csv")
print(ncol(df))
print(nrow(df))
data <- cor(df[sapply(df,is.numeric)])
data1 <- melt(data)
ggplot(data1, aes(x = Var1,y = Var2,fill = value))+geom_tile()
ggplot(data1,aes(x = Var1, y = Var2,fill = value))+geom_tile() + scale_fill_distiller(palette = "Spectral")
 
##################################################################
genes=df[,40:45]
n=cor(genes) #compute the correlation coefficient between variables
# The correlation plot to Check the correlation between variables
corrplot(n, method="shade", type = "upper",tl.cex=.6,tl.col="black", title="Correlation Plot",number.font = 2, mar=c(0,0,1,0),	)
###############################################
clinical = df[, c(1:20)] 
clinical_n <- select_if(clinical, is.numeric) 
hist(df$tumor_size , main = "", xlab = "Tumor size", freq = FALSE )
curve(dnorm(x, mean = mean(df$tumor_size), sd = sd(df$tumor_size)),col = "red", add = TRUE)

# The logarithm transformation
tumor_size_T <- log10((df$tumor_size));
hist(tumor_size_T , main = "logarithm transformation", xlab = "Tumor size ", freq = FALSE )
curve(dnorm(x, mean = mean(tumor_size_T), sd = sd(tumor_size_T)),col = "red", add = TRUE)

####################################################

#density of mutation count
hist(clinical$mutation_count , main = "mutation_count", xlab = "mutation_count", freq = FALSE )
curve(dnorm(x, mean = mean(df$mutation_count),sd = sd(df$mutation_count)),col = "green", add = TRUE)

#The logarithm transformation of mutation count

mutation_count_T <- log10((clinical$mutation_count));
hist(mutation_count_T , main = "logarithm mutation_count", xlab = "mutation_count", freq = FALSE )
curve(dnorm(x, mean = mean(mutation_count_T),sd = sd(mutation_count_T)),col = "blue", add = TRUE)

###################################################################3
#boxplot of tumor size
boxplot(df$tumor_size, xlab = "Tumor size ")
###################################
age=df$age_at_diagnosis
tumor_size=df$tumor_size
tumor_stage=df$tumor_stage
unique(tumor_stage)
stages <- ifelse(tumor_stage == 1, "stage 1", ifelse(tumor_stage ==2 , "stage 2"
                                                     ,ifelse(tumor_stage ==3 , "stage 3","stage 4")))
plot(age,tumor_size,pch = 19,col = factor(stages))
legend("topleft",legend = levels(factor(stages)),pch = 20,
       col = factor(levels(factor(stages))))

###################################################################

correlation_data=data[40:100,40:100]
heatmap(correlation_data,main="Correlation Heatmap")

#######################################################3
genes=as.matrix(df[,40:100])
#identify colors
my_colors <- colorRampPalette(c("cyan", "black")) 
heatmap(genes,main="genes Heatmap",col = my_colors(100))

#####################################
treated=subset(df,df["chemotherapy"]==1)
treated_matrix=as.matrix(treated[,40:100])
my_colors <- colorRampPalette(c("cyan", "deeppink")) 
heatmap(treated_matrix,col = my_colors(100),main = "treated samples with chemotherapy")

control=subset(df,df["chemotherapy"]==0)
control_matrix=as.matrix(treated[,40:100])
heatmap(control_matrix,main="control samples with chemotherapy")
##########################################
t1=treated[sapply(treated,is.numeric)]
res=t.test(t1)
res
res$p.value
#########################################################
 #classify our treated_data according to the stages of tumor
stage1=subset(treated,treated["tumor_stage"]==1)
s1=stage1[sapply(stage1,is.numeric)]
stage2=subset(treated,treated["tumor_stage"]==2)
s2=stage2[sapply(stage2,is.numeric)]
stage3=subset(treated,treated["tumor_stage"]==3)
s3=stage3[sapply(stage3,is.numeric)]
stage4=subset(treated,treated["tumor_stage"]==4)
s4=stage4[sapply(stage4,is.numeric)]
############################333
density1=density(as.matrix(s1))
plot(density1,col="red" )
density2=density(as.matrix(s2))
lines(density2,col="green")
density3=density(as.matrix(s3))
lines(density3,col="blue")
density4=density(as.matrix(s4))
lines(density4,col="pink")
##########################################
res_stage1=t.test(s1)
res_stage1
res_stage1$p.value

res_stage2=t.test(s2)
res_stage2
res_stage2$p.value

res_stage3=t.test(s3)
res_stage3
res_stage3$p.value

res_stage4=t.test(s4)
res_stage4
res_stage4$p.value
#####################################################
#qq-plot
age_of_survival=subset(df,df["overall_survival"]==1)

qqnorm(age_of_survival$age_at_diagnosis, pch = 1, frame = FALSE,main="overall survival =1")
qqline(age_of_survival$age_at_diagnosis, col = "red", lwd = 2)

####
age_of_nonsurvival=subset(df,df["overall_survival"]==0)
qqnorm(age_of_nonsurvival$age_at_diagnosis, pch = 1, frame = FALSE,main="overall survival =0")
qqline(age_of_nonsurvival$age_at_diagnosis, col = "steelblue", lwd = 2)
########################
#boxplots
boxplot(age_of_nonsurvival$age_at_diagnosis,age_of_survival$age_at_diagnosis,notch=TRUE,main="boxplot"
        ,names=c("age of non_survival","age of survival"), col=c("red","yellow"))
##########################################

plot(age_of_nonsurvival$age_at_diagnosis,type = "l",col = "red", xlab = "age", ylab = "density", 
     main = "relation between age of survival and non-survival")

lines(age_of_survival$age_at_diagnosis, type = "l", col = "blue")

