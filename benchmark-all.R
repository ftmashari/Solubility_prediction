library("ROCR")
library("caret")
library("UpSetR")
library("ComplexHeatmap")
library("circlize")
library("scales")
library("ggplot2")
library("scales")

Matt_Coef <- function (conf_matrix)
{
  TP <- conf_matrix$table[1,1]
  TN <- conf_matrix$table[2,2]
  FP <- conf_matrix$table[1,2]
  FN <- conf_matrix$table[2,1]
  
  mcc_num <- (TP*TN - FP*FN)
  mcc_den <- 
    as.double((TP+FP))*as.double((TP+FN))*as.double((TN+FP))*as.double((TN+FN))
  
  mcc_final <- mcc_num/sqrt(mcc_den)
  return(mcc_final)
}

colors <- c("#990000", "#FF9999", "#8c510a", "#d8b365", "#01665e", "#5ab4ac", "black", "#7b3294", "#c2a5cf", "#2166ac",
            "#9ecae1")
benchmarktab <- read.table("data/method-comparison.tsv", header = TRUE)
labels <- benchmarktab[,2]
names(labels) <- benchmarktab[,1]
#pdf("figures.pdf")

########################
### Probability methods
########################
# Data preparationsens
probtab <- benchmarktab[,seq(3,13)]
row.names(probtab) <- benchmarktab[,1]
colnames(probtab) <- c("SOLpro","ccSOL","RPSP","PaRSnIP", "PROSOII","PROSO", "Protein-Sol", "DeepSol1", "DeepSol2",
                       "DeepSol3", "SoluProt")
probtab <-probtab[,c(11,10,1,8,4,9,7,6,5,3,2)]
# ROC curve and AUC
pred <- prediction(probtab, matrix(labels, nrow = length(labels), ncol = ncol(probtab)) )
perf <- performance(pred,"tpr","fpr")
auc <- performance(pred,measure = "auc")@y.values
probaucs <- unlist(auc)
names(probaucs)<-colnames(probtab)
auc <- lapply(auc,round,3)
auc <- unlist(auc)*100
names(auc)<-colnames(probtab)
bitmap("results/ROC.tiff", height = 7, width = 7, units = 'in', res=600, family="sans")
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1, 0, 0), cex.axis=2, family='sans',mgp=c(3,1,0))
plot(perf,lty=1,lwd=4,cex.lab=2, main="ROC curves", cex.main=2, ylim=c(0,1),
     col=as.list(colors), xlab='Proportion of soluble proteins', ylab='Proportion of insoluble proteins')
legend("bottomright", inset=-0.01,names(probaucs), col=colors[1:length(probaucs)],pch = 15,cex=1.8, bty="n")
dev.off()
bitmap("results/AUC.tiff", height = 7, width = 7, units = 'in', res=600, family="sans")
par(mar = mar.default + c(0, 7.5, 0, 0), cex.axis=2, family='sans',mgp=c(4,1,0))
barplot(auc,col = colors,las=2, main="AUC comparison",cex.main=2,xlim=c(0,70),xlab="AUC (%)", cex.lab=2,
        horiz = TRUE, border=NA)
box()
abline(v=50,lwd=4)
dev.off()
# Maximise Accuracy
thresholds <- c(0.5,0.5,0.5,0.5,0.5,0.5,0.45,0.5,0.6,50,50)
bitmap("results/threshold.tiff", height = 7, width = 7, units = 'in', res=600, family="sans")
op <- par(mar = c(1, 2, 4, 3), cex.axis=2, family='sans', mfrow=c(4,3), oma = c(5,5,0,0))
bestcutoffs <- c()
pred <- prediction(probtab, matrix(labels, nrow = length(labels), ncol = ncol(probtab)) )
perf <- performance(pred,'sens','spec')
for (i in seq(1,ncol(probtab)))
{
  balac <- (perf@x.values[[i]]+perf@y.values[[i]])/2 # balanced accuracy
  maxacc = max(balac[!is.na(balac) & is.finite(perf@alpha.values[[i]])])
  cutoff = min(perf@alpha.values[[i]][balac==maxacc & !is.na(balac) & is.finite(perf@alpha.values[[i]])])
  bestcutoffs <-c(bestcutoffs,cutoff)
  plot(perf@alpha.values[[i]],(balac*100),type='l',lwd=4, xlab=NA, ylab=NA,cex.lab=2,
       main=colnames(probtab)[i],cex.main=2,las=1,ylim=c(30,70))
  abline(v=cutoff,col=colors[6],lwd=4)
  if(colnames(probtab)[i]=='Protein-Sol')
  {
    abline(v=thresholds[i],col=colors[3],lwd=4,lty=3)
  }
  else
  {
    abline(v=thresholds[i],col=colors[4],lwd=4,lty=2)
  }
}
par(mar=c(0,0,1,0))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend("topleft",legend=c('Best threshold', 'Defined threshold', 'Average score of\nsoluble proteins'),
       col=colors[c(6,4,3)], lwd = 2, lty = c(1,2,3),cex=2, bty="n")
title(xlab = "Threshold",
      ylab = "Balanced accuracy (%)",
      outer = TRUE, line = 3, cex.lab=2.5)
par(op)
dev.off()
names(bestcutoffs) <- colnames(probtab)
#thresholds <- thresholds[-11]
#probtab <- probtab[,-11]

########################
### Accuracy for cut-offs
########################
classification <- probtab[,-7]
newthresholds <- thresholds[-7]
newcolors <- colors[-7]
for (i in seq(1,ncol(classification)))
{
  classification[,i][classification[,i] >= newthresholds[i]] <- 'yes'
  classification[,i][classification[,i] != 'yes'] <- 'no'
}
accuracies <- c()
for (item in colnames(classification))
{
  cm <- confusionMatrix(as.factor(classification[,item]), labels, positive = "yes")
  cm$byClass <- c(cm$byClass,MCC=Matt_Coef(cm),cm$overall['Accuracy'])
  accuracies <- c(accuracies,list(cm$byClass))
}
names(accuracies) <- colnames(classification)

########################
### All measures
########################
# Calculate sensitivity
sens <- c()
for (item in names(accuracies))
{
  sens <-c(sens,accuracies[[item]][1])
}
names(sens) <- names(accuracies)
sens <- sens*100
# Calculate specificity
spec <- c()
for (item in names(accuracies))
{
  spec <-c(spec,accuracies[[item]][2])
}
names(spec) <- names(accuracies)
spec <- spec*100
# Calculate accuracies
accs <-c()
for (item in names(accuracies))
{
  accs <-c(accs,accuracies[[item]][11])
}
names(accs) <- names(accuracies)
accs <- accs*100
# Plot
shapes <- c(15,15,19,19,20,20,17,17,18,18)
bitmap("results/sensitivity-vs-specificity.tiff", height = 7, width = 7, units = 'in', res=600, family="sans")
par(mar = mar.default + c(0, 1, 0, 0), cex.axis=2, family='sans',mfrow=c(1,1),mgp=c(3,1,0), xaxs="i", yaxs="i")
plot(sens,spec,col=newcolors,cex=4,xlab ="Sensitivity (%)", ylab="Specificity (%)",
     cex.lab=2,xlim=c(0,100),ylim=c(0,100), main="Sensitivity vs. specificity", cex.main=2,
     cex.axis=2, pch=shapes)
abline(100,-1,lwd=4)
legend("bottomleft", inset=-0.01,names(accs)[order(accs, decreasing = TRUE)], col=newcolors[order(accs, decreasing = TRUE)],
       pch =shapes[order(accs, decreasing = TRUE)],cex=2, bty="n")
dev.off()

bitmap("results/accuracy.tiff", height = 7, width = 7, units = 'in', res=600, family="sans")
par(mar = mar.default + c(0, 6.5, 0, 0), cex.axis=2, family='sans',mgp=c(4,1,0),  xaxs="r", yaxs="r")
barplot(sort(accs,decreasing = TRUE),col = newcolors[order(accs, decreasing = TRUE)],las=2, main="Balanced accuracy comparison",cex.main=2,
        xlim=c(min(0,min(accs)),65),xlab="Balanced accuracy (%)", cex.lab=2,
        horiz = TRUE, border=NA)
box()
abline(v=50,lwd=4)
dev.off()

########################
### Simulation
########################
numsim = 10000
numdat = nrow(benchmarktab)
aucs = list(rep(0, numsim))
bacs = list(rep(0, numsim))
set.seed(123)
probs <- matrix(runif(numdat*numsim), numdat, numsim)
classification <- probs
for (i in seq(numsim))
{
  pred <- prediction(probs[,i], labels)
  perf <- performance(pred,"tpr","fpr")
  aucs[i] <- performance(pred,measure = "auc")@y.values
  
  classification[,i][classification[,i] >= 0.5] <- 'yes'
  classification[,i][classification[,i] != 'yes'] <- 'no'
  cm <- confusionMatrix(as.factor(classification[,i]), labels, positive = "yes")
  bacs[i] = cm$byClass['Balanced Accuracy']
}
aucs<-unlist(aucs)
bacs <- unlist(bacs)
# plot balanced accuracies
bacs <- bacs*100
dens <- density(bacs)
multiplier <- length(bacs)/sum(dens$y)
dens$y <- multiplier * dens$y
q95 <- quantile(bacs, .95)
q100 <- quantile(bacs, 1)
xmid <- mean(bacs)
x1 <- min(which(dens$x > q95))  
x2 <- max(which(dens$x <=  q100))

bitmap("results/accuracy-vs-random.tiff", height = 7, width = 7, units = 'in', res=600, family="sans")
par(mar = mar.default + c(0, 1, 0, 0), cex.axis=2, family='sans',mfrow=c(1,1),mgp=c(3,1,0))
plot(dens,lwd=4, xlab="Balanced accuracy (%)", ylab="Frequency", cex.lab=2, main = "Balanced accuracy of random predictions",
     cex.main=2)
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="gray"))
points(accs,c(5,8,5,8,5,5,5,5,5,8),col=newcolors,pch=shapes,cex=3)
abline(v=xmid,col="gray",lwd=4,lty=2)
legend("topleft", inset=-0.01,names(accs)[order(accs, decreasing = TRUE)], col=newcolors[order(accs, decreasing = TRUE)],
       pch = shapes[order(accs, decreasing = TRUE)],cex=2, bty="n")
legend("topright", c("Mean","Pval < 0.05"), col="gray", lwd=4,lty=c(2,NA),pch=c(NA,15),cex=2, bty = "n")
dev.off()
acc95 <- q95
# plot AUCs
aucs <- aucs*100
dens <- density(aucs)
multiplier <- length(aucs)/sum(dens$y)
dens$y <- multiplier * dens$y
q95 <- quantile(aucs, .95)
q100 <- quantile(aucs, 1)
xmid <- mean(aucs)
x1 <- min(which(dens$x > q95))  
x2 <- max(which(dens$x <=  q100))
bitmap("results/AUC-vs-random.tiff", height = 7, width = 7, units = 'in', res=600, family="sans")
par(mar = mar.default + c(0, 1, 0, 0), cex.axis=2, family='sans',mfrow=c(1,1),mgp=c(3,1,0))
plot(dens,lwd=4, xlab="AUC (%)", ylab="Frequency", cex.lab=2, main = "AUC of random predictions",
     cex.main=2)
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="gray"))
shapes_all <- c(15,15,19,19,20,20,10,17,17,18,18)
points(probaucs*100,c(5,8,5,8,5,5,5,5,8,11,5),col=colors,pch=shapes_all,cex=3)
abline(v=xmid,col="gray",lwd=4,lty=2)
legend("topleft", inset=-0.01,names(probaucs), col=colors[1:length(probaucs)],pch = shapes_all,cex=2, bty="n")
legend("topright", c("Mean","Pval < 0.05"), col="gray", lwd=4,lty=c(2,NA),pch=c(NA,15),cex=2, bty = "n")
dev.off()
auc95 <- q95

########################
### Majority voting
########################
maxauc <- max(probaucs)*100
maxacc <- max(accs)
# combmeasures <- list()
# 
# pdf('combinations.pdf')
# # SoluProt+DeepSol3+SOLpro
# bests <- probtab[,c(1,2,3)] # best predictors
# majority <- rowSums(bests)/3 # average the probabilities from the best predictors
# majorityyn <- majority >= 0.5
# majorityyn[majorityyn] <- "yes"
# majorityyn[majorityyn!="yes"] <- "no"

# pred <- prediction(majority, labels )
# perf <- performance(pred,"tpr","fpr")
# auc <- performance(pred,measure = "auc")@y.values
# par(mar = mar.default + c(0, 1, 0, 0), cex.axis=2, family='sans',mfrow=c(1,1),mgp=c(3,1,0))
# plot(perf,lty=1,lwd=4,cex.lab=2, main="Solpar3d", cex.main=2, ylim=c(0,1),
#      col='gray')
# cm <- confusionMatrix(as.factor(majorityyn), labels, positive = "yes")
# measures <- c(cm$byClass,AUC=auc[[1]][1])
# measures <- measures[c(1,2,11,12)]*100
# names(measures)[3]="Balanced accuracy"
# combmeasures[['SoluProt, DeepSol3,\nSOLpro']] <- measures
# par(mar = mar.default + c(11, 1, 0, 0), cex.axis=2, family='sans',mfrow=c(1,1),mgp=c(3,1,0))
# df.bar <- barplot(measures,col = 'gray', main="DeepSol3, SOLpro, PaRSnIP (Solpar3d)",las=2,border=NA, ylim=c(0,100),cex.main=2)
# box()
# text(df.bar[1],measures[1]+6, lapply(measures[1],round,0),cex=2)
# text(df.bar[2],measures[2]+6, lapply(measures[2],round,0),cex=2)
# text(df.bar[3],measures[3]+6, lapply(measures[3],round,0),cex=2)
# text(df.bar[4],measures[4]+6, lapply(measures[4],round,0),cex=2)
# segments(df.bar[3]-0.49,maxacc,df.bar[3]+0.49,maxacc,col=colors[3],lwd=4)
# segments(df.bar[3]-0.49,acc95,df.bar[3]+0.49,acc95,col=colors[5],lwd=4)
# segments(df.bar[4]-0.49,maxauc,df.bar[4]+0.49,maxauc,col=colors[3],lwd=4)
# segments(df.bar[4]-0.49,auc95,df.bar[4]+0.49,auc95,col=colors[5],lwd=4)
# legend("bottomright", inset=-0.01,c('Best predictor', 'Pval = 0.05'), col=colors[c(3,5)],lwd=4,cex=2, bty="n")

# # SoluProt+DeepSol3
# df <- probtab[,c(1,2)]
# majority <- rowSums(df)/2
# majorityyn <- majority >= 0.5
# majorityyn[majorityyn] <- "yes"
# majorityyn[majorityyn!="yes"] <- "no"
# pred <- prediction(majority, labels )
# perf <- performance(pred,"tpr","fpr")
# auc <- performance(pred,measure = "auc")@y.values
# cm <- confusionMatrix(as.factor(majorityyn), labels, positive = "yes")
# measures <- c(cm$byClass,AUC=auc[[1]][1])
# measures <- measures[c(1,2,11,12)]*100
# names(measures)[3]="Balanced accuracy"
# combmeasures[['SoluProt, DeepSol3']] <- measures
# par(mar = mar.default + c(11, 1, 0, 0), cex.axis=2, family='sans',mfrow=c(1,1),mgp=c(3,1,0))
# df.bar <- barplot(measures,col = 'gray', main="DeepSol3, SOLpro",las=2,border=NA, ylim=c(0,100),cex.main=2)
# box()
# text(df.bar[1],measures[1]+6, lapply(measures[1],round,0),cex=2)
# text(df.bar[2],measures[2]+6, lapply(measures[2],round,0),cex=2)
# text(df.bar[3],max(measures[3],acc95)+6, lapply(measures[3],round,0),cex=2)
# text(df.bar[4],max(measures[4],maxauc)+6, lapply(measures[4],round,0),cex=2)
# segments(df.bar[3]-0.49,maxacc,df.bar[3]+0.49,maxacc,col=colors[3],lwd=4)
# segments(df.bar[3]-0.49,acc95,df.bar[3]+0.49,acc95,col=colors[5],lwd=4)
# segments(df.bar[4]-0.49,maxauc,df.bar[4]+0.49,maxauc,col=colors[3],lwd=4)
# segments(df.bar[4]-0.49,auc95,df.bar[4]+0.49,auc95,col=colors[5],lwd=4)
# legend("bottomright", inset=-0.01,c('Best predictor', 'Pval = 0.05'), col=colors[c(3,5)],lwd=4,cex=2, bty="n")

# # SoluProt+SOLpro
# df <- probtab[,c(1,3)]
# majority <- rowSums(df)/2
# majorityyn <- majority >= 0.5
# majorityyn[majorityyn] <- "yes"
# majorityyn[majorityyn!="yes"] <- "no"
# pred <- prediction(majority, labels )
# perf <- performance(pred,"tpr","fpr")
# auc <- performance(pred,measure = "auc")@y.values
# cm <- confusionMatrix(as.factor(majorityyn), labels, positive = "yes")
# measures <- c(cm$byClass,AUC=auc[[1]][1])
# measures <- measures[c(1,2,11,12)]*100
# names(measures)[3]="Balanced accuracy"
# combmeasures[['SoluProt, SOLpro']] <- measures
# par(mar = mar.default + c(11, 1, 0, 0), cex.axis=2, family='sans',mfrow=c(1,1),mgp=c(3,1,0))
# df.bar <- barplot(measures,col = 'gray', main="DeepSol3, PaRSnIP",las=2,border=NA, ylim=c(0,100),cex.main=2)
# box()
# text(df.bar[1],measures[1]+6, lapply(measures[1],round,0),cex=2)
# text(df.bar[2],measures[2]+6, lapply(measures[2],round,0),cex=2)
# text(df.bar[3],max(measures[3],acc95)+6, lapply(measures[3],round,0),cex=2)
# text(df.bar[4],max(measures[4],maxauc)+6, lapply(measures[4],round,0),cex=2)
# segments(df.bar[3]-0.49,maxacc,df.bar[3]+0.49,maxacc,col=colors[3],lwd=4)
# segments(df.bar[3]-0.49,acc95,df.bar[3]+0.49,acc95,col=colors[5],lwd=4)
# segments(df.bar[4]-0.49,maxauc,df.bar[4]+0.49,maxauc,col=colors[3],lwd=4)
# segments(df.bar[4]-0.49,auc95,df.bar[4]+0.49,auc95,col=colors[5],lwd=4)
# legend("bottomright", inset=-0.01,c('Best predictor', 'Pval = 0.05'), col=colors[c(3,5)],lwd=4,cex=2, bty="n")

# # DeepSol3+PaRSnIP+PROSOII
# df <- probtab[,c(1,4,8)]
# deepsnipproso <- rowSums(df)/3
# deepsnipprosoyn <- deepsnipproso >= (0.5+0.5+0.6)/3
# deepsnipprosoyn[deepsnipprosoyn] <- "yes"
# deepsnipprosoyn[deepsnipprosoyn!="yes"] <- "no"
# pred <- prediction(deepsnipproso, labels )
# perf <- performance(pred,"tpr","fpr")
# auc <- performance(pred,measure = "auc")@y.values
# cm <- confusionMatrix(as.factor(deepsnipprosoyn), labels, positive = "yes")
# measures <- c(cm$byClass,AUC=auc[[1]][1])
# measures <- measures[c(1,2,11,12)]*100
# names(measures)[3]="Balanced accuracy"
# par(mar = mar.default + c(11, 1, 0, 0), cex.axis=2, family='sans',mfrow=c(1,1),mgp=c(3,1,0))
# df.bar <- barplot(measures,col = 'gray', main="DeepSol3, PaRSnIP, PROSOII",las=2,border=NA, ylim=c(0,100),cex.main=2)
# box()
# text(df.bar[1],measures[1]+6, lapply(measures[1],round,0),cex=2)
# text(df.bar[2],measures[2]+6, lapply(measures[2],round,0),cex=2)
# text(df.bar[3],max(measures[3],acc95)+6, lapply(measures[3],round,0),cex=2)
# text(df.bar[4],max(measures[4],maxauc)+6, lapply(measures[4],round,0),cex=2)
# segments(df.bar[3]-0.49,maxacc,df.bar[3]+0.49,maxacc,col=colors[3],lwd=4)
# segments(df.bar[3]-0.49,acc95,df.bar[3]+0.49,acc95,col=colors[5],lwd=4)
# segments(df.bar[4]-0.49,maxauc,df.bar[4]+0.49,maxauc,col=colors[3],lwd=4)
# segments(df.bar[4]-0.49,auc95,df.bar[4]+0.49,auc95,col=colors[5],lwd=4)
# legend("bottomright", inset=-0.01,c('Best predictor', 'Pval = 0.05'), col=colors[c(3,5)],lwd=4,cex=2, bty="n")
# 
# # DeepSol3+PaRSnIP+PROSOII+SOLpro
# df <- probtab[,c(1,4,8,2)]
# deepsnipprososolpro <- rowSums(df)/4
# deepsnipprososolproyn <- deepsnipprososolpro >= (0.5+0.5+0.5+0.6)/4
# deepsnipprososolproyn[deepsnipprososolproyn] <- "yes"
# deepsnipprososolproyn[deepsnipprososolproyn!="yes"] <- "no"
# pred <- prediction(deepsnipprososolpro, labels )
# perf <- performance(pred,"tpr","fpr")
# auc <- performance(pred,measure = "auc")@y.values
# cm <- confusionMatrix(as.factor(deepsnipprososolproyn), labels, positive = "yes")
# measures <- c(cm$byClass,AUC=auc[[1]][1])
# measures <- measures[c(1,2,11,12)]*100
# names(measures)[3]="Balanced accuracy"
# combmeasures[['SOLpro, DeepSol3,\nPaRSnIP, PROSOII']] <- measures
# par(mar = mar.default + c(11, 1, 0, 0), cex.axis=2, family='sans',mfrow=c(1,1),mgp=c(3,1,0))
# df.bar <- barplot(measures,col = 'gray', main="DeepSol3, SOLpro, PaRSnIP, PROSOII",las=2,border=NA, ylim=c(0,100),cex.main=2)
# box()
# text(df.bar[1],measures[1]+6, lapply(measures[1],round,0),cex=2)
# text(df.bar[2],measures[2]+6, lapply(measures[2],round,0),cex=2)
# text(df.bar[3],max(measures[3],acc95)+6, lapply(measures[3],round,0),cex=2)
# text(df.bar[4],max(measures[4],maxauc)+6, lapply(measures[4],round,0),cex=2)
# segments(df.bar[3]-0.49,maxacc,df.bar[3]+0.49,maxacc,col=colors[3],lwd=4)
# segments(df.bar[3]-0.49,acc95,df.bar[3]+0.49,acc95,col=colors[5],lwd=4)
# segments(df.bar[4]-0.49,maxauc,df.bar[4]+0.49,maxauc,col=colors[3],lwd=4)
# segments(df.bar[4]-0.49,auc95,df.bar[4]+0.49,auc95,col=colors[5],lwd=4)
# legend("bottomright", inset=-0.01,c('Best predictor', 'Pval = 0.05'), col=colors[c(3,5)],lwd=4,cex=2, bty="n")

# for (item in names(combmeasures))
# {
#   names(combmeasures[[item]]) <- c("Sensitivity","Specificity","Balanced\naccuracy","AUC")
# }
# measures2plot <- t(as.data.frame(matrix(unlist(strsplit(names(unlist(combmeasures)),split="[.]")),nrow=2)))
# rownames(measures2plot) <- c()
# measures2plot <- as.data.frame(cbind(measures2plot,as.numeric(unlist(combmeasures))))
# colnames(measures2plot) <- c('Methods','Measure', 'Value(%)')
# measures2plot$`Value(%)` <- as.numeric(as.matrix(measures2plot)[,3])
# ggplot(measures2plot, aes(Measure, `Value(%)`, fill = Methods)) +
#   geom_bar(stat="identity", position = "dodge2") + scale_y_continuous(limits=c(4,85),oob = rescale_none) +
#   scale_fill_manual(values=colors[c(4,6,9)]) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         panel.border = element_rect(colour = "black", fill=NA), axis.title.x = element_text(colour = "black"),
#         axis.title.y = element_text(colour = "black"), text = element_text(size=20),
#         axis.text.x = element_text(angle=90, hjust=1), legend.key = element_rect(colour = "transparent", fill = "white"),
#         legend.position="bottom", legend.direction = "vertical", legend.key.size = unit(2.7,"line")) +
#   geom_segment(aes(x = 0.5, y = maxauc, xend = 1.5, yend = maxauc,color=colors[1]),lwd=1,linetype=2) +
#   geom_segment(aes(x = 0.5, y = auc95, xend = 1.5, yend = auc95,color=colors[10]),lwd=1,linetype=2) +
#   geom_segment(aes(x = 1.5, y = maxacc, xend = 2.5, yend = maxacc),color=colors[10],lwd=1,linetype=2) +
#   geom_segment(aes(x = 1.5, y = acc95, xend = 2.5, yend = acc95),color=colors[1],lwd=1,linetype=2) +
#   scale_color_manual(values = colors[c(1,10)], labels = c("Pval = 0.05", "SoluProt"),name='Lines') +
#   coord_flip() + guides(fill = guide_legend(reverse = TRUE))
# dev.off()

########################
### Study top n (upset)
########################
#probtab$Solpar3d <- majority
#probtab <- probtab[,c(11,seq(1,10))]
#pdf("upsetr.pdf")
all <- probtab
for (item in colnames(all))
{
  all[,item] <- rownames(all)[order(all[,item], decreasing=TRUE)]
}
rownames(all) <- c()
# n=5
# upset(fromList(all[1:n,]),nsets = ncol(all),set_size.show=TRUE,nintersects = NA, text.scale = 2, point.size = 4,
#       mb.ratio = c(0.3, 0.7))
# dev.off()

########################
### Heatmap all
########################
bitmap("results/heatmap.tiff", height = 12, width = 11, units = 'in', res=600, family="sans")
#pdf("heatmap.pdf",height=15,width=15)
all.labs <- data.frame(matrix(ncol = ncol(all), nrow = nrow(all)))
colnames(all.labs) <- colnames(all)
for (i in seq(1,nrow(all)))
{
  for (j in seq(1,ncol(all)))
  {
    all.labs[i,j] <- as.character(labels[all[i,j]])
  }
}
colormat <- all
for (i in seq(1,nrow(colormat)))
{
  for (j in seq(1,ncol(colormat)))
  {
    if (probtab[colormat[i,j],j] < thresholds[j])
    {
      colormat[i,j] <- "white"
    }
    else
    {
      colormat[i,j] <- "black"
    }
  }
}
ha = rowAnnotation(foo = anno_text(c("Highest", "predicted", "solubility",rep(' ',nrow(all.labs)-6),"Lowest", "predicted",
                                     "solubility"),gp=gpar(fontsize=12)))
ht= Heatmap(as.matrix(all.labs),col=colors[c(3,6)], name = 'Actual_solubility', right_annotation = ha, column_names_side = "top",
            show_heatmap_legend = FALSE,
            cell_fun = function(j, i, x, y, width, height, fill)
        {
          grid.text(all[i,j], x = x, y = y, gp=gpar(col=colormat[i,j], fontsize=9))
        })
lgd_1 = Legend(at = c("yes","no"), title = "Actual_solubility", legend_gp = gpar(fill = colors[c(6,3)],fontsize=12))
lgd_2 = Legend(at = c("yes","no"), title = "Predicted_solubility", legend_gp = gpar(fill = c("black","white"),fontsize=12),border="black")
draw(ht,annotation_legend_list=c(lgd_1,lgd_2))
dev.off()
