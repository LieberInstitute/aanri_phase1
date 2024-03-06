dat1 <- read.csv("../out/dmr.csv")
dat2 <- read.table("./out/dmr.csv",header=T)
dat1$id <- paste0(dat1[,2],"_",dat1$start,"_",dat1$end)
dat2$id <- paste0(dat2[,1],"_",dat2$start,"_",dat2$end)
id <- intersect(dat1$id,dat2$id)
dat1 <- dat1[match(id,dat1$id),]
dat2 <- dat2[match(id,dat2$id),]
dat <- cbind(dat1,dat2)
write.csv(dat, "./out/dmr_global_local_combined.csv")

r=cor(dat1$beta,dat2$beta)
idx1 <- dat1$fdr <= 0.05 & dat2$fdr <= 0.05 
idx2 <- dat1$fdr <= 0.05 & dat2$fdr > 0.05
idx3 <- dat1$fdr > 0.05 & dat2$fdr <= 0.05
idx4 <- ! (idx1 | idx2 | idx3)

pdf("./out/dmr_global_local_effect_scatter_plot.pdf")
xlim=c(min(dat2$beta),max(dat2$beta))
ylim=c(min(dat1$beta),max(dat1$beta))
plot(dat2$beta[idx4],dat1$beta[idx4],xlim=xlim,ylim=ylim,xlab="Local Ancestry Effects",ylab="Global Ancestry Effects",main=paste0("DLPFC r=",round(r,2)))
points(dat2$beta[idx1],dat1$beta[idx1],pch=16,col="red")
points(dat2$beta[idx2],dat1$beta[idx2],pch=16,col="orange")
points(dat2$beta[idx3],dat1$beta[idx3],pch=16,col="blue")
legend <- c("Both FDR < 0.05", "Global FDR < 0.05", "Local FDR < 0.05","Both FDR > 0.05")
legend <- paste0(legend," (n=",c(sum(idx1),sum(idx2),sum(idx3),sum(idx4)),")")
legend("bottomright", legend = legend, col=c("red", "orange","blue","black"), pch=16)
abline(coef = c(0,1),col="yellow")
dev.off()