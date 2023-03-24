dat1 <- read.table("./out/dmr.csv",header=T) # local
dat2 <- read.table("../out/dmr.csv",header=T) # global
dat2$id <- paste(dat2[,1],dat2[,2],dat2[,3],sep="_")

id <- intersect(dat1[,1],dat2$id)
dat1 <- dat1[match(id,dat1[,1]),]
dat2 <- dat2[match(id,dat2$id),]

p <- 0.05
sum(dat2$p <= p)
sum(dat1$p <= p)
sum(dat2$p <= p & dat1$p <= p)

p <- 0.05/nrow(dat1)
sum(dat2$p <= p)
sum(dat1$p <= p)
sum(dat2$p <= p & dat1$p <= p)
