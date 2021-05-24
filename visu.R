setwd("E:/FengLab/glioma")
library(ggplot2)
umap_tx<-read.csv("umap_tx.csv",row.names=1)
umap_tx$subpop<-as.character(umap_tx$subpop)
pdf("umap.pdf",height=10,width=10)
umap_tx_cluster0<-subset(umap_tx,subpop=="0")
cluster0_x<-(min(umap_tx_cluster0$UMAP_1)+max(umap_tx_cluster0$UMAP_1))/2
cluster0_y<-(min(umap_tx_cluster0$UMAP_2)+max(umap_tx_cluster0$UMAP_2))/2
umap_tx_cluster1<-subset(umap_tx,subpop=="1")
cluster1_x<-(min(umap_tx_cluster1$UMAP_1)+max(umap_tx_cluster1$UMAP_1))/2
cluster1_y<-(min(umap_tx_cluster1$UMAP_2)+max(umap_tx_cluster1$UMAP_2))/2
umap_tx_cluster2<-subset(umap_tx,subpop=="2")
cluster2_x<-(min(umap_tx_cluster2$UMAP_1)+max(umap_tx_cluster2$UMAP_1))/2
cluster2_y<-(min(umap_tx_cluster2$UMAP_2)+max(umap_tx_cluster2$UMAP_2))/2
umap_tx_cluster3<-subset(umap_tx,subpop=="3")
cluster3_x<-(min(umap_tx_cluster3$UMAP_1)+max(umap_tx_cluster3$UMAP_1))/2
cluster3_y<-(min(umap_tx_cluster3$UMAP_2)+max(umap_tx_cluster3$UMAP_2))/2
umap_tx_cluster4<-subset(umap_tx,subpop=="4")
cluster4_x<-(min(umap_tx_cluster4$UMAP_1)+max(umap_tx_cluster4$UMAP_1))/2
cluster4_y<-(min(umap_tx_cluster4$UMAP_2)+max(umap_tx_cluster4$UMAP_2))/2
umap_tx_cluster5<-subset(umap_tx,subpop=="5")
cluster5_x<-(min(umap_tx_cluster5$UMAP_1)+max(umap_tx_cluster5$UMAP_1))/2
cluster5_y<-(min(umap_tx_cluster5$UMAP_2)+max(umap_tx_cluster5$UMAP_2))/2
umap_tx_cluster6<-subset(umap_tx,subpop=="6")
cluster6_x<-(min(umap_tx_cluster6$UMAP_1)+max(umap_tx_cluster6$UMAP_1))/2
cluster6_y<-(min(umap_tx_cluster6$UMAP_2)+max(umap_tx_cluster6$UMAP_2))/2
umap_tx_cluster7<-subset(umap_tx,subpop=="7")
cluster7_x<-(min(umap_tx_cluster7$UMAP_1)+max(umap_tx_cluster7$UMAP_1))/2
cluster7_y<-(min(umap_tx_cluster7$UMAP_2)+max(umap_tx_cluster7$UMAP_2))/2
umap_tx_cluster8<-subset(umap_tx,subpop=="8")
cluster8_x<-(min(umap_tx_cluster8$UMAP_1)+max(umap_tx_cluster8$UMAP_1))/2
cluster8_y<-(min(umap_tx_cluster8$UMAP_2)+max(umap_tx_cluster8$UMAP_2))/2
umap_tx_cluster9<-subset(umap_tx,subpop=="9")
cluster9_x<-(min(umap_tx_cluster9$UMAP_1)+max(umap_tx_cluster9$UMAP_1))/2
cluster9_y<-(min(umap_tx_cluster9$UMAP_2)+max(umap_tx_cluster9$UMAP_2))/2
umap_tx_cluster10<-subset(umap_tx,subpop=="10")
cluster10_x<-(min(umap_tx_cluster10$UMAP_1)+max(umap_tx_cluster10$UMAP_1))/2
cluster10_y<-(min(umap_tx_cluster10$UMAP_2)+max(umap_tx_cluster10$UMAP_2))/2
umap_tx_cluster11<-subset(umap_tx,subpop=="11")
cluster11_x<-(min(umap_tx_cluster11$UMAP_1)+max(umap_tx_cluster11$UMAP_1))/2
cluster11_y<-(min(umap_tx_cluster11$UMAP_2)+max(umap_tx_cluster11$UMAP_2))/2
umap_tx_cluster12<-subset(umap_tx,subpop=="12")
cluster12_x<-(min(umap_tx_cluster12$UMAP_1)+max(umap_tx_cluster12$UMAP_1))/2
cluster12_y<-(min(umap_tx_cluster12$UMAP_2)+max(umap_tx_cluster12$UMAP_2))/2
umap_tx_cluster13<-subset(umap_tx,subpop=="13")
cluster13_x<-(min(umap_tx_cluster13$UMAP_1)+max(umap_tx_cluster13$UMAP_1))/2
cluster13_y<-(min(umap_tx_cluster13$UMAP_2)+max(umap_tx_cluster13$UMAP_2))/2
umap_tx_cluster14<-subset(umap_tx,subpop=="14")
cluster14_x<-(min(umap_tx_cluster14$UMAP_1)+max(umap_tx_cluster14$UMAP_1))/2
cluster14_y<-(min(umap_tx_cluster14$UMAP_2)+max(umap_tx_cluster14$UMAP_2))/2
umap_tx_cluster15<-subset(umap_tx,subpop=="15")
cluster15_x<-(min(umap_tx_cluster15$UMAP_1)+max(umap_tx_cluster15$UMAP_1))/2
cluster15_y<-(min(umap_tx_cluster15$UMAP_2)+max(umap_tx_cluster15$UMAP_2))/2
umap_tx_cluster16<-subset(umap_tx,subpop=="16")
cluster16_x<-(min(umap_tx_cluster16$UMAP_1)+max(umap_tx_cluster16$UMAP_1))/2
cluster16_y<-(min(umap_tx_cluster16$UMAP_2)+max(umap_tx_cluster16$UMAP_2))/2
umap_tx_cluster17<-subset(umap_tx,subpop=="17")
cluster17_x<-(min(umap_tx_cluster17$UMAP_1)+max(umap_tx_cluster17$UMAP_1))/2
cluster17_y<-(min(umap_tx_cluster17$UMAP_2)+max(umap_tx_cluster17$UMAP_2))/2
umap_tx_cluster18<-subset(umap_tx,subpop=="18")
cluster18_x<-(min(umap_tx_cluster18$UMAP_1)+max(umap_tx_cluster18$UMAP_1))/2
cluster18_y<-(min(umap_tx_cluster18$UMAP_2)+max(umap_tx_cluster18$UMAP_2))/2
umap_tx_cluster19<-subset(umap_tx,subpop=="19")
cluster19_x<-(min(umap_tx_cluster19$UMAP_1)+max(umap_tx_cluster19$UMAP_1))/2
cluster19_y<-(min(umap_tx_cluster19$UMAP_2)+max(umap_tx_cluster19$UMAP_2))/2
umap_tx_cluster20<-subset(umap_tx,subpop=="20")
cluster20_x<-(min(umap_tx_cluster20$UMAP_1)+max(umap_tx_cluster20$UMAP_1))/2
cluster20_y<-(min(umap_tx_cluster20$UMAP_2)+max(umap_tx_cluster20$UMAP_2))/2
umap_tx_cluster21<-subset(umap_tx,subpop=="21")
cluster21_x<-(min(umap_tx_cluster21$UMAP_1)+max(umap_tx_cluster21$UMAP_1))/2
cluster21_y<-(min(umap_tx_cluster21$UMAP_2)+max(umap_tx_cluster21$UMAP_2))/2
umap_tx_cluster22<-subset(umap_tx,subpop=="22")
cluster22_x<-(min(umap_tx_cluster22$UMAP_1)+max(umap_tx_cluster22$UMAP_1))/2
cluster22_y<-(min(umap_tx_cluster22$UMAP_2)+max(umap_tx_cluster22$UMAP_2))/2
umap_tx_cluster23<-subset(umap_tx,subpop=="23")
cluster23_x<-(min(umap_tx_cluster23$UMAP_1)+max(umap_tx_cluster23$UMAP_1))/2
cluster23_y<-(min(umap_tx_cluster23$UMAP_2)+max(umap_tx_cluster23$UMAP_2))/2
umap_tx_cluster24<-subset(umap_tx,subpop=="24")
cluster24_x<-(min(umap_tx_cluster24$UMAP_1)+max(umap_tx_cluster24$UMAP_1))/2
cluster24_y<-(min(umap_tx_cluster24$UMAP_2)+max(umap_tx_cluster24$UMAP_2))/2
umap_tx_cluster25<-subset(umap_tx,subpop=="25")
cluster25_x<-(min(umap_tx_cluster25$UMAP_1)+max(umap_tx_cluster25$UMAP_1))/2
cluster25_y<-(min(umap_tx_cluster25$UMAP_2)+max(umap_tx_cluster25$UMAP_2))/2
umap_tx_cluster26<-subset(umap_tx,subpop=="26")
cluster26_x<-(min(umap_tx_cluster26$UMAP_1)+max(umap_tx_cluster26$UMAP_1))/2
cluster26_y<-(min(umap_tx_cluster26$UMAP_2)+max(umap_tx_cluster26$UMAP_2))/2
umap_tx_cluster27<-subset(umap_tx,subpop=="27")
cluster27_x<-(min(umap_tx_cluster27$UMAP_1)+max(umap_tx_cluster27$UMAP_1))/2
cluster27_y<-(min(umap_tx_cluster27$UMAP_2)+max(umap_tx_cluster27$UMAP_2))/2
umap_tx_cluster28<-subset(umap_tx,subpop=="28")
cluster28_x<-(min(umap_tx_cluster28$UMAP_1)+max(umap_tx_cluster28$UMAP_1))/2
cluster28_y<-(min(umap_tx_cluster28$UMAP_2)+max(umap_tx_cluster28$UMAP_2))/2
umap_tx_cluster29<-subset(umap_tx,subpop=="29")
cluster29_x<-(min(umap_tx_cluster29$UMAP_1)+max(umap_tx_cluster29$UMAP_1))/2
cluster29_y<-(min(umap_tx_cluster29$UMAP_2)+max(umap_tx_cluster29$UMAP_2))/2
umap_tx$x_val<-c(rep(cluster0_x,dim(umap_tx_cluster0)[1]), rep(cluster1_x,dim(umap_tx_cluster1)[1]), 
                 rep(cluster2_x,dim(umap_tx_cluster2)[1]), rep(cluster3_x,dim(umap_tx_cluster3)[1]), 
                 rep(cluster4_x,dim(umap_tx_cluster4)[1]), rep(cluster5_x,dim(umap_tx_cluster5)[1]), 
                 rep(cluster6_x,dim(umap_tx_cluster6)[1]), rep(cluster7_x,dim(umap_tx_cluster7)[1]), 
                 rep(cluster8_x,dim(umap_tx_cluster8)[1]), rep(cluster9_x,dim(umap_tx_cluster9)[1]), 
                 rep(cluster10_x,dim(umap_tx_cluster10)[1]), rep(cluster11_x,dim(umap_tx_cluster11)[1]), 
                 rep(cluster12_x,dim(umap_tx_cluster12)[1]), rep(cluster13_x,dim(umap_tx_cluster13)[1]), 
                 rep(cluster14_x,dim(umap_tx_cluster14)[1]), rep(cluster15_x,dim(umap_tx_cluster15)[1]), 
                 rep(cluster16_x,dim(umap_tx_cluster16)[1]), rep(cluster17_x,dim(umap_tx_cluster17)[1]), 
                 rep(cluster18_x,dim(umap_tx_cluster18)[1]), rep(cluster19_x,dim(umap_tx_cluster19)[1]), 
                 rep(cluster20_x,dim(umap_tx_cluster20)[1]), rep(cluster21_x,dim(umap_tx_cluster21)[1]),
                 rep(cluster22_x,dim(umap_tx_cluster22)[1]), rep(cluster23_x,dim(umap_tx_cluster23)[1]), 
                 rep(cluster24_x,dim(umap_tx_cluster24)[1]), rep(cluster25_x,dim(umap_tx_cluster25)[1]), 
                 rep(cluster26_x,dim(umap_tx_cluster26)[1]), rep(cluster27_x,dim(umap_tx_cluster27)[1]), 
                 rep(cluster28_x,dim(umap_tx_cluster28)[1]), rep(cluster29_x,dim(umap_tx_cluster29)[1]))
umap_tx$y_val<-c(rep(cluster0_y,dim(umap_tx_cluster0)[1]), rep(cluster1_y,dim(umap_tx_cluster1)[1]), 
                 rep(cluster2_y,dim(umap_tx_cluster2)[1]), rep(cluster3_y,dim(umap_tx_cluster3)[1]), 
                 rep(cluster4_y,dim(umap_tx_cluster4)[1]), rep(cluster5_y,dim(umap_tx_cluster5)[1]), 
                 rep(cluster6_y,dim(umap_tx_cluster6)[1]), rep(cluster7_y,dim(umap_tx_cluster7)[1]), 
                 rep(cluster8_y,dim(umap_tx_cluster8)[1]), rep(cluster9_y,dim(umap_tx_cluster9)[1]), 
                 rep(cluster10_y,dim(umap_tx_cluster10)[1]), rep(cluster11_y,dim(umap_tx_cluster11)[1]), 
                 rep(cluster12_y,dim(umap_tx_cluster12)[1]), rep(cluster13_y,dim(umap_tx_cluster13)[1]), 
                 rep(cluster14_y,dim(umap_tx_cluster14)[1]), rep(cluster15_y,dim(umap_tx_cluster15)[1]), 
                 rep(cluster16_y,dim(umap_tx_cluster16)[1]), rep(cluster17_y,dim(umap_tx_cluster17)[1]), 
                 rep(cluster18_y,dim(umap_tx_cluster18)[1]), rep(cluster19_y,dim(umap_tx_cluster19)[1]), 
                 rep(cluster20_y,dim(umap_tx_cluster20)[1]), rep(cluster21_y,dim(umap_tx_cluster21)[1]),
                 rep(cluster22_y,dim(umap_tx_cluster22)[1]), rep(cluster23_y,dim(umap_tx_cluster23)[1]), 
                 rep(cluster24_y,dim(umap_tx_cluster24)[1]), rep(cluster25_y,dim(umap_tx_cluster25)[1]), 
                 rep(cluster26_y,dim(umap_tx_cluster26)[1]), rep(cluster27_y,dim(umap_tx_cluster27)[1]), 
                 rep(cluster28_y,dim(umap_tx_cluster28)[1]), rep(cluster29_y,dim(umap_tx_cluster29)[1]))
umap_tx$name<-c(rep("C0",dim(umap_tx_cluster0)[1]), rep("C1",dim(umap_tx_cluster1)[1]), 
                rep("C2",dim(umap_tx_cluster2)[1]), rep("C3",dim(umap_tx_cluster3)[1]), 
                rep("C4",dim(umap_tx_cluster4)[1]), rep("C5",dim(umap_tx_cluster5)[1]), 
                rep("C6",dim(umap_tx_cluster6)[1]), rep("C7",dim(umap_tx_cluster7)[1]), 
                rep("C8",dim(umap_tx_cluster8)[1]), rep("C9",dim(umap_tx_cluster9)[1]), 
                rep("C10",dim(umap_tx_cluster10)[1]), rep("C11",dim(umap_tx_cluster11)[1]), 
                rep("C12",dim(umap_tx_cluster12)[1]), rep("C13",dim(umap_tx_cluster13)[1]), 
                rep("C14",dim(umap_tx_cluster14)[1]), rep("C15",dim(umap_tx_cluster15)[1]), 
                rep("C16",dim(umap_tx_cluster16)[1]), rep("C17",dim(umap_tx_cluster17)[1]), 
                rep("C18",dim(umap_tx_cluster18)[1]), rep("C19",dim(umap_tx_cluster19)[1]), 
                rep("C20",dim(umap_tx_cluster20)[1]), rep("C21",dim(umap_tx_cluster21)[1]), 
                rep("C22",dim(umap_tx_cluster22)[1]), rep("C23",dim(umap_tx_cluster23)[1]), 
                rep("C24",dim(umap_tx_cluster24)[1]), rep("C25",dim(umap_tx_cluster25)[1]), 
                rep("C26",dim(umap_tx_cluster26)[1]), rep("C27",dim(umap_tx_cluster27)[1]), 
                rep("C28",dim(umap_tx_cluster28)[1]), rep("C29",dim(umap_tx_cluster29)[1]))


labeldata<-data.frame(x_value=c(cluster0_x, cluster1_x, cluster2_x, cluster3_x, cluster4_x, cluster5_x, 
                                cluster6_x, cluster7_x, cluster8_x, cluster9_x, cluster10_x, cluster11_x, 
                                cluster12_x, cluster13_x, cluster14_x, cluster15_x, cluster16_x, cluster17_x, 
                                cluster18_x, cluster19_x, cluster20_x, cluster21_x, cluster22_x, cluster23_x, 
                                cluster24_x, cluster25_x, cluster26_x, cluster27_x, cluster28_x, cluster29_x),
                      y_value=c(cluster0_y, cluster1_y, cluster2_y, cluster3_y, cluster4_y, cluster5_y, 
                                cluster6_y, cluster7_y, cluster8_y, cluster9_y, cluster10_y, cluster11_y, 
                                cluster12_y, cluster13_y, cluster14_y, cluster15_y, cluster16_y, cluster17_y, 
                                cluster18_y, cluster19_y, cluster20_y, cluster21_y, cluster22_y, cluster23_y, 
                                cluster24_y, cluster25_y, cluster26_y, cluster27_y, cluster28_y, cluster29_y),
                     name=c("c0", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", 
                            "c13", "c14", "c15", "c16", "c17", "c18", "c19", "c20", "c21", "c22", "c23", "c24", 
                            "c25", "c26", "c27", "c28", "c29"))
ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2, color=subpop)) + 
  geom_point(stat= "identity",aes(),alpha=0.4,show.legend = TRUE)+
  geom_text(aes(x=x_val,y=y_val,label = name))
dev.off()
umap_tx_gliomaP<-umap_tx[1:38052,]
umap_tx_gliomaF<-umap_tx[38053:75941,]
pdf("umap.gliomaF.pdf",height=10,width=10)
ggplot(umap_tx_gliomaF, aes(x=UMAP_1, y=UMAP_2, color=subpop)) + 
  geom_point(stat= "identity",aes(),alpha=0.4,show.legend = TRUE)
dev.off()
pdf("umap.gliomaP.pdf",height=10,width=10)
ggplot(umap_tx_gliomaP, aes(x=UMAP_1, y=UMAP_2, color=subpop)) + 
  geom_point(stat= "identity",aes(),alpha=0.4,show.legend = TRUE)
dev.off()
pdf("umap.gliomaP1.pdf",height=10,width=10)
umap_tx_gliomaP1<-umap_tx[1:22489,]
ggplot(umap_tx_gliomaP1, aes(x=UMAP_1, y=UMAP_2, color=subpop)) + 
  geom_point(stat= "identity",aes(),alpha=0.4,show.legend = TRUE)
dev.off()
pdf("umap.gliomaP2.pdf",height=10,width=10)
umap_tx_gliomaP2<-umap_tx[22490:32395,]
ggplot(umap_tx_gliomaP2, aes(x=UMAP_1, y=UMAP_2, color=subpop)) + 
  geom_point(stat= "identity",aes(),alpha=0.4,show.legend = TRUE)
dev.off()
pdf("umap.gliomaP3.pdf",height=10,width=10)
umap_tx_gliomaP3<-umap_tx[32396:38052,]
ggplot(umap_tx_gliomaP3, aes(x=UMAP_1, y=UMAP_2, color=subpop)) + 
  geom_point(stat= "identity",aes(),alpha=0.4,show.legend = TRUE)
dev.off()
pdf("umap.gliomaF2.pdf",height=10,width=10)
umap_tx_gliomaF2<-umap_tx[38053:55141,]
ggplot(umap_tx_gliomaF2, aes(x=UMAP_1, y=UMAP_2, color=subpop)) + 
  geom_point(stat= "identity",aes(),alpha=0.4,show.legend = TRUE)
dev.off()
pdf("umap.gliomaF3.pdf",height=10,width=10)
umap_tx_gliomaF3<-umap_tx[55142:66082,]
ggplot(umap_tx_gliomaF3, aes(x=UMAP_1, y=UMAP_2, color=subpop)) + 
  geom_point(stat= "identity",aes(),alpha=0.4,show.legend = TRUE)
dev.off()
pdf("umap.gliomaF4.pdf",height=10,width=10)
umap_tx_gliomaF4<-umap_tx[66083:75941,]
ggplot(umap_tx_gliomaF4, aes(x=UMAP_1, y=UMAP_2, color=subpop)) + 
  geom_point(stat= "identity",aes(),alpha=0.4,show.legend = TRUE)
dev.off()



tsne_tx<-read.csv("tsne_tx.csv",row.names=1)
tsne_tx$subpop<-as.character(tsne_tx$subpop)
pdf("tsne.pdf",height=10,width=10)
ggplot(tsne_tx, aes(x=tSNE_1, y=tSNE_2, color=subpop)) + 
  geom_point(stat= "identity",aes(),alpha=0.4,show.legend = TRUE)
dev.off()
tsne_tx_gliomaP<-tsne_tx[1:38052,]
tsne_tx_gliomaF<-tsne_tx[38053:75941,]
pdf("tsne.gliomaF.pdf",height=10,width=10)
ggplot(tsne_tx_gliomaF, aes(x=tSNE_1, y=tSNE_2, color=subpop)) + 
  geom_point(stat= "identity",aes(),alpha=0.4,show.legend = TRUE)
dev.off()
pdf("tsne.gliomaP.pdf",height=10,width=10)
ggplot(tsne_tx_gliomaP, aes(x=tSNE_1, y=tSNE_2, color=subpop)) + 
  geom_point(stat= "identity",aes(),alpha=0.4,show.legend = TRUE)
dev.off()
pdf("tsne.gliomaP1.pdf",height=10,width=10)
tsne_tx_gliomaP1<-tsne_tx[1:22489,]
ggplot(tsne_tx_gliomaP1, aes(x=tSNE_1, y=tSNE_2, color=subpop)) + 
  geom_point(stat= "identity",aes(),alpha=0.4,show.legend = TRUE)
dev.off()
pdf("tsne.gliomaP2.pdf",height=10,width=10)
tsne_tx_gliomaP2<-tsne_tx[22490:32395,]
ggplot(tsne_tx_gliomaP2, aes(x=tSNE_1, y=tSNE_2, color=subpop)) + 
  geom_point(stat= "identity",aes(),alpha=0.4,show.legend = TRUE)
dev.off()
pdf("tsne.gliomaP3.pdf",height=10,width=10)
tsne_tx_gliomaP3<-tsne_tx[32396:38052,]
ggplot(tsne_tx_gliomaP3, aes(x=tSNE_1, y=tSNE_2, color=subpop)) + 
  geom_point(stat= "identity",aes(),alpha=0.4,show.legend = TRUE)
dev.off()
pdf("tsne.gliomaF2.pdf",height=10,width=10)
tsne_tx_gliomaF2<-tsne_tx[38053:55141,]
ggplot(tsne_tx_gliomaF2, aes(x=tSNE_1, y=tSNE_2, color=subpop)) + 
  geom_point(stat= "identity",aes(),alpha=0.4,show.legend = TRUE)
dev.off()
pdf("tsne.gliomaF3.pdf",height=10,width=10)
tsne_tx_gliomaF3<-tsne_tx[55142:66082,]
ggplot(tsne_tx_gliomaF3, aes(x=tSNE_1, y=tSNE_2, color=subpop)) + 
  geom_point(stat= "identity",aes(),alpha=0.4,show.legend = TRUE)
dev.off()
pdf("tsne.gliomaF4.pdf",height=10,width=10)
tsne_tx_gliomaF4<-tsne_tx[66083:75941,]
ggplot(tsne_tx_gliomaF4, aes(x=tSNE_1, y=tSNE_2, color=subpop)) + 
  geom_point(stat= "identity",aes(),alpha=0.4,show.legend = TRUE)
dev.off()