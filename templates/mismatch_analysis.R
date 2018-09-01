#!/usr/bin/env Rscript

library(dplyr)
library(xtable)
library(RColorBrewer)
library(ggplot2)

mismatch <- read.table("$txt",header=T,sep="\\t")

###pie ofsgRNA readcount
mapped_readcount <- mismatch %>% subset(mis=="0") %>% unique() %>% nrow()
mis_readcount<- mismatch %>% subset(mis=="1"|mis=="2"|mis=="3") %>%unique() %>% nrow()
unmapped_readcount <- mismatch %>% subset(mis=="unmapped") %>%unique() %>% nrow()
sgRNA_readcount <- data.frame("perfect_mapped"=mapped_readcount,"mismatch"=mis_readcount,"unmapped"=unmapped_readcount)
ratio <- sprintf("%.2f",100*sgRNA_readcount[1,]/sum(sgRNA_readcount[1,])) 
ratio <- paste(ratio,"%",sep="")
pie.sales <- c(mapped_readcount,mis_readcount,unmapped_readcount)
names(pie.sales) <- paste(c("perfect mapped","mismatch","unmapped"), ratio,sep="\\n")
png(file='01pie.png')
pie(pie.sales,
  col = brewer.pal(3,"Pastel1"),
  main='The readcounts of extracted sgRNA aligned to designed library',col.main=4)
dev.off()
###his
sgRNA_readcount2 <- sgRNA_readcount %>% t() %>%as.data.frame()
his <- ggplot(sgRNA_readcount2,aes(x=c("perfect mapped","mismatch","unmapped"),y=V1))+ geom_col(fill=brewer.pal(3,"Pastel1"))+
  xlab("sgRNA type")+ylab("sgNA readcount")+
  geom_text(label=sgRNA_readcount2\$V1,colour = "black", vjust=00)+
  theme(axis.text.x = element_text(size = 11, family = "myFont", vjust = 0.5, hjust = 0.5))
ggsave("02his.png")
  
perfect_map <- mismatch %>% subset(mis=="0")
final_map_num <- perfect_map %>% count(lib_id)
su <- summary(final_map_num\$n)##统计summary信息
su_value <- sub("Min.   :|1st Qu.:|Median :|Mean   :|3rd Qu.:|Max.   :","",su)

###distribution
q <- quantile(final_map_num\$n,c(.9,.1))##计算百分位数
ggplot(final_map_num,aes(x=final_map_num\$n))+ geom_density(alpha = 0.3,color="red",,fill="blue")+xlab("sgRNA readcount")+ylab("percentage of sgRNA" )+
  geom_vline(data=q[2],linetype=2,size=0.2,xintercept = q[2])+geom_vline(data=q[1],linetype=2,size=0.2,xintercept = q[1])+
  annotate(geom="text", x=3, y=0.001, label=paste0("10%:",q[2]), color="black",size = 3)+
  annotate(geom="text", x=q[1]+10, y=0.0016, label=paste0("90%:",q[1]), color="black",size = 3)+
  annotate(geom="text", x=3500, y=0.0027, label=paste0("Min:",su[1]), color="black",size = 4)+
  annotate(geom="text", x=3500, y=0.0026, label=paste0("Mean:",format(su[4],digits=2)), color="black",size = 4)+
  annotate(geom="text", x=3500, y=0.0025, label=paste0("Median:",su[3]), color="black",size = 4)+
  annotate(geom="text", x=3500, y=0.0024, label=paste0("Max:",su[6]), color="black",size = 4)
ggsave("03distribution.png")


###norm
final_map_num\$log_2 <- log(final_map_num\$n,2)
final_map_num <- final_map_num %>% arrange(log_2)
q_norm <- quantile(final_map_num\$log_2,c(.9,.1))
ggplot(final_map_num,aes(x=1:nrow(final_map_num),y=final_map_num\$log_2))+ geom_point(alpha = 0.3,color="black")+xlab("sgRNA rank")+ylab("log2(sgRNA readcount)")+
  geom_hline(data=q_norm[2],linetype=2,size=0.2,yintercept = q_norm[2])+
  geom_hline(data=q_norm[1],linetype=2,size=0.2,yintercept = q_norm[1])+
  annotate(geom="text", x=nrow(final_map_num)/4*2, y=q_norm[2],label=paste0("(",format(q_norm[2],digits=2),",0.1)"), color="black",size = 4)+
  annotate(geom="text", x=nrow(final_map_num)/4*1, y=q_norm[1], label=paste0("(",format(q_norm[1],digits=2),",0.9)"), color="black",size = 4)+
  annotate(geom="segment",x=nrow(final_map_num)/4*3.5,xend=nrow(final_map_num)/4*3.5, y=q_norm[1],yend =q_norm[2], color="black")+
  annotate(geom="text", x=nrow(final_map_num)/4*2.5, y=(q_norm[1]+q_norm[2])/2, label=paste0("均一性：",format(2^(q_norm[1]-q_norm[2]),digits=2)), color="red",size = 4)
ggsave("04norm.png")


#fail <- read.table("/f/mulinlab/ke/test/fail_align_to_vector.txt",header=T)
#pair_mis <- fail %>% subset(cigar1=="*"&cigar2=="*") %>% nrow()
#fail_3 <- fail %>% subset(cigar1=="*"|cigar2=="*") %>% nrow()
#single_mis <- fail_3-pair_mis

