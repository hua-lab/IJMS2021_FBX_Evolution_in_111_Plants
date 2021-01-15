library("tidyr")
library("ggplot2")


##################### Figure S2B
	
d<-read.csv("Data_S2_Number_of_sequences_processed_through_CTT_annotation_in_each_plant_genome.txt",head=T)

	d_rownames<-rownames(d)<-d[,1]
	d<-d[,-1]
	d<-apply(d,2,as.numeric)
	rownames(d)<-d_rownames
	
	cor<-cor(t(d),method="spearman")
	datasets<-c("4c","4d","5p","7eu","7p","prior","size_g")
	cor.df<-data.frame(datasets,cor)
	
	cor.long<-pivot_longer(data=cor.df,
						cols=-c(1),
						names_to="Class",
						values_to="cor")
	
	my_Theme = theme(
 	  	 		axis.title.x = element_text(size = 8),
 				axis.text.x = element_text(size = 7),
 	 			axis.title.y = element_text(size = 6))
		

	heatmap_s2b<-ggplot(data=cor.long, aes(x=datasets,y= Class,fill=cor) )+
						geom_tile()+
						xlab(label="Sample") +			
						scale_y_discrete(limits = as.factor(cor.long$Class)[1:7]) +		
						scale_fill_gradient(name = "rho",
			                		  		low = "#FFFFFF",
			                    			high = "#2196f3") +
						theme(strip.placement = "outside",plot.title = element_text(hjust = 0.5))+
						my_Theme 


######################################################################################################### Figure 1		

	sum<-rowSums(d[1:6,])
	annot<-d[4:6,]
	annot<-annot[rev(rownames(annot)),]
	fbx_peps<-colSums(d[5:6,])


	pdf("Figure 1_fbx_ctt_step_annotation_statistics.pdf",width=9,height=9)

			par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=2)
			
			layout(matrix(c(1,4,2,4,3,4),nrow=2,ncol=3))

			barplot(sum,lwd = 0.3,cex.axis=0.5, cex.names=0.5,beside=FALSE,legend = FALSE,ylab="Sequence Dataset Size")
			
			plot(d[6,],d[5,],cex=1,pch=1, xlab="prior",ylab="new")
			trend<-lm( d[5,] ~ d[6,])
			abline(trend, col="blue",lwd=2)
			text(2000, 500, "rho=0.54, p-value=6.4e-10")

			plot(d[4,],d[5,],cex=1,pch=1,xlab="pseudo",ylab="new")
			trend<-lm( d[5,] ~ d[4,])
			abline(trend, col="blue",lwd=2)
			text(2000, 500, "rho=0.79, p-value=2.2e-16")	
			
			barplot(annot,col=c("light blue","gray","red"),ylim=c(0,4000),lwd = 0.3,cex.axis=0.5, cex.names=0.5,beside=FALSE,
			legend = rownames(annot), args.legend = list(x = "top", bty = "n", inset=c(-0.15, 0)))
			abline(h=c("500","1000","1500","2000","2500","3000","3500","4000"), col = "gray", lty = 5)				

	dev.off()


##################### Figure S2A


	Sample.name<-c("4c","4d","5p","7eu","7p","prior")
	d.df<-data.frame(Sample.name,d[1:6,])

	d.long<-pivot_longer(data=d.df,
						cols=-c(1),
						names_to="Class",
						values_to="count")

	d.long$log.count <- log(d.long$count,2)

	my_Theme = theme(
 	  	 		axis.title.x = element_text(size = 8),
 				axis.text.x = element_text(size = 7),
 	 			axis.title.y = element_text(size = 6))

	heatmap_s2a<-ggplot(data=d.long, aes(x=Sample.name,y= Class,fill=log.count) )+
						geom_tile()+
						xlab(label="Sample") +			
						scale_y_discrete(limits = as.factor(d.long$Class)[1:112]) +		
						scale_fill_gradient(name = "og.count",
			                      low = "#FFFFFF",
			                      high = "#012345") +
						theme(strip.placement = "outside",plot.title = element_text(hjust = 0.5))+
						my_Theme
			
	pdf("Figure S2A_fbx_ctt_annotation_seq_counts.pdf",width=3,height=10)
			par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=2)
			heatmap_s2a
	dev.off()


##################### Figure S2B

	pdf("Figure S2B_ctt_annotation_seq_count_cor.pdf",width=6,height=3)
			par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=2)
			heatmap_s2b
	dev.off()



##################### Figure S2C

	pdf("Figure S2C_fbx_pep_count_vs_genome_size.pdf",width=3,height=9)
	
			par(mar=c(4.1,4.1,1.1,2.1), mgp=c(3,1,0), las=2)	
			layout(matrix(c(1,2,3),nrow=3,ncol=1))
	
			plot(d[7,],d[5,],xlab="Genome Size",ylab="New")
			trend<-lm(d[5,] ~ d[7,])
			abline(trend, col="blue",lwd=2)
			text(3e+9, 500, "rho=0.5, p-value=2.9e-8")

			plot(d[7,],d[6,],cex=0.5,pch=1,xlab="Genome Size", ylab="Prior")
			trend<-lm(d[6,] ~ d[7,])
			abline(trend, col="blue",lwd=2)
			text(3e+9, 2500, "rho=0.25, p-value=7.7e-3")
		
			plot(d[7,],colSums(d[5:6,]),cex=0.5,pch=1,xlab="Genome Size",ylab="Total_Pep")
			trend<-lm(colSums(d[5:6,]) ~ d[7,])
			abline(trend, col="blue",lwd=2)
			text(3e+9, 2500, "rho=0.31, p-value=8.3e-4")	
	
	dev.off()













