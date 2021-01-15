#libraries
library("ggplot2")
library("gplots")
library("phytools")
library("tidyr")
library("fitdistrplus")
library("KSgeneral")
library("scales")
library("dendextend")


####################### retrieve count file of subfamilie  ####################### 

	ortho_count<-read.table("Data_S4_Count_data_matrix_of_5,858_FBX_subfamilies_in_111_plant_genomes.txt",header=T)
	# re-oder the rows so that the species are grouped as their order as in the species tree matrix
	ids<-ortho_count$id
	ortho_count<-ortho_count[,-1]
	ortho_count<-apply(ortho_count,2,as.numeric)
	rownames(ortho_count)<-ids
	ortho_count<-ortho_count[order(ortho_count[,1]),]
	ortho_count<-ortho_count[,-1]

	#replace "NA" with 0
	ortho_count<-as.data.frame(ortho_count) # need to convert to dataframe to replace NA with 0
 	for(k in 1:ncol(ortho_count)){ 	
		ortho_count[[k]][is.na(ortho_count[[k]])] <- 0 		
		} 


	dim(ortho_count)
	#	[1]  111 5858 # in total 111 species and 5,858 OrthoMCL groups (FBX subfamilies)

	species_count<-as.matrix(colSums(ortho_count !=0)) # Count number of species in each OthoMCL group
	species_count_none_zero<-species_count[species_count[,1]>0,] # some ortho groups are specific to different accessions
		 														# within a species due to initial analysis
	length(species_count_none_zero) #	[1] 5858

	ortho_count<-ortho_count[,colnames(ortho_count)%in%names(species_count_none_zero)]
	dim(ortho_count)
	#	[1]  111 5858
	
	
####################### read in four clusters of subfamilies  ####################### 
 
	clusters<-read.table("Data_S6_Data_matrix_of_four_clusters_of_FBX_subfamilies.txt",header=T)

	table(clusters[,2])
			
	#	c_1  c_2  c_3  c_4 
	#	 95   45 5523  195 

	# Clusters are re ordered based on the mean number of species

	cluster1<-clusters[clusters[,1]==3,]
	cluster2<-clusters[clusters[,1]==4,]
	cluster3<-clusters[clusters[,1]==2,]
	cluster4<-clusters[clusters[,1]==1,]

			
############################################################################
# Figure 5 & S5														   #
#																		   #
# Size comparison across species in different groups				       #
#  			1) Complete set,Figure S5A					                   #
#			2) All set of subfamilies, Figure 5A                           #
#			3) Orphan genes, Figure 5B                                     #
#			4) Cluster 1, Figure 5C,                                       #
#			5) Cluster 2, Figure 5D,                                       #
#			6) Cluster 3, Figure 5E,                                       #
#			7) Cluster 4, Figure 5F,                                       #
#			8) Group correlation, Figure S5B,                               #
############################################################################	
			
############################################################################# complete set, Figure S5A

	all<-read.csv("Data_S2_Number_of_sequences_processed_through_CTT_annotation_in_each_plant_genome.txt",header=T)
	all<-all[all$X%in%c("New_FBX","Prior_FBX"),]
	all<-all[,-1]
	all<-all[,rev(colnames(all))]
	size<-colSums(all)

	size_com<-size	
	
	var(size_com, na.rm = TRUE, use= "all.obs") #[1] 277286.8

	species<-seq(1,length(size_com)[1],by=1)
	species_family_size_df<-as.data.frame(cbind(species,size_com))


	pdf("Figure S5A_complete_family_size_trend.pdf",height=5,width=5)
		par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=2)
		ggplot(species_family_size_df, aes(species,size_com)) + geom_point() + geom_smooth(method = "loess")+ 
		scale_x_continuous(name = "species",
						breaks = seq(1, dim(species_family_size_df)[1], 1),
						labels=rownames(species_family_size_df)) +						
		theme(axis.text.x = element_text(face="italic", color="black", 
					                    size=3, angle=45,vjust = 0.5)) +											
		theme( axis.line = element_line(colour = "darkblue", 
					        size = 0.5, linetype = "solid"))
	dev.off()


############################################################################# all homologues, Figure 5A

	d<-t(ortho_count)				
	dx<-d[rownames(d)%in%rownames(clusters),]
	size_all<-colSums(dx)
	
	var(size_all, na.rm = TRUE, use= "all.obs") #[1] 219489.5

	species<-seq(1,length(size_all)[1],by=1)
	species_family_size_df<-as.data.frame(cbind(species,size_all))

	pdf("Figure 5A_all_subfamilies_size_trend.pdf",height=5,width=5)
		par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=2)
		ggplot(species_family_size_df, aes(species,size_all)) + geom_point() + geom_smooth(method = "loess") + 
		scale_x_continuous(name = "species",
							breaks = seq(1, dim(species_family_size_df)[1], 1),
							labels=rownames(species_family_size_df)) +					
		theme(axis.text.x = element_text(face="italic", color="black", 
						                    size=2, angle=45,vjust = 0.5)) +											
		theme( axis.line = element_line(colour = "darkblue", 
						        size = 0.5, linetype = "solid"))

	dev.off()


############################################################################# Orphan fbxes, Figure 5B

	orphan<-size_com-size_all
	size_orphan<-orphan
	
	var(size_orphan, na.rm = TRUE, use= "all.obs") #[1] 6314.636

	species<-seq(1,length(size_orphan)[1],by=1)
	species_family_size_df<-as.data.frame(cbind(species,size_orphan))


	pdf("Figure 5B_orphan_fbx_size_trend.pdf",height=5,width=5)
		par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=2)
		ggplot(species_family_size_df, aes(species,size_orphan)) + geom_point() + geom_smooth(method = "loess")+ 

		scale_x_continuous(name = "species",
							breaks = seq(1, dim(species_family_size_df)[1], 1),
							labels=rownames(species_family_size_df)) +
						
		theme(axis.text.x = element_text(face="italic", color="black", 
						                    size=3, angle=45,vjust = 0.5)) +				
							
		theme( axis.line = element_line(colour = "darkblue", 
						        size = 0.5, linetype = "solid"))

	dev.off()



#############################################################################  cluster 1, Figure 5C


	d<-t(ortho_count)				
	dx<-d[rownames(d)%in%rownames(cluster1),]
	size_c1<-colSums(dx)

	var(size_c1, na.rm = TRUE, use= "all.obs") #[1] 94258.1

	species<-seq(1,length(size_c1)[1],by=1)
	species_family_size_df<-as.data.frame(cbind(species,size_c1))


	pdf("Figure 5C_cluster1_family_size_trend.pdf",height=5,width=5)
		par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=2)
		ggplot(species_family_size_df, aes(species,size_c1)) + geom_point() + geom_smooth(method = "loess") + 

		scale_x_continuous(name = "species",
							breaks = seq(1, dim(species_family_size_df)[1], 1),
							labels=rownames(species_family_size_df)) +
						
		theme(axis.text.x = element_text(face="italic", color="black", 
						                    size=2, angle=45,vjust = 0.5)) +				
							
		theme( axis.line = element_line(colour = "darkblue", 
						        size = 0.5, linetype = "solid"))

	dev.off()


			
#############################################################################  cluster 2, Figure 5D

	d<-t(ortho_count)				
	dx<-d[rownames(d)%in%rownames(cluster2),]
	size_c2<-colSums(dx)
	
	var(size_c2, na.rm = TRUE, use= "all.obs") #[1] 29443.69

	species<-seq(1,length(size_c2)[1],by=1)
	species_family_size_df<-as.data.frame(cbind(species,size_c2))


	pdf("Figure 5D_cluster2_family_size_trend.pdf",height=5,width=5)
		par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=2)
		ggplot(species_family_size_df, aes(species,size_c2)) + geom_point() + geom_smooth(method = "loess") + 

		scale_x_continuous(name = "species",
							breaks = seq(1, dim(species_family_size_df)[1], 1),
							labels=rownames(species_family_size_df)) +
						
		theme(axis.text.x = element_text(face="italic", color="black", 
						                    size=2, angle=45,vjust = 0.5)) +				
							
		theme( axis.line = element_line(colour = "darkblue", 
						        size = 0.5, linetype = "solid"))

	dev.off()



############################################################################# cluster 3, Figure 5E

	d<-t(ortho_count)				
	dx<-d[rownames(d)%in%rownames(cluster3),]
	size_c3<-colSums(dx)
	
	var(size_c3, na.rm = TRUE, use= "all.obs") #[1] 1218.677

	species<-seq(1,length(size_c3)[1],by=1)
	species_family_size_df<-as.data.frame(cbind(species,size_c3))


	pdf("Figure 5E_cluster3_family_size_trend.pdf",height=5,width=5)
		par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=2)

		ggplot(species_family_size_df, aes(species,size_c3)) + geom_point() + geom_smooth(method = "loess") + 

		scale_x_continuous(name = "species",
							breaks = seq(1, dim(species_family_size_df)[1], 1),
							labels=rownames(species_family_size_df)) +
						
		theme(axis.text.x = element_text(face="italic", color="black", 
						                    size=2, angle=45,vjust = 0.5)) +				
							
		theme( axis.line = element_line(colour = "darkblue", 
						        size = 0.5, linetype = "solid"))

	dev.off()



############################################################################# cluster 4, figure 5F

	d<-t(ortho_count)				
	dx<-d[rownames(d)%in%rownames(cluster4),]
	size_c4<-colSums(dx)
	
	var(size_c4, na.rm = TRUE, use= "all.obs") #[1] 6606.752

	species<-seq(1,length(size_c4)[1],by=1)
	species_family_size_df<-as.data.frame(cbind(species,size_c4))


	pdf("Figure 5F_cluster4_family_size_trend.pdf",height=5,width=5)
		par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=2)

		ggplot(species_family_size_df, aes(species,size_c4)) + geom_point() + geom_smooth(method = "loess") + 

		scale_x_continuous(name = "species",
							breaks = seq(1, dim(species_family_size_df)[1], 1),
							labels=rownames(species_family_size_df)) +
						
		theme(axis.text.x = element_text(face="italic", color="black", 
						                    size=2, angle=45,vjust = 0.5)) +				
							
		theme( axis.line = element_line(colour = "darkblue", 
						        size = 0.5, linetype = "solid"))

	dev.off()
	
	sd(size_com, na.rm = TRUE)/mean(size_com, na.rm = TRUE) #[1] 0.7448664
	sd(size_all, na.rm = TRUE)/mean(size_all, na.rm = TRUE) #[1] 0.7379687
	sd(size_orphan, na.rm = TRUE)/mean(size_orphan, na.rm = TRUE) #[1] 1.102159
	sd(size_c1, na.rm = TRUE)/mean(size_c1, na.rm = TRUE) # [1] 0.9359437
	sd(size_c2, na.rm = TRUE)/mean(size_c2, na.rm = TRUE) # [1] 1.58103
	sd(size_c3, na.rm = TRUE)/mean(size_c3, na.rm = TRUE) # [1] 0.7236153
	sd(size_c4, na.rm = TRUE)/mean(size_c4, na.rm = TRUE) # [1] 0.5417169
	
	
############################################################################# size correlation, Figure S5B

	fbx_sizes<-cbind(size_all,orphan,size_c1,size_c2,size_c3,size_c4)

	colnames(fbx_sizes)<-c("All","orphan","c1","c2","c3","c4")
		
	cor<-cor(fbx_sizes, method="spearman")

	rowDistance=dist(cor,method="manhattan")
	rowCluster = hclust(rowDistance,method="ward.D2")
	rowDend = as.dendrogram(rowCluster)
	rowDend = reorder(rowDend, rowSums(cor))

	colDistance=dist(t(cor),method="manhattan")
	colCluster = hclust(colDistance,method="ward.D2")
	colDend = as.dendrogram(colCluster)
	colDend = reorder(colDend, colSums(cor))

	mycol <- colorpanel(n=9,low="blue",mid="light yellow",high="red")

	pdf ("Figure S5B_fbx_size_all_orphan_1_2_3_4_correlation_among_clusters.pdf", family="Times", height=10, width=10)
		par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=0)
		heatmap.2(cor,Rowv=rowDend, Colv=colDend, col=mycol,  keysize=2,trace=c("none"),density.info=c("none"),cexRow=0.5,cexCol=0.5)
	dev.off()


############################################################################
# Figure S9 															   #
#																		   #
#  Cluster size ratio comparison across species in different groups		   #
#			1) Ratio of all subfamilies FBXes, Figure S9A, 				   #
#           2) Ratio of orphan FBXes, Figure S9B,  						   # 
#			3) Ratio of cluster 1 FBXes, Figure S9C,                       #
#			4) Ratio of cluster 2 FBXes, Figure S9D,                       #
#			5) Ratio of cluster 3 FBXes, Figure S9E,                       #
#			6) Ratio of cluster 4 FBXes, Figure S9F,                       #
############################################################################	


############################################################################# homologues_ratio distribution , Figure S9A
	size_all_ratio<-size_all/size_com
	size<-size_all_ratio

	species<-seq(1,length(size)[1],by=1)
	species_family_size_df<-as.data.frame(cbind(species,size))

	pdf("Figure S9A_size_all_ratio.pdf",height=5,width=5)
		par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=2)

		ggplot(species_family_size_df, aes(species,size)) + geom_point() + geom_smooth(method = "loess")+ 

		scale_x_continuous(name = "species",
							breaks = seq(1, dim(species_family_size_df)[1], 1),
							labels=rownames(species_family_size_df)) +
						
		theme(axis.text.x = element_text(face="italic", color="black", 
						                    size=3, angle=45,vjust = 0.5)) +				
							
		theme( axis.line = element_line(colour = "darkblue", 
						        size = 0.5, linetype = "solid"))

	dev.off()
	


############################################################################# orphan_ratio distribution , Figure S6B
	orphan_ratio<-orphan/size_com
	size<-orphan_ratio

	species<-seq(1,length(size)[1],by=1)
	species_family_size_df<-as.data.frame(cbind(species,size))

	pdf("Figure S9B_orphan_ratio.pdf",height=5,width=5)
		par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=2)

		ggplot(species_family_size_df, aes(species,size)) + geom_point() + geom_smooth(method = "loess")+ 

		scale_x_continuous(name = "species",
							breaks = seq(1, dim(species_family_size_df)[1], 1),
							labels=rownames(species_family_size_df)) +
						
		theme(axis.text.x = element_text(face="italic", color="black", 
						                    size=3, angle=45,vjust = 0.5)) +				
							
		theme( axis.line = element_line(colour = "darkblue", 
						        size = 0.5, linetype = "solid"))

	dev.off()



############################################################################# c1_ratio distribution, Figure S7C

	size_c1_ratio<-size_c1/size_com
	size<-size_c1_ratio

	species<-seq(1,length(size)[1],by=1)
	species_family_size_df<-as.data.frame(cbind(species,size))

	pdf("Figure S9C_size_c1_ratio.pdf",height=5,width=5)
	par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=2)

	ggplot(species_family_size_df, aes(species,size)) + geom_point() + geom_smooth(method = "loess")+ 

		scale_x_continuous(name = "species",
							breaks = seq(1, dim(species_family_size_df)[1], 1),
							labels=rownames(species_family_size_df)) +
						
		theme(axis.text.x = element_text(face="italic", color="black", 
						                    size=3, angle=45,vjust = 0.5)) +				
							
		theme( axis.line = element_line(colour = "darkblue", 
						        size = 0.5, linetype = "solid"))

	dev.off()
	


############################################################################# c2_ratio distribution , Figure S6D
	size_c2_ratio<-size_c2/size_com
	size<-size_c2_ratio

	species<-seq(1,length(size)[1],by=1)
	species_family_size_df<-as.data.frame(cbind(species,size))

	pdf("Figure S9D_size_c2_ratio.pdf",height=5,width=5)
	par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=2)

	ggplot(species_family_size_df, aes(species,size)) + geom_point() + geom_smooth(method = "loess")+ 

		scale_x_continuous(name = "species",
							breaks = seq(1, dim(species_family_size_df)[1], 1),
							labels=rownames(species_family_size_df)) +
						
		theme(axis.text.x = element_text(face="italic", color="black", 
						                    size=3, angle=45,vjust = 0.5)) +				
							
		theme( axis.line = element_line(colour = "darkblue", 
						        size = 0.5, linetype = "solid"))

	dev.off()

	
############################################################################# c3_ratio distribution , Figure S9E

	size_c3_ratio<-size_c3/size_com
	size<-size_c3_ratio

	species<-seq(1,length(size)[1],by=1)
	species_family_size_df<-as.data.frame(cbind(species,size))

	pdf("Figure S9E_size_c3_ratio.pdf",height=5,width=5)
		par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=2)
		
		ggplot(species_family_size_df, aes(species,size)) + geom_point() + geom_smooth(method = "loess")+ 

		scale_x_continuous(name = "species",
							breaks = seq(1, dim(species_family_size_df)[1], 1),
							labels=rownames(species_family_size_df)) +
						
		theme(axis.text.x = element_text(face="italic", color="black", 
						                    size=3, angle=45,vjust = 0.5)) +				
							
		theme( axis.line = element_line(colour = "darkblue", 
						        size = 0.5, linetype = "solid"))

	dev.off()
	

############################################################################# c4_ratio distribution , Figure S9F
	size_c4_ratio<-size_c4/size_com
	size<-size_c4_ratio

	species<-seq(1,length(size)[1],by=1)
	species_family_size_df<-as.data.frame(cbind(species,size))

	pdf("Figure S9F_size_c4_ratio.pdf",height=5,width=5)
		par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=2)

		ggplot(species_family_size_df, aes(species,size)) + geom_point() + geom_smooth(method = "loess")+ 

		scale_x_continuous(name = "species",
							breaks = seq(1, dim(species_family_size_df)[1], 1),
							labels=rownames(species_family_size_df)) +
						
		theme(axis.text.x = element_text(face="italic", color="black", 
						                    size=3, angle=45,vjust = 0.5)) +				
							
		theme( axis.line = element_line(colour = "darkblue", 
						        size = 0.5, linetype = "solid"))

	dev.off()


############################################################################
# Figures 8, S12, S10												  	   #
#																		   #
# Heatmap, SCP modeling, dendrogram cor				           		       #
#			1) Clusters in 6 groups of species, Figure 8A,                 #
#			2) heatmap Cluster 1, Figure 8C,                               #
#			3) heatmap Cluster 2, Figure 8D,                               #
#			4) heatmap Cluster 3, Figure 8E,                               #
#			5) heatmap Cluster 4, Figure 8F,   							   #
#			6) modeling SCPs in Cluster 4, Figure S11,                     #
#			6) modeling SCPs in Cluster 4, Figure 8G,                      #
#			7) modeling SCPs in Cluster 3, Figure S12,                     #
#			8) Cluster dendrogram cor, Figure S10,                         #
############################################################################	

		
############################################################################ Cluster size comparison in 6 groups of species, Figure 8A
	
		# Read in species tree file
		library("phytools")
		
		phytozome.tree<-read.tree("Data_S2_A_bifurcated_newick_species_tree_of_111_plant_genomes.txt")

		# convert species phylogeny to a data matrix
		tree_matrix<-compute.mr(phytozome.tree, type = "matrix")
		species_short<-substring(rownames(tree_matrix),1,5)
		species_short<-gsub("\\.","",species_short)
		
		# make a data frame for proportion of each group of FBX genes in each plant FBX family
		
		species_order<-seq(1,length(species_short),1)	
		fbx_sizes<-cbind(size_all,orphan,size_c1,size_c2,size_c3,size_c4)
		fbx_size_ratios<-fbx_sizes/size_com	
		fbx_sizes.df<-data.frame(species_order,fbx_size_ratios)
	
	
		library("tidyr")
		fbx_sizes.long<-pivot_longer(data=fbx_sizes.df,
								cols=-c(1,2),
								names_to="Class",
								values_to="size")

		my_Theme = theme(
			  			axis.title.x = element_text(size = 8),
		 			   	axis.text.x = element_text(size = 7),
						axis.title.y = element_text(size = 6))


		heatmap<-ggplot(data=fbx_sizes.long, aes(x=species_order,y= Class,fill=size) )+
												geom_tile()+
												xlab(label="Sample") +			
					scale_y_discrete(limits = as.factor(fbx_sizes.long$Class)[1:6]) +			
					scale_fill_gradient(name = "Count",
					                      low = "#FFFFFF",					
					                      high = "black") +
					theme(strip.placement = "outside",plot.title = element_text(hjust = 0.5))+
					my_Theme 
			
		pdf("Figure 8A_species_6_group_size_proportion_in_clusters_heatmap.pdf",height=3,width=15)	
			heatmap
		dev.off()
	
		
		rowDistance=dist(fbx_size_ratios[,2:6],method="manhattan")
		rowCluster = hclust(rowDistance,method="ward.D2")
	
		mycl <- cutree(rowCluster, 6)
	
	    clusterCols <- rainbow(length(unique(mycl)))
	    myClusterSideBar <- clusterCols[mycl]
	
	    myheatcol <- rev(redgreen(75))
	
		pdf("Figure 8A_Sidebar_species_6_group_size_proportion_in_clusters_heatmap2_w_rowsidecolors.pdf",height=5,width=10)
			
			heatmap.2(fbx_size_ratios[,2:6], main="Hierarchical Cluster", Rowv=FALSE, Colv=FALSE, dendrogram=c("none"),scale="row", 
			col=myheatcol, density.info="none", trace="none", RowSideColors= myClusterSideBar)
		dev.off()
	
	
		
############################################################################ heatmap_Cluster 1, Figure 8C

		#Color lable 6 groups of species in Figure 8A
		
		br<-species_short[1:29]
		s_r<-species_short[30:60]
		nr_dicots<-species_short[61:72]
		poaceae<-species_short[73:87]
		basal<-species_short[88:101]
		algae<-species_short[102:111]
		
		col_br<-rep(2,length(br))
		col_s_r<-rep(1,length(s_r))
		col_nr_dicots<-rep(4,length(nr_dicots))
		col_poaceae<-rep(5,length(poaceae))
		col_basal<-rep(3,length(basal))
		col_algae<-rep(6,length(algae))
		
		names(col_br)<-br
		names(col_s_r)<-s_r
		names(col_nr_dicots)<-nr_dicots
		names(col_poaceae)<-poaceae
		names(col_basal)<-basal
		names(col_algae)<-algae
		
		
		mylist<-c(col_br,col_s_r,col_nr_dicots,col_poaceae,col_basal,col_algae)
		
	
	    clusterCols <- rainbow(length(unique(mylist)))
	    myClusterSideBar <- clusterCols[mylist]
		
#################################		
		
		d<-t(ortho_count)				
		d1<-d[rownames(d)%in%rownames(cluster1),]
		
		d1<-t(d1)	
		
	 	dim(d1)[2]
				#[1] 5523

		rowDistance=dist(d1,method="manhattan")
		rowCluster = hclust(rowDistance,method="ward.D2")
		rowDend = as.dendrogram(rowCluster)
		rowDend = reorder(rowDend, rowSums(d1))
		dnd_cluster1<-rowDend

		colDistance=dist(t(d1),method="manhattan")
		colCluster = hclust(colDistance,method="ward.D2")
		colDend = as.dendrogram(colCluster)
		colDend = reorder(colDend, colSums(d1))

		my_palette<-colorRampPalette(c("light blue","yellow","red"))(n=19)
		col_breaks=c(seq(0,0,length=1),seq(0.1,3,length=14),seq(3.1,max(d1),length=5))
			
		hc1<-heatmap.2(d1,Rowv=rowDend, Colv=colDend,col=my_palette,breaks=col_breaks,key=FALSE,trace=c("none"),density.info=c("none"),
		cexRow=0.2,cexCol=0.2,RowSideColors= myClusterSideBar)

		pdf(file ="Figure 8C_ortho_spec_df1_heatmap.pdf", family="Times", height=10, width=10)
			par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=0)
			heatmap.2(d1,Rowv=rowDend, Colv=colDend,col=my_palette,breaks=col_breaks,key=FALSE,trace=c("none"),density.info=c("none"),
			cexRow=0.2,cexCol=0.2,RowSideColors= myClusterSideBar)
		dev.off()


###########################################################################	heatmap_Cluster 2, Figure 8D	

		d2<-d[rownames(d)%in%rownames(cluster2),]		
		d2<-t(d2)			
		dim(d2)
			#[1] 111 195

		rowDistance=dist(d2,method="manhattan")
		rowCluster = hclust(rowDistance,method="ward.D2")
		rowDend = as.dendrogram(rowCluster)
		rowDend = reorder(rowDend, rowSums(d2))
		dnd_cluster2<-rowDend

		colDistance=dist(t(d2),method="manhattan")
		colCluster = hclust(colDistance,method="ward.D2")
		colDend = as.dendrogram(colCluster)
		colDend = reorder(colDend, colSums(d2))

	
		my_palette<-colorRampPalette(c("light blue","yellow","red"))(n=19)
		col_breaks=c(seq(0,0,length=1),seq(0.1,3,length=14),seq(3.1,max(d2),length=5))
			
		h2=heatmap.2(d2,Rowv=rowDend, Colv=colDend,col=my_palette,breaks=col_breaks,keysize=2,margins=c(5,5),trace=c("none"),
			density.info=c("none"),cexRow=0.2,cexCol=0.2,RowSideColors= myClusterSideBar)
				
		pdf(file ="Figure 8D_ortho_spec_df2_heatmap.pdf", family="Times", height=10, width=10)
			par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=0)
			heatmap.2(d2,Rowv=rowDend, Colv=colDend,col=my_palette,breaks=col_breaks,key=FALSE,trace=c("none"),density.info=c("none"),
			cexRow=0.2,cexCol=0.2,RowSideColors= myClusterSideBar)
		dev.off()

		

############################################################################ heatmap_Cluster 3, Figure 8E	
					
		d3<-d[rownames(d)%in%rownames(cluster3),]		
		d3<-t(d3)		
		dim(d3)
				#[1] 111  45	

		rowDistance=dist(d3,method="manhattan")
		rowCluster = hclust(rowDistance,method="ward.D2")
		rowDend = as.dendrogram(rowCluster)
		rowDend = reorder(rowDend, rowSums(d3))
		dnd_cluster3<-rowDend
		

		colDistance=dist(t(d3),method="manhattan")
		colCluster = hclust(colDistance,method="ward.D2")
		colDend = as.dendrogram(colCluster)
		colDend = reorder(colDend, colSums(d3))


		my_palette<-colorRampPalette(c("light blue","yellow","red"))(n=19)
		col_breaks=c(seq(0,0,length=1),seq(0.1,3,length=14),seq(3.1,max(d3),length=5))

		h3<-heatmap.2(d3,Rowv=rowDend, Colv=colDend,col=my_palette,breaks=col_breaks,keysize=2,margins=c(5,5),
		trace=c("none"),density.info=c("none"),cexRow=0.2,cexCol=0.2,RowSideColors= myClusterSideBar)
			
		pdf(file ="Figure 8E_ortho_spec_df3_heatmap.pdf", family="Times", height=10, width=10)
			par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=0)
			heatmap.2(d3,Rowv=rowDend, Colv=colDend,col=my_palette,breaks=col_breaks,key=FALSE,trace=c("none"),density.info=c("none"),
			cexRow=0.2,cexCol=0.2,RowSideColors= myClusterSideBar)
		dev.off()


############################################################################ heatmap_Cluster 4, Figure 8F	

		d4<-d[rownames(d)%in%rownames(cluster4),]		
		d4<-t(d4)		
	    dim(d4)
	  	 	#[1] 111  95	

		rowDistance=dist(d4,method="manhattan")
		rowCluster = hclust(rowDistance,method="ward.D2")
		rowDend = as.dendrogram(rowCluster)
		rowDend = reorder(rowDend, rowSums(d4))
		dnd_cluster4<-rowDend	

		colDistance=dist(t(d4),method="manhattan")
		colCluster = hclust(colDistance,method="ward.D2")
		colDend = as.dendrogram(colCluster)
		colDend = reorder(colDend, colSums(d4))


		my_palette<-colorRampPalette(c("light blue","yellow","red"))(n=19)
		col_breaks=c(seq(0,0,length=1),seq(0.1,3,length=14),seq(3.1,max(d4),length=5))

		h4<-heatmap.2(d4,Rowv=rowDend, Colv=colDend,col=my_palette,breaks=col_breaks,keysize=2,margins=c(5,5),trace=c("none"),
		density.info=c("none"),cexRow=0.2,cexCol=0.2,RowSideColors= myClusterSideBar)
			
		pdf(file ="Figure 8F_ortho_spec_df4_heatmap.pdf", family="Times", height=10, width=10)
			par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=0)
			heatmap.2(d,Rowv=rowDend, Colv=colDend,col=my_palette,breaks=col_breaks,key=FALSE,trace=c("none"),density.info=c("none"),
			cexRow=0.2,cexCol=0.2,RowSideColors= myClusterSideBar)
		dev.off()

############################################################################ Modeling SCPs in Cluster 4, Figure 8G

		library("fitdistrplus")
		library("KSgeneral")
	
		hc<-as.matrix(rownames(h4$carpet))
		
		hc1<-hc[1:22,] #Multi-copy subfamilies
		hc2<-hc[23:95,] # Intermediate subfamlies
	
		hc1_count<-d4[,colnames(d4)%in%hc1]		 	
		for(k in 1:ncol(hc1_count)){ 	
			hc1_count[,k][hc1_count[,k]>1] <- 0 		
			} 
	
		scp_rate<-colSums(hc1_count)/111*100
		rate_c1<-as.data.frame(scp_rate)		
		group<-rep("v1",dim(rate_c1)[1])		
		rate_c1<-cbind(rate_c1,group)
	

		hc2_count<-d4[,colnames(d4)%in%hc2]		
	 	for(k in 1:ncol(hc2_count)){ 	
			hc2_count[,k][hc2_count[,k]>1] <- 0 		
			} 
	
		scp_rate<-colSums(hc2_count)/111*100
		rate_c2<-as.data.frame(scp_rate)		
		group<-rep("v2",dim(rate_c2)[1])		
		rate_c2<-cbind(rate_c2,group)
	
		rate_c1c2<-rbind(rate_c1,rate_c2)	
		group<-rep("v3",dim(rate_c1c2)[1])	
		rate_c1c2[,2]<-group
	
	
		rate<-rbind(rate_c1,rate_c2,rate_c1c2)
	
######################
		
		x<-rate_c1c2$scp_rate

		fw <- fitdist(x, "weibull")
		fg <- fitdist(x, "gamma")
		fln <- fitdist(x, "lnorm")
		fn <- fitdist(x, "norm")
		fo<-fitdist(x, "logis")
		fu<-fitdist(x, "unif")

		gofstat(list(fw, fg,fln,fn,fo,fu))

				# Goodness-of-fit statistics
				# Goodness-of-fit statistics
				#                              1-mle-weibull 2-mle-gamma 3-mle-lnorm 4-mle-norm
				# Kolmogorov-Smirnov statistic     0.1545166   0.2014643    0.219318  0.1254945
				# Cramer-von Mises statistic       0.5203379   0.9412989    1.396080  0.3400282
				# Anderson-Darling statistic       3.3436908   5.3036210    7.724317  2.1063107
				#                              5-mle-logis 6-mle-unif
				# Kolmogorov-Smirnov statistic  0.09585994  0.2640565
				# Cramer-von Mises statistic    0.22907469  2.2537483
				# Anderson-Darling statistic    1.83445033        Inf

				# Goodness-of-fit criteria
				#                                1-mle-weibull 2-mle-gamma 3-mle-lnorm 4-mle-norm
				# Akaike's Information Criterion      848.1928    874.7934    915.3326   836.0093
				# Bayesian Information Criterion      853.3006    879.9012    920.4404   841.1170
				#                                5-mle-logis 6-mle-unif
				# Akaike's Information Criterion    839.7339         NA
				# Bayesian Information Criterion    844.8416         NA


		pdf("Figure S11_cluster4_scp_fw_fn_fo_fit.pdf",width=10,height=10)
			par(mfrow = c(2, 2))
			plot.legend <- c("Weibull", "normal","logis")
			denscomp(list(fw, fn,fo), legendtext = plot.legend)
			qqcomp(list(fw, fn,fo), legendtext = plot.legend)
			cdfcomp(list(fw, fn,fo), legendtext = plot.legend)
			ppcomp(list(fw, fn,fo), legendtext = plot.legend)
		dev.off()
						
		fo
						#
						#Fitting of the distribution ' logis ' by maximum likelihood 
						#Parameters:
						#         estimate Std. Error
						#location 49.65587  2.0249386
						#scale    11.19333  0.9492182
					
		disc_ks_test(unique(x), "plogis", exact = TRUE,		location =49.65587  ,scale  = 11.19333, 
						tol = 1e-08, sim.size = 1e+06, num.sim = 10) #	D = 0.13126, p-value = 0.3261

###########################################################################Figure 8G
			
		library("scales")
	
		args=list(location =49.65587  ,scale  = 11.19333)
		pdf("Figure 8G_Cluster 4_Subgroup_single copy rate_density_w_dlogis.pdf",width=5,height=5)			
				ggplot(rate_c1c2, aes(x=scp_rate,fill=group)) +
				geom_density(alpha=0.25,aes(y = 5*..count../sum(..count..)))+
			    scale_y_continuous(labels = percent, name = "percent",limits=c(0,0.03) )+
			    theme_classic()+
				xlim(-10,100)+	
				stat_function(fun = dlogis, args = args, color = "red")							
		dev.off()


		pdf("Figure 8G_Cluster 4_Subgroup_single copy rate_density.pdf",width=5,height=5)
				ggplot(rate, aes(x=scp_rate,fill=group)) +
				geom_density(alpha=0.25,aes(y = 10*..count../sum(..count..)))+
			    scale_y_continuous(labels = percent, name = "percent",limits=c(0,0.03) )+
			    theme_classic()+
				xlim(-10,100)				
		dev.off()
			
############################################################################ Modeling SCPs in Cluster 3, Figure S12
		
		hc<-as.matrix(rownames(h3$carpet))					
		hc3<-hc
	
		hc3_count<-d3[,colnames(d3)%in%hc3]		
		 for(k in 1:ncol(hc3_count)){ 	
			hc3_count[,k][hc3_count[,k]>1] <- 0 		
				} 
	
		scp_rate<-colSums(hc3_count)/111*100
		rate_c3<-as.data.frame(scp_rate)		
		group<-rep("v1",dim(rate_c3)[1])		
		rate_c3<-cbind(rate_c3,group)	
			
	
		x<-rate_c3$scp_rate

		fw <- fitdist(x, "weibull")
		fg <- fitdist(x, "gamma")
		fln <- fitdist(x, "lnorm")
		fn <- fitdist(x, "norm")
		fo<-fitdist(x, "logis")
		fu<-fitdist(x, "unif")

		gofstat(list(fw, fg,fln,fn,fo,fu))

				# Goodness-of-fit statistics
					#
				#	Goodness-of-fit statistics
				#	                             1-mle-weibull 2-mle-gamma 3-mle-lnorm 4-mle-norm
				#	Kolmogorov-Smirnov statistic    0.06978292  0.10094033   0.1263824 0.07171368
				#	Cramer-von Mises statistic      0.03440195  0.06304752   0.1196631 0.03310618
				#	Anderson-Darling statistic      0.30213753  0.57201365   0.9657831 0.28368246
				#	                             5-mle-logis 6-mle-unif
				#	Kolmogorov-Smirnov statistic  0.07595956  0.1400000
				#	Cramer-von Mises statistic    0.03359728  0.1632889
				#	Anderson-Darling statistic    0.29003425        Inf
				#
				#	Goodness-of-fit criteria
				#	                               1-mle-weibull 2-mle-gamma 3-mle-lnorm 4-mle-norm
				#	Akaike's Information Criterion      350.8068    354.9814    359.5544   351.5569
				#	Bayesian Information Criterion      354.4201    358.5948    363.1677   355.1702
				#	                               5-mle-logis 6-mle-unif
				#	Akaike's Information Criterion    353.5330         NA
				#	Bayesian Information Criterion    357.1464         NA
				
		fn
						#
						#
					#	Fitting of the distribution ' norm ' by maximum likelihood 
					#	Parameters:
					#	     estimate Std. Error
					#	mean 34.69469   1.715103
					#	sd   11.50526   1.212761

		disc_ks_test(unique(x), "pnorm", exact = TRUE,		mean =34.69469  ,sd  = 11.50526,
			 			tol = 1e-08, sim.size = 1e+06, num.sim = 10) #	D = 0.10701, p-value = 0.7922
		
		
##################################################################### Figure S12
		
		args=list(mean =34.69469  ,sd  = 11.50526)

		pdf("Figure S12_Cluster 3_Subgroup_single copy rate_density_w_dlogis.pdf",width=5,height=5)
				ggplot(rate_c3, aes(x=scp_rate,fill=group)) +
				geom_density(alpha=0.25,aes(y = 5*..count../sum(..count..)))+
			    scale_y_continuous(labels = percent, name = "percent",limits=c(0,0.05) )+
			    theme_classic()+
				xlim(-10,100)+	
				stat_function(fun = dnorm, args = args, color = "red")						
		dev.off()

	
############################################################################ Dngrogram correlation, Figure S10
		
		library("dendextend")

		phylo<-read.tree("Data_S1_A_bifurcated_newick_species_tree_of_111_plant_genomes.txt")
	
		tree_matrix<-compute.mr(phylo, type = "matrix")
		species_short<-substring(rownames(tree_matrix),1,5)
		species_short<-gsub("\\.","",species_short)
		rownames(tree_matrix)<-species_short
		tree_matrix_dist<-as.matrix(dist(tree_matrix))
			
		tree_rowDistance=dist(tree_matrix_dist,method="manhattan")
		tree_hc = hclust(tree_rowDistance,method="ward.D2")		
		
		rowDistance=dist(tree_matrix_dist,method="manhattan")
		rowCluster = hclust(rowDistance,method="ward.D2")
		rowDend = as.dendrogram(rowCluster)
		rowDend = reorder(rowDend, rowSums(d))
		dnd_tree<-rowDend

		
		group_dendlist<-dendlist(dnd_tree,dnd_cluster1,dnd_cluster2,dnd_cluster3,dnd_cluster4)
	
		names(group_dendlist) <- c("tree","c1","c2","c3","c4")

		group_dendlist_cor_spearman <- cor.dendlist(group_dendlist, method_coef = "spearman")
	
		group_dendlist_cor_spearman
				#			#
				#           tree         c1         c2        c3         c4
				#  tree 1.0000000 0.18019523 0.34251257 0.5649266 0.10015909
				#  c1   0.1801952 1.00000000 0.27453994 0.1798856 0.08987023
				#  c2   0.3425126 0.27453994 1.00000000 0.3051712 0.08623967
				#  c3   0.5649266 0.17988563 0.30517120 1.0000000 0.16411817
				#  c4   0.1001591 0.08987023 0.08623967 0.1641182 1.00000000
	
		pdf("Figure S10_cluster_dendlist_cor_spearman_cophenetic_correlation.pdf",height=5,width=5)
			corrplot::corrplot(group_dendlist_cor_spearman, "pie", "lower")
		dev.off()





