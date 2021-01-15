#libraries

library("stringr")

	
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
	
	fbx_cluster1_4_size<-c(dim(cluster1)[1],dim(cluster4)[1])
	names(fbx_cluster1_4_size)<-c("Cluster 1","Cluster 4")

	angiosperm_lineage_specific_core_size<-as.numeric(c("69133","9178")) #see Li et al. (Plant Cell 2016(28):326-344) for the values
	names(angiosperm_lineage_specific_core_size)<-c("Lineage Specific","Core")
	
	pdf("Figure_S4_FBX_angiosperm_families_lineage_specificity_comparison.pdf",height=5,width=5);
			#
			par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=0)
			layout(matrix(c(1,2),nrow=1,ncol=2))
			barplot(fbx_cluster1_4_size,names.arg=names(fbx_cluster1_4_size),ylim=c(0,6000),ylab="Number of FBX Subfamilies",
						lwd = 0.3,cex.axis=0.5, cex.names=0.5,beside=TRUE)
			#
			barplot(angiosperm_lineage_specific_core_size,names.arg=names(angiosperm_lineage_specific_core_size),ylim=c(0,70000),
						ylab="Number of Angiosperm Gene Families",lwd = 0.3,cex.axis=0.5, cex.names=0.5,beside=TRUE)
	dev.off()

#########################################

 		#
		ortho<-read.delim("Data_S3_A_complete_list_of_5,858_orthoMCL_groups_of_FBX_genes.txt",header=F)
		library("stringr")

#
######################## Function for retrieving Arabidopsis FBX IDs in each subfamily ####################################
	
		ath_ortho_genes<-function(conserved_orthos){

			conserved_ath_ortho_genes<-c()

			for(i in 1:dim(ortho)[1]){
				b<-conserved_orthos[str_detect(ortho[i,],conserved_orthos)]	
				if(!identical(b, character(0)) && str_detect(ortho[i,],"Ath")){	
					conserved_ath_ortho_genes<-rbind(conserved_ath_ortho_genes,as.character(ortho[i,]))
									}						
								}		
			gene_array<-c()

			for(i in 1:dim(conserved_ath_ortho_genes)[1]){
					group<-conserved_ath_ortho_genes[i,]
					group<-unlist(strsplit(group," "))
					for(j in 1:length(group)){
						member<-gsub(".*\\|","",group[j])		
						if(str_detect(member,"Ath")){	
							gene<-cbind(group[1],member)
							gene_array<-rbind(gene_array,gene)		
									}
								}
						}	
			 gene_array	
			}


######################## RetrievE Arabidopsis FBX IDs in each cluster ####################################
		#
		ortho1<-rownames(cluster1)
		ortho2<-rownames(cluster2)
		ortho3<-rownames(cluster3)
		ortho4<-rownames(cluster4)
		
		ortho1_genes<-ath_ortho_genes(ortho1)
		ortho2_genes<-ath_ortho_genes(ortho2)
		ortho3_genes<-ath_ortho_genes(ortho3)
		ortho4_genes<-ath_ortho_genes(ortho4)
		
			
   		known_ids<-read.table("Data_S10_List_of_known_Arabidopsis_FBX_genes.txt", header=TRUE)
		
   		known_ortho1_genes<-ortho1_genes[ortho1_genes[,2]%in%known_ids[,3], ]
   		known_ortho2_genes<-ortho2_genes[ortho2_genes[,2]%in%known_ids[,3], ]
   		known_ortho3_genes<-ortho3_genes[ortho3_genes[,2]%in%known_ids[,3], ]
   		known_ortho4_genes<-ortho4_genes[ortho4_genes[,2]%in%known_ids[,3], ]
		
		
		ortho1_known_rate<-dim(known_ortho1_genes)[1]/dim(ortho1_genes)[1]*100
		ortho2_known_rate<-dim(known_ortho2_genes)[1]/dim(ortho2_genes)[1]*100
		ortho3_known_rate<-dim(known_ortho3_genes)[1]/dim(ortho3_genes)[1]*100
		ortho4_known_rate<-dim(known_ortho4_genes)[1]/dim(ortho4_genes)[1]*100
		
		
		ortho1_activities<-cbind(dim(known_ortho1_genes)[1],dim(ortho1_genes)[1])
		ortho2_activities<-cbind(dim(known_ortho2_genes)[1],dim(ortho2_genes)[1])
		ortho3_activities<-cbind(dim(known_ortho3_genes)[1],dim(ortho3_genes)[1])
		ortho4_activities<-cbind(dim(known_ortho4_genes)[1],dim(ortho4_genes)[1])
		   
		ath_fbx_activities<-rbind(ortho1_activities,ortho2_activities,ortho3_activities,ortho4_activities)
		   
		ath_fbx_unknown_known<-cbind((ath_fbx_activities[,2]-ath_fbx_activities[,1]),ath_fbx_activities[,1])
		rownames(ath_fbx_unknown_known)<-c("Cluster1","Cluster2","Cluster3","Cluster4")
		colnames(ath_fbx_unknown_known)<-c("unknown","known")
		
				
		barplot(t(ath_fbx_unknown_known),col=c("gray","white"),legend = colnames(ath_fbx_unknown_known))
	
		t(ath_fbx_unknown_known)
					
					#	        Cluster1 Cluster2 Cluster3 Cluster4
					#	unknown      270      247       34       78
					#	known         10       11        8       52
		
		pdf("Figure S14_ath_fbx_unknown_known_4_clusters.pdf",height=5,width=5)				
			barplot(t(ath_fbx_unknown_known),col=c("gray","white"),legend = colnames(ath_fbx_unknown_known))			
		dev.off()
		






