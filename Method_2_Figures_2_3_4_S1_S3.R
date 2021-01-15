#libraries

library("phytools")
library("ConsensusClusterPlus")
library("ggplot2")
library("ggbiplot")
library("fitdistrplus")
library("KSgeneral")



############################################################################
#  Order of analysis					   						           #
# 					   					 				                   #
#   Figure S1															   #
#   Figure 2, 3B-C														   #
#   Figure 4, 3A, S3													   #
#   Figure 6, S5														   #
#   Figure 5, S4														   #
#   Figure S6															   #
#   Figure 8 															   #
# 					   					 				                   #
############################################################################	


############################################################################
# 					   						        			           #
#   Figure S1															   #
# 					   					 				                   #
############################################################################	


############################################################################ Figure S1

	phytozome.tree<-read.tree("Data_S1_A_bifurcated_newick_species_tree_of_111_plant_genomes.txt")

	pdf("Figure S1_phytozome13_tree_single_accession.pdf",height=9.2,width=6.8)
		plotTree(phytozome.tree,ftype="i",fsize=0.5,lwd=1)
	dev.off()

	# convert species phylogeny to a data matrix
	tree_matrix<-compute.mr(phytozome.tree, type = "matrix")
	species_short<-substring(rownames(tree_matrix),1,5)
	species_short<-gsub("\\.","",species_short)
	rownames(tree_matrix)<-species_short
	tree_matrix_dist<-as.matrix(dist(tree_matrix))
	# rownames(tree_matrix_dist)

	id<-rownames(tree_matrix_dist)
	list<-seq(1,length(id),1) 

	species<-cbind(list,id) # order the species names in the "species" matrix
	                        # to allow the row arrangement of ortho_count matrix as follows

#############


## Count paralogues of each species in each OrthoMCL group

	ortho<-read.delim("Data_S3_A_complete_list_of_5,858_orthoMCL_groups_of_FBX_genes.txt",header=F)
	ortho_count<-species #

	for(i in 1:dim(ortho)[1]){

		group<-ortho[i,]
		group<-as.matrix(group)
		group<-unlist(strsplit(group," "))
	
		members<-matrix(group, ncol = length(group) , byrow = TRUE )

		group_array<-c()

		for(j in 2:dim(members)[2]){
			member<-gsub("\\|.*","",members[,j]) #retrieve the speciese name of each FBX id	
			group_array<-c(group_array,member)		
		}

		count<-as.matrix(table(group_array)) # Count the number of paralogues in each species
		id<-rownames(count) # the species names are ordered alphabetically
		count<-cbind(count,id)
		colnames(count)<-c(gsub(":","",group[1]),"id") # OrthoMCL group ID
	
		ortho_count<-merge(ortho_count, count, by="id", all=TRUE) # merge counts for each OrthoMCL group according 
																  # to "id" in "species" matrix.  "NA" is given to a species if it does not
																  # have members in the OrthoMCL group i.
				}

	ortho_count<-ortho_count[order(ortho_count$list),]

	write.table(ortho_count,"fbx_ortho_count_matrix.tab") # Counting takes time. Save the restult in a data matrix file for following analysis


####################### retrieve orthomcl cnt f-box genes  ####################### 


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

############################################################################
# Figure 2 & 3B-C													       #
#																		   #
# For each OrhoMCL Group i, calculate                                      #
#  			1) phylogenetic distance (dist), 					           #
#			2) numner of species having loss of members (gap),             #
#			3) number of species having members (cnt),                     #
#			4) maximum number of members within a species (max),           #
#			5) mean number of members of Group i (mean)                    #
############################################################################	
	
	ortho_spec<-c()

	for(i in 1:dim(ortho_count)[2]){
	
		ortho_i_spec<-c()
	
		ortho_count_i<-ortho_count[names(ortho_count)[i]]
		species<-rownames(ortho_count_i)[ortho_count_i[,1]>0]
	
		cnt<-length(species) # number of species are present in the ith OrhoMCL group

		max<-max(ortho_count_i[,1])
		mean<-mean(ortho_count_i[,1])
	
		dist<-c()
		gap<-c()
	
		if(cnt==1){	        # species specific duplications
				dist<-0
				gap<-0
				ortho_i_spec<-cbind(colnames(ortho_count_i),dist,gap,cnt,max,mean)	
			}
	    else{
			ortho_tree_i<-tree_matrix_dist[rownames(tree_matrix_dist)%in%species,]	
			ortho_tree_i<-ortho_tree_i[,colnames(tree_matrix_dist)%in%species]	
			dist<-sum(ortho_tree_i) #Use the sum of species phylogenetic distance to reflect the evolutionary conservation of an OrthoMCL group			
		
			# Find the coverage of species present OrthoMCL i in species tree matrix
			s<-rownames(ortho_tree_i)[1]
			e<-rownames(ortho_tree_i)[length(rownames(ortho_tree_i))]
			species_tree_i<-tree_matrix_dist[which(row.names(tree_matrix_dist) == s) : which(row.names(tree_matrix_dist) == e), , drop = FALSE]
			# rows containing species with missing members

			gap<-dim(species_tree_i)[1]-dim(ortho_tree_i)[1]			
			ortho_i_spec<-cbind(colnames(ortho_count_i),dist,gap,cnt,max,mean)		
		}

		ortho_spec<-rbind(ortho_spec,ortho_i_spec)					
	
	}

	
# convert ortho_spec to numeric data.matrix and data.frame
	rownames(ortho_spec)<-ortho_spec[,1]
	ortho_spec<-ortho_spec[,-1]

	ortho_spec_num<-apply(ortho_spec,2,as.numeric)
	rownames(ortho_spec_num)<-rownames(ortho_spec)
	
	ortho_spec<-ortho_spec_num
	
	ortho_spec_df<-as.data.frame(ortho_spec)
	
	write.csv(ortho_spec_df,"Data_matrix_of_multidimensional_features_of_FBX_subfamilies.tab")

############################################################################ Figure 2

	ortho_spec_df<-read.table("Data_S5_Data_matrix_of_multidimensional_features_of_FBX_subfamilies.txt",header=T)
	dist<-ortho_spec_df$dist
	gap<-ortho_spec_df$gap
	cnt<-ortho_spec_df$cnt

	h_dist<-hist(dist, breaks = seq(min(dist), max(dist), length.out = 108))
	h_gap<-hist(gap, breaks = seq(min(gap), max(gap), length.out = 108))
	h_cnt<-hist(cnt, breaks = seq(min(cnt), max(cnt), length.out = 108))

	dist_names<-c(h_dist$breaks)
	dist_names<-round(dist_names[-108], digits = 0)

	gap_names<-c(h_gap$breaks)
	gap_names<-round(gap_names[-108],digits = 0)

	cnt_names<-c(h_cnt$breaks)
	cnt_names<-round(cnt_names[-108],digits=0)
	


	pdf("Figure 2_fbx_ortho_cnt_dist_gap_hist.pdf",width=8,height=8) # Figure 2
			par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=0)
			layout(matrix(c(1,2,3,4),nrow=2,ncol=2))
			barplot(log(h_cnt$counts,2),names.arg=cnt_names,ylim=c(0,10),ylab="log(# of Subfamilies,2)",xlab="# of species",lwd = 0.3,cex.axis=0.5, cex.names=0.5,beside=TRUE,)
			barplot(log(h_dist$counts,2),names.arg=dist_names,ylim=c(0,10),ylab="log(# of Subfamilies,2)",xlab="dist",lwd = 0.3,cex.axis=0.5, cex.names=0.5,beside=TRUE,)
			plot(log(ortho_spec_df$cnt,2),log(ortho_spec_df$dist,2),xlab="log(# of Subfamilies,2)",ylab="log(dist)2",pch=1,cex=0.1,col="dark gray")
			text(2,14,"rho=0.93, p-value< 2.2e-16",cex=0.7)	
			barplot(log(h_gap$counts,2),names.arg=gap_names,ylim=c(0,10),ylab="log(# of Subfamilies,2)",xlab="gap",lwd = 0.3,cex.axis=0.5, cex.names=0.5,beside=TRUE,)
	dev.off()
	

	
############################################################################ ConsensusClusterPlus km clustering, Figure 3A
	
	d<-ortho_spec
	d<-t(d)

	title<-getwd()

	results = ConsensusClusterPlus(d,maxK=9,reps=1000,pItem=0.8,pFeature=1,title=title,innerLinkage="average",finalLinkage="average",clusterAlg="km",distance="euclidean",seed=1262118388.71279,plot="png")

	clusters<-as.matrix(results[[4]][["consensusClass"]])

	clusters<-cbind(clusters,paste("c",clusters[,1],sep="_"))

	write.table(clusters,"Data_matrix_of_four_clusters_of_FBX_subfamilies.tab") #ConsensusClusterPlus takes time.  Save the result for following analyses 


############################################################################ PCA, Figure 3B

	clusters<-read.table("Data_S6_Data_matrix_of_four_clusters_of_FBX_subfamilies.txt",header=T)

	table(clusters[,2])
			
	#	c_1  c_2  c_3  c_4 
	#	 95   45 5523  195
	
	# rearrange the order of clusters	
	clusters_arranged <- apply(clusters, 2, function(x) as.character(gsub("1", "four", x)))
	clusters_arranged<-apply(clusters_arranged, 2, function(x) as.character(gsub("2", "three", x)))
	clusters_arranged<-apply(clusters_arranged, 2, function(x) as.character(gsub("3", "1", x)))
	clusters_arranged<-apply(clusters_arranged, 2, function(x) as.character(gsub("4", "2", x)))
	clusters_arranged<-apply(clusters_arranged, 2, function(x) as.character(gsub("three", "3", x)))
	clusters_arranged<-apply(clusters_arranged, 2, function(x) as.character(gsub("four", "4", x)))
	
	table(clusters_arranged[,2])

		  #  c_1  c_2  c_3  c_4 
		  # 5523  195   45   95 
		
	ortho_spec.pca <- prcomp(ortho_spec, center = TRUE,scale. = TRUE)
	summary(ortho_spec.pca)

		#	                        PC1    PC2    PC3     PC4     PC5
		#	Standard deviation     1.6781 1.0374 0.9743 0.34186 0.20418
		#	Proportion of Variance 0.5632 0.2152 0.1898 0.02337 0.00834
		#	Cumulative Proportion  0.5632 0.7784 0.9683 0.99166 1.00000


	pdf("Figure 3B_ortho_spec_four_clusters.pca2.pdf",width=5,height=5)
		ggbiplot(ortho_spec.pca,ellipse=TRUE,obs.scale = 2, var.scale = 1, var.axes=FALSE,  labels=clusters_arranged[,1],groups= clusters_arranged[,2]) +
 		scale_colour_manual(name="Origin", values= c("purple", "forest green", "blue","magenta"))+
  		ggtitle("PCA of ortho_spec dataset")+
  		theme_minimal()+
  		theme(legend.position = "bottom")
	dev.off()



################################################################################################################################################ modeling


############################################################################
# Figure 3C, 4 & S3A												       #
#																		   #
# Modeling the conservation of subfamilies        				           #
#  			1) All set, Figure 3C					                       #
#			2) Cluster 1, Figure 4A, E                                     #
#			3) Cluster 2, Figure 4B, F                                     #
#			4) Cluster 3, Figure 4C, G                                     #
#			5) Cluster 4, Figure 4D, H                                     #
#			5) Cluster 1+2, Figure S3A                                     #
#			5) Cluster 3+4, Figure S3B                                     #
############################################################################	

	#functions

	getmode <- function(v) {
		   uniqv <- unique(v)
		   uniqv[which.max(tabulate(match(v, uniqv)))]
		}


# Clusters are re ordered based on the mean number of species

	cluster1<-clusters[clusters[,1]==3,]
	cluster2<-clusters[clusters[,1]==4,]
	cluster3<-clusters[clusters[,1]==2,]
	cluster4<-clusters[clusters[,1]==1,]

	ortho_spec_df1<-ortho_spec_df[rownames(ortho_spec_df)%in%rownames(cluster1),]
	ortho_spec_df2<-ortho_spec_df[rownames(ortho_spec_df)%in%rownames(cluster2),]
	ortho_spec_df3<-ortho_spec_df[rownames(ortho_spec_df)%in%rownames(cluster3),]
	ortho_spec_df4<-ortho_spec_df[rownames(ortho_spec_df)%in%rownames(cluster4),]



################################################################################  All set


	dat<-ortho_spec_df
	x<-dat$cnt
	
	fw <- fitdist(x, "weibull")
	fg <- fitdist(x, "gamma")
	fln <- fitdist(x, "lnorm")
	fe <- fitdist(x, "exp")

	gofstat(list(fw, fg,fln,fe))

		#	Goodness-of-fit statistics
		#	                             1-mle-weibull 2-mle-gamma 3-mle-lnorm   4-mle-exp
		#	Kolmogorov-Smirnov statistic     0.2416677   0.2407896   0.2003158   0.2794098
		#	Cramer-von Mises statistic      75.2643999  99.5508660  38.2696753 127.5468251
		#	Anderson-Darling statistic     428.7775567 521.0768497 225.2395296 633.6501583

		#	Goodness-of-fit criteria
		#	                               1-mle-weibull 2-mle-gamma 3-mle-lnorm 4-mle-exp
		#	Akaike's Information Criterion      31905.58    32550.28    28918.33  32728.24
		#	Bayesian Information Criterion      31918.93    32563.64    28931.68  32734.92


	fln  # choose the best fit to further examine its fitness

		# Fitting of the distribution ' lnorm ' by maximum likelihood 
		# Parameters:
		#         estimate  Std. Error
		# meanlog 1.0640034 0.012870869
		# sdlog   0.9851051 0.009101037

	b1 <- bootdist(fln, niter= 500)

	pdf("Figure 3C_Inset_all_fln_CIcdfplot.pdf",width=5,height=5)
		CIcdfplot(b1, CI.level= 95/100, CI.output = "probability",CI.fill = "light gray", CI.col = "red")
	dev.off()

	disc_ks_test(unique(x), "plnorm", exact = TRUE,meanlog = 1.0640034, sdlog =  0.9851051,
				 tol = 1e-08, sim.size = 1e+06, num.sim = 10) #	D = 0.80452, p-value < 2.2e-16

				
################## plot the distribution, Figure 3C

	args=list(meanlog = 1.0640034, sdlog =  0.9851051)

	pdf("Figure 3C_ggplot_cnt_lnormal_distribution_all.pdf",width=5,height=5)

		ggplot(dat, aes(x=cnt)) + 
				 geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
				binwidth=1,
				colour="black", fill="white") +
				geom_density(alpha=.2, fill="#FF6666")+  # Overlay with transparent density plot
				geom_vline(aes(xintercept=mean(cnt, na.rm=T)),   # Ignore NA values for mean
				               color="red", linetype="dashed", size=1) +
				geom_vline(aes(xintercept=getmode(cnt)),   # Ignore NA values for mean
							               color="green", linetype="dashed", size=1) +
				
				stat_function(fun = dlnorm, args = args, color = "red") 
							
	dev.off()


################################################################################ cluster 1

	dat<-ortho_spec_df1
	x<-dat$cnt

	fw <- fitdist(x, "weibull")
	fg <- fitdist(x, "gamma")
	fln <- fitdist(x, "lnorm")
	fn <- fitdist(x, "norm")

	gofstat(list(fw, fg,fln,fn))

		# Goodness-of-fit statistics
		#                              1-mle-weibull 2-mle-gamma 3-mle-lnorm 4-mle-norm
		# Kolmogorov-Smirnov statistic      0.210407   0.2327179   0.1984569   0.261106
		# Cramer-von Mises statistic       45.293135  45.5003887  30.2192643 101.159207
		# Anderson-Darling statistic      263.063818 258.5647895 185.1705589 549.077684

		# Goodness-of-fit criteria
		#                                1-mle-weibull 2-mle-gamma 3-mle-lnorm 4-mle-norm
		# Akaike's Information Criterion      24100.33    23681.32    22494.05   28642.25
		# Bayesian Information Criterion      24113.57    23694.55    22507.29   28655.48

	fln 

		# Fitting of the distribution ' lnorm ' by maximum likelihood 
		# Parameters:
		#          estimate  Std. Error
		# meanlog 0.9050703 0.010088984
		# sdlog   0.7497819 0.007133932

	b1 <- bootdist(fln, niter= 500)

	pdf("Figure 4E_cluster1_fln_CIcdfplot.pdf",width=5,height=5)
		CIcdfplot(b1, CI.level= 95/100, CI.output = "probability",CI.fill = "light gray", CI.col = "red")
	dev.off()

	disc_ks_test(unique(x), "plnorm", exact = TRUE,meanlog = 0.9050703, sdlog = 0.7497819, 
					tol = 1e-08, sim.size = 1e+06, num.sim = 10) 	# D = 0.61836, p-value = 1.575e-07

	
################## Figure 4A
	

	args=list(meanlog = 0.9050703, sdlog = 0.7497819)

	pdf("Figure 4A_ggplot_cnt_lnormal_distribution_cluster1.pdf",width=5,height=5)

		ggplot(dat, aes(x=cnt)) + 
				 geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
				binwidth=1,
				colour="black", fill="white") +
				geom_density(alpha=.2, fill="#FF6666")+  # Overlay with transparent density plot
				geom_vline(aes(xintercept=mean(cnt, na.rm=T)),   # Ignore NA values for mean
				               color="red", linetype="dashed", size=1) +
				geom_vline(aes(xintercept=getmode(cnt)),   # Ignore NA values for mean
							               color="green", linetype="dashed", size=1) +
				
				stat_function(fun = dlnorm, args = args, color = "red") 
							
	dev.off()



################################################################################ cluster 2

	dat<-ortho_spec_df2
	x<-dat$cnt

	fw <- fitdist(x, "weibull")
	fg <- fitdist(x, "gamma")
	fln <- fitdist(x, "lnorm")

	gofstat(list(fw, fg,fln))

		# Goodness-of-fit statistics
		#                              1-mle-weibull 2-mle-gamma 3-mle-lnorm
		# Kolmogorov-Smirnov statistic     0.2041602    0.148321   0.1383733
		# Cramer-von Mises statistic       1.9631340    0.875998   0.6883358
		# Anderson-Darling statistic      11.6674348    5.383939   4.2999221

		# Goodness-of-fit criteria
		#                                1-mle-weibull 2-mle-gamma 3-mle-lnorm
		# Akaike's Information Criterion      1223.309    1145.066    1130.022
		# Bayesian Information Criterion      1229.855    1151.612    1136.568

	fln
		
		# Fitting of the distribution ' lnorm ' by maximum likelihood 
		# Parameters:
		#          estimate  Std. Error
		# meanlog 3.1627097 0.013155542
		# sdlog   0.1837071 0.009301133
		
	b1 <- bootdist(fln, niter= 500)

	pdf("Figure 4F_cluster2_fln_CIcdfplot.pdf",width=5,height=5)
			CIcdfplot(b1, CI.level= 95/100, CI.output = "probability",CI.fill = "light gray", CI.col = "red")
	dev.off()

	disc_ks_test(unique(x), "plnorm", exact = TRUE,meanlog = 3.1627097, sdlog = 0.1837071, 
					tol = 1e-08, sim.size = 1e+06, num.sim = 10) # D = 0.36729, p-value = 0.002069


################## plot the distribution, Figure 4B

	args=list(meanlog = 3.1627097, sdlog = 0.1837071)

	pdf("Figure 4B_ggplot_cnt_lnormal_distribution_cluster2.pdf",width=5,height=5)

			ggplot(dat, aes(x=cnt)) + 
					 geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
					binwidth=1,
					colour="black", fill="white") +
					geom_density(alpha=.2, fill="#FF6666")+  # Overlay with transparent density plot
					geom_vline(aes(xintercept=mean(cnt, na.rm=T)),   # Ignore NA values for mean
					               color="red", linetype="dashed", size=1) +
					geom_vline(aes(xintercept=getmode(cnt)),   # Ignore NA values for mean
								               color="green", linetype="dashed", size=1) +
				
					stat_function(fun = dlnorm, args = args, color = "red") 
							
	dev.off()


################################################################################ cluster 3

	dat<-ortho_spec_df3
	x<-dat$cnt

	fu <- fitdist(x, "unif")
	fn <- fitdist(x, "norm")

	gofstat(list(fu,fn))

			# Goodness-of-fit statistics
			#                              1-mle-unif 2-mle-norm
			# Kolmogorov-Smirnov statistic 0.07741935  0.1238325
			# Cramer-von Mises statistic   0.02319343  0.1046098
			# Anderson-Darling statistic          Inf  0.6840456

			# Goodness-of-fit criteria
			#                                1-mle-unif 2-mle-norm
			# Akaike's Information Criterion         NA   331.0568
			# Bayesian Information Criterion         NA   334.6701

	fu 

		# Fitting of the distribution ' unif ' by maximum likelihood 
		# Parameters:
		#     estimate Std. Error
		# min       47         NA
		# max       78         NA

	b1 <- bootdist(fu, niter= 500)

	pdf("Figure 4G_cluster3_fu_CIcdfplot.pdf",width=5,height=5)
			CIcdfplot(b1, CI.level= 95/100, CI.output = "probability",CI.fill = "light gray", CI.col = "red")
	dev.off()

	#######

	disc_ks_test(unique(x), "punif", exact = TRUE,min = 47, max = 78, 
						tol = 1e-08, sim.size = 1e+06, num.sim = 10) #	D = 0.1067, p-value = 0.8983

	
################## plot the distribution, Figure 4C	

	args=list(min = 47, max = 78)

	pdf("Figure 4C_ggplot_cnt_unif_distribution_cluster3.pdf",width=5,height=5)

			ggplot(dat, aes(x=cnt)) + 
					 geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
					binwidth=1,
					colour="black", fill="white") +
					geom_density(alpha=.2, fill="#FF6666")+  # Overlay with transparent density plot
					geom_vline(aes(xintercept=mean(cnt, na.rm=T)),   # Ignore NA values for mean
					               color="red", linetype="dashed", size=1) +
					geom_vline(aes(xintercept=getmode(cnt)),   # Ignore NA values for mean
								               color="green", linetype="dashed", size=1) +
				
					stat_function(fun = dunif, args = args, color = "red") 
							
	dev.off()
		

################################################################################ cluster 4

	dat<-ortho_spec_df4
	x<-dat$cnt

	fw <- fitdist(x, "weibull")
	fg <- fitdist(x, "gamma")
	fln <- fitdist(x, "lnorm")
	fn <- fitdist(x, "norm")

	gofstat(list(fw, fg,fln,fn))

		# Goodness-of-fit statistics
		#                              1-mle-weibull 2-mle-gamma 3-mle-lnorm 4-mle-norm
		# Kolmogorov-Smirnov statistic    0.07723763   0.1079200   0.1104145  0.1027122
		# Cramer-von Mises statistic      0.06544069   0.1849941   0.1998752  0.1580878
		# Anderson-Darling statistic      0.53899157   1.0379249   1.1201434  0.8930617

		# Goodness-of-fit criteria
		#                                1-mle-weibull 2-mle-gamma 3-mle-lnorm 4-mle-norm
		# Akaike's Information Criterion      613.6026    614.4677    615.2638   613.1868
		# Bayesian Information Criterion      618.7104    619.5754    620.3716   618.2945

	fw 

		# Fitting of the distribution ' weibull ' by maximum likelihood 
		# Parameters:
		#       estimate Std. Error
		# shape 17.97388  1.3999859
		# scale 97.26900  0.5866933

	b1 <- bootdist(fw, niter= 500)

	pdf("Figure 4H_cluster4_fw_CIcdfplot.pdf",width=5,height=5)
			CIcdfplot(b1, CI.level= 95/100, CI.output = "probability",CI.fill = "light gray", CI.col = "red")
	dev.off()

	disc_ks_test(unique(x), "pweibull", exact = TRUE,shape = 17.97388, scale = 97.26900, 
						tol = 1e-08, sim.size = 1e+06, num.sim = 10) #	D = 0.13522, p-value = 0.6791

		
################## plot the distribution, Figure 4D

	args=list(shape = 17.97388, scale = 97.26900)
	
	pdf("Figure 4D_ggplot_cnt_weibull_distribution_cluster4.pdf",width=5,height=5)

			ggplot(dat, aes(x=cnt)) + 
					 geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
					binwidth=1,
					colour="black", fill="white") +
					geom_density(alpha=.2, fill="#FF6666")+  # Overlay with transparent density plot
			
					stat_function(fun = dweibull, args = args, color = "red") +
					geom_vline(aes(xintercept=mean(cnt, na.rm=T)),   # Ignore NA values for mean
						               color="red", linetype="dashed", size=1) +
					geom_vline(aes(xintercept=getmode(cnt)),   # Ignore NA values for mean
									              color="green", linetype="dashed", size=1) 
																		
	dev.off()

	

################################################################################ cluster 1 & 2


	dat<-rbind(ortho_spec_df1,ortho_spec_df2)
	x<-dat$cnt

	fw <- fitdist(x, "weibull")
	fg <- fitdist(x, "gamma")
	fln <- fitdist(x, "lnorm")
	fe <- fitdist(x, "exp")

	gofstat(list(fw, fg,fln,fe))
			# Goodness-of-fit statistics
			#                              1-mle-weibull 2-mle-gamma 3-mle-lnorm   4-mle-exp
			# Kolmogorov-Smirnov statistic     0.1993045    0.231958   0.2011263   0.2166955
			# Cramer-von Mises statistic      54.0470520   59.235539  32.2567213  52.6868949
			# Anderson-Darling statistic     311.1997307  326.126920 195.4613981 309.9334417

			# Goodness-of-fit criteria
			#                                1-mle-weibull 2-mle-gamma 3-mle-lnorm 4-mle-exp
			# Akaike's Information Criterion      27533.92    27320.37    25519.93  27558.52
			# Bayesian Information Criterion      27547.23    27333.67    25533.23  27565.18

	fln

			# Fitting of the distribution ' lnorm ' by maximum likelihood 
			# Parameters:
			#          estimate  Std. Error
			# meanlog 0.9820622 0.011159159
			# sdlog   0.8438272 0.007890667
		
	disc_ks_test(unique(x), "plnorm", exact = TRUE,meanlog = 0.9820622, sdlog = 0.8438272, 
						tol = 1e-08, sim.size = 1e+06, num.sim = 10) # D = 0.72828, p-value < 2.2e-16

			
################## plot the distribution, Figure #3A

	args=list(meanlog = 0.9820622, sdlog = 0.8438272)

	pdf("Figure S3A_ggplot_cnt_lnormal_distribution_cluster1&2.pdf",width=5,height=5)

			ggplot(dat, aes(x=cnt)) + 
					 geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
					binwidth=1,
					colour="black", fill="white") +
					geom_density(alpha=.2, fill="#FF6666")+  # Overlay with transparent density plot
					geom_vline(aes(xintercept=mean(cnt, na.rm=T)),   # Ignore NA values for mean
					               color="red", linetype="dashed", size=1) +
					geom_vline(aes(xintercept=getmode(cnt)),   # Ignore NA values for mean
								               color="green", linetype="dashed", size=1) +
				
					stat_function(fun = dlnorm, args = args, color = "red") 
							
	dev.off()
	

################################################################################ cluster 3 & 4
	
	dat<-rbind(ortho_spec_df3,ortho_spec_df4)
	x<-dat$cnt

	fw <- fitdist(x, "weibull")
	fg <- fitdist(x, "gamma")
	fln <- fitdist(x, "lnorm")
	fn <- fitdist(x, "norm")

	gofstat(list(fw, fg,fln,fn))

			# Goodness-of-fit statistics
			#                              1-mle-weibull 2-mle-gamma 3-mle-lnorm 4-mle-norm
			# Kolmogorov-Smirnov statistic     0.1503313   0.2031788   0.2138306  0.1785376
			# Cramer-von Mises statistic       0.9104529   1.3990938   1.5283071  1.1491011
			# Anderson-Darling statistic       5.6023830   7.7850078   8.4744962  6.4825938

			# Goodness-of-fit criteria
			#                                1-mle-weibull 2-mle-gamma 3-mle-lnorm 4-mle-norm
			# Akaike's Information Criterion      1168.133    1204.137    1214.127   1187.998
			# Bayesian Information Criterion      1174.016    1210.020    1220.010   1193.881
	
	fw

		# Fitting of the distribution ' weibull ' by maximum likelihood 
		# Parameters:
		#        estimate Std. Error
		# shape  6.726984  0.4907648
		# scale 90.722318  1.1897902
			
	disc_ks_test(unique(x), "pweibull", exact = TRUE,shape = 6.726984, scale = 90.722318, 
					tol = 1e-08, sim.size = 1e+06, num.sim = 10) #	D = 0.19963, p-value = 0.02718

	

################## plot the distribution, Figure #3B


	args=list(shape = 6.726984, scale = 90.7223180)

	pdf("Figure S3B_ggplot_cnt_Weibull_distribution_cluster3&4.pdf",width=5,height=5)
			ggplot(dat, aes(x=cnt)) + 
					 geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
					binwidth=1,
					colour="black", fill="white") +
					geom_density(alpha=.2, fill="#FF6666")+  # Overlay with transparent density plot
			
					stat_function(fun = dweibull, args = args, color = "red") +
					geom_vline(aes(xintercept=mean(cnt, na.rm=T)),   # Ignore NA values for mean
						               color="red", linetype="dashed", size=1) +
						geom_vline(aes(xintercept=getmode(cnt)),   # Ignore NA values for mean
									               color="green", linetype="dashed", size=1) 
																		
	dev.off()









