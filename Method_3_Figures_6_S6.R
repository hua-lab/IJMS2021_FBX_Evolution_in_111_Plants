#libraries

	library("ggplot2")
	library("fitdistrplus")
	library("KSgeneral")

#functions

	getmode <- function(v) {
		   uniqv <- unique(v)
		   uniqv[which.max(tabulate(match(v, uniqv)))]
		}



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
# Figure 6 & S6 													       #
#																		   #
# Modeling the distribution of number of FBX genes per plant               #
#  			1) Complete set, including orphan genes, Figure S6A            #
#			2) All subfamily set, Figure S6B                               #
#			3) Orphan genes, Figure S6C                                    #
#			4) Cluster 1 subfamilies, Figure 6A,E                          #
#			4) Cluster 2 subfamilies, Figure 6B,F                          #
#			4) Cluster 3 subfamilies, Figure 6C,G                          #
#			4) Cluster 4 subfamilies, Figure 6D,H                          #
############################################################################	



############################################################################  Complete set including orphan genes, Figure S6A
		
	all<-read.csv("Data_S2_Number_of_sequences_processed_through_CTT_annotation_in_each_plant_genome.txt",header=T)
	all<-all[all$X%in%c("New_FBX","Prior_FBX"),]
	all<-all[,-1]
	all<-all[,rev(colnames(all))]
	size<-colSums(all)
	
	size_complete<-as.data.frame(size)
	x<-size_complete$size

	fw <- fitdist(x, "weibull")
	fg <- fitdist(x, "gamma")
	fln <- fitdist(x, "lnorm")
	fn <- fitdist(x, "norm")
	fo<-fitdist(x, "logis")

	gofstat(list(fw, fg,fln,fn,fo))

		# Goodness-of-fit statistics
		#                              1-mle-weibull 2-mle-gamma 3-mle-lnorm 4-mle-norm
		# Kolmogorov-Smirnov statistic    0.06806925  0.08377607   0.1134408  0.1079775
		# Cramer-von Mises statistic      0.06948937  0.08879015   0.2914822  0.3622769
		# Anderson-Darling statistic      0.41165010  0.49290769   1.8017713  2.4743486
		#                              5-mle-logis
		# Kolmogorov-Smirnov statistic   0.0891432
		# Cramer-von Mises statistic     0.1639784
		# Anderson-Darling statistic     1.4652711

		# Goodness-of-fit criteria
		#                                1-mle-weibull 2-mle-gamma 3-mle-lnorm 4-mle-norm
		# Akaike's Information Criterion      1664.425    1664.352    1681.358   1709.141
		# Bayesian Information Criterion      1669.844    1669.771    1686.777   1714.561
		#                                5-mle-logis
		# Akaike's Information Criterion    1694.606
		# Bayesian Information Criterion    1700.025					

	fw
		# Fitting of the distribution ' weibull ' by maximum likelihood 
		# Parameters:
		#         estimate Std. Error
		# shape   1.400709  0.1012143
		# scale 776.555232 55.4391560
	
	disc_ks_test(unique(x), "pweibull", exact = TRUE,shape = 1.400709 ,scale = 776.555232, 
					tol = 1e-08, sim.size = 1e+06, num.sim = 10) #	D = 0.061015, p-value = 0.802


	################## plot the distribution, Figure S6A

	args=list(shape = 1.400709 ,scale = 776.555232)
	
	pdf("Figure S6A_ggplot_fbx_size_weibull_distribution_complete.pdf",width=5,height=5)

			ggplot(size_complete, aes(x=size)) + 
					 geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
					binwidth=50,
					colour="black", fill="white") +
					geom_density(alpha=.2, fill="#FF6666")+  # Overlay with transparent density plot
			
					stat_function(fun = dweibull, args = args, color = "red") +
					geom_vline(aes(xintercept=mean(size, na.rm=T)),   # Ignore NA values for mean
						               color="red", linetype="dashed", size=1) +
					geom_vline(aes(xintercept=getmode(size)),   # Ignore NA values for mean
									   color="green", linetype="dashed", size=1) 
																		
	dev.off()



############################################################################ All orthoMCL set, Figure S6B

	d<-t(ortho_count)				
	dx<-d[rownames(d)%in%rownames(clusters),]
	size_all_ortho<-colSums(dx)

	size_all<-as.data.frame(size_all_ortho)
	x<-size_all$size
	
	fw <- fitdist(x, "weibull")
	fg <- fitdist(x, "gamma")
	fln <- fitdist(x, "lnorm")
	fn <- fitdist(x, "norm")
	fo<-fitdist(x, "logis")

	gofstat(list(fw, fg,fln,fn,fo))

			# Goodness-of-fit statistics
			#                              1-mle-weibull 2-mle-gamma 3-mle-lnorm 4-mle-norm
			# Kolmogorov-Smirnov statistic    0.06271275  0.07324201   0.1144003  0.1001482
			# Cramer-von Mises statistic      0.08887766  0.14648753   0.5637286  0.3021593
			# Anderson-Darling statistic      0.65171741  0.97902058   3.5716756  2.1689300
			#                              5-mle-logis
			# Kolmogorov-Smirnov statistic  0.08258655
			# Cramer-von Mises statistic    0.13852124
			# Anderson-Darling statistic    1.22282315

			# Goodness-of-fit criteria
			#                                1-mle-weibull 2-mle-gamma 3-mle-lnorm 4-mle-norm
			# Akaike's Information Criterion      1643.614    1646.148    1679.802   1683.195
			# Bayesian Information Criterion      1649.033    1651.567    1685.221   1688.615
			#                                5-mle-logis
			# Akaike's Information Criterion    1669.367
			# Bayesian Information Criterion    1674.786
			
	fw

		# Fitting of the distribution ' weibull ' by maximum likelihood 
		# Parameters:
		#         estimate Std. Error
		# shape   1.363577  0.1003574
		# scale 691.243583 50.4935722

	disc_ks_test(unique(x), "pweibull", exact = TRUE,	shape=1.363577 ,scale = 691.243583, 
					tol = 1e-08, sim.size = 1e+06, num.sim = 10)	#		D = 0.062713, p-value = 0.7511
		
		
	################### plot the distribution, Figure S6B


	args=list(shape=1.363577 ,scale = 691.243583)

	pdf("Figure S6B_ggplot_fbx_size_weibull_distribution_all.pdf",width=5,height=5)

			ggplot(size_all, aes(x=size)) + 
					geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
					binwidth=50,
					colour="black", fill="white") +
					geom_density(alpha=.2, fill="#FF6666")+  # Overlay with transparent density plot
			
					stat_function(fun = dweibull, args = args, color = "red") +
					geom_vline(aes(xintercept=mean(size, na.rm=T)),   # Ignore NA values for mean
						               color="red", linetype="dashed", size=1) +
					geom_vline(aes(xintercept=getmode(size)),   # Ignore NA values for mean
									      color="green", linetype="dashed", size=1) 
																		
	dev.off()
	

############################################################################ orphan cluster  fbxes, Figure S6C

	orphan<-size_complete-size_all_ortho
	size<-orphan
			
	size_orphan<-as.data.frame(size)
	x<-size_orphan$size

	pdf("orphan_fbx_size_cullen_frey_graph.pdf",width=5,height=5)
		descdist(x, discrete = FALSE,boot = 1000)
	dev.off()

	fw <- fitdist(x, "weibull")
	fg <- fitdist(x, "gamma")
	fln <- fitdist(x, "lnorm")
	fn <- fitdist(x, "norm")
	fo<-fitdist(x, "logis")

	gofstat(list(fw, fg,fln,fn,fo))

			# Goodness-of-fit statistics
			#                              1-mle-weibull 2-mle-gamma 3-mle-lnorm 4-mle-norm
			# Kolmogorov-Smirnov statistic    0.07992166  0.09610606  0.06338173  0.2016964
			# Cramer-von Mises statistic      0.16703196  0.22263914  0.09723810  1.4217498
			# Anderson-Darling statistic      1.28940127  1.47381715  0.62427200  8.2201305
			#                              5-mle-logis
			# Kolmogorov-Smirnov statistic   0.1992266
			# Cramer-von Mises statistic     0.6652236
			# Anderson-Darling statistic     5.3562719

			# Goodness-of-fit criteria
			#                                1-mle-weibull 2-mle-gamma 3-mle-lnorm 4-mle-norm
			# Akaike's Information Criterion      1175.683    1174.898    1161.732   1289.319
			# Bayesian Information Criterion      1181.102    1180.317    1167.151   1294.738
			#                                5-mle-logis
			# Akaike's Information Criterion    1263.657
			# Bayesian Information Criterion    1269.076


	fln
		# Fitting of the distribution ' lnorm ' by maximum likelihood 
		# Parameters:
		#         estimate Std. Error
		# meanlog 3.767460 0.09767072
		# sdlog   1.029025 0.06906333

	disc_ks_test(unique(x), "plnorm", exact = TRUE,	meanlog = 3.767460 ,sdlog  = 1.029025, 
			tol = 1e-08, sim.size = 1e+06, num.sim = 10) #		D = 0.16133, p-value = 0.03396

				

	##################### plot the distribution, Figure S6C

	args=list(meanlog = 3.767460 ,sdlog  = 1.029025)
	
	pdf("Figure S6C_ggplot_fbx_size_lnorm_distribution_orphan.pdf",width=5,height=5)
			ggplot(size_orphan, aes(x=size)) + 
					 geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
					binwidth=50,
					colour="black", fill="white") +
					geom_density(alpha=.2, fill="#FF6666")+  # Overlay with transparent density plot
			
					stat_function(fun = dlnorm, args = args, color = "red") +
					geom_vline(aes(xintercept=mean(size, na.rm=T)),   # Ignore NA values for mean
						               color="red", linetype="dashed", size=1) +
						geom_vline(aes(xintercept=getmode(size)),   # Ignore NA values for mean
									               color="green", linetype="dashed", size=1) 
																		
	dev.off()



############################################################################ size cluster 1 fbxes, Figure 6A,E
		
	d<-t(ortho_count)				
	dx<-d[rownames(d)%in%rownames(cluster1),]

	size<-colSums(dx)
	size_c1<-as.data.frame(size)
	x<-size_c1$size

	pdf("cluster1_fbx_size_cullen_frey_graph.pdf",width=5,height=5)
		descdist(x, discrete = FALSE,boot = 1000)
	dev.off()

	fw <- fitdist(x, "weibull")
	fg <- fitdist(x, "gamma")
	fln <- fitdist(x, "lnorm")
	fn <- fitdist(x, "norm")
	fo<-fitdist(x, "logis")
	
	gofstat(list(fw, fg,fln,fn,fo))

			# Goodness-of-fit statistics
			#                              1-mle-weibull 2-mle-gamma 3-mle-lnorm 4-mle-norm
			# Kolmogorov-Smirnov statistic    0.05070428  0.04725203  0.08913252  0.1708628
			# Cramer-von Mises statistic      0.03651713  0.03428931  0.10983035  0.8927765
			# Anderson-Darling statistic      0.23423351  0.21353705  0.73451780  5.2133642
			#                              5-mle-logis
			# Kolmogorov-Smirnov statistic   0.1498547
			# Cramer-von Mises statistic     0.4534331
			# Anderson-Darling statistic     3.7352981

			# Goodness-of-fit criteria
			#                                1-mle-weibull 2-mle-gamma 3-mle-lnorm 4-mle-norm
			# Akaike's Information Criterion      1510.837    1510.525    1521.466   1589.371
			# Bayesian Information Criterion      1516.256    1515.944    1526.885   1594.790
			#                                5-mle-logis
			# Akaike's Information Criterion    1578.744
			# Bayesian Information Criterion    1584.163

	fg

		# Fitting of the distribution ' gamma ' by maximum likelihood 
		# Parameters:
		#          estimate   Std. Error
		# shape 1.163644255 0.1258372765
		# rate  0.003547612 0.0004339777
			
	b1 <- bootdist(fg, niter= 500)

	pdf("Figure 6E_cluster1_fbx_size_fg_CIcdfplot.pdf",width=5,height=5)
		CIcdfplot(b1, CI.level= 95/100, CI.output = "probability",CI.fill = "light gray", CI.col = "red")
	dev.off()

	disc_ks_test(unique(x), "pgamma", exact = TRUE,	shape=1.163644255 ,rate =0.003547612, 
					tol = 1e-08, sim.size = 1e+06, num.sim = 10) #		D = 0.062782, p-value = 0.8066
		
		
	####################### plot the distribution, Figure 6A


	args=list(shape=1.163644255 ,rate =0.003547612)
	
	pdf("Figure 6A_ggplot_fbx_size_gamma_distribution_cluster1.pdf",width=5,height=5)

			ggplot(size_c1, aes(x=size)) + 
					 geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
					binwidth=50,
					colour="black", fill="white") +
					geom_density(alpha=.2, fill="#FF6666")+  # Overlay with transparent density plot
			
					stat_function(fun = dgamma, args = args, color = "red") +
					geom_vline(aes(xintercept=mean(size, na.rm=T)),   # Ignore NA values for mean
						               color="red", linetype="dashed", size=1) +
						geom_vline(aes(xintercept=getmode(size)),   # Ignore NA values for mean
									               color="green", linetype="dashed", size=1) 
																		
	dev.off()



############################################################################ size cluster 2 fbxes, Figure 6B,F

	d<-t(ortho_count)				
	dx<-d[rownames(d)%in%rownames(cluster2),]

	size<-colSums(dx)
	size<-size[1:101]  # only embryophytes are analyzed because algal species do not have any

	size_c2<-as.data.frame(size)
	x<-size_c2$size
	
	fw <- fitdist(x, "weibull")
	fg <- fitdist(x, "gamma")
	fln <- fitdist(x, "lnorm")
	fn <- fitdist(x, "norm")
	fo<-fitdist(x, "logis")
	fu<-fitdist(x, "unif")

	gofstat(list(fw, fg,fln,fn,fo,fu))

			# Goodness-of-fit statistics
			#                              1-mle-weibull 2-mle-gamma 3-mle-lnorm 4-mle-norm
			# Kolmogorov-Smirnov statistic     0.1248769   0.1582106  0.08964688  0.2794755
			# Cramer-von Mises statistic       0.3404767   0.5606828  0.14592420  2.1785559
			# Anderson-Darling statistic       2.0816306   2.9839489  1.04806428 11.2145151
			#                              5-mle-logis 6-mle-unif
			# Kolmogorov-Smirnov statistic   0.2627781  0.6236198
			# Cramer-von Mises statistic     1.5456586 19.0670645
			# Anderson-Darling statistic     9.3446303        Inf

			# Goodness-of-fit criteria
			#                                1-mle-weibull 2-mle-gamma 3-mle-lnorm 4-mle-norm
			# Akaike's Information Criterion      1152.867    1158.415    1138.477   1334.441
			# Bayesian Information Criterion      1158.097    1163.645    1143.707   1339.671
			#                                5-mle-logis 6-mle-unif
			# Akaike's Information Criterion    1303.864         NA
			# Bayesian Information Criterion    1309.094         NA				


	fln
		# Fitting of the distribution ' lnorm ' by maximum likelihood 
		# Parameters:
		#         estimate Std. Error
		# meanlog 3.861544 0.13920314
		# sdlog   1.398974 0.09843126

	b1 <- bootdist(fln, niter= 500)

	pdf("Figure 6F_cluster4_fbx_size_fln_CIcdfplot.pdf",width=5,height=5)
		CIcdfplot(b1, CI.level= 95/100, CI.output = "probability",CI.fill = "light gray", CI.col = "red")
	dev.off()

	disc_ks_test(unique(x), "plnorm", exact = TRUE,	meanlog =3.861544  ,sdlog  = 1.398974, 
				tol = 1e-08, sim.size = 1e+06, num.sim = 10) #	D = 0.16849, p-value = 0.03314
		
	
	####################### plot the distribution, Figure 6B

	args=list(meanlog =3.861544  ,sdlog  = 1.398974)
	
	pdf("Figure 6B_ggplot_fbx_size_lnormal_distribution_cluster2.pdf",width=5,height=5)

			ggplot(size_c2, aes(x=size)) + 
					 geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
					binwidth=50,
					colour="black", fill="white") +
					geom_density(alpha=.2, fill="#FF6666")+  # Overlay with transparent density plot
			
					stat_function(fun = dlnorm, args = args, color = "red") +
					geom_vline(aes(xintercept=mean(size, na.rm=T)),   # Ignore NA values for mean
						               color="red", linetype="dashed", size=1) +
					geom_vline(aes(xintercept=getmode(size)),   # Ignore NA values for mean
									   color="green", linetype="dashed", size=1) 
																		
	dev.off()





############################################################################ size cluster 3 fbxes, Figure 6C, G

	d<-t(ortho_count)				
	dx<-d[rownames(d)%in%rownames(cluster3),]

	size<-colSums(dx)
	size<-size[1:101] # algal species are not included
	size<-size[-c(31:35)] # cottons are removed due to recent duplications

	size_c3<-as.data.frame(size)
	x<-size_c3$size

	fw <- fitdist(x, "weibull")
	fg <- fitdist(x, "gamma")
	fln <- fitdist(x, "lnorm")
	fn <- fitdist(x, "norm")
	fo<-fitdist(x, "logis")
	fu<-fitdist(x, "unif")


	gofstat(list(fw, fg,fln,fn,fo,fu))

			# Goodness-of-fit statistics
			#                              1-mle-weibull 2-mle-gamma 3-mle-lnorm 4-mle-norm
			# Kolmogorov-Smirnov statistic    0.09217573  0.07530657   0.1013495  0.1037846
			# Cramer-von Mises statistic      0.08617713  0.09325289   0.2293934  0.2233675
			# Anderson-Darling statistic      0.75715809  0.66699277   1.2902336  1.7961163
			#                              5-mle-logis 6-mle-unif
			# Kolmogorov-Smirnov statistic  0.08458266  0.5677885
			# Cramer-von Mises statistic    0.06890904 10.8481977
			# Anderson-Darling statistic    0.89488064        Inf

			# Goodness-of-fit criteria
			#                                1-mle-weibull 2-mle-gamma 3-mle-lnorm 4-mle-norm
			# Akaike's Information Criterion      902.9387    897.3694    898.8237   929.2837
			# Bayesian Information Criterion      908.0674    902.4981    903.9524   934.4124
			#                                5-mle-logis 6-mle-unif
			# Akaike's Information Criterion    915.7419         NA
			# Bayesian Information Criterion    920.8706         NA					

	fg

		# Fitting of the distribution ' gamma ' by maximum likelihood 
		# Parameters:
		#         estimate  Std. Error
		# shape 2.98394503 0.408509525
		# rate  0.06026557 0.008982213
		
	b1 <- bootdist(fg, niter= 500)

	pdf("Figure 6G_cluster3_fbx_size_fg_CIcdfplot.pdf",width=5,height=5)
			CIcdfplot(b1, CI.level= 95/100, CI.output = "probability",CI.fill = "light gray", CI.col = "red")
	dev.off()

	disc_ks_test(unique(x), "pgamma", exact = TRUE,	shape= 2.90585101 ,rate = 0.05806406, 
				tol = 1e-08, sim.size = 1e+06, num.sim = 10) #	D = 0.060198, p-value = 0.9763
		

	######################### plot the distribution, Figure 6C

	args=list(shape= 2.98394503 ,rate = 0.06026557)
	
	pdf("Figure 6C_ggplot_fbx_size_gamma_distribution_cluster3_adjusted.pdf",width=5,height=5)

			ggplot(size_c3, aes(x=size)) + 
					 geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
					binwidth=10,
					colour="black", fill="white") +
					geom_density(alpha=.2, fill="#FF6666")+  # Overlay with transparent density plot
			
					stat_function(fun = dgamma, args = args, color = "red") +
					geom_vline(aes(xintercept=mean(size, na.rm=T)),   # Ignore NA values for mean
						               color="red", linetype="dashed", size=1) +
						geom_vline(aes(xintercept=getmode(size)),   # Ignore NA values for mean
									    color="green", linetype="dashed", size=1) 
																		
	dev.off()



############################################################################ size cluster 4 fbxes, Figure 6D, H

	d<-t(ortho_count)				
	dx<-d[rownames(d)%in%rownames(cluster4),]
	size<-colSums(dx)

	size<-size[1:101] # algal species are not included
	size<-size[-c(31:35)] # cottons are removed due to recent duplications

	size_c4<-as.data.frame(size)
	x<-size_c4$size
	
	fw <- fitdist(x, "weibull")
	fg <- fitdist(x, "gamma")
	fln <- fitdist(x, "lnorm")
	fn <- fitdist(x, "norm")
	fo<-fitdist(x, "logis")
	fu<-fitdist(x, "unif")

	gofstat(list(fw, fg,fln,fn,fo,fu))

			# Goodness-of-fit statistics
			#                              1-mle-weibull 2-mle-gamma 3-mle-lnorm 4-mle-norm
			# Kolmogorov-Smirnov statistic     0.1842984   0.1645304   0.1441993  0.1992651
			# Cramer-von Mises statistic       0.5841905   0.3877616   0.2841403  0.6422419
			# Anderson-Darling statistic       3.2406280   1.9768939   1.4248178  3.4656684
			#                              5-mle-logis 6-mle-unif
			# Kolmogorov-Smirnov statistic   0.1540587  0.2963684
			# Cramer-von Mises statistic     0.4318158  3.3020032
			# Anderson-Darling statistic     2.8475747        Inf

			# Goodness-of-fit criteria
			#                                1-mle-weibull 2-mle-gamma 3-mle-lnorm 4-mle-norm
			# Akaike's Information Criterion      1030.109    1013.835    1009.561   1030.578
			# Bayesian Information Criterion      1035.237    1018.964    1014.689   1035.707
			#                                5-mle-logis 6-mle-unif
			# Akaike's Information Criterion    1028.573         NA
			# Bayesian Information Criterion    1033.702         NA
			
	fln
		# Fitting of the distribution ' lnorm ' by maximum likelihood 
		# Parameters:
		#          estimate Std. Error
		# meanlog 4.9791031 0.03197115
		# sdlog   0.3132521 0.02260598

	b1 <- bootdist(fln, niter= 500)

	pdf("Figure 6H_cluster1_fbx_size_fln_CIcdfplot.pdf",width=5,height=5)
			CIcdfplot(b1, CI.level= 95/100, CI.output = "probability",CI.fill = "light gray", CI.col = "red")
	dev.off()

	disc_ks_test(unique(x), "plnorm", exact = TRUE,meanlog= 4.9791031, sdlog=0.3132521,
		 				tol = 1e-08, sim.size = 1e+06, num.sim = 10) #	D = 0.12894, p-value = 0.1909

	
	########################## plot the distribution, Figure 6D

	args=list(meanlog= 4.9791031, sdlog=0.3132521)
	
	pdf("Figure 6D_ggplot_fbx_size_lognormal_distribution_cluster1.pdf",width=5,height=5)
			ggplot(size_c4, aes(x=size)) + 
					 geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
					binwidth=20,
					colour="black", fill="white") +
					geom_density(alpha=.2, fill="#FF6666")+  # Overlay with transparent density plot
			
					stat_function(fun = dlnorm, args = args, color = "red") +
					geom_vline(aes(xintercept=mean(size, na.rm=T)),   # Ignore NA values for mean
						               color="red", linetype="dashed", size=1) +
						geom_vline(aes(xintercept=getmode(size)),   # Ignore NA values for mean
									   color="green", linetype="dashed", size=1) 
																		
	dev.off()
		












