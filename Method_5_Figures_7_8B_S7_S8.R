#libraries

library("phytools")
library("ggplot2")
library("dplyr")
library("ggpubr")
library("dunn.test")

####################### retrieve species in phylogenetic order  ##################


		phytozome.tree<-read.tree("Data_S1_A_bifurcated_newick_species_tree_of_111_plant_genomes.txt")

		tree_matrix<-phytools::compute.mr(phytozome.tree, type = "matrix")
		species_short<-substring(rownames(tree_matrix),1,5)
		species_short<-gsub("\\.","",species_short)

		list<-seq(1,length(species_short),1)

		species<-cbind(list,species_short)


####################### retrieve count file of subfamilie  ##################

		#
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
		
		

		kaks<-read.delim("Data_S7_Evolutionary_selection_data_of_paralogous_FBX_genes_in_each_subfamily.txt",header=T)

		kaks<-kaks[kaks$Ks<5,] #only paralogues with ks<5 were kept
 
		cluster1_kaks<-kaks[kaks[,1]%in%rownames(cluster1),]
		cluster1_kaks_ks_fun<-cluster1_kaks[,c("Subfamily","Paralogous_1","Ka.Ks","Ks", "Neutrality_Test")]
 
		cluster2_kaks<-kaks[kaks[,1]%in%rownames(cluster2),]
		cluster2_kaks_ks_fun<-cluster2_kaks[,c("Subfamily","Paralogous_1","Ka.Ks","Ks", "Neutrality_Test")]
 
		cluster3_kaks<-kaks[kaks[,1]%in%rownames(cluster3),]
		cluster3_kaks_ks_fun<-cluster3_kaks[,c("Subfamily","Paralogous_1","Ka.Ks","Ks", "Neutrality_Test")]

		cluster4_kaks<-kaks[kaks[,1]%in%rownames(cluster4),]
		cluster4_kaks_ks_fun<-cluster4_kaks[,c("Subfamily","Paralogous_1","Ka.Ks","Ks", "Neutrality_Test")]
		  
			   
####################### read in Ks ranges for dups in 27 species  ####################### 


		ks_age<-read.csv("Data_S8_Ks-based_age_distribution_in_36_flowering_plants.txt",header=T)
		ks_age_species<-substr(ks_age[,1],1,4)
		ks_age<-cbind(ks_age_species,ks_age)

		ks_age_species<-unique(ks_age_species)
		length(ks_age_species) # 36

		kaks_species<-gsub("FBX.*","",kaks[,2])
		kaks<-cbind(kaks_species,kaks)
		dim(kaks)[1]
			# [1] 48151

		kaks_w_ages<-kaks[kaks[,1]%in%ks_age_species, ]
		dim(kaks_w_ages)[1]
			# [1] 9551

		kaks_w_age_species<-unique(kaks_w_ages[,1])
		length(kaks_w_age_species)
			#[1] 27

		ks_age_species<-ks_age_species[ks_age_species%in%kaks_w_age_species]
		

	
################################# Functions  #########################################

#retrieve counts of FBXes under SSD and 3 WGDs, functional and neutral evolution, ks and ka/ks means
		
		cluster_kaks_functions<-function(cluster_kaks){
				kaks_w_ages<-cluster_kaks[cluster_kaks[,1]%in%ks_age_species, ]			
				kaks_w_age_species<-unique(kaks_w_ages[,1]) # unique FBX ids, double check
				ks_age_species<-ks_age_species[ks_age_species%in%kaks_w_age_species]

				kaks_27species_dup_stat<-c()

				dups<-c("SSD","REC","KT","OLD");

				for(i in 1:length(ks_age_species)){
	
						species<-ks_age_species[i]
						species_ks_age<-ks_age[ks_age[,1]%in%species,]
	
					#	species_kaks<-kaks_w_ages[kaks_w_ages[,1]%in%species,]
				
						species_kaks<-kaks_w_ages[kaks_w_ages[,c("kaks_species")]%in%species,]
						species_ortho_count<-ortho_count[rownames(ortho_count)%in%species,]
						species_ortho_count<-sum(species_ortho_count)
						species_total_dups<-dim(species_kaks)[1]
		
						species_total_dup_stat<-c()
	
						for(j in 1:4){
		
								dup<-dups[j]		
								species_ks_dup<-species_ks_age[species_ks_age$WGD_type %in% dup,]		
								species_kaks_dup<-species_kaks[species_kaks$Ks >= species_ks_dup$L_bound,]
		
								species_kaks_dup<-species_kaks_dup[species_kaks_dup$Ks < species_ks_dup$H_bound,]
								dup_cnt<-dim(species_kaks_dup)[1]
						
								ks_mean<-mean(species_kaks_dup$Ks,na.rm=TRUE)
								kaks_mean<-mean(species_kaks_dup$Ka.Ks,,na.rm=TRUE)
		
								activity<-table(species_kaks_dup$Neutrality_Test)
								fun<-as.numeric(activity[2])
								neu<-as.numeric(activity[3])
		
								dup_stat<-cbind(species,species_ortho_count,species_total_dups,dup,dup_cnt,fun,neu,ks_mean,kaks_mean)	
								species_total_dup_stat<-rbind(species_total_dup_stat,dup_stat)
					
									}
	
							kaks_27species_dup_stat<-rbind(kaks_27species_dup_stat,species_total_dup_stat)	
		
						}
	
					 kaks_27species_dup_stat					
					}


# Assigne species into 6 plant groups as defiend in Figure 8A
				
		pair<-function(x){
		
				group_species<-species_group[species_group[,1]%in%x,]
				group<-group_species[,4]
				group
				}

# Calculate the mean Ks, three WGD proportions in each species. A character is assigned to indicate one of the six plant groups
					
		cluster_wgd_df<-function(cluster_dup_stat){
			
			
				cluster_wgd_stat<-cluster_dup_stat[!cluster_dup_stat[,c("dup")]%in%c("SSD"),]	
				cluster_wgd_stat<-cluster_wgd_stat[!cluster_wgd_stat[,c("ks_mean")]%in%c(NaN),]

				wgd_prop<-as.numeric(cluster_wgd_stat[,c("dup_cnt")])/as.numeric(cluster_wgd_stat[,c("species_total_dups")])
				ks<-cluster_wgd_stat[,c("ks_mean")]
				wgd<-cluster_wgd_stat[,c("dup")]
				species<-cluster_wgd_stat[,c("species")]
		
				groups<-c()
				group_characters<-c()
		
				for(i in 1:length(species)){
			
						group_character<-as.character(pair(species[i]))
						group_characters<-c(group_characters,group_character)
						group<-pair(species[i])
						groups<-c(groups,group)					
			
							}
				
					df<-cbind(species,ks,wgd_prop,wgd,groups,group_characters)
					df
				
				}

				
######################### use cluster_kaks_functions #################################  

		all_kaks<-kaks
		all_kaks_27species_dup_stat<-cluster_kaks_functions(all_kaks)
		
		cluster1_kaks<-kaks[kaks$Subfamily%in%rownames(cluster1),]		
		cluster1_kaks_27species_dup_stat<-cluster_kaks_functions(cluster1_kaks)
		
		cluster2_kaks<-kaks[kaks$Subfamily%in%rownames(cluster2),]		
		cluster2_kaks_27species_dup_stat<-cluster_kaks_functions(cluster2_kaks)
		
		cluster3_kaks<-kaks[kaks$Subfamily%in%rownames(cluster3),]		
		cluster3_kaks_27species_dup_stat<-cluster_kaks_functions(cluster3_kaks)
		
		cluster4_kaks<-kaks[kaks$Subfamily%in%rownames(cluster4),]		
		cluster4_kaks_27species_dup_stat<-cluster_kaks_functions(cluster4_kaks)



######################### use function cluster_wgd_df #################################  


 	    species_group<-read.table("Data_S9_List_of_five_groups_of_plants.txt",header=T)									

		all_wgd_df<-cluster_wgd_df(all_kaks_27species_dup_stat)
		cluster1_wgd_df<-cluster_wgd_df(cluster1_kaks_27species_dup_stat)
		cluster2_wgd_df<-cluster_wgd_df(cluster2_kaks_27species_dup_stat)
		cluster3_wgd_df<-cluster_wgd_df(cluster3_kaks_27species_dup_stat)
		cluster4_wgd_df<-cluster_wgd_df(cluster4_kaks_27species_dup_stat)
	
					
		cluster_wgd_list<-list(all_wgd_df,
					   cluster1_wgd_df,
					   cluster2_wgd_df,
					   cluster3_wgd_df,
					   cluster4_wgd_df)



######################### power-law goodness-fit #################################  
					   

############################################################################
# Figure 7A, S7, S8				  										   #
#																		   #
# #power-law goodness-fit analysis                                         #
#			1) Figure 7A                           						   #
#  			2) Figure S7 					        					   #
#			3) Figure S8   									               #
############################################################################	

		#

		pchs<-c("15","19","17")	
		mycol <- c("gray","purple","green","blue","magenta")
		dups<-c("REC","KT","OLD")

		x<-0
		y<-0

		plot(x,y,ylim=c(0,1),xlim=c(0,3.0),col="light gray",pch=0)

		s <- seq(from = 0, to = 3, length = 50)

		chi_pvalue<-c()

		for(i in 1:5){
	
			cluster_df<-cluster_wgd_list[[i]]
	
			for(j in 1:3){			
				dup<-dups[j]		
			 	dup_df<-cluster_df[cluster_df[,4]%in%dup,]
				points( as.numeric(dup_df[,2]), as.numeric(dup_df[,3]), pch=as.numeric(pchs[j]),col=mycol[i] )
					}
			
			x<-as.numeric(cluster_df[,2])
			y<-as.numeric(cluster_df[,3])
	
			df <- data.frame(x, y)
			
			m <- nls(y ~ I(x^(power)/lcp), data = df, start = list(power = -1,lcp=1), trace = F)	

			lines(s, predict(m, list(x = s)), col = mycol[i])
	
			predict<-predict(m, list(x<-as.numeric(cluster_df[,2])))
	
			chi2 = sum(((y- predict)^2)/predict)
	
			p_value<-pchisq(chi2,df=1,lower.tail=TRUE)
	
			chi<-cbind(i,chi2,p_value)
	
			chi_pvalue<-rbind(chi_pvalue,chi)
		
			}		
	

		 chi_pvalue
			
		 				#
				 #       i     chi2   p_value
				 #  [1,] 1 3.727026 0.9464614
				 #  [2,] 2 3.328405 0.9319071
				 #  [3,] 3 6.840038 0.9910863
				 #  [4,] 4 4.804251 0.9716104
				 #  [5,] 5 4.521524 0.9665290
							
	
######################################################################## draw 4 Clucters together, Figure 7A, S7, S8

		#
		pchs<-c("15","19","17")	
		mycol <- c("gray","gray","green","blue","magenta")
		dups<-c("REC","KT","OLD")

		x<-0
		y<-0

		s <- seq(from = 0, to = 3, length = 50)

		pdf("Figure 7A_cluster1_2_3_4_power_law_plot.pdf",height=5,width=5)
		
			plot(x,y,ylim=c(0,1),xlim=c(0,3.0),col="light gray",pch=0)
			for(i in 1){
				
			# replace "for(i in 2:5)"with "for(i in 1)", "for(i in 2)", "for(i in 3)", "for(i in 4)","for(i in 5)"
			# to get Figure S7 and Figure S8.
			
				cluster_df<-cluster_wgd_list[[i]]
				for(j in 1:3){			
					dup<-dups[j]		
				 	dup_df<-cluster_df[cluster_df[,4]%in%dup,]
					points( as.numeric(dup_df[,2]), as.numeric(dup_df[,3]), pch=as.numeric(pchs[j]),col=mycol[i] )
						}
			
				x<-as.numeric(cluster_df[,2])
				y<-as.numeric(cluster_df[,3])	
				df <- data.frame(x, y)			
				m <- nls(y ~ I(x^(power)/lcp), data = df, start = list(power = -1,lcp=1), trace = F)	
				lines(s, predict(m, list(x = s)), col = mycol[i])	
				}		
	
		dev.off()
		




#############################################################################
# calculate 1) Figure 7D	   												#
#			2) Figure 7E											        #
#			data for plot-ly graph	                       					#
#############################################################################	
	


		cluster_list<-list(all_kaks_27species_dup_stat,
				cluster1_kaks_27species_dup_stat,
				cluster2_kaks_27species_dup_stat,
				cluster3_kaks_27species_dup_stat,
				cluster4_kaks_27species_dup_stat)
						
	  	all_ssd<-c()
	  	all_rec<-c()
	  	all_kt<-c()
	  	all_old<-c()
						
	  	for(i in 1:5){
		
	  		cluster<-paste("Cluster",(i-1),sep="")
	  		cluster_kaks<-cluster_list[[i]]
			
	  		ssd<-cluster_kaks[cluster_kaks[,c("dup")]%in%c("SSD"),]
	  		rec<-cluster_kaks[cluster_kaks[,c("dup")]%in%c("REC"),]
	  		kt<-cluster_kaks[cluster_kaks[,c("dup")]%in%c("KT"),]
	  		old<-cluster_kaks[cluster_kaks[,c("dup")]%in%c("OLD"),]
		
	  		ssd<-cbind(ssd,cluster)
	  		rec<-cbind(rec,cluster)
	  		kt<-cbind(kt,cluster)
	  		old<-cbind(old,cluster)
		
	  		all_ssd<-rbind(all_ssd,ssd)
	  		all_rec<-rbind(all_rec,rec)
	  		all_kt<-rbind(all_kt,kt)
	  		all_old<-rbind(all_old,old)
		
	  			}
			
######################################################################################################### plot_ly Figure 7D-E
		
		clusters<-c("Cluster0","Cluster1","Cluster2","Cluster3","Cluster4" ) #use Cluster0 to indicate all subfamilies
		dup_names<-c("SSD","REC","KT","OLD")

		dup_list<-list(all_ssd,all_rec,all_kt,all_old)

		dup_cluster_fun_neus<-c()
		dup_cnts<-c()

		for(i in 1:4){
	
			dup_file<-dup_list[[i]]
			dup<-dup_names[i]
	
			for(j in 1:5){
		
				cluster<-clusters[j]
		
				cluster_dup_file<-dup_file[dup_file[,10]%in%cluster,]
		
				cluster_fun_neu<-colSums ( apply(cluster_dup_file[,c(6,7)],2,as.numeric)  )		
				cluster_dup_cnt<-colSums ( apply(cluster_dup_file[,c(3,5)],2,as.numeric)  )
		
				cluster_fun_neu<-t(as.matrix(cluster_fun_neu))
				cluster_dup_cnt<-t(as.matrix(cluster_dup_cnt))
		
				rownames(cluster_fun_neu)<-paste(dup,cluster,sep="_")
				rownames(cluster_dup_cnt)<-paste(dup,cluster,sep="_")
		
				dup_cluster_fun_neus<-rbind(dup_cluster_fun_neus,cluster_fun_neu)
				dup_cnts<-rbind(dup_cnts,cluster_dup_cnt)
						
							}	
					}

		dup_ratio<-round(dup_cnts[,2]/dup_cnts[,1]*100,1)	
		dup_cnts_ratio<-cbind(dup_cnts,dup_ratio)

		neutral_ratio<-round(dup_cluster_fun_neus[,2]/(rowSums(dup_cluster_fun_neus))*100,1)	
		dup_cluster_fun_neus<-cbind(dup_cluster_fun_neus,neutral_ratio)
	
		
######################################################################### Figure 7D
		
	   dup_cnts_ratio 

 				#    species_total_dups dup_cnt dup_ratio
				#   SSD_Cluster0               9551    3580      37.5
				#   SSD_Cluster1               5100    2690      52.7
				#   SSD_Cluster2               1569     516      32.9
				#   SSD_Cluster3                706     205      29.0
				#   SSD_Cluster4               2176     169       7.8
				#   REC_Cluster0               9551    1582      16.6
				#   REC_Cluster1               5100     456       8.9
				#   REC_Cluster2               1569     257      16.4
				#   REC_Cluster3                706     150      21.2
				#   REC_Cluster4               2176     719      33.0
				#   KT_Cluster0                9551    2974      31.1
				#   KT_Cluster1                5100    1477      29.0
				#   KT_Cluster2                1569     634      40.4
				#   KT_Cluster3                 706     154      21.8
				#   KT_Cluster4                2176     709      32.6
				#   OLD_Cluster0               9551    1243      13.0
				#   OLD_Cluster1               5100     434       8.5
				#   OLD_Cluster2               1569     148       9.4
				#   OLD_Cluster3                706     172      24.4
				#   OLD_Cluster4               2176     489      22.5	
		

			  				
		#Fisher Exact Test with Bonferroni correction
	  	dups<-c("SSD","REC","KT","OLD")	

	  	fe_dup_cnt_pvalues<-c()
	
	  	for(i in 1:4){
		
	  		dup<-dups[i];				
	  		dup_df<-dup_cnts_ratio[gsub("\\_.*","",rownames(dup_cnts_ratio) )%in%dup,]

	  		dup_cnt_all<-dup_df[1,2]
	  		other_dup_cnt_all<-dup_df[1,1]-dup_df[1,2]

	  		fisher_exact_test_pvalues<-c()

	  		for(j in 2:dim(dup_df)[1]){

	  				dup_cnt_cluster<-dup_df[j,2]
	  				other_dup_cnt_cluster<-dup_df[j,1]-dup_df[j,2]

	  				test<-matrix(c(dup_cnt_cluster,dup_cnt_all,other_dup_cnt_cluster,other_dup_cnt_all), 
	  							nr=2, dimnames=list( c("a","b"),c("w","wo") ) )

	  				fe_greater_test<-fisher.test(test,alternative="greater")
	  				fe_less_test<-fisher.test(test,alternative="less")					
	  				fe_greater_test_pvalue<-c(paste(dup,"cluster",(j-1),"greater",sep="_"),fe_greater_test$p.value)
	  				fe_less_test_pvalue<-c(paste(dup,"cluster",(j-1),"less",sep="_"),fe_less_test$p.value)

	  				fisher_exact_test_pvalues<-rbind(fisher_exact_test_pvalues,fe_greater_test_pvalue)	
	  				fisher_exact_test_pvalues<-rbind(fisher_exact_test_pvalues,fe_less_test_pvalue)

	  									}
				
	  		rownames(fisher_exact_test_pvalues)<-fisher_exact_test_pvalues[,1]

	  		pvalue<-p.adjust(fisher_exact_test_pvalues[,2], method = "bonferroni")			
	  		pvalue<-pvalue[(pvalue<0.05)]
	  		pvalue<-as.matrix(pvalue)

	  		fe_dup_cnt_pvalues<-rbind(fe_dup_cnt_pvalues,signif(pvalue,digits=3))					

	  			}						
	
	  	fe_dup_cnt_pvalues
							
				#
				#	                            [,1]
				#	 SSD_cluster_1_greater  9.08e-70
				#	 SSD_cluster_2_less     1.92e-03
				#	 SSD_cluster_3_less     2.52e-05
				#	 SSD_cluster_4_less    7.15e-191
				#	 REC_cluster_1_less     1.99e-38
				#	 REC_cluster_3_greater  8.45e-03
				#	 REC_cluster_4_greater  3.13e-61
				#	 KT_cluster_1_less      2.64e-02
				#	 KT_cluster_2_greater   3.91e-12
				#	 KT_cluster_3_less      4.69e-07
				#	 OLD_cluster_1_less     4.13e-16
				#	 OLD_cluster_2_less     2.01e-04
				#	 OLD_cluster_3_greater  3.19e-14
				#	 OLD_cluster_4_greater  4.04e-26

######################################################################### Figure 7E
		
		dup_cluster_fun_neus 
	
							#
			 #               fun  neu neutral_ratio
			   # SSD_Cluster0 1913 1667          46.6
			   # SSD_Cluster1 1427 1263          47.0
			   # SSD_Cluster2  259  257          49.8
			   # SSD_Cluster3  116   89          43.4
			   # SSD_Cluster4  111   58          34.3
			   # REC_Cluster0 1251  331          20.9
			   # REC_Cluster1  233  223          48.9
			   # REC_Cluster2  169   88          34.2
			   # REC_Cluster3  148    2           1.3
			   # REC_Cluster4  701   18           2.5
			   # KT_Cluster0  2732  242           8.1
			   # KT_Cluster1  1298  179          12.1
			   # KT_Cluster2   578   56           8.8
			   # KT_Cluster3   150    4           2.6
			   # KT_Cluster4   706    3           0.4
			   # OLD_Cluster0 1172   71           5.7
			   # OLD_Cluster1  383   51          11.8
			   # OLD_Cluster2  132   16          10.8
			   # OLD_Cluster3  168    4           2.3
			   # OLD_Cluster4  489    0           0.0

#Fisher Exact Test with Bonferroni correction
	
		dups<-c("SSD","REC","KT","OLD")	
	
		fe_neu_pvalues<-c()
			
		for(i in 1:4){
				
			dup<-dups[i];				
			dup_df<-dup_cluster_fun_neus[gsub("\\_.*","",rownames(dup_cluster_fun_neus) )%in%dup,]
		
			all_fun<-dup_df[1,1]
			all_neu<-dup_df[1,2]
		
			fisher_exact_test_pvalues<-c()
		
			for(j in 2:dim(dup_df)[1]){

					cluster_fun<-dup_df[j,1]
					cluster_neu<-dup_df[j,2]
				
					test<-matrix(c(cluster_neu,all_neu,cluster_fun,all_fun), 
								nr=2, dimnames=list( c("a","b"),c("w","wo") ) )

					fe_greater_test<-fisher.test(test,alternative="greater")
					fe_less_test<-fisher.test(test,alternative="less")					
					fe_greater_test_pvalue<-c(paste(dup,"cluster",(j-1),"greater",sep="_"),fe_greater_test$p.value)
					fe_less_test_pvalue<-c(paste(dup,"cluster",(j-1),"less",sep="_"),fe_less_test$p.value)

					fisher_exact_test_pvalues<-rbind(fisher_exact_test_pvalues,fe_greater_test_pvalue)	
					fisher_exact_test_pvalues<-rbind(fisher_exact_test_pvalues,fe_less_test_pvalue)

										}
						
			rownames(fisher_exact_test_pvalues)<-fisher_exact_test_pvalues[,1]

			pvalue<-p.adjust(fisher_exact_test_pvalues[,2], method = "bonferroni")			
			pvalue<-pvalue[(pvalue<0.05)]
			pvalue<-as.matrix(pvalue)
		
			fe_neu_pvalues<-rbind(fe_neu_pvalues,signif(pvalue,digits=3))					

				}						
			
		fe_neu_pvalues
	
				#
				####
                         
			   #                         [,1]
			   # SSD_cluster_4_less    8.56e-03
			   # REC_cluster_1_greater 3.14e-29
			   # REC_cluster_2_greater 3.11e-05
			   # REC_cluster_3_less    1.58e-11
			   # REC_cluster_4_less    3.99e-37
			   # KT_cluster_1_greater  1.32e-04
			   # KT_cluster_3_less     3.80e-02
			   # KT_cluster_4_less     8.20e-19
			   # OLD_cluster_1_greater 3.95e-04
			   # OLD_cluster_4_less    2.63e-10
				


#############################################################################
# calculate 1) Figure 7B			 										#
#			2) Figure 7C											        #
#			3) Figure 8B											        #
#																			#
#############################################################################	

################################# Function  #########################################

		#
		species_dup_ratio<-function(dup){
		
				dup_ratios<-c()
		
				for(j in 1:5){
		
						cluster<-clusters[j]
		
						cluster_dup<-dup[dup[,c("cluster")]%in%cluster,]
		
						dup_ratio<- as.numeric(cluster_dup[,c("dup_cnt")])/as.numeric(cluster_dup[,c("species_total_dups")])	
		
						species<-cluster_dup[,c("species")]
		
						groups<-c()
						group_letters<-c()
	
						for(i in 1:length(species)){
		
							group_letter<-as.character(pair(species[i]))
							group_letters<-c(group_letters,group_letter)
							group<-pair(species[i])
							groups<-c(groups,group)					
		
							}		
		
						cluster_name<-rep(cluster,length(dup_ratio))
		
						cluster_dup_ratio<-cbind(dup_ratio,cluster_name,species,groups,group_letters)
		
						dup_ratios<-rbind(dup_ratios,cluster_dup_ratio)
								
									}	
				dup_ratios
	
			}

	
########################	
	
		library("ggpubr")

		clusters<-c("Cluster0","Cluster1","Cluster2","Cluster3","Cluster4" ) #the cluster order is rearranged to fit the publication order
		dup_names<-c("SSD","REC","KT","OLD")

		dup_list<-list(all_ssd,all_rec,all_kt,all_old)


######################### use function species_dup_ratio #################################  
	
		ssd_species_dup_ratio<-species_dup_ratio(all_ssd)
		rec_species_dup_ratio<-species_dup_ratio(all_rec)
		kt_species_dup_ratio<-species_dup_ratio(all_kt)
		old_species_dup_ratio<-species_dup_ratio(all_old)
	
		wgd_species_ratio<-as.numeric(rec_species_dup_ratio[,1])+as.numeric(kt_species_dup_ratio[,1])+as.numeric(old_species_dup_ratio[,1])
	
		wgd_ssd_species_dup_ratio<-cbind(wgd_species_ratio,ssd_species_dup_ratio)

	
######################################################################### Figure 7B

	
		Measure<-as.numeric(wgd_ssd_species_dup_ratio[,1])
		Group<-as.character(wgd_ssd_species_dup_ratio[,3])
	
		wgd_dups<-data.frame(Measure,Group)	
			
		
   		pdf("Figure 7B_complete_&_4_clusters_wgd_proportion_boxplot_sum.pdf",height=5,width=5)
   		 	ggboxplot(wgd_dups, x = "Group", y = "Measure", 
   		          color = "Group", palette = c("gray","#6a0dad","#008000", "#0000FF","#FF00FF" ),
   		          order = c("Cluster0","Cluster1","Cluster2","Cluster3","Cluster4"),
   		          ylab = "Proportion", xlab = "Clusters")
   		dev.off()
				   			  
	
		library("dunn.test")
		
		dunn.test(Measure,Group, method="BH", list=TRUE)	
		
				#		
				#     		     Kruskal-Wallis rank sum test
				# List of pairwise comparisons: Z statistic (adjusted p-value)
				# -----------------------------------------
				# Cluster0 - Cluster1 :  2.003353 (0.0376)
				# Cluster0 - Cluster2 : -0.338545 (0.3675)
				# Cluster1 - Cluster2 : -2.341898 (0.0192)*
				# Cluster0 - Cluster3 : -1.486809 (0.0979)
				# Cluster1 - Cluster3 : -3.490162 (0.0012)*
				# Cluster2 - Cluster3 : -1.148263 (0.1394)
				# Cluster0 - Cluster4 : -2.882874 (0.0066)*
				# Cluster1 - Cluster4 : -4.886227 (0.0000)*
				# Cluster2 - Cluster4 : -2.544328 (0.0137)*
				# Cluster3 - Cluster4 : -1.396064 (0.1017)
				
			  
########################################################################## Figure 7C

		
		ssd_species_dup_ratio<-species_dup_ratio(all_ssd)

		Measure<-as.numeric(ssd_species_dup_ratio[,1])
		Group<-as.character(ssd_species_dup_ratio[,2])

		ssd_dups<-data.frame(Measure,Group)	
	   
  		pdf("Figure 7C_complete_&_4_clusters_ssd_proportion_boxplot_sum.pdf",height=5,width=5)

  			ggboxplot(ssd_dups, x = "Group", y = "Measure", 
  		          color = "Group", palette = c("gray","#6a0dad","#008000", "#0000FF","#FF00FF" ),
  		          order = c("Cluster0","Cluster1","Cluster2","Cluster3","Cluster4"),
  		          ylab = "Proportion", xlab = "Clusters")
  		dev.off()	   
		
		dunn.test(Measure,Group, method="BH", list=TRUE)	
		
				#
				#    Kruskal-Wallis rank sum test
				#	List of pairwise comparisons: Z statistic (adjusted p-value)
				#	-----------------------------------------
				#	Cluster0 - Cluster1 : -2.146647 (0.0265)
				#	Cluster0 - Cluster2 :  0.394282 (0.3467)
				#	Cluster1 - Cluster2 :  2.540929 (0.0111)*
				#	Cluster0 - Cluster3 :  1.706803 (0.0628)
				#	Cluster1 - Cluster3 :  3.853450 (0.0003)*
				#	Cluster2 - Cluster3 :  1.312521 (0.1052)
				#	Cluster0 - Cluster4 :  3.252389 (0.0019)*
				#	Cluster1 - Cluster4 :  5.399036 (0.0000)*
				#	Cluster2 - Cluster4 :  2.858107 (0.0053)*
				#	Cluster3 - Cluster4 :  1.545585 (0.0764)
					
	   
########################################################################## Figure 8B

	
	 	phytozome.tree<-read.tree("Data_S1_A_Bi-furcated_newick_species_tree_of_111_plant_genomes.txt")
		
	 	tree_matrix<-compute.mr(phytozome.tree, type = "matrix")
	 	species_short<-substring(rownames(tree_matrix),1,5)
	 	species_short<-gsub("\\.","",species_short)
		 
		
		brassicales<-species_short[1:29]	
		s_r<-species_short[30:62]				
		superasterids<-species_short[63:72]
		poaceae<-species_short[73:87]
		basal<-species_short[88:101]
		algae<-species_short[102:111]
	
		list_clusters<-list(cluster1_kaks_ks_fun,
							cluster2_kaks_ks_fun,
							cluster3_kaks_ks_fun,
							cluster4_kaks_ks_fun)
		list_groups<-list(brassicales,
							s_r,
							superasterids,
							poaceae,
							basal,
							algae)
	
		group_names<-c("brassicales","s_r","superasterids","poaceae","basal","algae")
		
		fun_neu_comparison<-c()
		
		for(i in 1:4){
			cluster<-list_clusters[[i]]
			cluster_name<-paste("cluster",i,sep="")
			
			for(j in 1:6){
				
				group<-list_groups[[j]]				
				cluster_group<-cluster[substring(cluster$Paralogous_1,1,4)%in%group,]				
				fun_neu<-table(cluster_group$Neutrality_Test)
								
				cluster_group_fun_neu<-cbind(cluster_name,group_names[j],fun_neu[c("Functional")],fun_neu[c("Neutral")])				
				fun_neu_comparison<-rbind(fun_neu_comparison,cluster_group_fun_neu)				
				
					}
			}
 
		colnames(fun_neu_comparison)<-c("Clusters","Species_Groups","Fun","Neu")
 	   
		fun_neu<-apply(fun_neu_comparison[,3:4],2,as.numeric)	
		neu_ratio<-round(fun_neu[,c("Neu")]/(rowSums(fun_neu))*100,2)	
		neu_ratio<-cbind(fun_neu_comparison,neu_ratio)
		
		neu_ratio 
					#            Clusters   Species_Groups  Fun    Neu    neu_ratio
				#  Functional "cluster1" "brassicales"   "3012" "4184" "58.14"  
				#  Functional "cluster1" "s_r"           "2579" "2426" "48.47"  
				#  Functional "cluster1" "superasterids" "2199" "1733" "44.07"  
				#  Functional "cluster1" "poaceae"       "4759" "1598" "25.14"  
				#  Functional "cluster1" "basal"         "720"  "516"  "41.75"  
				#  Functional "cluster1" "algae"         "194"  "94"   "32.64"  
				#  Functional "cluster2" "brassicales"   "4627" "3150" "40.5"   
				#  Functional "cluster2" "s_r"           "712"  "387"  "35.21"  
				#  Functional "cluster2" "superasterids" "424"  "246"  "36.72"  
				#  Functional "cluster2" "poaceae"       "75"   "4"    "5.06"   
				#  Functional "cluster2" "basal"         "60"   "30"   "33.33"  
				#  Functional "cluster2" "algae"         "0"    "0"    "NaN"    
				#  Functional "cluster3" "brassicales"   "721"  "260"  "26.5"   
				#  Functional "cluster3" "s_r"           "1270" "413"  "24.54"  
				#  Functional "cluster3" "superasterids" "432"  "107"  "19.85"  
				#  Functional "cluster3" "poaceae"       "141"  "13"   "8.44"   
				#  Functional "cluster3" "basal"         "105"  "42"   "28.57"  
				#  Functional "cluster3" "algae"         "0"    "0"    "NaN"    
				#  Functional "cluster4" "brassicales"   "2575" "234"  "8.33"   
				#  Functional "cluster4" "s_r"           "4129" "660"  "13.78"  
				#  Functional "cluster4" "superasterids" "868"  "48"   "5.24"   
				#  Functional "cluster4" "poaceae"       "1478" "100"  "6.34"   
				#  Functional "cluster4" "basal"         "746"  "67"   "8.24"   
				#  Functional "cluster4" "algae"         "0"    "0"    "NaN" 
	

							   	
	


