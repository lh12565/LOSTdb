array_pipe_func<-function(targets,meta,type="450k",input_type="array",is.normal=F,normal_samples=NULL,filter_sam=NULL,read_force=F,method="hc"){
	if(input_type=="array"){
		rgSet <- read.metharray.exp(targets=targets,force=read_force)
		sampleNames(rgSet) <- targets$Sample_Name
		save(rgSet,file="rgSet.RData")
		if(!is.null(filter_sam)){
			inter<-intersect(filter_sam,targets$Sample_Name)
			targets<-targets[ targets$Sample_Name%in%inter,]
			rgSet<-rgSet[,inter]
		}
		detP <- detectionP(rgSet)
		save(detP,file="detP.RData")
		
		##detP distribution
		pdf('detp.pdf')
		pal <- brewer.pal(8,"Dark2")
		par(mfrow=c(1,2))
		barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
				cex.names=0.8, ylab="Mean detection p-values")
		abline(h=0.05,col="red")
		legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal,
			   bg="white")

		barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
				cex.names=0.8, ylim=c(0,0.002), ylab="Mean detection p-values")
		abline(h=0.05,col="red")
		legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal, 
			   bg="white")
		dev.off()
		
		mSetSq <- preprocessFunnorm(rgSet)
		detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 
		keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
		table(keep)
		mSetSqFlt <- mSetSq[keep,]
		mSetSqFlt
		if(type=="450k"){
			ann450k= getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
			keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% 
															c("chrX","chrY")])
			table(keep)
			mSetSqFlt <- mSetSqFlt[keep,]
			mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
			dataDirectory <- system.file("extdata", package = "methylationArrayAnalysis")
			xReactiveProbes <- read.csv(file=paste(dataDirectory,
												   "48639-non-specific-probes-Illumina450k.csv",
												   sep="/"), stringsAsFactors=FALSE)
			keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
			table(keep)
			mSetSqFlt <- mSetSqFlt[keep,]
			id_gene_loc<-id_gene_loc_ann450
		}else{
			annEPIC= getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
			keep <- !(featureNames(mSetSqFlt) %in% annEPIC$Name[annEPIC$chr %in% c("chrX","chrY")])
			mSetSqFlt <- mSetSqFlt[keep,]
			mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
			id_gene_loc<-id_gene_loc_annEPIC
		}
		print(mSetSqFlt)
		save(mSetSq,file="mSetSq.RData")
		save(mSetSqFlt,file="mSetSqFlt.RData")
		
		bval_loc <- getBeta(mSetSqFlt)
	}else if(input_type=="bval"){
		targets<-na.omit(targets)
		if(type=="450k"){
			ann450k= getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
			keep <- !(rownames(targets) %in% ann450k$Name[ann450k$chr %in% 
															c("chrX","chrY")])
			table(keep)
			targets <- targets[keep,]
			dataDirectory <- system.file("extdata", package = "methylationArrayAnalysis")
			xReactiveProbes <- read.csv(file=paste(dataDirectory,
												   "48639-non-specific-probes-Illumina450k.csv",
												   sep="/"), stringsAsFactors=FALSE)
			keep <- !(rownames(targets) %in% xReactiveProbes$TargetID)
			table(keep)
			bval_loc <- targets[keep,]
			id_gene_loc<-id_gene_loc_ann450
		}else{
			annEPIC= getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
			keep <- !(rownames(targets) %in% annEPIC$Name[annEPIC$chr %in% c("chrX","chrY")])
			bval_loc <- targets[keep,]
			id_gene_loc<-id_gene_loc_annEPIC
		}
		print(paste0("After filter the probs: ",nrow(bval_loc)))
	}else{
		stop("input_type must be array or bval")
	}
	
	if(is.normal==F){
		bval_loc_t<-as.data.frame(bval_loc)
		inter<-intersect(rownames(bval_loc_t),id_gene_loc$id)
		bval_cpg_t<-bval_loc_t[inter,]
		bval_cpg_t$id<-rownames(bval_cpg_t)
		bval2<-merge(bval_cpg_t,id_gene_loc,by="id")
		bval2<-aggregate(bval2[,2:(ncol(bval2)-2)],by=list(bval2[,"gene"]),mean)
		rownames(bval2)<-bval2$Group.1
		bval_t<-bval2[,-1]
		save(bval_loc_t,file="bval_loc_t.RData")
		save(bval_cpg_t,file="bval_cpg_t.RData")
		save(bval_t,file="bval_t.RData")
		meta_t<-meta[match(colnames(bval_t),meta$Sample_ID),]
	}else{
		bval_loc_tn<-as.data.frame(bval_loc)
		inter<-intersect(rownames(bval_loc_tn),id_gene_loc$id)
		bval_cpg_tn<-bval_loc_tn[inter,]
		bval_cpg_tn$id<-rownames(bval_cpg_tn)
		bval2<-merge(bval_cpg_tn,id_gene_loc,by="id")
		bval2<-aggregate(bval2[,2:(ncol(bval2)-2)],by=list(bval2[,"gene"]),mean)
		rownames(bval2)<-bval2$Group.1
		bval_tn<-bval2[,-1]
		meta_tn<-meta[match(colnames(bval_tn),meta$Sample_ID),]
		save(bval_loc_tn,file="bval_loc_tn.RData")
		save(bval_cpg_tn,file="bval_cpg_tn.RData")
		save(bval_tn,file="bval_tn.RData")
		t_sam<-setdiff(colnames(bval_tn),normal_samples)
		bval_loc_t<-bval_loc_tn[,t_sam]
		bval_cpg_t<-bval_cpg_tn[,t_sam]
		bval_t<-bval_tn[,t_sam]
		meta_t<-meta[match(colnames(bval_t),meta$Sample_ID),]
		save(bval_loc_t,file="bval_loc_t.RData")
		save(bval_cpg_t,file="bval_cpg_t.RData")  #bval site cpgi+promoter
		save(bval_t,file="bval_t.RData")
	}
	
	bval2<-bval_t[rowMeans(bval_t)>0.2,]
	mval<-log2(bval2/(1-bval2))
	mval[mval==Inf]<-log2((1-0.000001)/(1-(1-0.000001)))
	mval[mval==-Inf]<-log2((0+0.000001)/(1-(0+0.000001)))
	print(paste0("NA number: ",sum(is.na(mval))))
	print(paste0("is Inf?: ",max(mval)))
	print(paste0("is -Inf?: ",min(mval)))
	mval.norm<-t(scale(t(mval)))
	if(method=="hc"){
		res = ConsensusClusterPlus(as.matrix(mval.norm),maxK=10,reps=1000,pItem=0.8,pFeature=1,
		 title="CC_png_tcga_100sample_bval0.2mean_mval_hc_scale",clusterAlg="hc",distance="pearson",seed=123,plot="png")
		save(res,file="res_cc_hc.RData")
	}else{
		res = ConsensusClusterPlus(as.matrix(mval.norm),maxK=10,reps=1000,pItem=0.8,pFeature=1,
		 title="CC_png_tcga_100sample_bval0.2mean_mval_pam_scale",clusterAlg="pam",distance="euclidean",seed=123,plot="png")
		save(res,file="res_cc_pam.RData")
	}
	my_t<-res[[3]][["consensusClass"]]
	c1<-cc_mean_mval_func(my_t,mval,1)
	c2<-cc_mean_mval_func(my_t,mval,2)
	c3<-cc_mean_mval_func(my_t,mval,3)
	myorder<-order(c(c1,c2,c3))
	my_t[my_t==myorder[1]]<-"CIMP-low"
	my_t[my_t==myorder[2]]<-"CIMP-int"
	my_t[my_t==myorder[3]]<-"CIMP-high"
	heatmap_meth_rna_func(bval2,my_t,paste0("heatmap_meth_genebval0.2_",method,"_scale_3cluster.pdf"))
	heatmap_meth_rna_func(mval,my_t,paste0("heatmap_meth_genebval0.2_",method,"_scale_3cluster_mval.pdf"))
	meta_t$Classical_subtype<-my_t
	save(meta_t,file="meta_t.RData")
	if(is.normal){
		meta_tn$Classical_subtype<-my_t[ meta_tn$Sample_ID]
		meta_tn$Classical_subtype[is.na(meta_tn$Classical_subtype)]<-"Normal"
		save(meta_tn,file="meta_tn.RData")
	}
	return(bval_t)
}

cc_mean_mval_func<-function(consensus_cluster,mval,cluster){
	sam1<-names(consensus_cluster[consensus_cluster==cluster])
	mval2<-mval[,sam1,drop=F]
	mean(rowMeans(mval2))
}

array_bval_pip_func<-function(bval_t,meta,method="pam"){
	bval2<-bval_t[rowMeans(bval_t)>0.2,]
	
	mysd<-apply(bval2,1,sd)
	index<-which(mysd!=0)
	bval2<-bval2[index,]
	
	mval<-log2(bval2/(1-bval2))
	mval[mval==Inf]<-log2((1-0.000001)/(1-(1-0.000001)))
	mval[mval==-Inf]<-log2((0+0.000001)/(1-(0+0.000001)))
	print(paste0("NA number: ",sum(is.na(mval))))
	print(paste0("is Inf?: ",max(mval)))
	print(paste0("is -Inf?: ",min(mval)))
	mval.norm<-t(scale(t(mval)))
	if(method=="hc"){
		res = ConsensusClusterPlus(as.matrix(mval.norm),maxK=10,reps=1000,pItem=0.8,pFeature=1,
		 title="CC_png_tcga_100sample_bval0.2mean_mval_hc_scale",clusterAlg="hc",distance="pearson",seed=123,plot="png")
		save(res,file="res_cc_hc.RData")
	}else if(method=="pam"){
		res = ConsensusClusterPlus(as.matrix(mval.norm),maxK=10,reps=1000,pItem=0.8,pFeature=1,
		 title="CC_png_tcga_100sample_bval0.2mean_mval_pam_scale",clusterAlg="pam",distance="euclidean",seed=123,plot="png")
		save(res,file="res_cc_pam.RData")
	}else if(method=="km"){
		res = ConsensusClusterPlus(as.matrix(mval.norm),maxK=10,reps=1000,pItem=0.8,pFeature=1,
		 title="CC_png_tcga_100sample_bval0.2mean_mval_kmeans_scale",clusterAlg="km",distance="euclidean",seed=123,plot="png")
		save(res,file="res_cc_km.RData")
	}else{
		stop("The method can only be one of hc, pam and km")
	}
	my_t<-res[[3]][["consensusClass"]]
	c1<-cc_mean_mval_func(my_t,mval,1)
	c2<-cc_mean_mval_func(my_t,mval,2)
	c3<-cc_mean_mval_func(my_t,mval,3)
	myorder<-order(c(c1,c2,c3))
	my_t[my_t==myorder[1]]<-"CIMP-low"
	my_t[my_t==myorder[2]]<-"CIMP-int"
	my_t[my_t==myorder[3]]<-"CIMP-high"
	heatmap_meth_rna_func(bval2,my_t,paste0("heatmap_meth_genebval0.2_",method,"_scale_3cluster.pdf"))
	heatmap_meth_rna_func(mval,my_t,paste0("heatmap_meth_genebval0.2_",method,"_scale_3cluster_mval.pdf"))
	meta_t<-meta[match(colnames(bval_t),meta$Sample_ID),]
	meta_t$Classical_subtype<-my_t
	save(meta_t,file="meta_t.RData")
	return(my_t)
}


cc_cluster_modify_func<-function(res,bval_t,cluster=4,is.normal=F,saves=T,meta_extra=NULL){
	bval2<-bval_t[rowMeans(bval_t)>0.2,]
	
	mysd<-apply(bval2,1,sd)
	index<-which(mysd!=0)
	bval2<-bval2[index,]
	
	mval<-log2(bval2/(1-bval2))
	mval[mval==Inf]<-log2((1-0.000001)/(1-(1-0.000001)))
	mval[mval==-Inf]<-log2((0+0.000001)/(1-(0+0.000001)))
	my_t<-res[[cluster]][["consensusClass"]]
	val<-list()
	for(i in 1:cluster){
		val[[i]]<-cc_mean_mval_func(my_t,mval,i)
	}
	myorder<-order(unlist(val),decreasing=T)
	my_t[my_t==myorder[1]]<-"CIMP-high"
	my_t[my_t==myorder[length(myorder)]]<-"CIMP-low"
	num<-length(myorder)-2
	if(num==1){
		my_t[my_t==myorder[2]]<-"CIMP-int"
	}else if(num==0){
		
	}else{
		for(i in 2:(2+num-1)){
			my_t[my_t==myorder[i]]<-paste0("CIMP-int-",(i-1))
		}
	}
	
	heatmap_meth_rna_func(mval,my_t,paste0("heatmap_meth_genebval0.2_modify_scale_",cluster,"cluster_mval.pdf"),meta_extra=meta_extra)
	
	if(saves){
		load("meta_t.RData")
		my_t<-my_t[meta_t$Sample_ID]
		meta_t$Classical_subtype<-my_t
		print(dim(meta_t))
		save(meta_t,file="meta_t.RData")
		if(is.normal){
			load("meta_tn.RData")
			meta_tn$Classical_subtype<-my_t[ meta_tn$Sample_ID]
			meta_tn$Classical_subtype[is.na(meta_tn$Classical_subtype)]<-"Normal"
			print(dim(meta_tn))
			save(meta_tn,file="meta_tn.RData")
		}
	}
	return(my_t)
}



gmt2list <- function(gmtfile){
 sets <- as.list(read_lines(gmtfile))
 for(i in 1:length(sets)){
    tmp = str_split(sets[[i]], '\t')
  n = length(tmp[[1]])
  names(sets)[i] = tmp[[1]][1]
  sets[[i]] = tmp[[1]][3:n]
  rm(tmp, n)
 }
 return(sets)
}


diff_path_func<-function(meta,bval,s.sets,type="array",arraytype="450K",genome="hg19",dmr_method="mCSEA",filter_probe=F,group="Classical_subtype",group_myname,dataset,cancertype){
	rownames(meta)<-meta$Sample_ID
	if(type=="array"){
		if(arraytype=="450K"){
			if(filter_probe){
				inter<-intersect(rownames(ann450k),rownames(bval))
				bval<-bval[inter,]
			}
			ann_tmp <- ann450k[match(rownames(bval),ann450k$Name),
								  c(1:3,19,24,26)]
			ann_tmp<-ann_tmp[!is.na(ann_tmp$pos),]
		}else{
			ann_tmp <- annEPIC[match(rownames(bval),annEPIC$Name),
					  c(1:3,19,22,24)]
		}
		gr<-GRanges(
			seqnames = Rle(ann_tmp$chr),
			ranges = IRanges(ann_tmp$pos, end = ann_tmp$pos, names = ann_tmp$Name),
			strand = Rle(strand(ann_tmp$strand)))
		
		meta<-del_info_func(meta,group)
		cluster_num<-table(meta[,group])
		subtypes<-names(cluster_num)[cluster_num>=5]
		if(length(subtypes)==0){
			print("After filter, there is no subtype. Terminated")
			return(list(dmp=NULL,dmr=NULL,pathway=NULL))
		}
		
		subtypes2<-c()
		for(i in subtypes){
			othersub<-setdiff(as.character(meta[,group]),i)
			meta_othersub<-meta[meta[,group]%in%othersub,]
			if(nrow(meta_othersub)>=5){
				subtypes2<-c(subtypes2,i)
			}
		}
		subtypes<-subtypes2
		if(length(subtypes)==1){
			print("After filter, there is only one subtype. Terminated")
			return(list(dmp=NULL,dmr=NULL,pathway=NULL))
		}
		
		print(subtypes)
		
		bval<-bval[,match(meta$Sample_ID,colnames(bval))]
		diff_list<-lapply(subtypes,function(i){
			meta$tmp_subtype<-as.character(meta[,group])
			meta$tmp_subtype<-ifelse(meta$tmp_subtype%in%i,i,"Control")
			obj<-SummarizedExperiment(assays=list(counts=as.matrix(bval)),rowRanges=gr,colData=meta)
			plot_i=gsub("\\/","--",i)
			tcgabio_res <- TCGAanalyze_DMC(data = obj,
								   groupCol = "tmp_subtype",
								   group1 = i,
								   group2 = 'Control',
								   plot.filename=paste0("methylation_volcano_",plot_i,".pdf"),
								   cores = 40
								   )
			print(table(tcgabio_res$status))
			colnames(tcgabio_res)[1:5]<-c("mean.Exp","mean.Control","FC","p.value","p.adj")
			tcgabio_res<-merge(tcgabio_res,ann_tmp,by='row.names')
			tcgabio_res$group<-i
			tcgabio_res<-tcgabio_res[order(tcgabio_res$p.value),]
			tcgabio_res$col_change<-group_myname
			tcgabio_res$dataset<-dataset
			tcgabio_res$cancertype<-cancertype
			tcgabio_res
		})
		names(diff_list)<-subtypes
		save(diff_list,file="diff_list.RData")
		diffs<-do.call(rbind,diff_list)
		##pathway
		path_list<-lapply(names(diff_list),function(x){
			mytype<-c("Hypermethylated","Hypomethylated")
			hyper_hypo<-lapply(mytype,function(hh){
				index<-grep(hh,diff_list[[x]]$status)
				sigCpGs<-diff_list[[x]]$Row.names[index]
				if(length(sigCpGs)>=50){
					if(length(sigCpGs)>10000){
						sigCpGs<-sigCpGs[1:10000]
					}
					all <- diff_list[[x]]$Row.names
					gset<-lapply(names(s.sets),function(i){
						sets<-gmt2list(s.sets[[i]])
						gsa <- gsameth(sig.cpg=as.character(sigCpGs), all.cpg=as.character(all), collection=sets,sig.genes = TRUE,array.type=arraytype)
						gsa$group<-x
						gsa$pathway<-i
						gsa<-gsa[order(gsa$"P.DE"),]
						gsa$col_change<-group_myname
						gsa$dataset<-dataset
						gsa$cancertype<-cancertype
						gsa$hyper_hypo<-hh
						gsa
					})
					names(gset)<-names(s.sets)
				}else{
					gset=NULL
				}
				gset
			})
			hyper_hypo<-do.call(rbind,hyper_hypo)
		})
		names(path_list)<-names(diff_list)
		nullclass<-sapply(path_list,class)
		index<-which(nullclass=="NULL")
		if(length(index)!=0){
			path_list<-path_list[-(index)]
		}
		save(path_list,file="path_list.RData")
		##DMR
		if(dmr_method=="mCSEA"){
			if(genome=='hg19'){
				gene.obj<-gene.obj_hg19
				gene_trans<-gene_trans_hg19
				cpg.obj<-cpg.obj_hg19
			}else if(genome=='hg38'){
				gene.obj<-gene.obj_hg38
				gene_trans<-gene_trans_hg38
				cpg.obj<-cpg.obj_hg38
			}
			dmr_list<-lapply(subtypes,function(x){
				meta$tmp_subtype<-as.character(meta[,group])
				meta$tmp_subtype<-ifelse(meta$tmp_subtype%in%x,x,"Control")
				index<-which(colnames(meta)=="tmp_subtype")
				pheno<-meta[,index,drop=F]
				pheno$tmp_subtype<-factor(pheno$tmp_subtype,levels=c(x,"Control"))
				if(sum(colnames(bval)==rownames(pheno))!=ncol(bval)){
					stop("DMR: colnames bval and rownames pheno are not equal!")
				}
				
				bval[bval==1]<-1-0.000001
				bval[bval==0]<-0.000001
				
				myRank <- rankProbes(bval, pheno, refGroup = "Control")
				if(arraytype=="450K"){
					arraytype2<-"450k"
				}else{
					arraytype2<-"EPIC"
				}
				
				
				myResults <- mCSEATest(myRank, bval, pheno,
                    regionsTypes = "CGI", platform = arraytype2)
				myres<-myResults[[1]]
				site<-strsplit(rownames(myres),split=":|-")
				site<-as.data.frame(do.call(rbind,site))
				colnames(site)<-c("chr","start","end")
				mydiff<-cbind(site,myres)
				anno1<-annotateWithGeneParts(as(site,"GRanges"),gene.obj)
				tss<-anno1@dist.to.TSS
				tss$geneid<-gene_trans[match(tss$feature.name,gene_trans$V5),"V6"]
				tss$gene<-gene_trans[match(tss$feature.name,gene_trans$V5),"V7"]
				anno2=annotateWithFeatureFlank(as(site,"GRanges"),
												cpg.obj$CpGi,cpg.obj$shores,
									 feature.name="CpGi",flank.name="shores")
				site2<-data.frame(promoter=anno1@members[,1],exon=anno1@members[,2],intron=anno1@members[,3],CpGi=anno2@members[,1],shores=anno2@members[,2])
				mydiff$group<-x
				mydiff$type<-ifelse(mydiff$NES>0,"Hyper","Hypo")
				tmp<-cbind(mydiff,site2)
				tmp$target.row<-1:nrow(tmp)
				mydiff<-merge(tmp,tss,by='target.row',all.x=T)
				mydiff<-mydiff[,-1]
				mydiff<-mydiff[order(mydiff$pval),]
				mydiff$col_change<-group_myname
				mydiff$dataset<-dataset
				mydiff$cancertype<-cancertype
				mydiff
			})
		}else if(dmr_method=='dmrcate'){
			grset=makeGenomicRatioSetFromMatrix(as.matrix(bval),what="Beta")
			mval = getM(grset)
			mval[mval==Inf]<-log2((1-0.000001)/(1-(1-0.000001)))
			mval[mval==-Inf]<-log2((0+0.000001)/(1-(0+0.000001)))
			print(max(mval))
			print(min(mval))
			dmr_list<-lapply(subtypes,function(x){
				x2<-gsub("-","_",x)
				
				tmp<-gsub("-","_",as.character(meta[,group]))
				meta$tmp_subtype<-tmp
				meta$tmp_subtype<-ifelse(meta$tmp_subtype%in%x2,"ex","con")
				print(table(meta$tmp_subtype))
				cla<-factor(meta$tmp_subtype)
				design <- model.matrix(~0+cla, data=meta)
				colnames(design)<-levels(cla)
				rownames(design)<-colnames(mval)
				contrast.matrix<-makeContrasts("ex-con",levels=design)
				myAnnotation <- cpg.annotate(object = mval, datatype = "array", what = "M", 
							 analysis.type = "differential", design = design, 
							 contrasts = TRUE, cont.matrix = contrast.matrix, 
							 coef = "ex-con", arraytype = arraytype)
				DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
				results.ranges <- extractRanges(DMRs,genome = genome)
				results.ranges$group<-x
				results.ranges
			})
		}
		names(dmr_list)<-subtypes
	}
	dmpall<-do.call(rbind,diff_list)
	dmrall<-do.call(rbind,dmr_list)
	path_tmp<-lapply(path_list,function(x){
		do.call(rbind,x)
	})
	pathall<-do.call(rbind,path_tmp)
	meth_enrich_all<-list(dmp=dmpall,dmr=dmrall,pathway=pathall)
	save(meth_enrich_all,file=paste0("meth_enrich_all_",group,".RData"))
	return(meth_enrich_all)
}


#seq
##methylKit
overlap_func<-function(mobj,ref_range){
	obj_list<-lapply(mobj,function(x){
		hits<-S4Vectors::subjectHits(IRanges::findOverlaps(ref_range,as(x,"GRanges")))
		hits<-unique(hits)
		x[hits]
	})
	obj_list<-new("methylRawList", obj_list, treatment=mobj@treatment)
}

seq_pipe_func<-function(myobj,meta,impute=F,is.normal=F,normal_pattern=NULL,min.per.group=NULL,method='pam',genome='hg19'){
	myobj.filt <- filterByCoverage(myobj,
						  lo.count=10,
						  lo.perc=NULL,
						  hi.count=NULL,
						  hi.perc=99.9)
	myobj.filt.norm <- normalizeCoverage(myobj.filt, method = "median")
	meth_all <- unite(myobj.filt.norm, destrand=FALSE,mc.cores=length(myobj),min.per.group=min.per.group)
	print(paste0("Number of meth_all site: ",nrow(meth_all)))
	save(meth_all,file="meth_all.RData")
	if(genome=='hg19'){
		cpgi_promoter<-cpgi_promoter_hg19
		cpgi_promoter_file<-cpgi_promoter_hg19_file
	}else if(genome=='hg38'){
		cpgi_promoter<-cpgi_promoter_hg38
		cpgi_promoter_file<-cpgi_promoter_hg38_file
	}
	myobj_reg<-regionCounts(myobj.filt.norm,unique(cpgi_promoter))
	meth_reg<-unite(myobj_reg, destrand=FALSE,mc.cores=length(myobj),min.per.group=min.per.group)
	save(meth_reg,file="meth_reg.RData")
	print(paste0("Number of loc_cpgi+promoter site: ",nrow(meth_reg)))
	if(sum(is.na(meth_reg)!=0)&impute==F){
		print("Note: your beta value have NA, please check it! or use impute=T!")
	}
	if(is.normal==F){
		bval_loc_t<-bval_mval_func(meth_all,NULL,type="loc",impute=impute)
		bval_t<-bval_mval_func(meth_reg,cpgi_promoter_file,impute=impute)
		print(paste0("Number of bval_t gene: ",nrow(bval_t)))
		bval_cpg_t<-bval_mval_func(meth_reg,NULL,type="loc",impute=impute)
		save(bval_cpg_t,file="bval_cpg_t.RData")  #bval site cpgi+promoter
		save(bval_t,file="bval_t.RData")
		save(bval_loc_t,file="bval_loc_t.RData")
		meta_t<-meta[match(colnames(bval_t),meta$Sample_ID),]
	}else{
		bval_loc_tn<-bval_mval_func(meth_all,NULL,type="loc",is.normal=T,normal_pattern=normal_pattern,impute=impute)
		bval_tn<-bval_mval_func(meth_reg,cpgi_promoter_file,is.normal=T,normal_pattern=normal_pattern,impute=impute)
		print(paste0("Number of bval_tn gene: ",nrow(bval_tn)))
		bval_cpg_tn<-bval_mval_func(meth_reg,NULL,type="loc",is.normal=T,normal_pattern=normal_pattern,impute=impute)
		meta_tn<-meta[match(colnames(bval_tn),meta$Sample_ID),]
		save(bval_cpg_tn,file="bval_cpg_tn.RData")  #bval site cpgi+promoter
		save(bval_tn,file="bval_tn.RData")
		save(bval_loc_tn,file="bval_loc_tn.RData")
		normal_samples<-grep(normal_pattern,colnames(bval_tn),value=T)
		t_sam<-setdiff(colnames(bval_tn),normal_samples)
		bval_cpg_t<-bval_cpg_tn[,t_sam]
		bval_t<-bval_tn[,t_sam]
		meta_t<-meta[match(colnames(bval_t),meta$Sample_ID),]
		save(bval_cpg_t,file="bval_cpg_t.RData")  #bval site cpgi+promoter
		save(bval_t,file="bval_t.RData")
	}
	
	bval2<-bval_t[rowMeans(bval_t)>0.2,]
	mval<-log2(bval2/(1-bval2))
	mval[mval==Inf]<-log2((1-0.000001)/(1-(1-0.000001)))
	mval[mval==-Inf]<-log2((0+0.000001)/(1-(0+0.000001)))
	
	print(paste0("NA number: ",sum(is.na(mval))))
	print(paste0("is Inf?: ",max(mval)))
	print(paste0("is -Inf?: ",min(mval)))
	mval.norm<-t(scale(t(mval)))
	if(method=="hc"){
		res = ConsensusClusterPlus(as.matrix(mval.norm),maxK=10,reps=1000,pItem=0.8,pFeature=1,
		 title="CC_png_tcga_100sample_bval0.2mean_mval_hc_scale",clusterAlg="hc",distance="pearson",seed=123,plot="png")
		save(res,file="res_cc_hc.RData")
	}else{
		res = ConsensusClusterPlus(as.matrix(mval.norm),maxK=10,reps=1000,pItem=0.8,pFeature=1,
		 title="CC_png_tcga_100sample_bval0.2mean_mval_pam_scale",clusterAlg="pam",distance="euclidean",seed=123,plot="png")
		save(res,file="res_cc_pam.RData")
	}
	my_t<-res[[3]][["consensusClass"]]
	c1<-cc_mean_mval_func(my_t,mval,1)
	c2<-cc_mean_mval_func(my_t,mval,2)
	c3<-cc_mean_mval_func(my_t,mval,3)
	myorder<-order(c(c1,c2,c3))
	my_t[my_t==myorder[1]]<-"CIMP-low"
	my_t[my_t==myorder[2]]<-"CIMP-int"
	my_t[my_t==myorder[3]]<-"CIMP-high"
	heatmap_meth_rna_func(bval2,my_t,paste0("heatmap_meth_genebval0.2_",method,"_scale_3cluster.pdf"))
	heatmap_meth_rna_func(mval,my_t,paste0("heatmap_meth_genebval0.2_",method,"_scale_3cluster_mval.pdf"))
	meta_t$Classical_subtype<-my_t
	save(meta_t,file="meta_t.RData")
	if(is.normal){
		meta_tn$Classical_subtype<-my_t[ meta_tn$Sample_ID]
		meta_tn$Classical_subtype[is.na(meta_tn$Classical_subtype)]<-"Normal"
		save(meta_tn,file="meta_tn.RData")
	}
	return(bval_t)
}


##DSS-bsseq
seq_dss_pipe_func<-function(bsseq,meta,is.normal=F,normal_pattern=NULL,method='pam',genome='hg19'){
	if(genome=='hg19'){
		cpgi_promoter<-cpgi_promoter_hg19
		cpgi_promoter_file<-cpgi_promoter_hg19_file
	}else if(genome=='hg38'){
		cpgi_promoter<-cpgi_promoter_hg38
		cpgi_promoter_file<-cpgi_promoter_hg38_file
	}
	
	if(is.normal==F){
		bval_loc_t<-bval_dss_func(bsseq,cpgi_promoter,type="loc")
		bval_cpg_t<-bval_dss_func(bsseq,cpgi_promoter,type="cpg")
		bval_t<-bval_dss_func(bsseq,cpgi_promoter,type="gene")
		print(paste0("Number of bval_t gene: ",nrow(bval_t)))
		save(bval_cpg_t,file="bval_cpg_t.RData")  #bval site cpgi+promoter
		save(bval_t,file="bval_t.RData")
		save(bval_loc_t,file="bval_loc_t.RData")
		meta_t<-meta[match(colnames(bval_t),meta$Sample_ID),]
	}else{
		bval_loc_tn<-bval_dss_func(bsseq,cpgi_promoter,type="loc")
		bval_cpg_tn<-bval_dss_func(bsseq,cpgi_promoter,type="cpg")
		bval_tn<-bval_dss_func(bsseq,cpgi_promoter,type="gene")
		print(paste0("Number of bval_tn gene: ",nrow(bval_tn)))
		meta_tn<-meta[match(colnames(bval_tn),meta$Sample_ID),]
		save(bval_cpg_tn,file="bval_cpg_tn.RData")  #bval site cpgi+promoter
		save(bval_tn,file="bval_tn.RData")
		save(bval_loc_tn,file="bval_loc_tn.RData")
		normal_samples<-grep(normal_pattern,colnames(bval_tn),value=T)
		t_sam<-setdiff(colnames(bval_tn),normal_samples)
		bval_loc_t<-bval_loc_tn[,t_sam]
		bval_cpg_t<-bval_cpg_tn[,t_sam]
		bval_t<-bval_tn[,t_sam]
		meta_t<-meta[match(colnames(bval_t),meta$Sample_ID),]
		save(bval_loc_t,file="bval_loc_t.RData")
		save(bval_cpg_t,file="bval_cpg_t.RData")  #bval site cpgi+promoter
		save(bval_t,file="bval_t.RData")
	}
	
	bval_t2<-unique(bval_t)
	bval2<-bval_t2[rowMeans(bval_t2)>0.2,]
	mval<-log2(bval2/(1-bval2))
	mval[mval==Inf]<-log2((1-0.000001)/(1-(1-0.000001)))
	mval[mval==-Inf]<-log2((0+0.000001)/(1-(0+0.000001)))
	
	print(paste0("NA number: ",sum(is.na(mval))))
	print(paste0("is Inf?: ",max(mval)))
	print(paste0("is -Inf?: ",min(mval)))
	mval.norm<-t(scale(t(mval)))
	if(method=="hc"){
		res = ConsensusClusterPlus(as.matrix(mval.norm),maxK=10,reps=1000,pItem=0.8,pFeature=1,
		 title="CC_png_tcga_100sample_bval0.2mean_mval_hc_scale",clusterAlg="hc",distance="pearson",seed=123,plot="png")
		save(res,file="res_cc_hc.RData")
	}else{
		res = ConsensusClusterPlus(as.matrix(mval.norm),maxK=10,reps=1000,pItem=0.8,pFeature=1,
		 title="CC_png_tcga_100sample_bval0.2mean_mval_pam_scale",clusterAlg="pam",distance="euclidean",seed=123,plot="png")
		save(res,file="res_cc_pam.RData")
	}
	my_t<-res[[3]][["consensusClass"]]
	c1<-cc_mean_mval_func(my_t,mval,1)
	c2<-cc_mean_mval_func(my_t,mval,2)
	c3<-cc_mean_mval_func(my_t,mval,3)
	myorder<-order(c(c1,c2,c3))
	my_t[my_t==myorder[1]]<-"CIMP-low"
	my_t[my_t==myorder[2]]<-"CIMP-int"
	my_t[my_t==myorder[3]]<-"CIMP-high"
	heatmap_meth_rna_func(bval2,my_t,paste0("heatmap_meth_genebval0.2_",method,"_scale_3cluster.pdf"))
	heatmap_meth_rna_func(mval,my_t,paste0("heatmap_meth_genebval0.2_",method,"_scale_3cluster_mval.pdf"))
	meta_t$Classical_subtype<-my_t
	save(meta_t,file="meta_t.RData")
	if(is.normal){
		meta_tn$Classical_subtype<-my_t[ meta_tn$Sample_ID]
		meta_tn$Classical_subtype[is.na(meta_tn$Classical_subtype)]<-"Normal"
		save(meta_tn,file="meta_tn.RData")
	}
	return(bval_t)
}

bval_dss_func<-function(bsseq,cpgi_promoter,type="gene"){
	if(type=="loc"){
		pm<-getMeth(bsseq)
		site<-as.data.frame(granges(bsseq))
		rownames(pm)<-paste0(site$seqnames,"__",site$start,"__",site$end)
		return(pm)
	}else{
		pm<-getMeth(bsseq,region=cpgi_promoter,what="perRegion")
		pm2<-as.data.frame(pm)
		na.num<-apply(pm,1,function(x){sum(is.na(x))})
		print(table(na.num))
		if(type=="cpg"){
			cpg<-as.data.frame(cpgi_promoter)
			pm2$site<-paste0(cpg$seqnames,"__",cpg$start,"__",cpg$end)
			pm2<-pm2[na.num!=ncol(pm),]
			if(sum(is.na(pm2))!=0){
				stop("There are some NA in bval!")
			}
			pm2<-pm2[!duplicated(pm2$site),]
			rownames(pm2)<-pm2$site
			pm2<-pm2[,-ncol(pm2)]
			return(pm2)
		}else if(type=="gene"){
			pm2$id<-as.data.frame(cpgi_promoter)$V10
			pm2<-pm2[na.num!=ncol(pm),]
			if(sum(is.na(pm2))!=0){
				stop("There are some NA in bval!")
			}
			pm2<-aggregate(pm2[,1:(ncol(pm2)-1)],by=list(pm2[,"id"]),mean)
			rownames(pm2)<-pm2$Group.1
			pm2<-pm2[,-1]
			print(dim(pm2))
			return(pm2)
		}
	}
}

knn_func<-function(bval_tcga,meta_tcga,bval_my,k=10,meta_cc=T,saves=T,probs=0.75){
	bval_tcga2<-bval_tcga[rowMeans(bval_tcga)>0.2,]
	
	mysd<-apply(bval_my,1,sd)
	index<-which(mysd!=0)
	bval_my<-bval_my[index,]
	
	bval_my2<-bval_my[rowMeans(bval_my)>0.2,]
	mval<-log2(bval_my2/(1-bval_my2))
	mval[mval==Inf]<-log2((1-0.000001)/(1-(1-0.000001)))
	mval[mval==-Inf]<-log2((0+0.000001)/(1-(0+0.000001)))

	inter<-intersect(rownames(bval_my2),rownames(bval_tcga2))
	print(paste0("intersect genes: ",length(inter)))
	bval_tcga3<-bval_tcga2[inter,]
	bval_my3<-bval_my2[inter,]
	bval_tcga3<-t(scale(t(bval_tcga3)))
	bval_my3<-t(scale(t(bval_my3)))
	print(sum(colnames(bval_tcga3)==meta_tcga$Sample_ID))
	set.seed(123)
	myknn<-knn(train=t(bval_tcga3),test=t(bval_my3),cl=meta_tcga$Classical_subtype,k=k,prob = TRUE)
	print(table(myknn))
	names(myknn)<-colnames(bval_my3)
	heatmap_meth_rna_func(unique(bval_my2),myknn,'heatmap_meth_knn_bval.pdf')
	heatmap_meth_rna_func(unique(mval),myknn,'heatmap_meth_knn_mval.pdf')
	prob<-attr(myknn,"prob")
	myknn2<-ifelse(prob<=probs,"Unknown",as.character(myknn))
	print(table(myknn2))
	if(meta_cc){
		load("res_cc_pam.RData")
		cc<-res[[3]][["consensusClass"]]
		print(table(myknn,cc))
		print(table(myknn2,cc))
	}
	names(myknn2)<-colnames(bval_my3)
	heatmap_meth_rna_func(unique(mval),myknn2,'heatmap_meth_knn_mval_probunknown.pdf')
	heatmap_meth_rna_func(unique(bval_my2),myknn2,'heatmap_meth_knn_bval_probunknown.pdf')
	if(saves){
		load("meta_t.RData")
		meta_t$Classical_subtype<-myknn2[meta_t$Sample_ID]
		save(meta_t,file="meta_t.RData")
		if(file.exists("meta_tn.RData")){
			load("meta_tn.RData")
			meta_tn$Classical_subtype<-myknn2[ meta_tn$Sample_ID]
			meta_tn$Classical_subtype[is.na(meta_tn$Classical_subtype)]<-"Normal"
			save(meta_tn,file="meta_tn.RData")
		}
	}
	return(myknn2)
}


knn_optimalk_func<-function(bval_tcga,meta_tcga,bval_my=NULL,is.probs=T,probs=0.75){
	bval_tcga2<-bval_tcga[rowMeans(bval_tcga)>0.2,]
	if(!is.null(bval_my)){
		bval_my2<-bval_my[rowMeans(bval_my)>0.2,]
		inter<-intersect(rownames(bval_my2),rownames(bval_tcga2))
		bval_tcga2<-bval_tcga2[inter,]
	}
	print(dim(bval_tcga2))
	bval_tcga2<-t(scale(t(bval_tcga2)))
	set.seed(123)
	ind <- sample(1:ncol(bval_tcga2), size = 0.7*ncol(bval_tcga2))
	train<-bval_tcga2[,ind]
	test<-bval_tcga2[,-ind]
	meta_tcga_train<-meta_tcga[ind,]
	meta_tcga_test<-meta_tcga[-ind,]
	ac<-c()
	for(i in 1:50){
	 set.seed(123)
	 myknn<-knn(train=t(train),test=t(test),cl=meta_tcga_train$Classical_subtype,k=i,prob = TRUE)
	 if(is.probs){
		prob<-attr(myknn,"prob")
		index<-which(prob>probs)
		myknn2<-myknn[index]
		true<-meta_tcga_test$Classical_subtype[index]
		ac[i]<-mean(myknn2 ==true)
	 }else{
		ac[i]<-mean(myknn ==meta_tcga_test$Classical_subtype)
	 }
	 cat("k=", i, " accuracy=", ac[i], "\n")
	}
	pdf('optimal_k.pdf')
	print(plot(ac, type="b", xlab="K",ylab="Accuracy"))
	dev.off()
	return(ac)
}


###nmf MPs
top_sd_func<-function(mt,top=5000){
	mt<-mt[grep("^ENSG",rownames(mt),invert=T),]
	mt<-mt[rowSums(mt)>0,]
	mysd<-apply(mt,1,sd)
	mysd<-sort(mysd,decreasing=T)
	mt<-mt[names(mysd)[1:top],]
}

nmf_input_meth_func<-function(mt){
	mt<-top_sd_func(mt)
	mt<-mt-rowMeans(mt)   #github CCLE_heterogeneity nmf_programs.R
	mt[mt<0]<-0
	return(mt)
}

nmf_run_func<-function(expr_tumor,dataset){
	w_basis_tumor <- list()
	h_coef_tumor <- list()
		   
	w <- NULL
	h <- NULL
	fit<-list()
	for(j in 2:9) {
		print(j)
		res<-nmf(expr_tumor,rank=j,nrun=50,seed=123)
		n<-list(w_basis=basis(res),h_coef=t(coef(res)))
		fit[[j]]<-res
		colnames(n$w_basis) <- paste0(dataset, "_", j, ".", 1:j)
		colnames(n$h_coef) <- paste0(dataset, "_", j, ".", 1:j)
		w <- cbind(w, n$w_basis)
		h <- cbind(h, n$h_coef)
	}
	nmf_wh<-list(w_basis=w,h_coef=h,fit=fit)
	save(nmf_wh,file="nmf_wh.RData")
	return(nmf_wh)
}


heatmap_meth_rna_func<-function(mt,meth_cluster,filename,meth=NULL,meta_extra=NULL,mt_filter=F){
	if(!is.null(meth)){
		inter<-intersect(rownames(mt),rownames(meth))
		mt<-mt[inter,]
	}
	if(mt_filter){
		mt<-mt[rowMeans(mt)>0.2,]
		
		mysd<-apply(mt,1,sd)
		index<-which(mysd!=0)
		mt<-mt[index,]
		
		mt<-log2(mt/(1-mt))
		mt[mt==Inf]<-log2((1-0.000001)/(1-(1-0.000001)))
		mt[mt==-Inf]<-log2((0+0.000001)/(1-(0+0.000001)))
	}
	
	inter<-intersect(colnames(mt),names(meth_cluster))
	mt<-mt[,inter]
	meth_cluster<-meth_cluster[inter]
	print(sum(colnames(mt)==names(meth_cluster)))
	meth_cluster<-sort(meth_cluster)
	mt<-mt[,match(names(meth_cluster),colnames(mt))]
	
	mysd<-apply(mt,1,sd)
	index<-which(mysd!=0)
	mt<-mt[index,]
	
	clust<-as.factor(meth_cluster)
	if(is.null(meta_extra)){
		annotation_col<-data.frame(Subtype=clust)
	}else{
		meta_extra<-meta_extra[match(names(meth_cluster),meta_extra$Sample_ID),]
		index<-which(colnames(meta_extra)=="Sample_ID")
		meta_extra<-meta_extra[,-index,drop=F]
		annotation_col<-data.frame(Subtype=clust)
		annotation_col<-cbind(annotation_col,meta_extra)
		colnames(annotation_col)[2:ncol(annotation_col)]<-colnames(meta_extra)
	}
	rownames(annotation_col) = colnames(mt)
	print(head(annotation_col))
	
	paletteLength <- 50
	myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
	
	tmp_colors=brewer.pal(n = 7, name = "Set1")
	
	
	ann_colors<-lapply(annotation_col,function(x){
		x<-as.factor(x)
		tmp_colors<-tmp_colors[1:length(levels(x))]
		names(tmp_colors)<-levels(x)
		tmp_colors
	})
	print(head(ann_colors))
	
	print(dim(mt))
	pdf(filename,height=6)
	print(pheatmap(mt,annotation_col=annotation_col,scale="row",
	  cluster_rows = T,treeheight_row=0,cluster_cols = F,annotation_colors = ann_colors,
	  show_rownames = F,show_colnames=T,color=myColor))
	dev.off()
}


del_info_func<-function(meta,group){
	meta[is.na(meta)]<-""
	meta<-meta[meta[,group]!="",]
	meta<-meta[meta[,group]!="N/A",]
	meta<-meta[meta[,group]!="NA",]
	meta<-meta[meta[,group]!="na",]
	meta<-meta[meta[,group]!="Unknown",]
	meta<-meta[meta[,group]!="unknown",]
	meta<-meta[meta[,group]!="Not Reported",]
	meta<-meta[meta[,group]!="not reported",]
	meta<-meta[meta[,group]!="[Not Available]",]
	meta<-meta[meta[,group]!="not profiled",]
	meta<-meta[meta[,group]!=".",]
	meta
}



