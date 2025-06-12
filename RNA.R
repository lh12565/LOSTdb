mt_meta_create_func<-function(mt,type="LUAD",is.normal=F,normal_pattern=NULL,extra_meta=NULL,tpm_filter=F,normal_index=NULL,gene_ratio=0.9){
	meta<-data.frame(Sample_ID=colnames(mt),stringsAsFactors=F)
	if(!is.normal){
		meta$TN<-"T"
		meta$type<-type
	}else{
		if(is.null(normal_index)){
			meta[grep(normal_pattern,meta$Sample_ID),"TN"]<-"N"
		}else{
			meta[normal_index,"TN"]<-"N"
		}
		meta[is.na(meta$TN),"TN"]<-"T"
		meta$type<-ifelse(meta$TN=="N","Normal",type)
		print(table(meta$type))
		print(table(meta$TN))
	}
	if(is.null(extra_meta)){
		meta<-meta
	}else{
		meta<-meta[!is.na(meta$Sample_ID),]
		meta<-merge(meta,extra_meta,by=c("Sample_ID","type","TN"))
	}
	
	if(is.normal){
		meta_tn<-meta[ meta$type%in%c(type,"Normal"),]
		save(meta_tn,file="meta_tn.RData")
		meta_nona<-na_trans_func(meta_tn,'meta_tn.txt')
		mt_tn<-mt[,match(meta_tn$Sample_ID,colnames(mt))]
		if(tpm_filter){
			mt_tn_nolog<-2^mt_tn-1
			mt_tn<-gene_tpm_filter_func(mt_tn_nolog,cutoff=1,ratio=gene_ratio)
		}
		mt_tn<-gene_filter_func(mt_tn)
		index<-apply(mt_tn,1,function(x){length(unique(x))!=1})  #delete constant values
		mt_tn<-mt_tn[index,]
		save(mt_tn,file="mt_tn.RData")
	}
	meta_t<-meta[ meta$type==type,]
	save(meta_t,file="meta_t.RData")
	meta_nona<-na_trans_func(meta_t,'meta_t.txt')
	mt_t<-mt[,match(meta_t$Sample_ID,colnames(mt))]
	if(tpm_filter){
		mt_t_nolog<-2^mt_t-1
		mt_t<-gene_tpm_filter_func(mt_t_nolog,cutoff=1,ratio=gene_ratio)
	}
	mt_t<-gene_filter_func(mt_t)
	index<-apply(mt_t,1,function(x){length(unique(x))!=1})  #delete constant values
	mt_t<-mt_t[index,]
	save(mt_t,file="mt_t.RData")
	return(mt_t)
}

na_trans_func<-function(meta,filename){
	meta[is.na(meta)]<-""
	write.table(meta,filename,row.names=F,sep="\t",quote=F)
	return(meta)
}


meta_create_func<-function(mt,type="LUAD",normal=NULL){   
	meta<-data.frame(Sample_ID=colnames(mt),stringsAsFactors=F)
	if(is.null(normal)){
		meta$tn<-"T"
		meta$type<-type
	}else{
		meta[grep(normal,meta$Sample_ID),"tn"]<-"N"
		meta[is.na(meta$tn),"tn"]<-"T"
		meta$type<-ifelse(meta$tn=="N","Normal",type)
	}
	save(meta,file="meta.RData")
	return(meta)
}

mt_meta_func<-function(mt,meta,type="LUAD",is.normal=F,normal=NULL){
	inter<-intersect(meta$Sample_ID,colnames(mt))
	meta<-meta[meta$Sample_ID%in%inter,]
	mt<-mt[,inter]
	if(is.normal){
		meta_tn<-meta[ meta$type%in%c(type,normal),]
		save(meta_tn,file="meta_tn.RData")
		mt_tn<-mt[,match(meta_tn$Sample_ID,colnames(mt))]
		mt_tn<-gene_filter_func(mt_tn)
		save(mt_tn,file="mt_tn.RData")
	}
	meta_t<-meta[ meta$type==type,]
	save(meta_t,file="meta_t.RData")
	mt_t<-mt[,match(meta_t$Sample_ID,colnames(mt))]
	mt_t<-gene_filter_func(mt_t)
	save(mt_t,file="mt_t.RData")
	return(mt_t)
}

tcga_mt_func<-function(obj,type="LUAD",expr_type="count",is.normal=T,gene_type="map"){
	if(exists("mt")){
		rm(mt)
	}
	if(expr_type=="count"){
		mt<-assays(obj)$unstrand
	}else if(expr_type=="tpm"){
		mt<-assays(obj)$tpm_unstrand
	}
	
	#delete duplicate samples TCGA-44-2666-01A-01R-A278-07  TCGA-44-2666-01A-01R-0946-07
	mt<-mt[,order(colnames(mt),decreasing=T)]
	tmp<-substr(colnames(mt),1,16)
	mt<-mt[,!duplicated(tmp)]
	
	meta<-as.data.frame(colData(obj))
	print(table(meta$sample_type))  #Solid Tissue Normal
	if(!is.normal){
		meta_t_id<-meta[meta$sample_type!="Solid Tissue Normal","barcode"]
		inter<-intersect(colnames(mt),meta_t_id)
		mt<-mt[,inter]
	}
	
	#delete duplicate patient
	if(is.normal){
		meta_n<-meta[meta$sample_type=="Solid Tissue Normal",]
		meta_n<-meta_n[order(meta_n$barcode),]
		tmp<-substr(meta_n$barcode,1,12)
		meta_n<-meta_n[!duplicated(tmp),]
		meta_n_id<-meta_n$barcode
		mt<-mt[,order(colnames(mt))]
		tmp<-substr(colnames(mt),1,12)
		mt_t<-mt[,!duplicated(tmp)]
		samples_filter<-union(colnames(mt_t),meta_n_id)
		mt<-mt[,samples_filter]
	}else{
		mt<-mt[,order(colnames(mt))]
		tmp<-substr(colnames(mt),1,12)
		mt<-mt[,!duplicated(tmp)]
		samples_filter<-colnames(mt)
	}
	save(samples_filter,file="samples_filter.RData")
	
	tmp<-substr(colnames(mt),1,16)
	colnames(mt)<-gsub("-",".",tmp)
	mt<-as.data.frame(mt)
	if(gene_type=="map"){
		gene_map<-as.data.frame(rowData(obj))
		print(sum(gene_map$gene_id==rownames(mt)))
		mt<-gene_dup_filter(mt,gene_map,5,7)
	}
	
	#meta
	tmp<-substr(meta$barcode,1,16)
	meta$Sample_ID<-gsub("-",".",tmp)
	print(length(intersect(colnames(mt),meta$Sample_ID)))
	meta<-meta[match(colnames(mt),meta$Sample_ID),]
	meta$TN<-ifelse(meta$sample_type=="Solid Tissue Normal","N","T")
	meta$type<-ifelse(meta$TN=="T",type,"Normal")
	
	if(expr_type=="count"){
		count<-mt
		save(count,file="count.RData")
	}else{
		save(mt,file="mt.RData")
	}
	meta_ori<-meta
	save(meta_ori,file="meta_ori.RData")
	return(mt)
}


gene_id_trans_name_func<-function(mt,hg="hg38"){
	if(hg=="hg38"){
		gtf<-read.table("~/ref/allGene.hg38.position")
	}else{
		gtf<-read.table("~/ref/hg19/allGene.v19.position")
	}
	gtf$id<-gsub("\\..*","",gtf$V5)
	inter<-intersect(gtf$id,rownames(mt))
	gene_id<-gtf[match(inter,gtf$id),7:6]
	colnames(gene_id)<-c("ENSEMBL","SYMBOL")
	diffs<-setdiff(rownames(mt),gtf$id)
	genes<-AnnotationDbi::select(org.Hs.eg.db,keys=diffs,keytype="ENSEMBL",column="SYMBOL")
	genes<-na.omit(genes)
	gene_id<-rbind(gene_id,genes)
	mt<-mt[match(gene_id$ENSEMBL,rownames(mt)),]
	mt<-gene_dup_filter(mt,gene_id,1,2)
}

gene_filter_func<-function(mt){
	sam0<-apply(mt,1,function(x){sum(x!=0)})
	index<-sam0>=ncol(mt)*0.3
	mt<-mt[index,]
	return(mt)
}

gene_tpm_filter_func<-function(mt,cutoff=1,ratio=0.9){
	a<-apply(mt,1,function(x){sum(x<cutoff)/ncol(mt)})
	mt<-mt[a<ratio,]
	mt<-log2(mt+1)
	return(mt)
}

counts_norm_func<-function(counts,group=NULL,genes=NULL){
	if(is.null(group)){
		y <- DGEList(counts=counts)
	}else{
		y <- DGEList(counts=counts,group=group)
	}
	##Remove rows conssitently have zero or very low counts
	keep <- filterByExpr(y)
	if(!is.null(genes)){
		genes<-intersect(genes,names(keep))
		keep[genes]<-TRUE
	}
	y <- y[keep,keep.lib.sizes=FALSE]
	##Perform TMM normalization and transfer to CPM (Counts Per Million)
	y <- calcNormFactors(y,method="TMM")
	count_norm=cpm(y)
	count_norm<-as.data.frame(count_norm)
	count_norm<-log2(count_norm+1)
	return(count_norm)
}

top_sd_func<-function(mt,top=5000){
	mt<-mt[rowSums(mt)>0,]
	mysd<-apply(mt,1,sd)
	mysd<-sort(mysd,decreasing=T)
	mt<-mt[names(mysd)[1:top],]
}

heatmap_subtype<-function(mt,centroids=NULL,cluster,filename,prefix="centroids"){
	if(!is.null(centroids)){
		inter<-intersect(rownames(mt),rownames(centroids))
		mt<-mt[inter,]
	}
	
	cluster<-sort(cluster)
	mt<-mt[,match(names(cluster),colnames(mt))]
	clust<-as.factor(cluster)
	annotation_col<-data.frame(Subtype=clust)
	rownames(annotation_col) = colnames(mt)
	
	save(mt,file=paste0("mt_heatmap_",prefix,".RData"))
	
	paletteLength <- 50
	myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
	
	tmp_colors=brewer.pal(n = 7, name = "Set1")
	
	tmp_colors<-tmp_colors[1:length(levels(clust))]
	names(tmp_colors)<-levels(clust)
	
	ann_colors = list(Subtype = tmp_colors)
	
	pdf(filename,height=6)
	print(pheatmap(mt,annotation_col=annotation_col,scale="row",
	  cluster_rows = T,treeheight_row=0,cluster_cols = F,annotation_colors = ann_colors,
	  show_rownames = F,show_colnames=F,color=myColor))
	dev.off()
}

centroids_func<-function(mt,centroids){
	medval<-apply(mt,1,median)
	mt2<-mt-medval 
	inBoth<- Reduce(intersect, list(rownames(mt2),rownames(centroids)))
	centroids_cluster_ori<- classify(mt2[inBoth,],centroids[inBoth,])
	if(ncol(centroids)==3){
		centroids_cluster_ori<-ifelse(centroids_cluster_ori==1,colnames(centroids)[1],ifelse(centroids_cluster_ori==2,colnames(centroids)[2],colnames(centroids)[3]))
	}else{
		centroids_cluster_ori<-ifelse(centroids_cluster_ori==1,colnames(centroids)[1],ifelse(centroids_cluster_ori==2,colnames(centroids)[2],ifelse(centroids_cluster_ori==3,colnames(centroids)[3],colnames(centroids)[4])))
	}
	save(centroids_cluster_ori,file="centroids_cluster_ori.RData")
	print(table(centroids_cluster_ori))
	unk<-sapply(colnames(centroids),function(x){
		sams<-names(centroids_cluster_ori[centroids_cluster_ori==x])
		myp<-sapply(sams,function(y){
				a<-cor.test(mt2[inBoth,y],centroids[inBoth,x])
				p<-a$p.value
			})
		names(myp[myp>=0.05])
	})
	unk<-unique(unlist(unk))
	centroids_cluster<-centroids_cluster_ori
	centroids_cluster[unk]<-"Unknown"
	print(table(centroids_cluster))
	save(centroids_cluster,file="centroids_cluster.RData")
	mt3<-mt[inBoth,]
	res<-heatmap_subtype(mt3,centroids,centroids_cluster,"heatmap_centroids.pdf")
	return(centroids_cluster)
}

cluster_assign_func<-function(meta,clusters,cluster_name,is.normal=F){
	if(!is.normal){
		num<-sum(meta$Sample_ID==names(clusters))
		if(num!=nrow(meta)){
			return(print("Error,meta$Sample_ID!=cluster name"))
		}else{
			meta[,cluster_name]<-clusters
		}
	}else{
		rownames(meta)<-meta$Sample_ID
		meta[ names(clusters),cluster_name]<-clusters
		meta[is.na(meta[,cluster_name]),cluster_name]<-"Normal"
	}
	print(table(meta$type,meta[,cluster_name]))
	return(meta)
}

mp_assign_func<-function(meta,study,mpscore,mpsub,score=F){
	if(score){
		mpscore_f<-mpscore[ mpscore$Study==study,]
		mpscore_f<-mpscore_f[,c(1,5,6)]  #Sample_ID,Meta_program,score
		mpscore_f<-reshape2::dcast(mpscore_f,Sample_ID~Meta_program,value.var="score",fun.aggregate = sum)
		mpscore_f<-mpscore_f[match(meta$Sample_ID,mpscore_f$Sample_ID),]
		if(sum(mpscore_f$Sample_ID==meta$Sample_ID)!=nrow(meta)){
			return(print(paste0("score-Samples not unique! ",study)))
		}
		if(sum(is.na(mpscore_f))!=0){
			return(print(paste0("NA score found! ",study)))
		}
		meta<-cbind(meta,mpscore_f[,-1])
	}
	mpsub_f<-mpsub[mpsub$Study==study,]
	mpsub_f<-mpsub_f[match(meta$Sample_ID,mpsub_f$Sample_ID),]
	if(sum(mpsub_f$Sample_ID==meta$Sample_ID)!=nrow(meta)){
		return(print(paste0("subtype-Samples not unique! ",study)))
	}
	if(sum(is.na(mpsub_f[,8]))!=0){
		return(print(paste0("NA subtype found! ",study)))
	}
	meta$MP_r_subtype<-mpsub_f[,8]
	return(meta)
}

tcga_subtype_func<-function(mt2,centroids_cluster,luad_centroids){
	unk<-sapply(colnames(luad_centroids),function(x){
		sams<-names(centroids_cluster[centroids_cluster==x])
		myp<-sapply(sams,function(y){
				a<-cor.test(mt2[inBoth,y],luad_centroids[inBoth,x])
				p<-a$p.value
			})
		names(myp[myp>=0.05])
	})
	unk<-unique(unlist(unk))
	return(unk)
}

marker_cluster_heatmap_func<-function(mt,genes,filename,scales="none",width=7,height=4,method=NULL){
	inter<-intersect(genes,rownames(mt))
	mt<-mt[inter,]
	
	pdf(filename,width=width,height=height)
	paletteLength <- 50
	myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
	if(is.null(method)){
		print(pheatmap(mt,scale=scales,color=myColor))
	}else{
		print(pheatmap(mt,scale=scales,color=myColor,clustering_distance_cols=method))
	}
	dev.off()
	if(is.null(method)){
		return(pheatmap(mt,scale=scales,color=myColor))
	}else{
		return(pheatmap(mt,scale=scales,color=myColor,clustering_distance_cols=method))
	}
}


marker_heatmap_func<-function(mt,meta,genes,filename,scales="none",add_meta_index=NULL,width=12,height=4){
	inter<-intersect(genes,rownames(mt))
	mt<-mt[inter,]
	
	subtype<-meta[match(colnames(mt),meta$Sample_ID),"Classical_subtype"]
	subtype[is.na(subtype)]<-"Unknown"
	type<-ifelse(subtype%in%c("SCLC-A","SCLC-N"),"NE",ifelse(subtype%in%c("SCLC-P","SCLC-Y"),"Non_NE",subtype))

	ann_colors = list(
		Subtype = c(`SCLC-A`="#000CCC",`SCLC-N`="#666FFF",`SCLC-P`="#CC9909",`SCLC-Y`="#FFCC00",Unknown="#CCCCCC"),
		Type = c(NE="#B0ABFF",Non_NE="#00DAE0",Unknown = "#CCCCCC")
	)
	
	if(!is.null(add_meta_index)){
		annotation_col<-data.frame(Subtype=factor(subtype),Type=factor(type),info=factor(meta[,add_meta_index]))
	}else{
		annotation_col<-data.frame(Subtype=factor(subtype),Type=factor(type))
	}
	rownames(annotation_col) = colnames(mt)
	
	pdf(filename,width=width,height=height)
	if(scales!="none"){
		paletteLength <- 50
		myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
		myplot<-pheatmap(mt,scale=scales,clustering_distance_cols   = 'correlation',
		annotation_colors = ann_colors,annotation_col=annotation_col,color=myColor)
	}else{
		myplot<-pheatmap(mt,clustering_distance_cols   = 'correlation',
		annotation_col=annotation_col,annotation_colors = ann_colors)
	}
	dev.off()
	return(myplot)
}


###nmf MP
nmf_input_func<-function(mt){
	mt<-top_sd_func(mt)
	mt<-mt-rowMeans(mt)   #github CCLE_heterogeneity nmf_programs.R
	mt[mt<0]<-0
	return(mt)
}

nmf_run_func<-function(expr_tumor,dataset){
	w_basis_tumor <- list() # nmf gene scores
	h_coef_tumor <- list() # nmf cell scores
		   
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

cpm_func<-function(count,sam,group_index){
	conditions<-factor(t(sam[,group_index]))
	y <- DGEList(counts=count,group=conditions)
	if(nrow(count)>25000){
		keep2<-rowSums(cpm(y)>1)>=ncol(count)*0.3  #Retain genes with CPM > 1 in at least 30% of samples
		y <- y[keep2, , keep.lib.sizes=FALSE]
	}
	if(nrow(y)<10000){
		stop("After cpm filter, the number of gene was less than 10000")
	}
	if(nrow(y)>25000){
		stop("Too many genes after cpm filter")
	}
	y <- calcNormFactors(y,method="TMM")
	count_norm=cpm(y)
	count_norm<-as.data.frame(count_norm)
	return(count_norm) #no log
}

gene_dup_filter<-function(mt,gene_anno,geneid_index,gene_name_index){
	mt$sum<-rowSums(mt,na.rm=T)
	mt$gene<-gene_anno[match(rownames(mt),gene_anno[,geneid_index]),gene_name_index]
	mt2<-mt[!is.na(mt$gene),]
	mt2<-mt2[order(mt2$sum,decreasing=T),]
	mt2<-mt2[!duplicated(mt2$gene),]    
	rownames(mt2)<-mt2$gene
	index<-which(colnames(mt2)%in%c("sum","gene"))
	mt2<-mt2[,-index]
	
	return(mt2)
}

par_wilcox_func<-function(gene_index,mt,sam,group_index){
	data<-cbind.data.frame(gene=as.numeric(t(mt[gene_index,])),sam)
	p=wilcox.test(gene~sam[,group_index], data)$p.value
}

wilcox_diff<-function(ex2,sam,group_index=5,count=F,con_exp=c("NE","non-NE")){
	cl <- makeCluster(60)
	pvalues<-parLapply(cl,1:nrow(ex2),par_wilcox_func,ex2,sam,group_index)
	stopCluster(cl)
	pvalues<-unlist(pvalues)
	pvalues<-as.numeric(sprintf("%.2e",pvalues))
	fdr=p.adjust(pvalues,method = "fdr")
	fdr<-as.numeric(sprintf("%.2e",fdr))

	dataCon1=ex2[,c(which(sam[,group_index]==con_exp[1]))]
	dataCon2=ex2[,c(which(sam[,group_index]==con_exp[2]))]
	mean_c=rowMeans(dataCon1)
	mean_e=rowMeans(dataCon2)
	
	mean_e[mean_e==0]<-0.0001
	mean_c[mean_c==0]<-0.0001
	foldChanges=log2(mean_e/mean_c)

	# Output results based on FDR threshold
	
	genes<-gsub(".*\\|","",rownames(ex2))
	diff_ex<-data.frame(mean_c=mean_c,mean_e=mean_e,log2FoldChange=foldChanges, pvalue=pvalues, padj=fdr,gene_name=genes,stringsAsFactors=F)
	rownames(diff_ex)=rownames(ex2)

	diff_ex<-diff_ex[order(diff_ex$padj),]
	return(diff_ex)
}

sig_gene_meta<-function(diff,cancer,dataset,group,subtype,cutoff_sig=0.05,cutoff_log=1,sig_type="q"){
	print(dim(diff))
	diff<-diff[!is.na(diff$padj),]
	diff<-diff[!is.na(diff$log2FoldChange),]
	if(sig_type=="q"){
		diff<-diff[diff$padj<cutoff_sig&abs(diff$log2FoldChange)>=cutoff_log,]
	}else{
		diff<-diff[diff$pvalue<cutoff_sig&abs(diff$log2FoldChange)>=cutoff_log,]
	}
	print(dim(diff))
	if(nrow(diff)==0){
		return(NULL)
	}

	diff$cancer<-cancer
	diff$dataset<-dataset
	diff$group<-group
	diff$subtye<-subtype
	diff$logp<-round((-log10(diff$pvalue)),2)
	diff$logfc_abs<-abs(diff$log2FoldChange)
	diff$direct<-ifelse(diff$log2FoldChange<0,"down","up")
	gene1<-rownames(diff[order(diff$logfc_abs,decreasing=T),][1:5,])
	gene2<-rownames(diff[order(diff$logp,decreasing=T),][1:5,])
	diff$highlight<-ifelse(rownames(diff)%in%union(gene1,gene2),"y","n")
	diff$cutoff<-sig_type

	return(diff)
}

gsea_func<-function(res_diff,s.sets,filename_prefix,saves=T){
	res_diff2<-res_diff[ order(res_diff$log2FoldChange,decreasing=T),]
	genes2<-gsub(".*\\|","",res_diff2$gene_name)
	gene_mt<-data.frame(gene=genes2,log2fc=res_diff2$log2FoldChange)
	gene_mt<-gene_mt[!duplicated(genes2),]
	geneList<-gene_mt$log2fc
	names(geneList)<-gene_mt$gene
	gse_msig <- GSEA(geneList, TERM2GENE = s.sets,pvalueCutoff = 1,nPermSimple = 1000)
	if(saves==T){
		save(gse_msig,file=paste0(filename_prefix,".RData"))
	}

	return(gse_msig)
}

gsea_mysql<-function(obj,cancer,dataset,group,subtype,pathtype){
	g<-obj@result
	g<-g[g$pvalue<0.05,]
	if(nrow(g)==0){
		return(NULL)
	}
	g2<-g[,1:9]
	g2$cancer<-cancer
	g2$dataset<-dataset
	g2$group<-group
	g2$subtype<-subtype
	g2$pathtype<-pathtype
	return(g2)
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

limma_func<-function(mt,group){
	group<-ifelse(group=="rest","con","ex")
	group<-factor(group)
	design<-model.matrix(~0+group)
	colnames(design)<-levels(group)
	rownames(design)<-colnames(mt)

	fit<-lmFit(mt,design)
	contrast.matrix<-makeContrasts(ex-con,levels=design)
	print(contrast.matrix)
	fit2 <- contrasts.fit(fit, contrast.matrix)
	fit2 <- eBayes(fit2)
	res<-topTable(fit2, coef=1, adjust="BH",n=Inf)
	return(list(fit2,res))
}

gsva_func<-function(ex,group,s.sets,cancer,dataset,group_data,subtype,pathtype,top=30,p_cutoff=0.01){
	gsva_mat <- gsva(expr=as.matrix(ex), 
				   gset.idx.list=s.sets, 
				   kcdf="Gaussian" ,#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
				   verbose=T,
				   parallel.sz = 120)
	lima_cbct_mt<-gsva_mat

	group<-group
	cbct_gs<-limma_func(lima_cbct_mt,group)

	degs<-cbct_gs[[2]]
	diff_gsva<-degs
	save(diff_gsva,file=paste0("gsva_",dataset,"_",group_data,"_",subtype,"_",pathtype,".RData"))
	top1<-min(nrow(subset(degs,logFC>0)),top)
	top2<-min(nrow(subset(degs,logFC<0)),top)
	Diff <- rbind(subset(degs,logFC>0)[1:top1,], subset(degs,logFC<0)[1:top2,])
	dat_plot <- data.frame(id  = row.names(Diff),
						   p   = Diff$P.Value,
						   lgfc= Diff$logFC)
	dat_plot$neg_pos <- ifelse(dat_plot$lgfc>0 ,1,-1)
	dat_plot$lg_p <- -log10(dat_plot$p)*dat_plot$neg_pos
	dat_plot$id[1:10]
	dat_plot$threshold <- factor(ifelse(abs(dat_plot$p) <= p_cutoff,
									   ifelse(dat_plot$lgfc >0 ,'Up','Down'),'Not'),
								levels=c('Up','Down','Not'))

	dat_plot <- dat_plot %>% arrange(lg_p)
	dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
	
	low1 <- dat_plot %>% filter(lg_p < log10(p_cutoff)) %>% nrow()
	low0 <- dat_plot %>% filter(lg_p < 0) %>% nrow()
	high0 <- dat_plot %>% filter(lg_p < -log10(p_cutoff)) %>% nrow()
	high1 <- nrow(dat_plot)
	
	library(ggthemes)
	library(ggprism)
	p <- ggplot(data = dat_plot,aes(x = id, y = lg_p, 
									fill = threshold)) +
	  geom_col()+
	  coord_flip() + 
	  scale_fill_manual(values = c('Up'= '#36638a','Not'='#cccccc','Down'='#7bcd7b')) +
	  geom_hline(yintercept = c(-log10(p_cutoff),log10(p_cutoff)),color = 'white',size = 0.5,lty='dashed') +
	  xlab('') + 
	  ylab('-log10(P.Value) of GSVA score') + 
	guides(fill="none")+
	  theme_prism(border = T) +
	  theme(
		axis.text.y = element_blank(),
		axis.ticks.y = element_blank()
	  ) +
	  geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
				  hjust = 0,color = 'black')
	a<-low1+high0
	if(low1!=top){
		p<-p+geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
				 hjust = 0,color = 'grey')
	}
	if(high0!=top){
		p<-p+geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
				 hjust = 1,color = 'grey')
	}
	p<-p+geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
				hjust = 1,color = 'black')
	ggsave(paste0("gsva_",dataset,"_",group_data,"_",subtype,"_",pathtype,".pdf"),p,width = 15,height  = 15)
	dat_plot$cancer<-cancer
	dat_plot$dataset<-dataset
	dat_plot$group<-group_data
	dat_plot$subtype<-subtype
	dat_plot$pathtype<-pathtype
	dat_plot<-dat_plot[!is.na(dat_plot$id),]
	dat_plot<-dat_plot[dat_plot$id!="NA",]
	return(dat_plot)
}

merge_list_func<-function(mylist,index,prefix){
	tmp<-lapply(mylist,function(x){
		x[[index]]
	})
	tmp<-do.call(rbind,tmp)
}

rna_run_func<-function(mt,meta,s.sets,cancer,dataset,group,method="wilcox",count=F,cutoff_sig=0.05,cutoff_log=1,sig_type="q",gsva_mt=NULL,sample_cutoff=5){
	meta<-meta[meta[,group]!="Unknown",]
	cluster_num<-table(meta[,group])
	subtypes<-names(cluster_num)[cluster_num>=sample_cutoff]
	if(length(subtypes)==0){
		print("After filter v1[5 sample of each subtype], there is no subtype. Terminated")
		return(NULL)
	}
	
	subtypes2<-c()
	for(i in subtypes){
		othersub<-setdiff(meta[,group],i)
		meta_othersub<-meta[meta[,group]%in%othersub,]
		if(nrow(meta_othersub)>=sample_cutoff){
			subtypes2<-c(subtypes2,i)
		}
	}
	if(length(subtypes2)<1){
		print("After filter v2, there is no subtype. Terminated")
		return(NULL)
	}
	
	mt<-mt[,match(meta$Sample_ID,colnames(mt))]
	meta$group<-meta[,group]
	print(dim(mt))
	print(dim(meta))
	index<-which(colnames(meta)=="group")
	clusters<-as.character(subtypes2)
	print(clusters)
	
	
	if(count==T&method=="wilcox"){
		mt<-cpm_func(mt,meta,index)
		print(paste0("after cpm filter: ",nrow(mt)))
	}else if(count==F&method=="wilcox"){
		mt<-2^mt-1
	}
	if(!is.null(gsva_mt)){
		gsva_mt<-gsva_mt[,match(meta$Sample_ID,colnames(gsva_mt))]
		if(sum(is.na(gsva_mt))!=0){
			stop("Error, NA!!")
		}
	}
	if(method=="wilcox"){
		run_all<-lapply(clusters,function(x){
			print(x)
			meta$group<-ifelse(meta[,group]==x,x,"rest")
			print("diff")
			diff_all<-wilcox_diff(mt,meta,index,count=count,con_exp=c("rest",as.character(x)))
			diff_sig<-sig_gene_meta(diff_all,cancer,dataset,group,x,cutoff_sig=cutoff_sig,cutoff_log=cutoff_log,sig_type=sig_type)
			save(diff_all,file=paste0("diff_",group,"_",x,"_all.RData"))
			gsea_alldb<-c()
			gsva_alldb<-c()
			for(i in names(s.sets)){
				print(paste0("gsea:",i))
				gsea_sets = read.gmt(s.sets[[i]])
				gse_msig<-gsea_func(diff_all,gsea_sets,paste0("gsea_",dataset,"_",group,"_",x,"_",i))
				gse_db<-gsea_mysql(gse_msig,cancer,dataset,group,x,i)
				gsea_alldb<-rbind(gsea_alldb,gse_db)
				print(paste0("gsva:",i))
				gsva_sets = gmt2list(s.sets[[i]])
				if(is.null(gsva_mt)){
					mt<-log2(mt+1)
					dat_plot<-gsva_func(mt,meta$group,gsva_sets,cancer,dataset,group,x,i,p_cutoff=cutoff_sig)
					gsva_alldb<-rbind(gsva_alldb,dat_plot)
				}else{
					dat_plot<-gsva_func(gsva_mt,meta$group,gsva_sets,cancer,dataset,group,x,i,p_cutoff=cutoff_sig)
					gsva_alldb<-rbind(gsva_alldb,dat_plot)
				}
			}
			list(diff_all=diff_all,diff_sig=diff_sig,gsea_db=gsea_alldb,gsva_db=gsva_alldb)
		})
	}else{
		run_all<-lapply(clusters,function(x){
			meta$group<-ifelse(meta[,group]==x,x,"rest")
			print("diff")
			diff_all<-mydeseq(mt,meta,index,count=count,con_exp=c("rest",as.character(x)))
			diff_sig<-sig_gene_meta(diff_all,cancer,dataset,group,x,cutoff_sig=cutoff_sig,cutoff_log=cutoff_log,sig_type=sig_type)
			save(diff_all,file=paste0("diff_",group,"_",x,"_all.RData"))
			gsea_alldb<-c()
			gsva_alldb<-c()
			for(i in names(s.sets)){
				print(paste0("gsea:",i))
				gsea_sets = read.gmt(s.sets[[i]])
				gse_msig<-gsea_func(diff_all,gsea_sets,paste0("gsea_",dataset,"_",group,"_",x,"_",i))
				gse_db<-gsea_mysql(gse_msig,cancer,dataset,group,x,i)
				gsea_alldb<-rbind(gsea_alldb,gse_db)
				print(paste0("gsva:",i))
				gsva_sets = gmt2list(s.sets[[i]])
				dat_plot<-gsva_func(gsva_mt,meta$group,gsva_sets,cancer,dataset,group,x,i,p_cutoff=cutoff_sig)
				gsva_alldb<-rbind(gsva_alldb,dat_plot)
			}
			list(diff_all=diff_all,diff_sig=diff_sig,gsea_db=gsea_alldb,gsva_db=gsva_alldb)
		})
	}
	names(run_all)<-clusters
	diff_sig_all<-merge_list_func(run_all,2)
	gsea_sig_all<-merge_list_func(run_all,3)
	gsva_sig_all<-merge_list_func(run_all,4)
	db_tmp<-list(diff_sig_all=diff_sig_all,gsea_sig_all=gsea_sig_all,gsva_sig_all=gsva_sig_all)
	return(list(mysplit=run_all,db=db_tmp))
}



rna_run_nogsva_func<-function(mt,meta,s.sets,cancer,dataset,group,method="wilcox",count=F,cutoff_sig=0.05,cutoff_log=0,sig_type="p",sample_cutoff=5,group_myname){
	meta<-meta[!is.na(meta[,group]),]
	cluster_num<-table(meta[,group])
	subtypes<-names(cluster_num)[cluster_num>=sample_cutoff]
	if(length(subtypes)==0){
		print("After filter v1[5 sample of each subtype], there is no subtype. Terminated")
		return(NULL)
	}
	
	subtypes2<-c()
	for(i in subtypes){
		othersub<-setdiff(meta[,group],i)
		meta_othersub<-meta[meta[,group]%in%othersub,]
		if(nrow(meta_othersub)>=sample_cutoff){
			subtypes2<-c(subtypes2,i)
		}
	}
	if(length(subtypes2)<1){
		print("After filter v2, there is no subtype. Terminated")
		return(NULL)
	}
	
	mt<-mt[,match(meta$Sample_ID,colnames(mt))]
	meta$group<-meta[,group]
	print(dim(mt))
	print(dim(meta))
	index<-which(colnames(meta)=="group")
	clusters<-as.character(subtypes2)
	print(clusters)
	
	
	if(count==T&method=="wilcox"){
		mt<-cpm_func(mt,meta,index)
		print(paste0("after cpm filter: ",nrow(mt)))
	}else if(count==F&method=="wilcox"){
		mt<-2^mt-1
	}
	if(method=="wilcox"){
		run_all<-lapply(clusters,function(x){
			print(x)
			meta$group<-ifelse(meta[,group]==x,x,"rest")
			print("diff")
			diff_all<-wilcox_diff(mt,meta,index,count=count,con_exp=c("rest",as.character(x)))
			diff_sig<-sig_gene_meta(diff_all,cancer,dataset,group_myname,x,cutoff_sig=cutoff_sig,cutoff_log=cutoff_log,sig_type=sig_type)
			gsea_alldb<-c()
			for(i in names(s.sets)){
				print(paste0("gsea:",i))
				gsea_sets = read.gmt(s.sets[[i]])
				gse_msig<-gsea_func(diff_all,gsea_sets,paste0("gsea_",dataset,"_",group_myname,"_",x,"_",i),saves=F)
				gse_db<-gsea_mysql(gse_msig,cancer,dataset,group_myname,x,i)
				gsea_alldb<-rbind(gsea_alldb,gse_db)
			}
			list(diff_all=diff_all,diff_sig=diff_sig,gsea_db=gsea_alldb)
		})
	}else{
		stop('only wilcox!')
	}
	names(run_all)<-clusters
	diff_sig_all<-merge_list_func(run_all,2)
	gsea_sig_all<-merge_list_func(run_all,3)
	db_tmp<-list(diff_sig_all=diff_sig_all,gsea_sig_all=gsea_sig_all)
	return(db_tmp)
}


#https://github.com/tiroshlab/3ca
robust_nmf_programs <- function(nmf_programs, intra_min = 35, intra_max = 10, inter_filter=T, inter_min = 10) {
  
  # Select NMF programs based on the minimum overlap with other NMF programs from the same cell line
  intra_intersect <- lapply(nmf_programs, function(z) apply(z, 2, function(x) apply(z, 2, function(y) length(intersect(x,y))))) 
  intra_intersect_max <- lapply(intra_intersect, function(x) apply(x, 2, function(y) sort(y, decreasing = T)[2]))             
  nmf_sel <- lapply(names(nmf_programs), function(x) nmf_programs[[x]][,intra_intersect_max[[x]]>=intra_min]) 
  names(nmf_sel) <- names(nmf_programs)
  
  # Select NMF programs based on i) the maximum overlap with other NMF programs from the same cell line and
  # ii) the minimum overlap with programs from another cell line
  nmf_sel_unlist <- do.call(cbind, nmf_sel)
  inter_intersect <- apply(nmf_sel_unlist , 2, function(x) apply(nmf_sel_unlist , 2, function(y) length(intersect(x,y))))
  
  final_filter <- NULL 
  for(i in names(nmf_sel)) {
    a <- inter_intersect[grep(i, colnames(inter_intersect), invert = T),grep(i, colnames(inter_intersect))]
    b <- sort(apply(a, 2, max), decreasing = T)
    if(inter_filter==T) b <- b[b>=inter_min]
    if(length(b) > 1) {
      c <- names(b[1]) 
      for(y in 2:length(b)) {
        if(max(inter_intersect[c,names(b[y])]) <= intra_max) c <- c(c,names(b[y])) # selects programs iteratively from top-down. Only selects programs that have a intersection smaller than 10 with a previously selected programs
      }
      final_filter <- c(final_filter, c)
    } else {
      final_filter <- c(final_filter, names(b))
    }
  }
  return(final_filter)                                                      
}

custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))



meta_mut_func<-function(meta,index){
	for(i in index){
		name<-colnames(meta)[i]
		name_f<-gsub("_type","",name)
		meta[,name_f]<-ifelse(is.na(meta[,i]),"n","y")
	}
	return(meta)
}


pca_func<-function(mt,batch){
	pca<-prcomp(t(mt))
	ggplot(,aes(x=pca$x[,1],y=pca$x[,2]))+geom_point(aes(color=factor(batch)))+labs(x="PC1",y="PC2",title="PCA",color="batch")+theme(plot.title = element_text(hjust = 0.5))
	ggsave("pca.pdf")
	return(pca)
}

check_gene_func<-function(mt){
	pattern <- "^[A-Za-z0-9_.-]+$"
	grep(pattern,rownames(mt),value=T,invert=T)
}

