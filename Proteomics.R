filter_impute_func<-function(mt,is.normal=F,normal_pattern=NULL,normal_index=NULL,impute=F,missing_per=0.3){
	if(is.normal){
		if(!is.null(normal_pattern)){
			mt_t<-mt[,grep(normal_pattern,colnames(mt),invert=T)]
			mt_n<-mt[,grep(normal_pattern,colnames(mt))]
		}else{
			mt_n<-mt[,normal_index]
			tcol<-setdiff(colnames(mt),colnames(mt_n))
			mt_t<-mt[,tcol]
		}
		na_n<-apply(mt_n,1,function(x){sum(is.na(x))})
		n_id<-na_n<=ncol(mt_n)*missing_per
		mt_n<-mt_n[n_id,]
		na_t<-apply(mt_t,1,function(x){sum(is.na(x))})
		t_id<-na_t<=ncol(mt_t)*missing_per
		mt_t<-mt_t[t_id,]
		inter_gene<-intersect(rownames(mt_t),rownames(mt_n))
		print(paste0("normal: ",nrow(mt_n)))
		print(paste0("tumor: ",nrow(mt_t)))
		print(paste0("intersect: ",length(inter_gene)))
		mt_t<-mt_t[match(inter_gene,rownames(mt_t)),]
		mt_n<-mt_n[match(inter_gene,rownames(mt_n)),]
		if(impute){
			mt_t<-impute.knn(as.matrix(mt_t))$data
			mt_n<-impute.knn(as.matrix(mt_n))$data
		}
		mt<-as.data.frame(cbind(mt_t,mt_n))
	}else{
		na_t<-apply(mt,1,function(x){sum(is.na(x))})
		t_id<-na_t<=ncol(mt)*missing_per
		mt<-mt[t_id,]
		print(paste0("tumor: ",nrow(mt)))
		if(impute){
			mt<-as.data.frame(impute.knn(as.matrix(mt))$data)
		}else{
			mt<-as.data.frame(mt)
		}
	}
	return(mt)
}


top_sd_pro_func<-function(mt,top=5000){
	mysd<-apply(mt,1,sd)
	mysd<-sort(mysd,decreasing=T)
	mt<-mt[names(mysd)[1:top],]
}


#ConsensusClusterPlus >50
class_run_func<-function(mt){
	mt<-top_sd_pro_func(mt)
	num<-sum(mt<0)
	if(num==0){
		mt <-  sweep(mt,1, apply(mt,1,median,na.rm=T))
	}else{
		min_val <- min(mt, na.rm = TRUE)
		mt <- mt - min_val
		mt <-  sweep(mt,1, apply(mt,1,median,na.rm=T))
	}
	res <-  ConsensusClusterPlus(as.matrix(mt),
									 maxK = 9,
									 reps = 1000,
									 pItem = 0.8,
									 clusterAlg = "km",
									 seed = 123,
									 title="ConsensusCluster/",
									 distance="euclidean",
									 plot="png")
	save(res,file="ConsensusCluster/res.RData")
	return(list(mt=mt,res=res))
}

#knn
knn_optimalk_func<-function(mt_t,meta_tcga,mt_t2=NULL,group,is.probs=T,probs=0.75){
	mt_t<-top_sd_pro_func(mt_t)
	if(!is.null(mt_t2)){
		inter<-intersect(rownames(mt_t2),rownames(mt_t))
		mt_t<-mt_t[inter,]
	}
	print(dim(mt_t))
	mt_t<-t(scale(t(mt_t)))
	set.seed(123)
	ind <- sample(1:ncol(mt_t), size = 0.7*ncol(mt_t))
	train<-mt_t[,ind]
	test<-mt_t[,-ind]
	meta_tcga_train<-meta_tcga[ind,]
	meta_tcga_test<-meta_tcga[-ind,]
	ac<-c()
	for(i in 1:50){
	 set.seed(123)
	 myknn<-knn(train=t(train),test=t(test),cl=meta_tcga_train[,group],k=i,prob = TRUE)
	 if(is.probs){
		prob<-attr(myknn,"prob")
		index<-which(prob>probs)
		myknn2<-myknn[index]
		true<-meta_tcga_test[,group][index]
		ac[i]<-mean(myknn2 ==true)
	 }else{
		ac[i]<-mean(myknn ==meta_tcga_test[,group])
	 }
	 cat("k=", i, " accuracy=", ac[i], "\n")
	}
	pdf('optimal_k.pdf')
	print(plot(ac, type="b", xlab="K",ylab="Accuracy"))
	dev.off()
	return(ac)
}

knn_func<-function(mt_t,meta_tcga,mt_t2,group,k=11,meta_cc=T,saves=T,probs=0.75){
	mt_t<-top_sd_pro_func(mt_t)

	inter<-intersect(rownames(mt_t),rownames(mt_t2))
	print(paste0("intersect genes: ",length(inter)))
	mt_t<-mt_t[inter,]
	mt_t2<-mt_t2[inter,]
	mt_t<-t(scale(t(mt_t)))
	mt_t2<-t(scale(t(mt_t2)))
	print(sum(colnames(mt_t)==meta_tcga$Sample_ID))
	print(nrow(mt_t))
	print(nrow(mt_t2))
	set.seed(123)
	myknn<-knn(train=t(mt_t),test=t(mt_t2),cl=meta_tcga[,group],k=k,prob = TRUE)
	prob<-attr(myknn,"prob")
	myknn2<-ifelse(prob<=probs,"Unknown",as.character(myknn))
	print(table(myknn2))
	if(saves){
		load("meta_t.RData")
		meta_t$Classical_subtype<-as.character(myknn)
		save(meta_t,file="meta_t.RData")
		if(file.exists("meta_tn.RData")){
			load("meta_tn.RData")
			meta_tn$Classical_subtype<-meta_t[match(meta_tn$Sample_ID,meta_t$Sample_ID),"Classical_subtype"]
			meta_tn$Classical_subtype[is.na(meta_tn$Classical_subtype)]<-"Normal"
			save(meta_tn,file="meta_tn.RData")
		}
	}
	return(list(myknn,myknn2))
}


#mp
nmf_input_pro_func<-function(mt,norm="zscore"){
	mt<-top_sd_pro_func(mt)
	if(norm=="zscore"){
		mt<-t(scale(t(mt)))
	}else{
		mt<-mt-rowMeans(mt)
		mt[mt<0]<-0
	}
	return(mt)
}

nmf_run_pro_func<-function(expr_tumor,dataset,is.wsigned=F){
	w_basis_tumor <- list()
	h_coef_tumor <- list()
		   
	w <- NULL
	h <- NULL
	fit<-list()
	tmp<-c()
	for(j in 2:9) {
		print(j)
		if(is.wsigned){
			expr_tumor<-posneg(expr_tumor)
			res<-nmf(expr_tumor,rank=j,nrun=50,seed=123)
			myw<-basis(res)
			wsigned<-apply(myw,2,function(x){
				a<-sapply(rownames(myw)[1:5000],function(y){
					index<-which(rownames(myw)%in%y)
					if(x[index[1]]>x[index[2]]){
						a<-x[index[1]]
					}else if(x[index[1]]<x[index[2]]){
						a<-x[index[2]]*(-1)
					}else{
						a<-NA
					}
				})
				a
			})
			rownames(wsigned)<-gsub("\\..*","",rownames(wsigned))
			print(sum(rownames(wsigned)==tmp))
			tmp<-rownames(wsigned)
			n<-list(w_basis=wsigned,h_coef=t(coef(res)))
		}else{
			res<-nmf(expr_tumor,rank=j,nrun=50,seed=123)
			n<-list(w_basis=basis(res),h_coef=t(coef(res)))
		}
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
	
	meta$MP_p_subtype<-mpsub_f[match(meta$Sample_ID,mpsub_f$Sample_ID),"Assign"]
	meta$MP_p_subtype[is.na(meta$MP_p_subtype)]<-"Unknown"
	return(meta)
}



#TMT normalization:https://pwilmart.github.io/TMT_analysis_examples/multiple_TMT_MQ.html
tmt_norm_func<-function(mt,index_list){
	mycolsum<-sapply(index_list,function(x){
		tmp<-mt[x]
		colSums(tmp)
	})
	target<-mean(unlist(mycolsum))
	if(class(mycolsum)!="list"){
		norm_facs<-target/mycolsum
		exp_sl_list<-lapply(1:length(index_list),function(x){
			sweep(mt[index_list[[x]]], 2, norm_facs[,x], FUN = "*")
		})
	}else{
		norm_facs<-lapply(mycolsum,function(x){
			target/x
		})
		exp_sl_list<-lapply(1:length(index_list),function(x){
			sweep(mt[index_list[[x]]], 2, norm_facs[[x]], FUN = "*")
		})
	}
	data_sl<-do.call(cbind,exp_sl_list)
	sl_tmm <- calcNormFactors(data_sl)
	data_sl_tmm <- sweep(data_sl, 2, sl_tmm, FUN = "/")
	exp_sl_sum<-lapply(exp_sl_list,function(x){
		rowSums(x)
	})
	irs<-do.call(cbind,exp_sl_sum)
	colnames(irs) <- c("sum1", "sum2", "sum3")
	irs.average <- apply(irs, 1, function(x) exp(mean(log(x))))
	irs_fac<-apply(irs,2,function(x){
		irs.average/x
	})
	data_irs_list<-lapply(1:length(exp_sl_list),function(x){
		exp_sl_list[[x]]*irs_fac[,x]
	})
	data_irs<-do.call(cbind,data_irs_list)
	irs_tmm <- calcNormFactors(data_irs)
	data_irs_tmm <- sweep(data_irs, 2, irs_tmm, FUN = "/")
	print(paste0("min: ",min(data_irs_tmm)))
	mt<-log2(data_irs_tmm)
}

MedianNorm<-function(x){
  x/median(x, na.rm=T);
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

merge_list_func<-function(mylist,index,prefix){
	tmp<-lapply(mylist,function(x){
		x[[index]]
	})
	tmp<-do.call(rbind,tmp)
}


pro_run_nogsva_func<-function(mt,meta,s.sets,cancer,dataset,group,method="wilcox",count=F,cutoff_sig=0.05,cutoff_log=0,sig_type="p",sample_cutoff=5,group_myname){
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
		stop("count found!")
		mt<-cpm_func(mt,meta,index)
		print(paste0("after cpm filter: ",nrow(mt)))
	}else if(count==F&method=="wilcox"){
		mt<-2^mt
		if(sum(mt<0,na.rm=T)!=0){
			stop("negtive!!")
		}
	}
	if(method=="wilcox"){
		run_all<-lapply(clusters,function(x){
			print(x)
			meta$group<-ifelse(meta[,group]==x,x,"rest")
			print("diff")
			diff_all<-wilcox_diff(mt,meta,index,count=count,con_exp=c("rest",as.character(x)))
			diff_all<-diff_all[ !is.na(diff_all$log2FoldChange),]
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



ngchm_func<-function(mt,meta,mycol,subtype,genes,dataname,filename,scale="none",display_col,continuous_col,extra_col){
	mycol<-intersect(colnames(meta),mycol)
	inter<-intersect(genes,rownames(mt))
	mt<-mt[match(inter,rownames(mt)),]
	mt<-mt[rowSums(mt)!=0,]
	unadjustedLayer <- chmNewDataLayer('Unadjusted', as.matrix(mt))
	if(scale=='row'){
			mt <- t(scale(t(mt)))
	}
	scaleDataLayer<-chmNewDataLayer('Scaled', mt)
	mychm <- chmNew(dataname,scaleDataLayer,unadjustedLayer,colOrder=colnames(mt))
	rownames(meta)<-meta$Sample_ID
	inter_meta<-intersect(colnames(meta),mycol)
	meta<-meta[,inter_meta,drop=F]
	no_continuous_col<-setdiff(inter_meta,continuous_col)
	hidden_col<-setdiff(inter_meta,display_col)
	index<-which(no_continuous_col%in%subtype)
	other_index<-setdiff(1:length(no_continuous_col),index)
	no_continuous_col<-no_continuous_col[c(index,other_index)]
	for(i in no_continuous_col){
			print("no_continuous_col")
			print(i)
			mycol<-meta[,i]
			names(mycol)<-rownames(meta)
			if(i=="MP_p_subtype"){
				mpname<-as.character(unique(meta$MP_p_subtype))
				mycolors<-c()
				mycolor=c(brewer.pal(n = 12, "Paired"),brewer.pal(n = 12, "Set3")[c(1:8,10:12)])[1:length(mpname)]
				mutationColorMap <- chmNewColorMap(mpname,mycolor)
				covariateBar <- chmNewCovariate(ori_new[[i]],mycol,mutationColorMap)
			}else{
				covariateBar <- chmNewCovariate(ori_new[[i]],mycol,type="discrete")
			}
			if(i%in%hidden_col){
				mychm <- chmAddCovariateBar(mychm, 'column', covariateBar,display="hidden")
			}else{
				mychm <- chmAddCovariateBar(mychm, 'column', covariateBar)
			}
	}
	for(i in continuous_col){
			print("continuous_col")
			print(i)
			mycol<-meta[,i]
			names(mycol)<-rownames(meta)
			meta_min<-min(mycol,na.rm=T)
			meta_max<-max(mycol,na.rm=T)
			meta_median<-median(mycol,na.rm=T)
			if(meta_median==meta_min|meta_median==meta_max){
				meta_median<-mean(mycol,na.rm=T)
			}
			mycolor<-unlist(sample(ngchm_continuous_color,1))
			ColorMap <- chmNewColorMap(c(meta_min,meta_median,meta_max), mycolor)
			covariateBar <- chmNewCovariate(ori_new[[i]],mycol,ColorMap)
			if(i%in%hidden_col){
				mychm <- chmAddCovariateBar(mychm, 'column', covariateBar,display="hidden")
			}else{
				mychm <- chmAddCovariateBar(mychm, 'column', covariateBar)
			}
	}
	chmExportToFile(mychm,filename)
}
ngchm_data_func<-function(meta_t,mycol,display_col,continuous_col,type="LUAD",subtype="Classical_p_subtype",datasets,cancertype,omics,extra_col=NULL){
	myorder<-meta_t[order(meta_t[,subtype]),"Sample_ID"]
	load("../mt_t.RData")
	mt<-mt_t
	mt<-mt[,match(myorder,colnames(mt))]
	if(subtype=="Classical_p_subtype"){
		if(type=="LUAD"){
			genes<-ngchm_pro_gene$luad
		}else if(type=="LUSC"){
			genes<-ngchm_pro_gene$lusc
			
		}else if(type=="SCLC"){
			genes<-ngchm_pro_gene$sclc
		}
	}else{
		mp<-unique(meta_t$MP_p_subtype)
		load("mps.RData")
		inter<-intersect(names(mps),mp)
		genes<-unlist(mps[inter])
	}
	
	ngchm_func(mt,meta_t,mycol,subtype,genes,paste0(datasets," ",cancertype," ",omics),paste0("ngchm_",subtype,".ngchm"),scale='row',
	  display_col=display_col,continuous_col=continuous_col,extra_col)
}










