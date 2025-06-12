tcga_maf_func<-function(obj,type="LUAD",is.normal=T,gene_type="map"){
	if(exists("mt")){
		rm(mt)
	}
	if(expr_type=="count"){
		mt<-assays(obj)$unstrand
	}else if(expr_type=="tpm"){
		mt<-assays(obj)$tpm_unstrand
	}
	
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


tcga_meta_merge<-function(meta,meta_t,type){
	inter<-intersect(colnames(meta),colnames(meta_t))
	meta_o<-meta[,inter]
	meta_o$barcode<-meta$Tumor_Sample_Barcode
	meta_o$sample<-substr(meta$Tumor_Sample_Barcode,1,16)
	meta_o$TN<-"T"
	meta_o$type<-type
	meta_o$Age<-meta$age_at_initial_pathologic_diagnosis
	meta_o$Stage<-gsub("Stage ","",meta$pathologic_stage)
	meta_o$ajcc_pathologic_t<-gsub("Stage ","",meta$pathologic_T)
	meta_o$ajcc_pathologic_n<-gsub("Stage ","",meta$pathologic_stage)
	meta_o$ajcc_pathologic_m<-gsub("Stage ","",meta$pathologic_stage)
	meta_o$Race<-meta$race
	meta_o$Gender<-meta$gender
	meta_o$Pack_years<-meta$number_pack_years_smoked
	meta_o$Smoking_status<-meta$tobacco_smoking_history
	meta_o$OS_status<-meta$vital_status
	clin<-meta[,c("vital_status", "days_to_death", "days_to_last_followup")]
	notDead <- is.na(clin$days_to_death)
	if (any(notDead == TRUE)) {
		clin[notDead, "days_to_death"] <- clin[notDead, "days_to_last_followup"]
	}
	meta_o$days_to_last_follow_up<-meta$days_to_last_followup
	meta_o$OS<-clin$days_to_death
	meta_o$Early_metastasis<-"Primary"
	meta_o$TNM<-paste0(meta_o$ajcc_pathologic_t,meta_o$ajcc_pathologic_n,meta_o$ajcc_pathologic_m)

	meta_o[meta_o=="[Not Available]"]<-""
	meta_o[meta_o=="[Not Applicable]"]<-""
	meta_o[meta_o=="[Unknown]"]<-""
	mut_sam<-meta_o$barcode
	mut_sam<-substr(mut_sam,1,16)
	mut_sam<-gsub("-",".",mut_sam)
	meta_o$Sample_ID<-mut_sam
	meta_o$Patient<-gsub("\\.","-",substr(mut_sam,1,12))

	xena_meta<-as.data.frame(fread(paste0("../RNA-seq/xena/TCGA-",type,".GDC_phenotype_XENA.tsv")))
	if(type=="LUAD"){
		xena_meta_extra<-xena_meta[,c(2:5,7,9,18,23,28:29,33:36,39,46:53,55,111,1,90,94)]
		colnames(xena_meta_extra)[which(colnames(xena_meta_extra)=="primary_diagnosis.diagnoses")]<-"Histology_1"
		colnames(xena_meta_extra)[which(colnames(xena_meta_extra)=="site_of_resection_or_biopsy.diagnoses")]<-"Site_Of_Finding"
	}else{
		xena_meta_extra<-xena_meta[,c(2:5,7,8,15:16,18,20:21,23,28:29,31,33:36,39,41,46:53,55,58:60,111,1,90,94)] #lusc
		colnames(xena_meta_extra)[which(colnames(xena_meta_extra)=="primary_diagnosis.diagnoses")]<-"Histology_1"
		colnames(xena_meta_extra)[which(colnames(xena_meta_extra)=="site_of_resection_or_biopsy.diagnoses")]<-"Site_Of_Finding"
	}
	cols<-setdiff(colnames(xena_meta_extra),colnames(meta_o))
	xena_meta_f<-xena_meta_extra[,cols]
	meta_o<-merge(meta_o,xena_meta_f,by.x="sample",by.y="submitter_id.samples",all.x=T)
	xena_clin<-as.data.frame(fread(paste0("../RNA-seq/xena/survival_",type,"_survival.txt")))
	xena_clin<-xena_clin[,-c(1,3,4,11)]
	colnames(xena_clin)<-c("patient","DSS_status","DSS","DFS_status","DFS","PFS_status","PFS")
	xena_clin2<-xena_clin[ !duplicated(xena_clin$patient),]
	meta_o<-merge(meta_o,xena_clin2,by.x="Patient",by.y="patient",all.x=T)

	print(setdiff(colnames(meta_t),colnames(meta_o)))
	
	other_sam<-setdiff(meta_o$Sample_ID,meta_t$Sample_ID)
	other_col<-setdiff(colnames(meta_t),colnames(meta_o))
	for(i in other_col){
		meta_o[,i]<-meta_t[match(meta_o$Sample_ID,meta_t$Sample_ID),i]
	}
	meta_o<-meta_o[,match(colnames(meta_t),colnames(meta_o))]
	
	meta_rnamut<-rbind(meta_t,meta_o[meta_o$Sample_ID%in%other_sam,])
	save(meta_rnamut,file="../meta_rnamut.RData")
	write.table(meta_rnamut,"../meta_rnamut.txt",row.names=F,quote=F,sep="\t")
	
	meta_t<-meta_o
	save(meta_t,file="meta_t.RData")
	write.table(meta_t,"meta_t.txt",row.names=F,quote=F,sep="\t")
	
	return(list(meta_mut=meta_o,meta_rnamut=meta_rnamut))
}


mut_func<-function(mut,meta,dataset,cancertype,onco_mutmap_index=NULL,incomplete=F,filter_meta=T){
	if(cancertype=="LUAD"){
		genes<-c("TP53","KRAS","KEAP1","EGFR","STK11","SMARCA4",
		  "RBM10","RB1","NF1","ARID1A","BRAF","ERBB2","SETD2","MGA",
		  "FTSJD1","MET","ATM","CDKN2A","U2AF1","RIT1","DOT1L","ARID2",
		  "SMAD4","PTPRU","CTNNB1","ARHGEF12","APC","KLHL5","PIK3CA","PPP3CA",
		  "ATF7IP","KARS","RAF1","MLL3","FANCM","STIM1","NRAS","MAP2K1")
	}else if(cancertype=="LUSC"){
		genes<-c("TP53","CDKN2A","NFE2L2","PTEN","MLL2","RB1",
		  "FAT1","NOTCH1","RASA1","NF1","ARID1A","KDM6A","PIK3CA","CUL3",
		  "HRAS","IRF6","FBXW7","ARHGAP35","PASK","NSD1")
	}else{
		genes<-c("TP53","RB1","CRACD","KIAA1211","COL22A1","RGS7","FPR1",
		  "EP300","CREBBP","ASPM","ALMS1","PDE4DIP","XRN1","PTGFRN","TP73","FMN2",
		  "NOTCH1","NOTCH2","NOTCH2NL","NOTCH3","NOTCH4","EPHB1","KMT2D","MLL2",
		  "CNTNAP2","ZFHX3",  "PTEN","MYCL","MYC","PIK3CA","RICTOR","CCNE1",
		  "FGF10","CDKN2A","FGFR1","ARID1A","EGFR","SOX2","KRAS","NF1","ZNF703")
	}
	mycol<-c("Sample_ID","Hugo_Symbol","HGVSp","Variant_Classification","Chromosome","Start_Position",
		  "End_Position","Reference_Allele","Tumor_Seq_Allele2","Variant_Type")
	if(incomplete==F){
		onco<-mut[,onco_mutmap_index]
	}else{
		inter<-intersect(colnames(mut),mycol)
		onco<-as.data.frame(matrix(ncol=length(mycol),nrow=nrow(mut)))
		colnames(onco)<-mycol
		onco[,inter]<-mut[,inter]
		onco[is.na(onco)]<-"-"
	}
	colnames(onco)<-mycol
	onco$HGVSp[onco$HGVSp=="."]<-""
	
	#meta
	sam_gene<-onco[,1:2]
	sam_gene<-sam_gene[sam_gene$Hugo_Symbol%in%genes,]
	sam_gene%>%group_by(Sample_ID)%>%mutate(mutgene=paste(sort(unique(Hugo_Symbol)),collapse=","))%>%as.data.frame()->sam_gene
	sam_gene$Classical_subtype<-sam_gene$mutgene
	sam_gene$Classical_subtype[grep(",",sam_gene$Classical_subtype)]<-"Multiple"
	gene_freq<-sort(table(sam_gene$mutgene),decreasing=T)
	gene_sel<-names(gene_freq)
	print(gene_sel[grep(",",gene_sel,invert=T)][1:20])
	gene_sel<-as.character(na.omit(gene_sel[grep(",",gene_sel,invert=T)][1:10]))
	others<-setdiff(genes,gene_sel)
	sam_gene$Classical_subtype[sam_gene$Classical_subtype%in%others]<-"Others"
	inter<-intersect(meta$Sample_ID,onco$Sample_ID)
	if(filter_meta){
		meta_t<-meta[meta$Sample_ID%in%inter,]
	}else{
		meta_t<-meta
	}
	onco<-onco[onco$Sample_ID%in%inter,]
	meta_t$mutgene<-sam_gene[match(meta_t$Sample_ID,sam_gene$Sample_ID),"mutgene"]
	meta_t$mutgene[is.na(meta_t$mutgene)]<-"Others"
	meta_t$Classical_subtype<-sam_gene[match(meta_t$Sample_ID,sam_gene$Sample_ID),"Classical_subtype"]
	meta_t$Classical_subtype[is.na(meta_t$Classical_subtype)]<-"Others"
	save(meta_t,file="meta_t.RData")
	write.table(meta_t,"meta_t.txt",row.names=F,sep="\t",quote=F)
	
	#mut
	onco$datasets<-dataset
	onco$cancertype<-cancertype
	mut_maf<-onco
	mut_maf$oncoprint_type<-mut_maf$Variant_Classification
	print(table(mut_maf$Variant_Classification))
		del_class<-c("RNA","Silent","3'Flank","3'UTR","5'UTR","5'Flank","Intron","IGR")  
		myclass<-setdiff(unique(mut_maf$Variant_Classification),del_class)
		mut_maf<-mut_maf[mut_maf$Variant_Classification%in%myclass,]
		mut_maf[ mut_maf$oncoprint_type=="Missense_Mutation","oncoprint_type"]<-"MISSENSE"
		mut_maf[ mut_maf$oncoprint_type=="Nonsense_Mutation","oncoprint_type"]<-"TRUNC"
		mut_maf[ grep("Frame_Shift",mut_maf$oncoprint_type),"oncoprint_type"]<-"TRUNC"
		mut_maf[ mut_maf$oncoprint_type=="In_Frame_Del","oncoprint_type"]<-"INFRAME"
		mut_maf[ mut_maf$oncoprint_type=="In_Frame_Ins","oncoprint_type"]<-"INFRAME"
		mut_maf[ mut_maf$oncoprint_type=="In_Frame_Mut","oncoprint_type"]<-"INFRAME"
		mut_maf[ grep("In_Frame|In-Frame|InFrame",mut_maf$oncoprint_type,ignore.case=T),"oncoprint_type"]<-"INFRAME"
		mut_maf[ grep("splice",mut_maf$oncoprint_type,ignore.case=T),"oncoprint_type"]<-"SPLICE"
		mut_maf[ grep("fusion",mut_maf$oncoprint_type,ignore.case=T),"oncoprint_type"]<-"FUSION"
		othermuts<-setdiff(unique(mut_maf$Variant_Classification),c("MISSENSE","TRUNC","INFRAME","SPLICE","PROMOTER","FUSION","Frame_Shift_Ins","Frame_Shift_Del"))
		mut_maf[ mut_maf$oncoprint_type%in%othermuts,"oncoprint_type"]<-"OTHER"
	colnames(mut_maf)[1]<-"Tumor_Sample_Barcode"
	save(mut_maf,file="mut_maf.RData")
	write.table(mut_maf,"mut_maf.txt",row.names=F,sep="\t",quote=F)
	return(list(onco=mut_maf,meta=meta_t))
}

mut_count_func<-function(mut,del_extra=NULL,type=1){
	if(type==1){
		print(table(mut$Variant_Classification))
		del_class<-c("UTR","Flank","Intron","Silent","RNA","IGR")
		if(!is.null(del_extra)){
			del_class<-union(del_class,del_extra)
		}
		for(i in del_class){
			mut<-mut[grep(i,mut$Variant_Classification,invert=T),]
		}
		mut2<-dcast(mut,Tumor_Sample_Barcode~Hugo_Symbol)
		mut2[,2:ncol(mut2)][mut2[,2:ncol(mut2)]>1]<-1
		mut2[is.na(mut2)]<-0
		load("meta_t.RData")
		other_sam<-setdiff(meta_t$Sample_ID,mut$Tumor_Sample_Barcode)
		if(length(other_sam)!=0){
			for(i in other_sam){
				mut2[nrow(mut2)+1,]<-0
				mut2$Tumor_Sample_Barcode[nrow(mut2)]<-i
			}
		}
		write.table(mut2,file="mut_count.txt",quote=F,sep="\t",row.names=F)
		rownames(mut2)<-mut2$Tumor_Sample_Barcode
		mut2<-mut2[,-1]
	}else{
		mut2<-mut
		sam<-mut2[,1]
		mut2[!is.na(mut2)]<-1
		mut2[is.na(mut2)]<-0
		mut2[,1]<-sam
		write.table(mut2,file="mut_count.txt",quote=F,sep="\t",row.names=F)
		rownames(mut2)<-mut2[,1]
		mut2<-mut2[,-1]
	}
	return(mut2)
}



mut_pathway_func<-function(mut,pathdb,dataset,cancertype){
	all_num<-as.data.frame(table(pathdb$Pathway))
	mut_list<-lapply(unique(mut$Tumor_Sample_Barcode),function(x){
		mut2<-mut[ mut$Tumor_Sample_Barcode==x,]
		inter<-intersect(pathdb$Gene,mut2$Hugo_Symbol)
		tmp<-pathdb
		tmp[ tmp$Gene%in%inter,"type"]<-"y"
		tmp2<-tmp[tmp$type%in%"y",]
		num<-as.data.frame(table(tmp2$Pathway))
		all_num$mut_freq<-num[match(all_num$Var1,num$Var1),"Freq"]
		all_num$ratio<-all_num$mut_freq/all_num$Freq
		all_num$ratio[is.na(all_num$ratio)]<-0
		tmp2%>%group_by(Pathway)%>%mutate(genes=paste(unique(Gene),collapse=","))->mygenes
		mygenes<-as.data.frame(unique(mygenes[,c("Pathway","genes")]))
		all_num$mut_gene<-mygenes[match(all_num$Var1,mygenes$Pathway),"genes"]
		all_num$sam<-x
		all_num
	})
	mut_path<-do.call(rbind,mut_list)
	colnames(mut_path)<-c("pathway","freq","mut_freq","ratio","mutgene","samples")
	mut_path$datasets<-dataset
	mut_path$cancertype<-cancertype
	mut_path[is.na(mut_path)]<-0
	load("meta_t.RData")
	other_sam<-setdiff(meta_t$Sample_ID,mut$Tumor_Sample_Barcode)
	print(other_sam)
	if(length(other_sam)!=0){
		other_mt<-c()
		for(i in other_sam){
			tmp<-mut_path[1:10,]
			tmp$mut_freq<-0
			tmp$ratio<-0
			tmp$mutgene<-NA
			tmp$samples<-i
			other_mt<-rbind(other_mt,tmp)
		}
		other_mt<-as.data.frame(other_mt)
		mut_path<-rbind(mut_path,other_mt)
	}
	print(head(mut_path))
	write.table(mut_path,file="mut_path.txt",row.names=F,sep="\t",quote=F)
	return(mut_path)
}


mut_nmf_func<-function(mut_maf,study,ref="BSgenome.Hsapiens.UCSC.hg19"){
	mut_maf$Chromosome<-gsub("chr","",mut_maf$Chromosome)
	mut<-read_maf(mut_maf)
	set.seed(123)
	mt_tally <- sig_tally(
		mut,
		ref_genome = ref,
		use_syn = TRUE
	)
	save(mt_tally,file="mt_tally.RData")
	mt_est <- sig_estimate(mt_tally$nmf_matrix,
	  range = 2:9,
	  nrun = 200, # increase this value if you wana a more stable estimation
	  use_random = FALSE, # if TRUE, add results from randomized input
	  cores = 40,
	  verbose = TRUE,
	  keep_nmfObj =T
	)
	save(mt_est,file="mt_est.RData")
    w <- NULL
    h <- NULL
	for(j in names(mt_est$nmfEstimate$fit)){
		w_basis=basis(mt_est$nmfEstimate$fit[[j]])
		h_coef=t(coef(mt_est$nmfEstimate$fit[[j]]))
		colnames(w_basis)<-paste0(study,"_",j,".",1:j)
		colnames(h_coef)<-paste0(study,"_",j,".",1:j)
		w<-cbind(w,w_basis)
		h<-cbind(h,h_coef)
	}
	save(w,file="w.RData")
	save(h,file="h.RData")
}













