#RNA example
diff_path_freq<-function(diff_rna,path_rna){
	#DGP
	diff_f<-diff_rna[ diff_rna$group%in%c("MP_r_subtype","Classical_r_subtype","Tumor/Normal"),]
	path_f<-path_rna[ path_rna$group%in%c("MP_r_subtype","Classical_r_subtype","Tumor/Normal"),]
	diff_sub_sig<-diff_f[ diff_f$log2FoldChange>=1&diff_f$padj<0.05,]

	diff_sub_sig$label<-paste0(diff_sub_sig$groupid,"__",diff_sub_sig$variables,"__",diff_sub_sig$gene)
	diff_sub_sig<-diff_sub_sig[ !duplicated(diff_sub_sig$label),]
	
	diff_sub_sig%>%group_by(cancer,group,variables,gene) %>%
	  summarize(count = n(), .groups = 'drop')%>% arrange(cancer,group,variables, desc(count))%>%as.data.frame()->diff_freq

	#pathway gene
	path_sub_sig<-path_f[ path_f$NES>=1&path_f$padj<0.05,]
	path_sub_sig%>%group_by(cancer,group,variables,pathway) %>%
	  summarize(count = n(), .groups = 'drop')%>% arrange(cancer,group,variables, desc(count))%>%as.data.frame()->path_freq_tmp
	path_freq_tmp%>%group_by(cancer,group,variables) %>%top_n(n=50,wt=count)%>%as.data.frame()->path_freq
	path_freq<-path_freq[ path_freq$count>1,]

	path<-unique(path_freq$pathway)
	path_gene<-c()
	for(i in mysets){
		tmp<-read.gmt(i)
		path_gene<-rbind(path_gene,tmp)
	}
	path_gene<-path_gene[ path_gene$term%in%path,]

	path_gene<-path_gene[ path_gene$gene%in%diff_freq$gene,]
	path_res<-merge(path_freq,path_gene,by.x="pathway",by.y="term")
	path_res$id<-paste0(path_res$cancer,"__",path_res$group,"__",path_res$variables,"__",path_res$gene)
	diff_freq$id<-paste0(diff_freq$cancer,"__",diff_freq$group,"__",diff_freq$variables,"__",diff_freq$gene)
	path_res<-path_res[ path_res$id%in%diff_freq$id,]
	path_res<-path_res[,-which(colnames(path_res)=="id")]
	path_res%>%group_by(cancer,group,variables,gene) %>%
	  summarize(count_path = n(), .groups = 'drop')%>% arrange(cancer,group,variables, desc(count_path))%>%as.data.frame()->path_res_gene

	diff_freq<-diff_freq[,-which(colnames(diff_freq)=="id")]
	target_rna<-merge(diff_freq,path_res_gene,by.x=c("cancer","group","variables","gene"),by.y=c("cancer","group","variables","gene"),all.x=T)
	target_rna$score<-target_rna$count+target_rna$count_path
	target_rna$score[is.na(target_rna$score)]<-target_rna$count[is.na(target_rna$score)]
	target_rna%>%arrange(cancer,group,variables, desc(score))%>%as.data.frame()->target_rna
	return(target_rna)
}


omics_target_score_func<-function(target_rna,gene_id_pdb,drugai_f,file_prefix){
	target_all<-merge(target_rna,gene_id_pdb,by="gene",all.x=T)
	target_all<-merge(target_all,drugai_f,by.x="gene",by.y="Gene Name",all.x=T)

	target_all$count_path[ is.na(target_all$count_path)]<-0
	target_all$Targeted_score[ is.na(target_all$Targeted_score)]<-0

	tmp<-target_all$count
	target_all$count_scale<-(tmp-min(tmp))/(max(tmp)-min(tmp))
	tmp<-target_all$count_path
	target_all$count_path_scale<-(tmp-min(tmp))/(max(tmp)-min(tmp))

	target_all$score_scale<-target_all$count_scale*0.4+target_all$count_path_scale*0.2+target_all$Targeted_score*0.4
	target_all%>%arrange(cancer,group,variables, desc(score_scale))%>%as.data.frame()->target_all

	target_rna_all<-target_all[,-which(colnames(target_all)=="score")]
	target_rna_all$Targeted_score[target_rna_all$Targeted_score==0]<-NA
	target_rna_all[is.na(target_rna_all)]<-""
	colnames(target_rna_all)<-c("Target","Cancer type","Group","Subtype","DEG_count","Pathway_count","Uniprot","PDB","Targeted_score","DEG_count_scale","Pathway_count_scale","Score")
	target_db<-target_rna_all[target_rna_all$Score!=0,]
	write.csv(target_db,paste0(file_prefix,".csv"),row.names=F)
	return(target_db)
}

mut_target_score_func<-function(target_mut,gene_id_pdb,drugai_f,file_prefix){
	target_all<-merge(target_mut,gene_id_pdb,by="gene",all.x=T)
	target_all<-merge(target_all,drugai_f,by.x="gene",by.y="Gene Name",all.x=T)

	target_all$Targeted_score[ is.na(target_all$Targeted_score)]<-0

	target_all$score_scale<-target_all$mean_val*0.5+target_all$Targeted_score*0.5
	target_all%>%arrange(cancer, desc(score_scale))%>%as.data.frame()->target_all

	target_all$Targeted_score[target_all$Targeted_score==0]<-NA
	target_all[is.na(target_all)]<-""
	colnames(target_all)<-c("gene","mean_ratio","cancer","uniprot","pdb","Targeted_score","score_scale")
	target_db<-target_all[target_all$score_scale!=0,]
	target_db$groups<-"Driver_subtype"
	target_db$subtypes<-"Driver_subtype"
	write.csv(target_db,paste0(file_prefix,".csv"),row.names=F)
	return(target_db)
}


