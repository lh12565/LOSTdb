seu_list_input<-function(seu_t,seqtype="10X",prefix,filter_50=T,cpm=T,is.mylist=F){
	if(!is.mylist){
		seu_t_list<-SplitObject(seu_t,split.by = "Sample_ID")
	}else{
		seu_t_list<-seu_t
	}
	seu_t_list<-lapply(seu_t_list,function(x){
		meta<-x@meta.data
		if(seqtype=="10X"){
			cells<-rownames(meta[meta$nFeature_RNA>=1000,])
		}else{
			cells<-rownames(meta[meta$nFeature_RNA>=2000,])
		}
		x<-subset(x,cells=cells)
	})
	subc_list2<-list()
	if(filter_50){
		for(i in names(seu_t_list)){
			if(ncol(seu_t_list[[i]])>=50){
				subc_list2[[i]]<-seu_t_list[[i]]
			}
		}
	}else{
		subc_list2<-seu_t_list
	}
	if(cpm){
		expr_tumor<-lapply(subc_list2,function(x){
			expr<-as.matrix(x@assays$RNA@counts)
			cpm <- t(t(expr)/colSums(expr))*1000000
		})
	}else{
		expr_tumor<-lapply(subc_list2,function(x){
			expr<-as.matrix(x@assays$RNA@counts)
		})
	}
	names(expr_tumor)<-paste0(prefix,"_",names(expr_tumor))
	return(expr_tumor)
}

r3ca_func<-function(counts,cells,genes){
	mt<-readMM(counts)
	cells<-read.csv(cells)
	genes<-read.table(genes)
	rownames(mt)<-genes$V1
	colnames(mt)<-cells$cell_name
	seu<-CreateSeuratObject(mt)
	meta<-seu@meta.data
	print(sum(rownames(meta)==cells$cell_name))
	meta<-cbind(meta,cells)
	seu@meta.data<-meta
	return(seu)
}

seurat_marker_png<-function(obj,genes,filename,type="feature",height=5,order=F){
	DefaultAssay(obj)<-"RNA"
	if(type=="feature"){
		png(filename)
		print(FeaturePlot(obj, features = genes,reduction="umap",order=order))
		dev.off()
	}else if(type=="vln"){
		png(filename,height=height)
		print(VlnPlot(obj,feature=genes,pt.size=0,ncol=1,legend="none"))
		dev.off()
	}
}

read10x_func<-function(name){
	mt_list<-list()
	for(i in name){
		setwd(i)
		sam_name<-gsub(".*_","",i)
		mt_list[[sam_name]]<-Read10X(".")
		setwd("..")
	}
	return(mt_list)
}

seurat_pip_single_func<-function(mt,name,dims=20,max_nfeature=7500,is.seurat=F,rerun_seurat=F,read_filter=T,anno=NULL,type="umi",meta_add=NULL,add_feature=NULL){
	if(is.seurat){
		seu<-mt
		if(rerun_seurat==T){
			meta<-mt@meta.data
			tmp<-mt@assays$RNA@counts
			seu<-CreateSeuratObject(counts=tmp)
			meta2<-seu@meta.data
			meta<-cbind(meta2,meta)
			seu@meta.data<-meta
		}
	}else{
		if(read_filter){
			seu<-CreateSeuratObject(counts = mt, project = name, min.cells = 3, min.features = 200)
		}else{
			seu<-CreateSeuratObject(counts=mt)
		}
	}
	if(!is.null(meta_add)){
		tmp<-seu@meta.data
		meta_add2<-meta_add[match(rownames(tmp),rownames(meta_add)),]
		tmp<-cbind(tmp,meta_add2)
		seu@meta.data<-tmp
	}
	seu$Sample_ID<-name
	seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
	pdf(paste0('vlnplot_qc_',name,'.pdf'))
	print(VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
	dev.off()
	pdf(paste0('dotplot_',name,'qc.pdf'))
	plot1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
	plot2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
	print(plot1 + plot2)
	dev.off()
	seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < max_nfeature & percent.mt < 20)
	if(type=="umi"){
		seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
		seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
	}else if(type=="tpm"){
		seu@assays$RNA@data<-as(as.matrix(log2((seu@assays$RNA@counts/10) + 1)), "dgCMatrix")
		seu <- FindVariableFeatures(seu, selection.method = "mvp", nfeatures = 2000)
	}else{
		stop(paste0("Unknown parameter: ",type))
	}
	if(!is.null(add_feature)){
		genes<-rownames(seu@assays$RNA@counts)
		genes<-intersect(genes,add_feature)
		VariableFeatures(object = seu)<-union(VariableFeatures(object = seu),genes)
	}
	all.genes <- rownames(seu)
	seu <- ScaleData(seu, features = all.genes)
	seu <- RunPCA(seu, features = VariableFeatures(object = seu))
	pdf(paste0('ElbowPlot_',name,'.pdf'))
	print(ElbowPlot(seu,ndims=50))
	dev.off()
	seu <- FindNeighbors(seu, dims = 1:dims)
	seu <- FindClusters(seu, resolution = 0.3)
	seu <- RunUMAP(seu, dims = 1:dims)
	pdf(paste0('umap_',name,'.pdf'))
	print(DimPlot(seu, reduction = "umap"))
	dev.off()
	if(!is.null(anno)){
		for(items in anno){
			Idents(seu)<-items
			pdf(paste0('umap_anno_',name,"_",items,'.pdf'))
			print(DimPlot(seu, reduction = "umap"))
			dev.off()
		}
	}
	markers<-c("CD3D","PTPRC","EPCAM","ASCL1","NEUROD1","POU2F3","YAP1","NCAM1","INSM1","SYP","CHGA")
	genes<-intersect(rownames(seu),markers)
	for(i in genes){
		seurat_marker_png(seu,i, paste0("umap_marker_",name,"_",i,".png"))
	}
	return(seu)
}


seurat_interg_func<-function(seu,dims=30,anno=NULL,extra_marker=NULL,type="umi"){
	if(type=="umi"){
		seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
		VariableFeatures(seu) <- split(row.names(seu@meta.data), seu@meta.data[,"Sample_ID"]) %>% lapply(function(cells_use) {
			seu[,cells_use] %>%
				FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
				VariableFeatures()
		}) %>% unlist %>% unique
	}else if(type=="tpm"){
		seu@assays$RNA@data<-as(as.matrix(log2(seu@assays$RNA@counts + 1)), "dgCMatrix")
		VariableFeatures(seu) <- split(row.names(seu@meta.data), seu@meta.data[,"Sample_ID"]) %>% lapply(function(cells_use) {
			seu[,cells_use] %>%
				FindVariableFeatures(selection.method = "mvp", nfeatures = 2000) %>% 
				VariableFeatures()
		}) %>% unlist %>% unique
	}else{
		stop(paste0("Unknown parameter: ",type))
	}
	seu <- seu %>% 
		ScaleData(verbose = FALSE) %>% 
		RunPCA(features = VariableFeatures(seu), npcs = 50, verbose = FALSE)
	seu <- RunHarmony(seu, group.by.vars = "Sample_ID",plot_convergence=T)
	seu <- RunUMAP(seu, reduction = "harmony", dims = 1:dims)
	seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:dims)
	seu <- FindClusters(seu,resolution=0.1)
	if(!is.null(anno)){
		for(items in anno){
			Idents(seu)<-items
			pdf(paste0('umap_anno__',items,'.pdf'))
			print(DimPlot(seu, reduction = "umap"))
			dev.off()
		}
	}
	markers<-c("CD3D","PTPRC","EPCAM","ASCL1","NEUROD1","POU2F3","YAP1","NCAM1","INSM1","SYP","CHGA")
	markers<-union(markers,extra_marker)
	genes<-intersect(rownames(seu),markers)
	for(i in genes){
		seurat_marker_png(seu,i, paste0("umap_marker_",i,".png"))
	}
	return(seu)
}


seurat_pip_single_nomt_func<-function(mt,name,dims=20,max_nfeature=7500,is.seurat=F){
	if(is.seurat){
		seu<-mt
	}else{
		seu<-CreateSeuratObject(counts = mt, project = name, min.cells = 3, min.features = 200)
	}
	DefaultAssay(seu)<-"RNA"
	seu$Sample_ID<-name
	Idents(seu)<-"Sample_ID"
	pdf(paste0('vlnplot_qc_',name,'.pdf'))
	print(VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2))
	dev.off()
	pdf(paste0('dotplot_',name,'qc.pdf'))
	plot2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
	print(plot2)
	dev.off()
	seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < max_nfeature)
	seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
	seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
	all.genes <- rownames(seu)
	seu <- ScaleData(seu, features = all.genes)
	seu <- RunPCA(seu, features = VariableFeatures(object = seu))
	pdf(paste0('ElbowPlot_',name,'.pdf'))
	print(ElbowPlot(seu,ndims=50))
	dev.off()
	seu <- FindNeighbors(seu, dims = 1:dims)
	seu <- FindClusters(seu, resolution = 0.2)
	seu <- RunUMAP(seu, dims = 1:dims)
	pdf(paste0('umap_',name,'.pdf'))
	print(DimPlot(seu, reduction = "umap"))
	dev.off()
	Idents(seu)<-"assigned_cell_type"
	pdf(paste0('umap_',name,'_celltype.pdf'))
	print(DimPlot(seu, reduction = "umap"))
	dev.off()
	markers<-c("CD3D","PTPRC","EPCAM","ASCL1","NEUROD1","POU2F3","YAP1","NCAM1","INSM1","SYP","CHGA")
	genes<-intersect(rownames(seu),markers)
	for(i in genes){
		seurat_marker_png(seu,i, paste0("umap_marker_",name,"_",i,".png"))
	}
	return(seu)
}


seurat_pip_subset_func<-function(seu,malig_col,malig_label,dims=20,type="umi"){
	Idents(seu)<-malig_col
	subc<-subset(seu,idents=malig_label)
	sam<-unique(subc$Sample_ID)
	subc_list<-lapply(sam,function(x){
		Idents(subc)<-"Sample_ID"
		subc_t<-subset(subc,idents=x)
		if(type=="umi"){
			subc_t <- NormalizeData(subc_t, normalization.method = "LogNormalize", scale.factor = 10000)
			subc_t <- FindVariableFeatures(subc_t, selection.method = "vst", nfeatures = 2000)
		}else if(type=="tpm"){
			subc_t@assays$RNA@data<-as(as.matrix(log2(subc_t@assays$RNA@counts + 1)), "dgCMatrix")
			subc_t <- FindVariableFeatures(subc_t, selection.method = "mvp", nfeatures = 2000)
		}else{
			stop(paste0("Unknown parameter: ",type))
		}
		all.genes <- rownames(subc_t)
		subc_t <- ScaleData(subc_t, features = all.genes)
		subc_t <- RunPCA(subc_t, features = VariableFeatures(object = subc_t))
		pdf(paste0('split_tumor/ElbowPlot_',x,'.pdf'))
		print(ElbowPlot(subc_t,ndims=50))
		dev.off()
		subc_t <- FindNeighbors(subc_t, dims = 1:dims)
		subc_t <- FindClusters(subc_t, resolution = 0.2)
		subc_t <- RunUMAP(subc_t, dims = 1:dims)
		pdf(paste0('split_tumor/umap_',x,'.pdf'))
		print(DimPlot(subc_t, reduction = "umap"))
		dev.off()
		markers<-c("CD3E","PTPRC","AGER","CLDN18","CAV1","HOPX","SFTPC","SFTPD","SFTPB","SFTPA1","SFTPA2","ABCA3","PGC","EMP2","LAMP3","SCGB1A1","SCGB3A2",
					  "SCGB3A1","WFDC2","NAPSA","FOXJ1","RFX2","TMEM190","CETN2","PIFO","HYDIN","CFAP299","TPPP3","EFHC1","TM4SF1","MUC5B","KRT17","KRT5","KRT6A",
					  "TP63","NTRK2","SOX2","CHGA","SYP","INSM1","NCAM1","TUBA1A","VIM","SERPINE1","CDH1","CDH2","MALAT1","NEAT1","TOP2A","MKI67","TACSTD2","AGR2",
					  "CD24","MUC1","NKX2-1","KRT7","MSLN","EPCAM","B2M","NPAS3","TENM1","NDRG1","SLC2A1","VEGFA","KLF6","ZFP36L1","RHOB","HSPA6","DUSP1","FOSB")
		markers<-intersect(markers,rownames(subc_t))
		pdf(paste0('split_tumor/dotplot_markers_',x,'.pdf'),width=50)
		print(DotPlot(subc_t,features=markers,group.by="seurat_clusters"))
		dev.off()
		for(i in markers){
			png(paste0('split_tumor/umap_marker_',i,"_",x,".png"))
			print(FeaturePlot(subc_t,i))
			dev.off()
		}
		subc_t
	})
	return(subc_list)
}


copykat_func<-function(seu,tcell_clusters,index="seurat_clusters"){
	library(copykat);library(Seurat)
	meta<-seu@meta.data
	sam<-as.character(seu$Sample_ID[1])
	print(sam)
	tcell_c<-tcell_clusters[[sam]]
	print(tcell_c)
	tcells<-rownames(meta[ meta[,index]%in%tcell_c,])
	counts<-as.matrix(seu@assays$RNA@counts)
	if(length(tcells)<50){
		sc_cnv = copykat(rawmat = counts,ngene.chr = 5,sam.name = sam,
		  output.seg=T,n.cores = 70)
	}else{
		sc_cnv = copykat(rawmat = counts,ngene.chr = 5,sam.name = sam,
	      norm.cell.names=tcells,output.seg=T,n.cores = 70)
	}
	seu$copykat.pred<-sc_cnv$prediction[ match(colnames(seu),sc_cnv$prediction$cell.names),]$copykat.pred
	DefaultAssay(seu)<-"RNA"
	Idents(seu)<-"copykat.pred"
	pdf(paste0("umap_",sam,"copykat.pdf"))
	print(DimPlot(seu))
	dev.off()
	return(seu)
}


copykat_extrat_cell_func<-function(seu_list,epi_list,anno=F){
	subc_t_list<-list()
	for(i in names(seu_list)){
		seu<-seu_list[[i]]
		if(anno){
			seu$seurat_clusters<-seu$assigned_cell_type
		}
		meta<-seu@meta.data
		sam<-as.character(seu$Sample_ID[1])
		print(sam)
		if(sam%in%names(epi_list)){
			if(epi_list[[sam]][1]!="all"){
				epi_c<-epi_list[[sam]]
				print(epi_c)
				meta$seurat_clusters<-as.character(meta$seurat_clusters)
				meta$copykat.pred<-as.character(meta$copykat.pred)
				meta[meta$seurat_clusters%in%epi_c&meta$copykat.pred=="aneuploid","copykat.pred"]<-"aneuploid_t"
				meta[meta$copykat.pred=="aneuploid","copykat.pred"]<-"aneuploid_imm"
				meta[meta$copykat.pred=="aneuploid_t","copykat.pred"]<-"aneuploid"
				print(sum(is.na(meta)))
				if("aneuploid"%in%meta$copykat.pred){
					seu@meta.data<-meta
					Idents(seu)<-"copykat.pred"
					seu<-subset(seu,idents="aneuploid")
					subc_t_list[[sam]]<-seu
				}
			}else if(epi_list[[sam]][1]=="all"){
				Idents(seu)<-"copykat.pred"
				seu<-subset(seu,idents="aneuploid")
				subc_t_list[[sam]]<-seu
			}
		}
	}
	return(subc_t_list)
}

seurat_pip_func<-function(mt,name,dims=20,max_nfeature=7500){
	seu<-CreateSeuratObject(counts = mt, project = "sclc", min.cells = 3, min.features = 200)
	seu_list<-SplitObject(seu,split.by = "Sample_ID")
	seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
	pdf('vlnplot_qc.pdf')
	VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	dev.off()
	pdf('dotplot_qc.pdf')
	plot1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
	plot2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
	plot1 + plot2
	dev.off()
	seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20)
	seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
	seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
	all.genes <- rownames(seu)
	seu <- ScaleData(seu, features = all.genes)
	seu <- RunPCA(seu, features = VariableFeatures(object = seu))
	pdf('ElbowPlot.pdf')
	ElbowPlot(seu,ndims=50)
	dev.off()
	seu <- FindNeighbors(seu, dims = 1:dims)
	seu <- FindClusters(seu, resolution = 0.2)
	seu <- RunUMAP(seu, dims = 1:dims)
	pdf('umap.pdf')
	DimPlot(seu, reduction = "umap")
	dev.off()
	markers<-c("CD3D","PTPRC","EPCAM","ASCL1","NEUROD1","POU2F3","YAP1")
	genes<-intersect(rownames(seu),markers)
	for(i in genes){
		seurat_marker_png(seu,i, paste0("umap_marker_",i,".png"))
	}
	return(seu)
}



mpscore_func<-function(study_map,subc_list,mps,Cluster_list,malig_col,malig_label){
	study_map<-as.data.table(study_map)
	cancer_type<-unique(study_map$cancer_type)
	setkey(study_map, cancer_type)

	study_contrib<-lapply(Cluster_list,function(x){
		gsub("(.*)_.*","\\1",x)
	})
	
	library(matkot)
	mpscore<-c()
	for(ct in cancer_type) {
		
		cat(ct, '\n')
		ctscore<-c()
		mps_ct<-mps
		study_map2<-as.data.frame(study_map)
		studys<-study_map2[study_map2$cancer_type==ct,"study"]
		for(id in studys){
			dataset<-study_map2[study_map2$cancer_type==ct&study_map2$study==id,"dataset"]
			mysam<-gsub(paste0(dataset,"_"),"",id)
			if(!mysam%in%names(subc_list)){
				stop(paste0(mysam," is not in list!"))
			}
			subc<-subc_list[[mysam]]
			if(!malig_label%in%subc@meta.data[,malig_col]){
				next
			}
			subc$Sample_ID<-rownames(subc@meta.data)
			print(id)
			Idents(subc)<-malig_col
			subc<-subset(subc,idents=malig_label)
			if(ncol(subc)<50){
				next
			}
			meta<-as.data.table(subc@meta.data)
			mt_t<-subc@assays$RNA@data
			
			genes_ori<-rownames(mt_t)
			
			mps_ct_genes<-unique(unlist(mps_ct))
			mp50alias_f<-mp50alias[mp50alias$study%in%id&mp50alias$new%in%mps_ct_genes,]
			if(nrow(mp50alias_f)!=0){
				for(i in 1:nrow(mp50alias_f)){
					genes_ori[which(genes_ori==mp50alias_f$ori[i])]<-mp50alias_f$new[i]
				}
				dup_genes<-genes_ori[duplicated(genes_ori)]
				if(length(dup_genes)!=0){
					print(paste0("duplicated ",dup_genes,"_",dataset))
				}
			}
			rownames(mt_t)<-genes_ori
			mps_ct2 <- slapply(mps_ct, function(x) x[x %in% rownames(mt_t)])
			mps_ct2<-mps_ct2[sapply(mps_ct2,function(x) length(x)>=10)]
			
			set.seed(123)
			scores <- lapply(
				names(mps_ct2),
				function(mp) meta[, .(Sample_ID = Sample_ID, Dataset=dataset, Study=id, Omics=ct, Meta_program = mp, score = sig_score(mt_t, mps_ct2[[mp]], nbin = 50, n = 50))]
			) %>% rbindlist
			ctscore<-rbind(ctscore,scores)
		}
		mpscore<-rbind(mpscore,ctscore)
	}
	save(mpscore,file="mpscore.RData")
	mpscore %>% arrange(Study,Sample_ID,desc(score))%>% group_by(Study,Sample_ID) %>% slice_max(n=2,score) %>% mutate(cutoff=max(score)*0.9>min(score)) %>% slice_max(n=1,score)->mpsub
	mpsub<-as.data.frame(mpsub)
	print(table(mpsub$Study,mpsub$cutoff))
	
	mpsub$Assign<-mpsub$Meta_program
	mpsub[mpsub$cutoff==F,"Assign"]<-"Unknown"
	save(mpsub,file="mpsub.RData")
}


classical_anno_func<-function(seu,malig_col,malig_label,integ=T,type="umi",add_feature=NULL){
	Idents(seu)<-malig_col
	subc_t<-subset(seu,idents=malig_label)
	if(type=="umi"){
		subc_t <- NormalizeData(subc_t, normalization.method = "LogNormalize", scale.factor = 10000)
		VariableFeatures(subc_t) <- split(row.names(subc_t@meta.data), subc_t@meta.data[,"Sample_ID"]) %>% lapply(function(cells_use) {
			subc_t[,cells_use] %>%
				FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
				VariableFeatures()
		}) %>% unlist %>% unique
	}else if(type=="tpm"){
		subc_t@assays$RNA@data<-as(as.matrix(log2(subc_t@assays$RNA@counts + 1)), "dgCMatrix")
		VariableFeatures(subc_t) <- split(row.names(subc_t@meta.data), subc_t@meta.data[,"Sample_ID"]) %>% lapply(function(cells_use) {
			subc_t[,cells_use] %>%
				FindVariableFeatures(selection.method = "mvp", nfeatures = 2000) %>% 
				VariableFeatures()
		}) %>% unlist %>% unique
	}else{
		stop(paste0("Unknown parameter: ",type))
	}
	if(!is.null(add_feature)){
		VariableFeatures(subc_t)<-union(VariableFeatures(subc_t),add_feature)
	}
	subc_t <- subc_t %>% 
		ScaleData(verbose = FALSE) %>% 
		RunPCA(features = VariableFeatures(subc_t), npcs = 50, verbose = FALSE)
	if(integ){
		subc_t <- RunHarmony(subc_t, group.by.vars = "Sample_ID",plot_convergence=T)
		subc_t <- RunUMAP(subc_t, reduction = "harmony", dims = 1:30)
		subc_t <- FindNeighbors(subc_t, reduction = "harmony", dims = 1:30)
	}else{
		subc_t <- RunUMAP(subc_t, dims = 1:30)
		subc_t <- FindNeighbors(subc_t, dims = 1:30)
	}
	dir.create("classical_anno")
	for(i in c(0.1,0.2,0.3,0.4,0.5)){
		subc_t<-FindClusters(subc_t,resolution = i)
		pdf(paste0('classical_anno/umap_alltumor_integ_res',i,'.pdf'))
		print(DimPlot(subc_t, group.by = c("Sample_ID", "seurat_clusters"), ncol = 1))
		dev.off()
	}
	markers<-c("CD3E","PTPRC","AGER","CLDN18","CAV1","HOPX","SFTPC","SFTPD","SFTPB","SFTPA1","SFTPA2","ABCA3","PGC","EMP2","LAMP3","SCGB1A1","SCGB3A2",
				  "SCGB3A1","WFDC2","NAPSA","FOXJ1","RFX2","TMEM190","CETN2","PIFO","HYDIN","CFAP299","TPPP3","EFHC1","TM4SF1","MUC5B","KRT17","KRT5","KRT6A",
				  "TP63","NTRK2","SOX2","CHGA","SYP","INSM1","NCAM1","TUBA1A","VIM","SERPINE1","CDH1","CDH2","MALAT1","NEAT1","TOP2A","MKI67","TACSTD2","AGR2",
				  "CD24","MUC1","NKX2-1","KRT7","MSLN","EPCAM","B2M","NPAS3","TENM1","NDRG1","SLC2A1","VEGFA","KLF6","ZFP36L1","RHOB","HSPA6","DUSP1","FOSB",
				  "ASCL1","NEUROD1","POU2F3","YAP1","IFITM3")
	markers<-intersect(markers,rownames(subc_t))
	pdf(paste0('classical_anno/dotplot_markers.pdf'),width=50)
	print(DotPlot(subc_t,features=markers,group.by="seurat_clusters"))
	dev.off()
	for(i in markers){
		png(paste0('classical_anno/umap_marker_',i,".png"))
		print(FeaturePlot(subc_t,i))
		dev.off()
	}
	return(subc_t)
}


diff_enricher_func<-function(diffs,cluster=0,type="c2",top=50){
	diffs<-diffs[ diffs$cluster==cluster&diffs$avg_log2FC>0.25&diffs$p_val_adj<0.05,]
	print(head(diffs,n=20))
	if(nrow(diffs)>=top){
		genes<-diffs[diffs$cluster==cluster,"gene"][1:top]
	}else{
		stop(paste0("no ",top," genes!"))
	}
	s.sets<-read.gmt(mysets2[[type]])
	ms<-enricher(genes,TERM2GENE=s.sets)
}


subc_list_all_func<-function(subc_list,subc_t=NULL,replace_col,dataset=NULL,Classical_col=NULL,Sample_ID_v2_col=NULL,Cell_type_col=NULL){
	sam<-names(subc_list)
	mp_sam<-unique(mpsub$Study)
	mp_sam<-gsub(paste0(dataset,"_"),"",mp_sam)
	sam<-intersect(sam,mp_sam)
	dir.create("allcells")
	subc_list_all<-lapply(sam,function(x){
		subc<-subc_list[[x]]
		meta<-subc@meta.data
		meta$Sample_ID_v2<-meta[,Sample_ID_v2_col]
		meta$Sample_ID<-rownames(meta)
		if(!is.null(Classical_col)){
			classical_index<-which(colnames(meta)==Classical_col)
			colnames(meta)[classical_index]<-"Classical_subtype"
		}
		meta$MP_r_subtype<-as.character(meta[,replace_col])
		meta$Classical_subtype<-as.character(meta[,replace_col])
		colnames(meta)[which(colnames(meta)==Cell_type_col)]<-"Cell_type"
		id<-paste0(dataset,"_",x)
		mpsub_f<-mpsub[mpsub$Study==id,]
		if(nrow(mpsub_f)==0){
			stop("No cell map!")
		}
		inter<-intersect(mpsub_f$Sample_ID,rownames(meta))
		mpsub_f<-mpsub_f[match(inter,mpsub_f$Sample_ID),]
		meta[match(mpsub_f$Sample_ID,rownames(meta)),"MP_r_subtype"]<-mpsub_f$Assign
		if(!is.null(subc_t)){
			inter<-intersect(colnames(subc_t),rownames(meta))
			meta[inter,"Classical_subtype"]<-as.character(subc_t@meta.data[inter,"Classical_subtype"])
		}
		subc@meta.data<-meta
		Idents(subc)<-"MP_r_subtype"
		celltype<-levels(Idents(subc))
		mycolors=hue_pal()(70)
		index<-which(celltype=="MP_r_MT")
		mycolors[index]<-"#D3D3D3"
		pdf(paste0('allcells/umap_',x,"_mp.pdf"),width=12)
		print(DimPlot(subc,reduction='umap',cols=mycolors))
		dev.off()
		pdf(paste0('allcells/umap_',x,"_classical.pdf"),width=10)
		print(DimPlot(subc,reduction='umap',cols=mycolors,group.by="Classical_subtype"))
		dev.off()
		subc
	})
	names(subc_list_all)<-sam
	return(subc_list_all)
}


subc_list_tumor_color_func<-function(subc_list_t,Mincell_per=NULL,score_filter=NULL){
	sam<-names(subc_list_t)
	dir.create("tumor_mp_filter")
	subc_list_t<-lapply(sam,function(x){
		subc_t<-subc_list_t[[x]]
		subc_t$MP_r_subtype_filter<-subc_t$MP_r_subtype
		if(!is.null(score_filter)){
			cells<-mpsub[mpsub$score<score_filter,"Sample_ID"]
			subc_t$MP_score<-mpsub[match(colnames(subc_t),mpsub$Sample_ID),"score"]
			inter<-intersect(cells,colnames(subc_t))
			subc_t@meta.data[inter,"MP_r_subtype_filter"]<-"Low_score"
		}
		Idents(subc_t)<-"MP_r_subtype_filter"
		celltype<-levels(Idents(subc_t))
		index<-which(celltype=="MP_r_MT")
		index2<-which(celltype=="Unknown")
		index3<-which(celltype=="Low_score")
		num<-table(subc_t$MP_r_subtype_filter)
		if(is.null(Mincell_per)){
			less_mp<-names(num[num<10])
		}else{
			less_mp<-names(num[num<ncol(subc_t)*Mincell_per])
		}
		less_index<-which(celltype%in%less_mp)
		
		other_index<-setdiff(1:length(celltype),c(index,index2,index3,less_index))
		
		mycolors<-c()
		mycolors[other_index]=c(brewer.pal(n = 12, "Paired"),brewer.pal(n = 12, "Set3")[c(1:8,10:12)])[1:length(other_index)]
		mycolors[index]<-"#D3D3D3"
		mycolors[index2]<-"#999999"
		if(!is.null(score_filter)){
			mycolors[index3]<-"#DCDCDC"
		}
		mycolors[less_index]<-"#F5F5F5"
		if(is.null(Mincell_per)&is.null(score_filter)){
			pdf(paste0('tumor_mp_filter/umap_',x,'_tumor_mp_mincell10.pdf'),width=12)
		}else{
			pdf(paste0('tumor_mp_filter/umap_',x,'_tumor_mp_',Mincell_per,'_score',score_filter,'.pdf'),width=12)
		}
		if(is.null(score_filter)){
			print(DimPlot(subc_t,reduction='umap',cols=mycolors))
		}else{
			print(DimPlot(subc_t,reduction='umap',cols=mycolors))
		}
		dev.off()
		subc_t
	})
	names(subc_list_t)<-sam
	return(subc_list_t)
}


meta_list_func<-function(subc_list_all,subc_list_t,meta_tmp,index_all,index_t,classical_anno=T){
	sam<-names(subc_list_all)
	meta_list<-lapply(sam,function(x){
		subc<-subc_list_all[[x]]
		meta<-subc@meta.data[,index_all]
		umap<-subc@reductions$umap@cell.embeddings
		meta<-cbind(meta,umap)
		meta<-merge(meta,meta_tmp,by=c("Sample_ID_v2","Patient"))
		meta$TN<-"T"
		if(classical_anno==F){
			meta$Classical_subtype<-NA
		}
		meta
	})
	names(meta_list)<-sam
	meta_all<-do.call(rbind,meta_list)
	save(meta_all,file="meta_all.RData")
	
	meta_t_list<-lapply(sam,function(x){
		subc<-subc_list_t[[x]]
		meta<-subc@meta.data[,index_t]
		umap<-subc@reductions$umap@cell.embeddings
		meta<-cbind(meta,umap)
		meta<-merge(meta,meta_tmp,by=c("Sample_ID_v2","Patient"))
		meta$TN<-"T"
		if(classical_anno==F){
			meta$Classical_subtype<-NA
		}
		meta
	})
	names(meta_t_list)<-sam
	meta_t<-do.call(rbind,meta_t_list)
	save(meta_t,file="meta_t.RData")
}


mt_t_list_func<-function(subc_list_t){
	sam<-names(subc_list_t)
	mt_t_list<-lapply(sam,function(x){
		subc<-subc_list_t[[x]]
		mt<-subc@assays$RNA@data
		mt<-reshape2::melt(as.matrix(mt))
		mt<-mt[ mt$value!=0,]
		colnames(mt)<-c("gene","cells","expression")
		mt
	})
	names(mt_t_list)<-sam
	mt_t<-do.call(rbind,mt_t_list)
	save(mt_t,file="mt_t.RData")
	return(mt_t)
}


heatmap_list_func<-function(subc_list_t,mps=NULL,group="MP_r_subtype",type="nsclc",mincell=10,filters=F){
	sam<-names(subc_list_t)
	heat_list<-lapply(sam,function(x){
		subc<-subc_list_t[[x]]
		if(group=="MP_r_subtype"){
			num<-table(subc@meta.data[,group])
			mymps<-names(num[num>=mincell])
			mymps<-mymps[mymps!="Unknown"]
			if(length(mymps)<=1){
				mymps<-names(num[num>=5])
				mymps<-mymps[mymps!="Unknown"]
			}
			tmps<-mps[mymps]
			inter_list<-lapply(tmps,function(x){intersect(x,rownames(subc))})
			mps_g<-unique(unlist(inter_list))
		}else if(group=="Classical_subtype"){
			if(type=="nsclc"){
				mps_g<-c("AGER","CLDN18","CAV1","HOPX","SFTPC","SFTPD","SFTPB","SFTPA1","SFTPA2","ABCA3","PGC","EMP2","LAMP3","SCGB1A1","SCGB3A2",
				  "SCGB3A1","WFDC2","NAPSA","FOXJ1","RFX2","TMEM190","CFAP299","HYDIN","PIFO","CETN2","TPPP3","EFHC1","TM4SF1","MUC5B","KRT17","KRT5","KRT6A",
				  "TP63","NTRK2","SOX2","CHGA","SYP","INSM1","NCAM1","TUBA1A","VIM","SERPINE1","CDH1","CDH2","MALAT1","NEAT1","TOP2A","MKI67","TACSTD2","AGR2",
				  "CD24","MUC1","NKX2-1","KRT7","MSLN","EPCAM","B2M","NDRG1","SLC2A1","VEGFA","KLF6","ZFP36L1","RHOB","CD3E","PTPRC")  #Classical_subtype
			}else if(type=="sclc"){
				ne1<-c('ASCL1','SYP','UCHL1','CHGA','CHGB','INSM1','TAGLN3','CRMP1','SCG3',
					  'MYT1L','SYN1','SYT11','RUNDC3A','SEZ6','CELF3','ROBO1','PIEZO2')
				none1<-c('HES1','SCGB1A1','TGFBR2','S100A10','LGALS3','EPHA2','CCN1',
				  'MYOF','EMP1','CAV2','AHNAK','IFITM3','CCND1','S100A16')

				ne2<-c('BEX1','ASCL1','INSM1','CHGA','TAGLN3','KIF5C','CRMP1','SCG3','SYT4','RTN1',
				  'MYT1','SYP','KIF1A','TMSB15A','SYN1','SYT11','RUNDC3A','TFF3','CHGB','TLCD3B',
				  'SH3GL2','BSN','SEZ6','TMSB15B','CELF3','FAM57B')  #FAM57B: TLCD3B
				none2<-c('RAB27B','TGFBR2','SLC16A5','S100A10','ITGB4','YAP1','LGALS3','EPHA2','S100A16',
				  'PLAU','ABCC3','ARHGDIB','CCN1','PTGES','CCND1','IFITM2','IFITM3','AHNAK','CAV2',
				  'TACSTD2','TGFBI','EMP1','CAV1','ANXA1','MYOF','CYR61')  #CYR61: CCN1
				ne<-union(ne1,ne2)
				ne<-c(ne,"GRP","NCAM1","CEACAM5","SCG2","DLK1","MS4A8B","NEUROD1","SLC17A6","MYC","NEUROD2","DLL3")
				none<-union(none1,none2)
				none<-union(none,c("POU2F3","AZGP1","ASCL2","POU2AF2","C11orf53","SERPINB1","HLA-A","HLA-C","SOX2",
				  "GBP1","GBP2","GBP3","GBP4","GBP5","HLA-DMB","HLA-DQB1","HLA-DMA","HLA-DRB1","HLA-DPB1","HLA-DPA1",
				  "IFI44L","ICAM1","CD74","HLA-DRA","TRIM22","SERPING1","XAF1","EPST11","IFIT3"))
				mps_g<-union(ne,none)
				luad<-c("AGER","CLDN18","CAV1","HOPX","SFTPC","SFTPD","SFTPB","SFTPA1","SFTPA2","ABCA3","PGC","EMP2","LAMP3","SCGB1A1","SCGB3A2",
				  "SCGB3A1","WFDC2","NAPSA","FOXJ1","RFX2","TMEM190","CFAP299","HYDIN","PIFO","CETN2","TPPP3","EFHC1","TM4SF1","MUC5B","KRT17","KRT5","KRT6A",
				  "TP63","NTRK2","SOX2","CHGA","SYP","INSM1","NCAM1","TUBA1A","VIM","SERPINE1","CDH1","CDH2","MALAT1","NEAT1","TOP2A","MKI67","TACSTD2","AGR2",
				  "CD24","MUC1","NKX2-1","KRT7","MSLN","EPCAM","B2M","NDRG1","SLC2A1","VEGFA","KLF6","ZFP36L1","RHOB","CD3E","PTPRC")
				gene_add1<-c("ACTA2","ELN","MGP","EDIL3","FBLN5","THY1") #kinker_2020 sclc EMT
				gene_add2<-c("SOX4","NEUROD4","SOX9","YBX3","HOXC10","MYCN","MYCL","NOTCH1","HES6","TCF4","MARCKS",
				  "ADCYAP1","NRXN1","SEMA6A","EFNB1","EPHB2","NRP2","HIF1A","FOXO3","COL1A2","TWIST1","ZEB1","ZEB2",
				  "TGFB1","TGFBR1","BMP2","BMPR2","STAT3","IL13RA1","PHLDA1","SMAD3","SMAD2","NFAT5")
				gene_add<-union(gene_add1,gene_add2)
				mps_g<-union(mps_g,gene_add)
				mps_g_ori<-union(mps_g,luad)
				
				inter<-intersect(rownames(subc),mps_g_ori)
				if(filters){
					mt<-subc@assays$RNA@data[inter,]
					num<-table(subc@meta.data[,group])
					mymps<-names(num[num>=mincell])
					mps_g<-c()
					for(i in mymps){
						cells<-rownames(subc@meta.data[ subc$Classical_subtype==i,])
						mt2<-mt[,cells,drop=F]
						num_0<-apply(mt2,1,function(x){sum(x==0)})
						
						tmp_genes<-names(num_0[num_0<=ncol(mt2)*0.7])
						mps_g<-union(mps_g,tmp_genes)
					}
				}else{
					mps_g<-inter
				}
				
			}
		}
		mt<-DoHeatmap(subc,mps_g,group.by=group)+
			  scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 11, name =
				   "RdBu")))(50))+theme(axis.text.y=element_blank())
		mt2<-mt$data
		mt2<-mt2[!is.na(mt2$Expression),]
		mt2<-reshape2::dcast(mt2[,1:3],Feature~Cell)
		rownames(mt2)<-mt2$Feature
		mt2<-mt2[,-1]
		mt2
	})
	names(heat_list)<-sam
	if(group=="MP_r_subtype"){
		heat_list_mp<-heat_list
		save(heat_list_mp,file="heat_list_mp.RData")
	}else if(group=="Classical_subtype"){
		heat_list_class<-heat_list
		save(heat_list_class,file="heat_list_class.RData")
	}
	return(heat_list)
}


ngchm_func<-function(mt,meta,dataname,filename,subtype,scale="none",display_col=NULL,continuous_col=NULL,extra_col=NULL,mincell=10){
	mycol<-c("MP_r_subtype",
	  "Classical_subtype","Histology_1","Histology_2")
	if(sum(is.na(meta$Classical_subtype))==nrow(meta)){
		mycol<-c("MP_r_subtype",
		  "Histology_1","Histology_2")
	}
	if(!is.null(extra_col)){
		mycol<-union(mycol,extra_col)
	}
	mt<-mt[rowSums(mt)!=0,]
	scaleDataLayer <- chmNewDataLayer('Scaled', as.matrix(mt))
	if(subtype=="MP_r_subtype"){
		mychm <- chmNew(dataname,scaleDataLayer,colOrder=colnames(mt),rowOrder=rownames(mt))
	}else if(subtype=="Classical_subtype"){
		mychm <- chmNew(dataname,scaleDataLayer,colOrder=colnames(mt))
	}
	rownames(meta)<-meta$Sample_ID
	inter_meta<-intersect(colnames(meta),mycol)
	meta<-meta[,inter_meta,drop=F]
	no_continuous_col<-setdiff(inter_meta,continuous_col)
	hidden_col<-setdiff(inter_meta,display_col)
	print(paste0("inter_meta: ",inter_meta))
	print(paste0("hidden_col: ",hidden_col))
	print(sum(is.na(hidden_col)))
	print(sum(is.na(inter_meta)))
	print(sum(is.na(no_continuous_col)))
	print(sum(is.na(continuous_col)))
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
			covariateBar <- chmNewCovariate(i,mycol,ColorMap)
			if(i%in%hidden_col){
				mychm <- chmAddCovariateBar(mychm, 'column', covariateBar,display="hidden")
			}else{
				mychm <- chmAddCovariateBar(mychm, 'column', covariateBar)
			}
	}
	for(i in no_continuous_col){
			print("no_continuous_col")
			print(i)
			mycol<-as.character(meta[,i])
			names(mycol)<-rownames(meta)
			if(i=="MP_r_subtype"|i=="Classical_subtype"){
				mpname<-as.character(unique(meta[,i]))
				num<-table(meta[,i])
				print(num)
				less_mp<-names(num[num<mincell])
				print(less_mp)
				setdiff_mp<-setdiff(mpname,less_mp)
				setdiff_mp<-setdiff_mp[setdiff_mp!="Unknown"]
				if(length(setdiff_mp)<=1){
					less_mp<-names(num[num<5])
				}
				less_index<-which(mpname%in%less_mp)
				print(less_index)
				more_index<-setdiff(1:length(mpname),less_index)
				print(more_index)
				mycolors<-c()
				mycolors[more_index]=c(brewer.pal(n = 12, "Paired"),brewer.pal(n = 12, "Set3")[c(1:8,10:12)],"#D01C8B","#F1B6DA")[1:length(more_index)]
				mycolors[less_index]<-"#D3D3D3"
				mutationColorMap <- chmNewColorMap(mpname,mycolors)
				covariateBar <- chmNewCovariate(i,mycol,mutationColorMap)
			}else{
				covariateBar <- chmNewCovariate(i,mycol,type="discrete")
			}
			if(i%in%hidden_col){
				print(i)
				mychm <- chmAddCovariateBar(mychm, 'column', covariateBar,display="hidden")
			}else{
				print(i)
				mychm <- chmAddCovariateBar(mychm, 'column', covariateBar)
			}
	}
	chmExportToFile(mychm,filename)
}

ngchm_data_func<-function(meta_t,mt,sam,genes=NULL,display_col=NULL,continuous_col=NULL,type="nsclc",subtype="Classical_subtype",extra_col=NULL,mincell=10){
	if(!is.null(genes)){
		inter<-intersect(genes,rownames(mt))
		mt<-mt[inter,]
	}
	myorder<-meta_t[order(meta_t[,subtype]),"Sample_ID"]
	mt<-mt[,match(myorder,colnames(mt))]
	index<-which(colnames(meta_t)==subtype)
	other_index<-setdiff(1:ncol(meta_t),index)
	meta_t<-meta_t[,c(index,other_index)]
		
	ngchm_func(mt,meta_t,'NGCHM',paste0("ngchm_",subtype,"_",sam,".ngchm"),subtype=subtype,scale='none',
	  display_col=display_col,continuous_col=continuous_col,extra_col=extra_col,mincell=mincell)
}

ngchm_scrna_run_func<-function(heat_list,meta_t,subtype="MP_r_subtype",type="nsclc",display_col=c("Histology_1","Classical_subtype","MP_r_subtype"),continuous_col=NULL,mincell=10,extra_col=NULL){
	sam<-names(heat_list)
	sapply(sam,function(x){
		mt<-heat_list[[x]]
		meta_t<-meta_t[match(colnames(mt),meta_t$Sample_ID),]
		ngchm_data_func(meta_t,as.matrix(mt),x,display_col=display_col,subtype=subtype,
		  continuous_col=continuous_col,type=type,mincell=mincell,extra_col=extra_col)
	})
}


diff_path_func<-function(subc_list_t,s.sets,cancer,dataset,group,sample_cutoff=20){
	sam<-names(subc_list_t)
	diff_path_all<-lapply(sam,function(x){
		subc<-subc_list_t[[x]]
		meta<-subc@meta.data
		cell_num<-table(subc@meta.data[,group])
		if(length(cell_num)==1){
			return(NULL)
		}
		clusters<-names(cell_num[cell_num>=sample_cutoff])
		if(length(clusters)==0){
			print("After filter v1[5 sample of each subtype], there is no subtype. Terminated")
			return(NULL)
		}
		
		print(clusters)
		
		diffs <- wilcoxauc(subc, group)
		diffs$cancer<-cancer
		diffs$dataset<-dataset
		diffs$study<-paste0(dataset,"_",x)
		
		cl <- makeCluster(20)
		path<-parLapply(cl,clusters,path_func,diffs,s.sets,x,cancer,dataset)
		stopCluster(cl)
		path<-do.call(rbind,path)
		diff_path<-list(diff_all=diffs,gsea_all=path)
		diff_path
	})
	names(diff_path_all)<-sam
	if(group=="MP_r_subtype"){
		diff_path_mp<-diff_path_all
		save(diff_path_mp,file="diff_path_mp.RData")
	}else if(group=="Classical_subtype"){
		diff_path_class<-diff_path_all
		save(diff_path_class,file="diff_path_class.RData")
	}
	return(diff_path_all)
}


path_func<-function(cluster,diffs,s.sets,sam,cancer,dataset){
	library(clusterProfiler);library(tidyverse)
	diff_f<- diffs %>%
	  dplyr::filter(group == cluster) %>%
	  arrange(desc(auc)) %>%
	  dplyr::select(feature, auc)
	ranks<- deframe(diff_f)
	gsea_alldb<-c()
	for(i in names(s.sets)){
		print(paste0("gsea:",i))
		gsea_sets = read.gmt(s.sets[[i]])
		gse <- GSEA(ranks, TERM2GENE = gsea_sets,pvalueCutoff = 1)
		gse<-gse@result
		gse$group<-cluster
		gse$cancer<-cancer
		gse$dataset<-dataset
		gse$study<-paste0(dataset,"_",sam)
		gse$pathtype<-i
		gsea_alldb<-rbind(gsea_alldb,gse)
	}
	gsea_alldb
}


RenameGenesSeurat_v2 <- function(obj,newnames,gene.use=NULL,de.assay="RNA") {
	print("Run this before integration. It only changes obj@assays$*@counts, @data and @scale.data, @var.features,@reductions$pca@feature.loadings")
	lassays <- Assays(obj)
	assay.use <- obj@reductions$pca@assay.used
	DefaultAssay(obj) <- de.assay
	if (is.null(gene.use)) {
		all_genenames <- rownames(obj)
	}else{
		all_genenames <- gene.use
		obj <- subset(obj,features=gene.use)
	}

	order_name <- function(v1,v2,ref){
		v2 <- make.names(v2,unique=T)
		df1 <- data.frame(v1,v2)
		rownames(df1) <- df1$v1
		df1 <- df1[ref,]
		return(df1)
	}

	df1 <- order_name(v1=all_genenames,v2=newnames,ref=rownames(obj))
	all_genenames <- df1$v1
	newnames <- df1$v2

	if ('SCT' %in% lassays) {
		if ('SCTModel.list' %in%  slotNames(obj@assays$SCT)) {
		obj@assays$SCT@SCTModel.list$model1@feature.attributes <- obj@assays$SCT@SCTModel.list$model1@feature.attributes[all_genenames,]
		rownames(obj@assays$SCT@SCTModel.list$model1@feature.attributes) <- newnames
		}
	}
	change_assay <- function(a1=de.assay,obj,newnames=NULL,all_genenames=NULL){
	  RNA <- obj@assays[a1][[1]]
	  if (nrow(RNA) == length(newnames)) {
		if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
		if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
		if (length(RNA@var.features)) {
			df1 <- order_name(v1=all_genenames,v2=newnames,ref=RNA@var.features)
			all_genenames1 <- df1$v1
			newnames1 <- df1$v2
			RNA@var.features            <- newnames1
		}
		if (length(RNA@scale.data)){
			df1 <- order_name(v1=all_genenames,v2=newnames,ref=rownames(RNA@scale.data))
			all_genenames1 <- df1$v1
			newnames1 <- df1$v2
			rownames(RNA@scale.data)    <- newnames1
		}

	  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
	  obj@assays[a1][[1]] <- RNA
	  return(obj)
	}

	for (a in lassays) {
		DefaultAssay(obj) <- a
		df1 <- order_name(v1=all_genenames,v2=newnames,ref=rownames(obj))
		all_genenames1 <- df1$v1
		newnames1 <- df1$v2
		obj <- change_assay(obj=obj,a1=a,newnames=newnames1,all_genenames=all_genenames1)
	}

	hvg <- VariableFeatures(obj,assay=assay.use)
	if (length(obj@reductions$pca)){
		df1 <- order_name(v1=all_genenames,v2=newnames,ref=hvg)
		all_genenames1 <- df1$v1
		newnames1 <- df1$v2
		rownames(obj@reductions$pca@feature.loadings) <- newnames1
	}
	return(obj)
}




