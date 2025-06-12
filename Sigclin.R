library(dunn.test)
library(dplyr)
library(ggplot2)


# one-sided hypergeometric test
hypergeom_test <- function(k, K, n, N, alternative = "greater") {
  # Fold Enrichment
  fold_enrichment <- (k / n) / (K / N)
  logFC <- log2(fold_enrichment)
  
  # P value
  if (alternative == "greater") {
    p_value <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
  } else if (alternative == "less") {
    p_value <- phyper(k, K, N - K, n, lower.tail = TRUE)
  } else {
    stop("Alternative must be 'greater' or 'less'")
  }
  
  return(list(
    LogFC = logFC,
    P_value = p_value
  ))
}

# Iterate
analyze_enrichment <- function(data,group1,group2,dataset,cancertype,omics,omics_id) {
  print(paste0("enrich","--",group1,"--",group2))
  results <- c()
  
	data<-del_info_func(data,group1)
	mygroup_num<-table(data[,group1])
	mygroup<-names(mygroup_num[mygroup_num>=3])
	data<-data[data[,group1]%in%mygroup,]
	
	data<-del_info_func(data,group2)
	mygroup_num<-table(data[,group2])
	mygroup<-names(mygroup_num[mygroup_num>=3])
	data<-data[data[,group2]%in%mygroup,]
	if(length(unique(data[,group1]))<=1|length(unique(data[,group2]))<=1){
		return(NULL)
	}
	
  protein_subtypes <- unique(data[,group1])
  histologies <- unique(data[,group2])
  
  
  for (subtype in protein_subtypes) {
    for (histology in histologies) {
      k <- sum(data[,group1] == subtype & data[,group2] == histology)  # a
      K <- sum(data[,group1] == subtype)  # a + b
      n <- sum(data[,group2] == histology)  # a + c
      N <- nrow(data)  # a + b + c + d
      
      if (K > 0 && n > 0) {
        enrichment_result <- hypergeom_test(k, K, n, N, alternative = "greater")
		enrichment_result<-as.data.frame(enrichment_result)
		enrichment_result$comparisons<-paste(subtype, histology, sep = "--")
		enrichment_result$group1_var<-subtype
		enrichment_result$group2_var<-histology
		enrichment_result<-add_info_func(enrichment_result,omics_id,dataset,cancertype,omics,group1,group2)
		results<-rbind(results,enrichment_result)
      }
    }
  }
  
  return(results)
}


# Kruskal-Wallis   Dunn's
kruskal_dunn_test <- function(data, group_col, value_col,dataset,cancertype,omics,omics_id) {
  #filter
  print(paste0("krus","--",group_col,"--",value_col))
	data[,value_col]<-as.numeric(data[,value_col])
	if(sum(is.na(data[,value_col]))==nrow(data)){
		return(list(
			krus=NULL,
			dunn=NULL
		))
	}
	data<-del_info_func(data,group_col)
	data<-del_info_func(data,value_col)
	mygroup_num<-table(data[,group_col])
	mygroup<-names(mygroup_num[mygroup_num>=3])
	data<-data[data[,group_col]%in%mygroup,]
	if(length(unique(data[,group_col]))<=1){
		return(list(
			krus=NULL,
			dunn=NULL
		))
	}


  # Kruskal-Wallis
  kruskal_test <- kruskal.test(data[[value_col]] ~ data[[group_col]])
  krus<-kruskal_func(kruskal_test)
  krus<-add_info_func(krus,omics_id,dataset,cancertype,omics,group_col,value_col)
  
  # Dunn's
  if (krus$p.value!="NaN"&krus$p.value <=1 ) {
    dunn_test <- dunn.test(data[[value_col]], data[[group_col]], method = "bonferroni")
    dunn<-as.data.frame(dunn_test)
    dunn<-add_info_func(dunn,omics_id,dataset,cancertype,omics,group_col,value_col)
  } else {
    dunn <- NULL
  }
  
  return(list(
	krus=krus,
	dunn=dunn
  ))
}

# ANOVA
anova_test <- function(data, group_col, value_col,dataset,cancertype,omics,omics_id) {
  #filter
  print(paste0("aov","--",group_col,"--",value_col))
	data[,value_col]<-as.numeric(data[,value_col])
	if(sum(is.na(data[,value_col]))==nrow(data)){
		return(list(
			krus=NULL,
			dunn=NULL
		))
	}
	data<-del_info_func(data,group_col)
	data<-del_info_func(data,value_col)
	mygroup_num<-table(data[,group_col])
	mygroup<-names(mygroup_num[mygroup_num>=3])
	data<-data[data[,group_col]%in%mygroup,]
	if(length(unique(data[,group_col]))<=1){
		return(list(
			krus=NULL,
			dunn=NULL
		))
	}

  mygroup<-as.factor(data[[group_col]])
  anova_result <- aov(data[[value_col]] ~ mygroup, data = data)
  myaov<-aov_func(anova_result)
  myaov<-add_info_func(myaov,omics_id,dataset,cancertype,omics,group_col,value_col)
  anova_p_value <- summary(anova_result)[[1]]$`Pr(>F)`[1]
  
  # Tukey's HSD
  if (anova_p_value!="NaN"&anova_p_value <= 1) {
    tukey_result <- TukeyHSD(anova_result)
    tukey_results <- tukey_result$`mygroup`
	tukey_results<-as.data.frame(tukey_results)
	tukey_results$comparisons<-rownames(tukey_results)
	tukey_results<-add_info_func(tukey_results,omics_id,dataset,cancertype,omics,group_col,value_col)
  } else {
    tukey_results <- NULL
  }
  
  return(list(
    myaov = myaov,
    Tukey_Results = tukey_results
  ))
}

# Spearman
continuous_correlation_test <- function(data, var1, var2, method = "pearson",dataset,cancertype,omics,omics_id) {
  #filter
  print(paste0("correlation","--",var1,"--",var2))
	data[,var1]<-as.numeric(data[,var1])
	data[,var2]<-as.numeric(data[,var2])
	if(sum(is.na(data[,var1]))==nrow(data)|sum(is.na(data[,var2]))==nrow(data)){
		return(NULL)
	}
	
	data_tmp1<-data[!is.na(data[,var1]),]
	data_tmp2<-data_tmp1[!is.na(data_tmp1[,var2]),]
	if(nrow(data_tmp2)<3){
		return(NULL)
	}

  if (method == "pearson") {
    cor_test <- cor.test(data[[var1]], data[[var2]], method = "pearson")
  } else if (method == "spearman") {
    cor_test <- cor.test(data[[var1]], data[[var2]], method = "spearman")
  } else {
    stop("Method must be 'pearson' or 'spearman'")
  }
  
  res<-data.frame(Correlation=cor_test$estimate,P_value=cor_test$p.value,stringsAsFactors=F)
  res<-add_info_func(res,omics_id,dataset,cancertype,omics,var1,var2)
  return(res)
}

kruskal_func<-function(kruskal_result){
	kruskal_df <- data.frame(
	  Statistic = kruskal_result$statistic,
	  df = kruskal_result$parameter,
	  p.value = kruskal_result$p.value
	)
}

aov_func<-function(aov){
	myaov<-as.data.frame(summary(aov)[[1]][1,])
}

add_info_func<-function(meta,omics_id,dataset,cancertype,omics,group1,group2){
	meta$dataset<-dataset
	meta$cancertype<-cancertype
	meta$omics<-omics
	meta$group1<-group1
	meta$group2<-group2
	meta$group1_change<-omics_id[ omics_id$datasets==dataset&omics_id$omics==omics&omics_id$cancertype==cancertype&omics_id$col_names==group1,"col_change"]
	meta$group2_change<-omics_id[ omics_id$datasets==dataset&omics_id$omics==omics&omics_id$cancertype==cancertype&omics_id$col_names==group2,"col_change"]
	meta
}

##cell
add_info_func<-function(meta,scrna_id_display,studyid,cancertype,omics,group1,group2){
	meta$dataset<-studyid
	meta$cancertype<-cancertype
	meta$omics<-omics
	meta$group1<-group1
	meta$group2<-group2
	meta$group1_change<-scrna_id_display[ scrna_id_display$study==studyid&scrna_id_display$omics==omics&scrna_id_display$cancertype==cancertype&scrna_id_display$col_names==group1,"col_change"]
	meta$group2_change<-scrna_id_display[ scrna_id_display$study==studyid&scrna_id_display$omics==omics&scrna_id_display$cancertype==cancertype&scrna_id_display$col_names==group2,"col_change"]
	meta
}


del_info_func<-function(meta,group){
	meta<-meta[!is.na(meta[,group]),]
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
	meta<-meta[meta[,group]!="NA (lost-to-follow-up)",]
	meta<-meta[meta[,group]!="NA (resection)",]
	meta<-meta[meta[,group]!="NA (therapy declined; lost-to-follow-up)",]
	meta<-meta[meta[,group]!=".",]
	meta
}


# main
analyze_clinical_data <- function(data,omics_id,dataset,cancer,omics) {
  
  ctype<-omics_id[ omics_id$datasets==dataset&omics_id$omics==omics&omics_id$cancertype==cancer,]
  categorical_vars<-ctype[ ctype$con_cat=="cat"&ctype$display=="y","col_names"]
  continuous_vars<-ctype[ ctype$con_cat=="con"&ctype$display=="y","col_names"]
  
  # 1. cat vs con
  krus<-c()
  dunn<-c()
  myaov<-c()
  tukey<-c()
  for (cat_var in categorical_vars) {
    for (cont_var in continuous_vars) {
      anova_result <- anova_test(data, cat_var, cont_var,dataset,cancer,omics,omics_id)
	  myaov<-rbind(myaov,anova_result[[1]])
	  tukey<-rbind(tukey,anova_result[[2]])
      
      # Kruskal-Wallis  Dunn's
      kruskal_dunn_result <- kruskal_dunn_test(data, cat_var, cont_var,dataset,cancer,omics,omics_id )
	  krus<-rbind(krus,kruskal_dunn_result[[1]])
	  dunn<-rbind(dunn,kruskal_dunn_result[[2]])
    }
  }
  
  # 2. cat vs cat
  enrich_phyper<-c()
  for (cat_var1 in categorical_vars) {
    for (cat_var2 in categorical_vars) {
      if (cat_var1 != cat_var2) {
        enrichment<-analyze_enrichment(data,cat_var1,cat_var2,dataset,cancer,omics,omics_id)
		enrich_phyper<-rbind(enrich_phyper,enrichment)
      }
    }
  }
  
  # 3. con vs con
  spearman<-c()
  for (cont_var1 in continuous_vars) {
    for (cont_var2 in continuous_vars) {
      if (cont_var1 != cont_var2) {
        spearman_result <- continuous_correlation_test(data, cont_var1, cont_var2, method = "spearman",dataset,cancer,omics,omics_id)
		spearman<-rbind(spearman,spearman_result)
      }
    }
  }
  
  results<-list(myaov=myaov,tukey=tukey,krus=krus,dunn=dunn,enrich=enrich_phyper,spearman=spearman)
  
  return(results)
}


