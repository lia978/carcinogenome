library(Biobase)

clean_dose<-function(x){
	i<-x$pert_dose
	res<-sapply(i, function(k){
		if(k == -666) return(NA)
		if(k < 1) return(round(k, digits = 2))
		else return(round(k, digits = 0))
		})
	return(res)
}

reduce_pdata<-function(x){
	sig_id<-as.character(x$sig_id)
	pdat<-data.frame(sig_id = sig_id)
	rownames(pdat)<-sig_id
	pData(x)<-pdat
	colnames(x)<-sig_id
	return(x)
}


annot_chem<-readRDS("~/Desktop/git_projects/environcology/CRCGN_breast/data/CRCGN_breast_chemical_list_by_ss.RDS")

ge<-readRDS("~/Desktop/git_projects/environcology/CRCGN_breast/data/eset/CRCGN_015_017_cp_modz_gene_symbol_median.RDS")
ge
#ge<-add_tas(ge)

x<-pData(ge)

dose_new<-clean_dose(x)

tab<-read.csv("~/Desktop/git_projects/environcology/breast_lung_carcinogens/data/sherr_dose.csv")

well_num<- as.numeric(substr(x$rna_well, 2, 3))


my_round<-function(i){
	if(i > 1) return(round(i, digits = 0))
	else return(round(i, digits = 2))
}

dose_new2<-unlist(lapply(1:nrow(x), function(i){
	if(is.na(dose_new[i])){
		ind<-match(as.character(x$Broad_external_Id[i]), as.character(tab$Broad_external_Id))
		print(ind)
		curr_dose<-tab$dose_um[ind]
		curr_level<-well_num[i] %% 3
		if(curr_level == 0) return(my_round(curr_dose / 9))
		else if (curr_level == 1) return(my_round(curr_dose))
		else if (curr_level == 2) return(my_round(curr_dose/3))
		return(curr_dose)
	} else return(dose_new[i]) 
	}))

dose_new3<-unlist(lapply(1:nrow(x), function(i){
	if(is.na(dose_new[i])){
		ind<-match(as.character(x$Broad_external_Id[i]), as.character(tab$Broad_external_Id))
		curr_dose<-tab$dose_um2[ind]
		if(is.na(curr_dose)) return("")
		else {
			curr_level<-well_num[i] %% 3
			if(curr_level == 0) return(my_round(curr_dose / 9))
			else if (curr_level == 1) return(my_round(curr_dose))
			else if (curr_level == 2) return(my_round(curr_dose/3))
			return(curr_dose)
		}
	} else return("") 
	}))

dose_combined<-sapply(1:nrow(x), function(i) {
	if(dose_new3[i] != "")
		return(paste0(dose_new2[i], "+",dose_new3[i]))
	else return(as.character(dose_new2[i]))
	})
#df<-data.frame(x$"Broad_external_Id", dose1=dose_new, dose2=dose_new2, dose3 = dose_new3, dose4  = dose_combined)
#df[is.na(df$dose1),]
#table(well_num %%3, dose_new)

x$"Chemical Name"<-x$Chemical.name
x$"dose (uM)"<-dose_combined

carc<-as.character(x$breast_carcinogen_any)
carc[carc %in% "POSITIVE"]<-"+"
carc[carc %in% "NEGATIVE"]<-"-"
carc[!(carc %in% c("+", "-"))]<-"N/A"

gtx<-as.character(x$Genotoxicity)
gtx[gtx %in% "POSITIVE"]<-"+"
gtx[gtx %in% "NEGATIVE"]<-"-"
gtx[!(gtx %in% c("+", "-"))]<-"N/A"

x$"Carcinogenicity"<-carc
x$"Genotoxicity"<-gtx
x$"TAS"<-round(x$tas, 3)
x$"Genotype"<-sapply(1:nrow(x), function(i){
	if(grepl("WT", x$sig_id[i])) return("WT")
	else return("TP53")
	})
x$"unique_ID_by_chem"<-paste0(x$"Genotype", "_", x$"dose (uM)")
cols<-c("sig_id", "BUID", "TAS", "Chemical Name", "CAS", "Carcinogenicity", "Genotoxicity", "dose (uM)", "Genotype", "unique_ID_by_chem")
x<-x[, cols]
rownames(x)<-x$sig_id


dat<-list()
dat[["Profile Annotation"]]<-x

x2<-do.call(rbind, lapply(unique(x$"Chemical Name"), function(i){
	ind<-which(x$"Chemical Name" %in% i)
	x[ind[1],c("BUID", "Chemical Name", "CAS", "Carcinogenicity", "Genotoxicity")]
	}))
x2$"TAS (mean)"<-do.call(c, lapply(unique(x$"Chemical Name"), function(i){
	ind<-which(x$"Chemical Name" %in% i)
	mean(x[ind,"TAS"])
	}))
x2$"TAS (mean)"<-round(x2$"TAS (mean)", 3)

dat[["Chemical Annotation"]]<-x2


ge<-reduce_pdata(ge)
lm<-as.character(fData(ge)$pr_is_lmark)
lm[lm %in% "1"]<-"Yes"
lm[lm %in% "0"]<-"No"
fData(ge)$pr_is_lmark <-lm
#fData(ge)<-fData(ge)[, c("id", "pr_gene_symbol", "pr_is_lmark")]
#colnames(fData(ge))<-c("Affy ID", "Gene Symbol", "Landmark Gene")
fData(ge)<-fData(ge)[, c( "pr_gene_symbol", "pr_is_lmark")]
colnames(fData(ge))<-c("Gene Symbol", "Landmark Gene")

dat[["Gene Expression"]]<-ge


gsscores.dir<-"~/Desktop/git_projects/environcology/CRCGN_breast/results/gsscores_cp"
gsscores.files<-list.files(gsscores.dir, pattern = "^gsscores")
gsscores<-lapply(gsscores.files, function(i)
	readRDS(file = paste0(gsscores.dir, "/", i))
)

gsscores<-lapply(gsscores, function(i){
	i<-reduce_pdata(i)
	return(i)
	})

names(gsscores)<-gsscores.files

dat[["Gene Set Enrichment"]]<-gsscores

connectivitypcl<-readRDS("~/Desktop/git_projects/environcology/cmap_sync/data/connectivity_breast/ps_pcl_summary.RDS")
connectivitypert<-readRDS("~/Desktop/git_projects/environcology/cmap_sync/data/connectivity_breast/ps_pert_summary.RDS")

connectivitypcl<-reduce_pdata(connectivitypcl)
connectivitypert<-reduce_pdata(connectivitypert)

con<-list(pcl = connectivitypcl, pert = connectivitypert)
dat[["Connectivity"]]<-con


dat[["title"]]<-"MCF10A Portal"
dat[["about page"]]<-"introduction_MCF10A.Rmd"
dir.create("./data")
dir.create("./data/MCF10A")
saveRDS(dat, file = "./data/MCF10A/data.RDS")

##save gct files

library(cmapR)
dirout<-"data/MCF10A/gct"
dir.create(dirout, recursive = TRUE)
dat<-readRDS("data/MCF10A/data.RDS")
eset2gct<-function(eset, annot, cols = NA, match_id = "sig_id"){

	if(is.na(cols)) cols <-colnames(annot)

	pdat<-annot[match(pData(eset)[, match_id], annot[, match_id]), cols]
	for(j in colnames(pdat)){
		pdat[,j]<-as.character(pdat[,j])
	}
	cas<-pdat$CAS
	pdat$CAS<-gsub("\n", "_", cas)

	fdat<-fData(eset)
	cid<- as.character(pData(eset)[, match_id])
	mat<-exprs(eset)
	ds<-new("GCT", mat = mat, 
	rid = rownames(eset), cid = cid,
	rdesc = fdat, 
	cdesc = pdat, #version = 3, 
	src = "")
}

res<-eset2gct(eset = dat[["Gene Expression"]], annot = dat[["Profile Annotation"]])
fout<-paste0(dirout, "/", "ge.gct")
cmapR::write.gct(ds = res, ofile = fout, appenddim = FALSE)

lm<-dat[["Gene Expression"]]
lm<-lm[fData(lm)[, "Landmark Gene"] %in% "Yes",]
fout<-paste0(dirout, "/", "ge_lm.gct")
res<-eset2gct(eset = lm, annot = dat[["Profile Annotation"]])
cmapR::write.gct(ds = res, ofile = fout, appenddim = FALSE)

gs<-dat[["Gene Set Enrichment"]]
for(i in names(gs)){
	eset<-gs[[i]]
	for(j in colnames(fData(eset))){
		fData(eset)[, j]<-as.character(fData(eset)[, j])
	}
	for(j in colnames(pData(eset))){
		pData(eset)[, j]<-as.character(pData(eset)[, j])
	}
	res<-eset2gct(eset = eset, annot = dat[["Profile Annotation"]])
	fout<-gsub(".RDS", ".gct", paste0(dirout, "/", i))
	cmapR::write.gct(ds = res, ofile = fout, appenddim = FALSE)
	print(fout)
}

conn<-dat[["Connectivity"]]
for(i in names(conn)){
	eset<-conn[[i]]
	res<-eset2gct(eset = eset, annot = dat[["Profile Annotation"]])
	fout<- paste0(dirout, "/conn_", i, ".gct")
	cmapR::write.gct(ds = res, ofile = fout, appenddim = FALSE)
	print(fout)
}
