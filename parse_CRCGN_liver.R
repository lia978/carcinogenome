library(Biobase)
add_tas<-function(dat){
	dat2<-read.csv("~/Desktop/git_projects/environcology/cmap_sync/data/annotation/tas_annotation.csv")
	dat$tas<-dat2$tas[match(dat$sig_id, dat2$sig_id)]
	return(dat)
}

reduce_pdata<-function(x){
	sig_id<-as.character(x$sig_id)
	pdat<-data.frame(sig_id = sig_id)
	rownames(pdat)<-sig_id
	pData(x)<-pdat
	colnames(x)<-sig_id
	return(x)
}
clean_dose<-function(x){
	i<-x$pert_dose
	j<-x$pert_dose_unit
	i2<-sapply(1:length(i), function(k){
		if (j[k] == "nm")
		return(i[k] * 0.001)
		else return(i[k])
		})
	bins<-c(40, 20, 10, 5, 2.5, 1.25, 0.625)
	get_nearest<-function(i, bins = bins){
		if (i < 0.5) return(i)
		else return(bins[which.min(abs(bins - i))])
	}
	i3<-sapply(i2, function(k){
		get_nearest(k, bins = bins)
		})
	#x[which(i3 < 0.0625), c("Chemical.name", "pert_idose")]
	i3
}
##make profile annotation
annot_chem<-readRDS("~/Desktop/git_projects/environcology/cmap_sync/data/annotation/CRCGN008-013_chemical_list_by_ss.RDS")
ge<-readRDS("~/Desktop/git_projects/environcology/cmap_sync/data/eset/subset/CRCGN008-013_full_median_cp_scq75.RDS")
ge
ge<-add_tas(ge)

x<-pData(ge)
x$"Chemical Name"<-x$Chemical.name
x$"dose (uM)"<-clean_dose(x)
carc<-as.character(x$carc_liv_final)
carc[carc %in% "POSITIVE"]<-"+"
carc[carc %in% "NEGATIVE"]<-"-"
carc[!(carc %in% c("+", "-"))]<-"N/A"
x$"Carcinogenicity"<-carc

gtx<-as.character(x$Genotoxicity)
gtx[gtx %in% "POSITIVE"]<-"+"
gtx[gtx %in% "NEGATIVE"]<-"-"
gtx[!(gtx %in% c("+", "-"))]<-"N/A"

x$"Carcinogenicity"<-carc
x$"Genotoxicity"<-gtx
x$"TAS"<-round(x$tas, 3)

cols<-c("sig_id", "BUID", "TAS", "Chemical Name", "CAS", "Carcinogenicity", "Genotoxicity", "dose (uM)")
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
lm[lm %in% "Y"]<-"Yes"
lm[lm %in% "N"]<-"No"
fData(ge)$pr_is_lmark <-lm
#fData(ge)<-fData(ge)[, c("id", "pr_gene_symbol", "pr_is_lmark")]
#colnames(fData(ge))<-c("Affy ID", "Gene Symbol", "Landmark Gene")
fData(ge)<-fData(ge)[, c( "pr_gene_symbol", "pr_is_lmark")]
colnames(fData(ge))<-c("Gene Symbol", "Landmark Gene")

ind<-which(fData(ge)[, "Gene Symbol"] %in% "-666")
ge<-ge[setdiff(1:nrow(ge), ind),]

dat[["Gene Expression"]]<-ge

gsscores.dir<-"~/Desktop/git_projects/environcology/cmap_sync/results/gsscores"
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

connectivitypcl<-readRDS("~/Desktop/git_projects/environcology/cmap_sync/data/connectivity/ps_pcl_summary.RDS")
connectivitypert<-readRDS("~/Desktop/git_projects/environcology/cmap_sync/data/connectivity/ps_pert_summary.RDS")

connectivitypcl<-reduce_pdata(connectivitypcl)
connectivitypert<-reduce_pdata(connectivitypert)

con<-list(pcl = connectivitypcl, pert = connectivitypert)
dat[["Connectivity"]]<-con


dat[["title"]]<-"HEPG2 Portal"
dat[["about page"]]<-"introduction_HEPG2.Rmd"
dir.create("./data")
dir.create("./data/HEPG2")
saveRDS(dat, file = "./data/HEPG2/data.RDS")


##write to gct files

library(cmapR)
dirout<-"data/HEPG2/gct"
dir.create(dirout, recursive = TRUE)
dat<-readRDS("data/HEPG2/data.RDS")

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


