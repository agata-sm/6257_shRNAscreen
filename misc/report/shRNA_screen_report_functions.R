# crispr_pipeline_report_functions.R

#save tables in tab-delimted format
save.table <- function(df=df, file=file, dir=dir) {
	write.table(df, file.path(dir, file), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, fileEncoding = "")
}


#numbering for html widgets (i.e. interactive tables)
#https://github.com/rstudio/bookdown/issues/372
test_table <- function(caption, label) {
  require(magrittr)
  lab <- paste0("(#tab:", label, ")")
  out <- c(
    "<table style='width: 99%;'>",
    paste0("<caption>", lab, caption, "</caption>"),
    "</table>"
  )
  cat(out, sep="\n")
}

# capt="Summary of gene level analysis of replicate 1, negative selection (MAGeCK). Shown are top 100 genes (by \"neg.rank\")."
# test_table(capt, label="tabsummary-neg-interact")
# DT::datatable(rep1g_neg[,1:9] %>% filter(neg.rank <=100), rownames = FALSE)



# save results of GSEA GO
save_GSEA_GO <- function(GSEAres=GSEAres, file=file, dir=dir){

	to_save=as.data.frame(GSEAres@result)
	save.table(df=to_save, file=file,dir=dir)

}


# save results of reactme GO
save_GSEA_react <- function(GSEAres=GSEAres, file=file, dir=dir){

	to_save=as.data.frame(GSEAres)
	save.table(df=to_save, file=file,dir=dir)

}

#generate GSEA plots: dotplot, treeplot, heatplot, GSEA ES plots of 6 top categories
GSEA_plots <- function(GSEAres=GSEAres, lfc_vect=lfc_vect){
	
	plot_list=list()

	gse.summary=as.data.frame(cbind(GSEAres@result$Description,GSEAres@result$ID,GSEAres@result$NES,GSEAres@result$p.adjust))
	colnames(gse.summary)=c("Description","ID","NES","p.adjust")
	gse.summary=gse.summary[order(gse.summary$NES, decreasing=TRUE)[1:20],]
	categories_dot=gse.summary$Description
	  
	dotplot.gsea=dotplot(GSEAres, showCategory=categories_dot, orderBy = "NES", label_format=50) +
	theme(text = element_text(size = 10)) + theme(axis.text.y = element_text(size=10)) +  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE), name = "FDR")

	heatplot.gsea <- heatplot(GSEAres, foldChange=lfc_vect, showCategory=10) +
	theme(text = element_text(size = 10)) + theme(axis.text.y = element_text(size=10)) +  theme(axis.text.x=element_blank()) +
	scale_fill_viridis_c(guide=guide_colorbar(reverse=TRUE))


	treeplot.gsea=treeplot(GSEAres) +  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE), name = "FDR")

  	
  	gseaplot.1=gseaplot2(GSEAres, geneSetID = 1, title = GSEAres$Description[1], base_size=7, rel_heights=c(3, .5, 2)) 
  	gseaplot.2=gseaplot2(GSEAres, geneSetID = 2, title = GSEAres$Description[2], base_size=7, rel_heights=c(3, .5, 2)) 
  	gseaplot.3=gseaplot2(GSEAres, geneSetID = 3, title = GSEAres$Description[3], base_size=7, rel_heights=c(3, .5, 2)) 
  	gseaplot.4=gseaplot2(GSEAres, geneSetID = 4, title = GSEAres$Description[4], base_size=7, rel_heights=c(3, .5, 2)) 
  	gseaplot.5=gseaplot2(GSEAres, geneSetID = 5, title = GSEAres$Description[5], base_size=7, rel_heights=c(3, .5, 2)) 
  	gseaplot.6=gseaplot2(GSEAres, geneSetID = 6, title = GSEAres$Description[6], base_size=7, rel_heights=c(3, .5, 2)) 

	gseaplot.gsea=cowplot::plot_grid(gseaplot.1,gseaplot.2,gseaplot.3,gseaplot.4,gseaplot.5,gseaplot.6, ncol=3, labels=LETTERS[1:6])

	plot_list[[1]]=dotplot.gsea
	plot_list[[2]]=heatplot.gsea
	plot_list[[3]]=treeplot.gsea
	plot_list[[4]]=gseaplot.gsea

	return(plot_list)
}



#GSEA ES plots for top n categories

GSEA_ES_plots <- function(GSEAres=GSEAres, dir=dir, n_plots=n_plots,sel=sel){

	for (j in (1:n_plots)){
		gseaplot.i=gseaplot2(GSEAres, geneSetID = j, title = GSEAres$Description[j]) + theme(text = element_text(size = 7))
		fname=paste(i,"GSEAplot",sel,j,"pdf",sep=".")

		pdf(file.path(dir,fname),paper = "a4r")
		print(gseaplot.i)
		dev.off()

	}

}


#for plotting correlation heatmaps
reorder_cormat <- function(cormat){
# Use correlation between variables as distance
dd <- as.dist((1-cormat)/2)
hc <- hclust(dd)
cormat <-cormat[hc$order, hc$order]
}


# plot heatmap of correlation coefficients

plot_corr_hm<-function(corr_matirx=corr_matirx){

	mincor=min(corr_matirx)
	maxcor=max(corr_matirx)
	mediancor=median(corr_matirx)

	cormat <- reorder_cormat(corr_matirx)
	melted_cormat <- melt(cormat)

	cor_hm=ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  		geom_tile() +
  		scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
  		 midpoint = mediancor, limit = c(mincor,maxcor), space = "Lab", 
    	name="Spearman\nCorrelation") + 
    	coord_fixed() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
		theme(text = element_text(size = 10)) +
		theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

		return(cor_hm)
}






#plot PCA on TMM normalised counts
#input is a matrix

#plot PCA on normalised counts
#input is pca object from prcomp

plot_pca_2<-function(pca=pca,sample_annot=sample_annot){ #pca is pca object from prcomp

	require(ggplot2)
	require(viridis)
	  
	  df_pca=as.data.frame(pca$x)
	  df_pca$library=rownames(df_pca)
	  df_pca=merge(df_pca,sample_annot,by="library")
		
		cond1=sample_annot$condition
		cond2=sample_annot$library


	legendname=colnames(sample_annot)[3]

	pca_plot=ggplot(df_pca, aes(x=PC1,y=PC2, colour=cond1, label=cond2)) + geom_point(size=6)
	pca_plot= pca_plot + theme(text = element_text(size = 6)) +
  		theme_bw() + scale_colour_viridis(discrete=TRUE,option="turbo",alpha=0.85,name = legendname)
	pca_plot=pca_plot+ geom_text_repel(size=4, max.overlaps=10) 

	return(pca_plot)

}



#take a vector of strings 1234 (0.45), return 1,234 (0.45)
format_parenth<-function(x=x){
  x.number=as.numeric(gsub("\\s*\\([^\\)]+\\)","",as.character(x) ))
  x.number.format=format(round(as.numeric(x.number), 0), nsmall=0, big.mark=",")

  x.parenth=sub("\\).*", "", as.character(x))
  x.parenth=sub(".*\\(", "", as.character(x.parenth))

  x.formatted=paste0(x.number.format," (",x.parenth,")")
  return(x.formatted)
}




#useful for reduced output of sessionInfo
#https://stackoverflow.com/questions/49660742/r-reduce-sessioninfo-output
mySIprint <- function(x, locale = TRUE, ...) {
    mkLabel <- function(L, n) {
        vers <- sapply(L[[n]], function(x) x[["Version"]])
        pkg <- sapply(L[[n]], function(x) x[["Package"]])
        paste(pkg, vers, sep = "_")
    }
    # ... don't forget this code that I'm snipping out for brevity
    if (!is.null(x$otherPkgs)) {
        cat("\nother attached packages:\n")
        print(mkLabel(x, "otherPkgs"), quote = FALSE, ...)
    }
    #if (!is.null(x$loadedOnly)) {
    #    cat("\nloaded via a namespace (and not attached):\n")
    #    print(mkLabel(x, "loadedOnly"), quote = FALSE, ...)
    #}
    invisible(x)
}



# for density coloured scatter plots
# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  require(MASS)
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}



