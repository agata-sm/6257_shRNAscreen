#!/usr/bin/env Rscript

# code for shRNA screen pipeline report
# use with 6257_screen_report_v0.1.Rmd

#### IMPORTANT!

#### VARIABLES TO DEFINE IN CLI
#proj.dir
#proj.name.pref
#data.dir
#mageck.dir
#metadata.dir
#setup="all"### change!!! all 0rm 0rm_noAlt
#library


#Rscript report_launcher.R  proj.dir proj.name.pref data.dir mageck.dir metadata.dir setup library &> ovcar3.report.stderrout
#OBS! positional arguments



## ---- dirs


#metadata
comparisons.file=file.path(metadata.dir,"comparisons.txt")
samples.file=file.path(metadata.dir,"metadata.txt")

#indata dirs
#mageck_datadir=file.path(mageck.dir,setup)
ctable_datadir=file.path(data.dir,"count_table_processed")

#infile
ctable_fname=paste(proj.name.pref,"counts_processed",setup,"tsv",sep=".")


#results dir
resdir=file.path(proj.dir,"results",setup,"rra_annotation")
#resdir=file.path("results","rra_annotation")
dir.create(resdir, recursive = TRUE)

plotdir=file.path(resdir,"plots")
dir.create(plotdir)



## ---- prep_environment
options(connectionObserver = NULL) # https://github.com/rstudio/rstudio/issues/9219


library(tidyverse)
library(ggplot2)
library(DT)
library(plotly)
library(viridis)
library(ggrepel)
library(dplyr)
library(reshape2)
library(htmltools)
library(magrittr)
library(MASS)
library(cowplot)

library(MAGeCKFlute)
#library(idr)
library(edgeR)

library(fgsea)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(biomaRt)
library(ReactomePA)
library(reactome.db)

require(org.Hs.eg.db)


# for plots
theme_pca=theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    strip.background=element_blank(),axis.text.x=element_text(colour="black"),
    axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"),
    legend.text=element_text(size=10) )

fontsize =theme(axis.text=element_text(size=16), axis.title=element_text(size=20))

fig_n=0
tab_n=0
tab_n_sv=0

#set here - significance cutoff for RRA (mainly used for plot labeling)
FDR.CO=0.05
mycolour = c('ns'="gray80", 'dn'= "#377eb8", 'up'="#e41a1c")


#######################################################
########## descriptive stats for read mapping and summarisation

## ---- data_mapstats

samples.tab=read.table(samples.file, sep="\t", header=TRUE, blank.lines.skip=TRUE)
samples.tab$library=samples.tab$sample

file.seqstats=file.path(data.dir,"read_logs","log_stats.txt")
seqstats=read.delim(file.seqstats, header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, row.names=NULL)

library_desc=read.delim(library, header = FALSE, sep = "\t", quote = "\"", fill = TRUE, row.names=NULL,colClasses = "character",blank.lines.skip=TRUE)
colnames(library_desc)=c("cloneId","targetSeq","oligoSeq","cloneName","region","id","NCBI_geneId","symbol","geneDesc","refseqSameGene","refseqAltGeneSameSpecies","refseqAltSpecies")

all_probes=nrow(library_desc)
all_genes=length(unique(library_desc$id))

seqstats$sample=factor(seqstats$sample, levels=samples.tab$library)


seqstats$frac_aln_passed=format(round(as.numeric(seqstats$alignments_mapq255_mismatch/seqstats$alignments),3), nsmall=0, big.mark=",")
seqstats$frac_tot_passed=format(round(as.numeric(seqstats$alignments_mapq255_mismatch/seqstats$all_read_pairs),3), nsmall=0, big.mark=",")
seqstats$frac_trimmed=format(round(as.numeric(seqstats$trimmed_read_pairs/seqstats$all_read_pairs),3), nsmall=0, big.mark=",")
seqstats$frac_mapq255=format(round(as.numeric(seqstats$alignments_mapq255/seqstats$all_read_pairs),3), nsmall=0, big.mark=",")

#formatted numbers for table
seqstats.tab=seqstats[,c(1:6)]

seqstats.tab$condition=samples.tab$condition

seqstats.tab$all_read_pairs=format(round(as.numeric(seqstats.tab$all_read_pairs), 0), nsmall=0, big.mark=",")
seqstats.tab$trimmed_read_pairs=format(round(as.numeric(seqstats.tab$trimmed_read_pairs), 0), nsmall=0, big.mark=",")
seqstats.tab$alignments=format(round(as.numeric(seqstats.tab$alignments), 0), nsmall=0, big.mark=",")
seqstats.tab$alignments_mapq255=format(round(as.numeric(seqstats.tab$alignments_mapq255), 0), nsmall=0, big.mark=",")
seqstats.tab$alignments_mapq255_mismatch=format(round(as.numeric(seqstats.tab$alignments_mapq255_mismatch), 0), nsmall=0, big.mark=",")

seqstats.tab$trimmed_read_pairs=paste0(seqstats.tab$trimmed_read_pairs," (",seqstats$frac_trimmed,")")
seqstats.tab$alignments_mapq255=paste0(seqstats.tab$alignments_mapq255," (",seqstats$frac_mapq255,")")
seqstats.tab$alignments_mapq255_mismatch=paste0(seqstats.tab$alignments_mapq255_mismatch," (",seqstats$frac_tot_passed,")")


colnames(seqstats.tab)=c("sample" ,"all read pairs","trimmed read pairs","alignments","high score alignments","high score alignments passing mismatch filter", "condition")


## ---- data-mapstats-table
knitr::kable(seqstats.tab, row.names = FALSE, caption = "Summary statistics of reads processing. Fractions of total sequenced read pairs given in parentheses.")




## ---- map-stats-plot
seqstats_plot=ggplot(data=seqstats, aes(x=sample, y=alignments_mapq255_mismatch, fill=sample)) + geom_bar(stat="identity") +
    scale_fill_viridis(discrete=TRUE, alpha=0.35,option="turbo") +
    theme_bw() +
    geom_text(aes(label=frac_tot_passed), vjust=1, hjust=1, angle=90, color="black", size=4) +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1), 
      axis.text.y = element_text(size=10)) +
    theme(axis.title.x = element_blank()) +
    labs(y="Read pairs passing filters")


seqstats_plot + 
  theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  + 
  theme(legend.text=element_text(size=10),legend.title=element_text(size=0)) 


## ---- map-stats-plot-save
fig_n=fig_n+1
fname=paste("Figure",fig_n,"shRNAprocessingStatistics.pdf",sep=".")
pdf(file.path(plotdir,fname))
seqstats_plot + 
  theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))  + 
  theme(legend.text=element_text(size=10),legend.title=element_text(size=0)) 
dev.off()



## ---- count-table
ctable.pth=file.path(ctable_datadir,ctable_fname)
reads.ctable=read.delim(ctable.pth, header = TRUE, sep = "\t", quote = "\"",dec = ".", fill = TRUE, row.names=NULL)

colnames(reads.ctable)=gsub("-","." , as.character(colnames(reads.ctable)) )
colnames(reads.ctable)=gsub("_","." , as.character(colnames(reads.ctable)) )

dat_reads_h = gather(reads.ctable[,-2], variable, value,-ProbeID)

#dat_reads_h_ord=dat_reads_h
#dat_reads_h_ord$variable=factor(dat_reads_h_ord$variable, levels=samples.tab$library[order(samples.tab$library)])


## ---- shRNA-stats

# all shRNA in the count table
total_shRNA_library=nrow(library_desc)
total_shRNA=nrow(reads.ctable) #all shRNA (will == total_expr_shRNA if 0-only rows have been filtered)
total_expr_shRNA=sum(rowSums(reads.ctable[,-c(1,2)] != 0) !=0) #shRNA with at least 1 count in at least 1 sample



#shRNA stats per sample
cnts=reads.ctable[,-2]
rownames(cnts)=cnts$ProbeID
cnts=as.matrix(cnts[,-1])

#https://rdrr.io/bioc/edgeR/man/gini.html
gini_idx=gini(cnts)


# expr per sample
# detected guides table
detect_guides=colSums(reads.ctable[,-c(1,2)] != 0)
detguides_table=as.data.frame(cbind(names(detect_guides),detect_guides))
detguides_table=cbind(detguides_table,gini_idx)
colnames(detguides_table)=c("sample","detected_shRNA","Gini_idx")

samples.tab.sub=samples.tab

samples.tab.sub$library=gsub("-","." , as.character(samples.tab.sub$library) )
samples.tab.sub$library=gsub("_","." , as.character(samples.tab.sub$library) )


detguides_table_ord=detguides_table %>% arrange(factor(detguides_table$sample, levels=samples.tab.sub$library))
colnames(detguides_table_ord)=c("sample","detected shRNAs","Gini index")



## ---- shRNA-stats-table
knitr::kable(detguides_table_ord, row.names = FALSE, caption = "Summary of shRNA read coverage. Detected shRNAs - shRNAs with at least 1 read pair assigned.")



## ---- shRNAboxplot
box_reads_shRNA=ggplot(dat_reads_h, aes(x=variable, y=value, color=variable))  + geom_boxplot()+
  theme_bw()+theme(legend.position="none")+
  scale_color_viridis(discrete=TRUE, alpha=0.35,option="turbo") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
    labs(y="Read pairs per shRNA",x="") +
    coord_flip()


fig_n=fig_n+1
fname=paste("Figure",fig_n,"reads_per_shRNA_boxplot.pdf",sep=".")
pdf(file.path(plotdir,fname))
box_reads_shRNA + 
  theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))  + 
  theme(legend.text=element_text(size=10),legend.title=element_text(size=0)) 
dev.off()



## ---- shRNA-stats-plot
box_reads_shRNA +
  theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  + 
  theme(legend.text=element_text(size=10),legend.title=element_text(size=0)) 





## ---- hist-reads-shRNA

x_lim=5000
y_lim=1000

reads_per_shRNA_hist=ggplot(dat_reads_h, aes(x=value)) + geom_histogram(binwidth=100, aes(fill=variable, color=variable) )+
  facet_wrap(~variable)+
  theme_bw() +
    coord_cartesian(xlim=c(0,x_lim), ylim=c(0,y_lim) ) +
  scale_fill_viridis(discrete=TRUE,option="turbo")+
    scale_colour_viridis(discrete=TRUE,option="turbo")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.x = element_text(margin=margin(t=20)) ) + #add margin to x-axis title
  labs(y="reads per shRNA",x="count")


fig_n=fig_n+1
fname=paste("Figure",fig_n,"reads_per_shRNA_histogram.pdf",sep=".")
pdf(file.path(plotdir,fname))
reads_per_shRNA_hist + 
  theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))  + 
  theme(legend.text=element_text(size=10),legend.title=element_text(size=0)) 
dev.off()



## ---- hist-reads-shRNA-plot
reads_per_shRNA_hist +
  theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  + 
  theme(legend.text=element_text(size=10),legend.title=element_text(size=0)) 




############ extra: overlap of genes with counts among samples

## ---- expr_overlap_cond

data_reads_bygene= reads.ctable %>%
  group_by(GeneID) %>%
  summarise(across(where(is.numeric), sum))

s.treated=samples.tab$library[samples.tab$condition=="treated"]
s.untreated=samples.tab$library[samples.tab$condition=="untreated"]
s.ctrl=samples.tab$library[samples.tab$condition=="ctrl"]
s.ctrl=gsub("-",".",s.ctrl)
s.ctrl=gsub("_",".",s.ctrl)


ctrl=as.tibble(data_reads_bygene)%>%rowwise(GeneID)%>%dplyr::select(all_of(s.ctrl))
ctrl_df=data.frame(ctrl[,-1])
rownames(ctrl_df)=ctrl$GeneID
ctrl_expr=ctrl_df[ rowSums( ctrl_df >=1) >=2,]
exprs_in_ctrl=rownames(ctrl_expr)

treated=as.tibble(data_reads_bygene)%>%rowwise(GeneID)%>%dplyr::select(all_of(s.treated))
treated_df=data.frame(treated[,-1])
rownames(treated_df)=treated$GeneID
treat_expr=treated_df[ rowSums( treated_df >=1) >=2,]
exprs_in_treat=rownames(treat_expr)

untreated=as.tibble(data_reads_bygene)%>%rowwise(GeneID)%>%dplyr::select(all_of(s.untreated))
untreated_df=data.frame(untreated[,-1])
rownames(untreated_df)=untreated$GeneID
utreat_expr=untreated_df[ rowSums( untreated_df >=1) >=2,]
exprs_in_utreat=rownames(utreat_expr)

library(ggVennDiagram)
exprs.venn=list(ctrl=exprs_in_ctrl, untreated=exprs_in_utreat, treated=exprs_in_treat)

exprs.venn_plot=ggVennDiagram(exprs.venn) +
 scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
 theme(legend.position="bottom")

fig_n=fig_n+1
fname=paste("Figure",fig_n,"detected_genes_venn.pdf",sep=".")
pdf(file.path(plotdir,fname))
exprs.venn_plot + 
  theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))  + 
  theme(legend.text=element_text(size=10),legend.title=element_text(size=0)) 
dev.off()



## ---- exprs-venn
exprs.venn_plot




#######################################################
####################### RRA


## ---- read-contrasts
contrasts.tab=read.table(comparisons.file, sep="\t", header=TRUE, blank.lines.skip=TRUE)
my.contrasts=contrasts.tab$name
n.cont=nrow(contrasts.tab)



## ---- contrasts-table
knitr::kable(contrasts.tab, row.names = FALSE, caption = "Comparisons analysed in this report.")


## ---- collect-rra-res
all.res.rra=vector(mode = "list", length = nrow(contrasts.tab))
all.res=vector(mode = "list", length = nrow(contrasts.tab))


## ---- contrasts-loop ## not in the report; instead this code is given directly in the chunk
res <- vector(mode = "list", length = n.cont)
options(knitr.duplicate.label = "allow")

for (i in my.contrasts) {
  res[[i]] <- knitr::knit_child("report_rra.Rmd", quiet = TRUE, envir = environment())
}

cat(unlist(res), sep = '\n')




## ---- pca_prep

## PCA on scaled data  after log2 transformation + pseudocount

## TMM scaling

cnts_mat=reads.ctable
rownames(cnts_mat)=cnts_mat$ProbeID

#shRNA
  dat_pca=as.matrix(cnts_mat[,-c(1,2)])
  dge = DGEList(counts=dat_pca) # TMM normalisation only
  dge = calcNormFactors(dge)
  dat_pca_lognorm_sgRNA <- cpm(dge, normalized.lib.sizes = TRUE,log = TRUE, prior.count = 1)

  pca=prcomp(t(dat_pca_lognorm_sgRNA),center = TRUE)


  pca_plot_shRNA=plot_pca_2(pca=pca,sample_annot=samples.tab.sub)

  df_pca_shRNA=as.data.frame(pca$x)
  df_pca_shRNA$library=rownames(df_pca_shRNA)


#genes
  df=cnts_mat

  data_reads_bygene= df %>%
  group_by(GeneID) %>%
  summarise(across(where(is.numeric), sum))

  dat_pca=as.data.frame(data_reads_bygene)
  dat_pca$GeneID[is.na(dat_pca$GeneID)]="na"
  rownames(dat_pca)=dat_pca$GeneID
  dat_pca=dat_pca[,-c(1)]

  dge = DGEList(counts=dat_pca) # TMM normalisation only
  dge = calcNormFactors(dge)
  dat_pca_lognorm_gene <- cpm(dge, normalized.lib.sizes = TRUE,log = TRUE, prior.count = 1)

  pca=prcomp(t(dat_pca_lognorm_gene),center = TRUE)

  pca_plot_gene=plot_pca_2(pca=pca,sample_annot=samples.tab.sub)



#get legend
mock_pca=pca_plot_gene+theme(legend.position="bottom")
legend <- cowplot::get_legend(mock_pca)


#pca.plot=cowplot::plot_grid(pca1,pca2,mock_pca, ncol=3, labels=LETTERS[1:2])

pca1=pca_plot_shRNA+
    theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  + 
    theme(legend.text=element_text(size=10),legend.title=element_text(size=0)) +
    theme(legend.position="none")

pca2=pca_plot_gene+
    theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  + 
    theme(legend.text=element_text(size=10),legend.title=element_text(size=0)) +
    theme(legend.position="none")

plot_pca_combined=plot_grid(plot_grid(pca1, pca2, ncol=2, labels=LETTERS[1:2]),
                 plot_grid(legend, ncol=1), ncol=1,
                  rel_widths=c(1, 0.2))


fig_n=fig_n+1
ggsave(filename=paste("Figure",fig_n,"PCA.pdf",sep="."),path=plotdir,device="pdf", plot=plot_pca_combined)


## ---- pca-plot
plot_pca_combined


## ---- pca_no_ctrls


## PCA on scaled data  after log2 transformation + pseudocount

## TMM scaling

cnts_mat=reads.ctable
rownames(cnts_mat)=cnts_mat$ProbeID

ctrl_smpls=samples.tab.sub$library[samples.tab.sub$condition=="ctrl"]
cnts_mat_noctrl <- cnts_mat [, !names(cnts_mat)%in%ctrl_smpls]
samples.tab.sub_noctrl=samples.tab.sub[!samples.tab.sub$condition=="ctrl",]


#shRNA
  dat_pca=as.matrix(cnts_mat_noctrl[,-c(1,2)])
  dge = DGEList(counts=dat_pca) # TMM normalisation only
  dge = calcNormFactors(dge)
  dat_pca_lognorm_sgRNA_subset <- cpm(dge, normalized.lib.sizes = TRUE,log = TRUE, prior.count = 1)

  pca=prcomp(t(dat_pca_lognorm_sgRNA_subset),center = TRUE)


  pca_plot_shRNA=plot_pca_2(pca=pca,sample_annot=samples.tab.sub_noctrl)

  df_pca_shRNA=as.data.frame(pca$x)
  df_pca_shRNA$library=rownames(df_pca_shRNA)


#genes
  df=cnts_mat_noctrl

  data_reads_bygene= df %>%
  group_by(GeneID) %>%
  summarise(across(where(is.numeric), sum))

  dat_pca=as.data.frame(data_reads_bygene)
  dat_pca$GeneID[is.na(dat_pca$GeneID)]="na"
  rownames(dat_pca)=dat_pca$GeneID
  dat_pca=dat_pca[,-c(1)]

  dge = DGEList(counts=dat_pca) # TMM normalisation only
  dge = calcNormFactors(dge)
  dat_pca_lognorm_gene_subset <- cpm(dge, normalized.lib.sizes = TRUE,log = TRUE, prior.count = 1)

  pca=prcomp(t(dat_pca_lognorm_gene_subset),center = TRUE)

  pca_plot_gene=plot_pca_2(pca=pca,sample_annot=samples.tab.sub_noctrl)



#get legend
mock_pca=pca_plot_gene+theme(legend.position="bottom")
legend <- cowplot::get_legend(mock_pca)


#pca.plot=cowplot::plot_grid(pca1,pca2,mock_pca, ncol=3, labels=LETTERS[1:2])

pca1=pca_plot_shRNA+
    theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  + 
    theme(legend.text=element_text(size=10),legend.title=element_text(size=0)) +
    theme(legend.position="none")

pca2=pca_plot_gene+
    theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  + 
    theme(legend.text=element_text(size=10),legend.title=element_text(size=0)) +
    theme(legend.position="none")

plot_pca_combined_noctrl=plot_grid(plot_grid(pca1, pca2, ncol=2, labels=LETTERS[1:2]),
                 plot_grid(legend, ncol=1), ncol=1,
                  rel_widths=c(1, 0.2))


fig_n=fig_n+1
ggsave(filename=paste("Figure",fig_n,"PCA_noctrl.pdf",sep="."),path=plotdir,device="pdf", plot=plot_pca_combined_noctrl)


## ---- pca-plot-noctrls
plot_pca_combined_noctrl



## ---- smpl_correlations
smpl_corr_sgRNA=cor(dat_pca_lognorm_sgRNA,use="pairwise.complete.obs", method="spearman")
smpl_corr_gene=cor(dat_pca_lognorm_gene,use="pairwise.complete.obs", method="spearman")


cor_hm_sgRNA=plot_corr_hm(corr_matirx=smpl_corr_sgRNA)
cor_hm_gene=plot_corr_hm(corr_matirx=smpl_corr_gene)


fig_n=fig_n+1
fname=paste("Figure",fig_n,"SpearmanCorrHeatmaps.pdf",sep=".")
pdf(file.path(plotdir,fname))

cor_hm_gene +theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  + 
  theme(legend.text=element_text(size=8),legend.title=element_text(size=0))  + ggtitle("Spearman Correlation of TMM normalised shRNA quantification at a gene level")

cor_hm_sgRNA +theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  + 
  theme(legend.text=element_text(size=8),legend.title=element_text(size=0)) + ggtitle("Spearman Correlation of TMM normalised shRNA quantification at a shRNA level")

dev.off()


## ---- smpl-correlations-gene
cor_hm_gene +theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  + 
  theme(legend.text=element_text(size=8),legend.title=element_text(size=0))


## ---- smpl-correlations-sgRNA
cor_hm_sgRNA +theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  + 
  theme(legend.text=element_text(size=8),legend.title=element_text(size=0)) 




## ---- contrast_scatters

my.contrasts=contrasts.tab$name


#contrasts combinations
contrasts.pairs.mtx=combn(my.contrasts,2)

contrasts.pairs.number=ncol(contrasts.pairs.mtx)

#produce the plots
scatters=vector(mode = "list", length = contrasts.pairs.number)


for (i in c(1:contrasts.pairs.number)){

  contr.pair.i=contrasts.pairs.mtx[,i]

  res.df=data.frame()
  for (j in c( contr.pair.i[1], contr.pair.i[2])){
 
      res.df.i=as.data.frame(cbind(all.res[[j]]$neg.lfc,all.res[[j]]$id ))
      colnames(res.df.i)=c("lfc","id")
      res.df.i$comparison=j
      res.df=rbind(res.df,res.df.i)

  }

  res.df$lfc=as.numeric(res.df$lfc)
  df_scatter=spread(res.df,comparison,lfc)

  #add this to avoid error in density distribution
  #https://stackoverflow.com/questions/53075331/error-using-geom-density-2d-in-r-computation-failed-in-stat-density2d-b
  #Error in MASS::kde2d(x, y, ...) : 
  #  missing or infinite values in the data are not allowed
  pseudocount=0.01
  df_scatter[,4]=df_scatter[,2]+pseudocount
  df_scatter[,5]=df_scatter[,3]+pseudocount

  df_scatter$density <- get_density(df_scatter[,2], df_scatter[,3], n = 100)

  pl2=ggplot(df_scatter, aes(x=df_scatter[,4],y=df_scatter[,5], text=paste(id,"; lfc replicate1",round(df_scatter[,2],digits=3),"; lfc replicate2",round(df_scatter[,3],digits=3)) ))+
    geom_point() +
    theme_bw() +
     xlab(paste0("log2 Fold Change in ",colnames(df_scatter)[2])) + ylab(paste0("log2 Fold Change in ",colnames(df_scatter)[3]))


  pl2.d=pl2+geom_point(aes(df_scatter[,4], df_scatter[,5], color = density)) + scale_color_viridis()

  fig_n=fig_n+1
  ggsave(filename=paste("Figure",fig_n,"log2FCscatterplot",colnames(df_scatter)[2],colnames(df_scatter)[3],"comparison",i,"pdf",sep="."),path=plotdir,device="pdf")


  pl2_int=ggplotly(pl2.d, tooltip=c("text"))

  scatters[[i]]=pl2_int

}

## ---- contrast-scatters-plot

for (i in c(1:contrasts.pairs.number)){
  print(scatters[[i]])
}


## ---- nn

