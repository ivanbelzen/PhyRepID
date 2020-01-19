### PhyRepID analysis and figures
### ivanbelzen
### Last update: 15 december 2019

### Analysis and Figures
## General dataframe and clan statistics
## For clans vs non-clan, and Schaper, Selectome, ExAC, Segmental duplication regions:
### Enrichment tests with Fisher's Exact Test
### Distribution comparison Wilcoxon
## Plots of main and supplementary figures

#install.packages("ggpubr")
#install.packages("ggplot2")
#dependancy RGL needs libraries 
#sudo apt install libftgl2 libcgal-dev libglu1-mesa-dev libglu1-mesa-dev libx11-dev libfreetype6-dev

#Source: https://stackoverflow.com/questions/14613355/how-to-get-something-like-matplotlibs-symlog-scale-in-ggplot-or-lattice
symlog_trans <- function(base = 10, thr = 1, scale = 1){
  trans <- function(x)
    ifelse(abs(x) < thr, x, sign(x) * 
             (thr + scale * suppressWarnings(log(sign(x) * x / thr, base))))
  
  inv <- function(x)
    ifelse(abs(x) < thr, x, sign(x) * 
             base^((sign(x) * x - thr) / scale) * thr)
  
  breaks <- function(x){
    sgn <- sign(x[which.max(abs(x))])
    if(all(abs(x) < thr))
      pretty_breaks()(x)
    else if(prod(x) >= 0){
      if(min(abs(x)) < thr)
        sgn * unique(c(pretty_breaks()(c(min(abs(x)), thr)),
                       log_breaks(base)(c(max(abs(x)), thr))))
      else
        sgn * log_breaks(base)(sgn * x)
    } else {
      if(min(abs(x)) < thr)
        unique(c(sgn * log_breaks()(c(max(abs(x)), thr)),
                 pretty_breaks()(c(sgn * thr, x[which.min(abs(x))]))))
      else
        unique(c(-log_breaks(base)(c(thr, -x[1])),
                 pretty_breaks()(c(-thr, thr)),
                 log_breaks(base)(c(thr, x[2]))))
    }
  }
  trans_new(paste("symlog", thr, base, scale, sep = "-"), trans, inv, breaks)
}

library(ggplot2)
library(ggpubr)
library(scales)
path=""

### Settings
gene_symbol_list = c("AHNAK","CDC20","VCAM1", "KNL1", "PRDM9", "SPATA31A3", "PRX", "TEP1", "CNTRL") #examples to display
top_clans_cnt=10 #Number of clans include in the top X most frequent list

# Output files
sink(paste0(path,"PhyRepID_Routput.txt"))
pdf(paste0(path,"PhyRepID_plots.pdf"),paper="a4r")
clan_summary_file=paste0(path,"SupTable_4_Pfam_clan_summary.tsv") #dataframe with clan statistics
SupFig_top_clans_file=paste0(path,"SupFig_top_clans_comparision.pdf")
SupFig_file=paste0(path,"Supplementary_figures.pdf")
Fig1_prd_score_file=paste0(path,"Fig1_prd_score_landscape.svg")
Fig2_topclan_file = paste0(path,"Fig2_top5_clans_prdscore_all_zoom.svg")
  
# Input files

scores=read.table(paste0(path,"SupTable_1_PRD_score_ranking.tsv"),sep="\t",na.strings = NULL,header=T)
human_anno=read.table(paste0(path,"SupTable_2_Human_lineage_full.tsv"),sep="\t",na.strings = NULL,header=T)
mscores = merge(x=scores, y=human_anno[,c("identifier","human_dup","human_frac","human_protein_cnt")], by="identifier",all=T)
scores=mscores


### Dataframes

## Process clans and construct dataframe 
scores$clan = factor(scores$clan)#, levels=clans[order(-clans$identifier),c("clan")]) 
clan_summary = aggregate(scores[,c("prd_score","netto_dup","gt_dup")], by=list(clan=scores$clan), FUN=mean)
clan_summary$count = aggregate(scores$identifier, by=list(clan=scores$clan), FUN=length)$x
clan_summary$frac = clan_summary$count/nrow(scores)
top_clans = head(clan_summary[order(-clan_summary$count),],n=top_clans_cnt)
#clan_list = top_clans$clan #short list for display purposes

## Genetree dataframe
gt_genesymbol = aggregate(gene_symbol ~ genetree_id, data = scores, paste, collapse=",")
gt_dataframe = aggregate(scores[,c("prd_score","netto_dup","gt_dup")],by=list(genetree_id=scores$genetree_id), data = scores, FUN=mean)
mgt_dataframe = merge(gt_genesymbol,gt_dataframe,by="genetree_id")
gt_dataframe = mgt_dataframe

## Scores dataframe 
scores$schaper_pos = factor(scores$schaper_pos, levels=c("True","False", ""))
scores$selectome = factor(scores$selectome, levels=c("True","False"))
scores$segdup = factor(scores$segdup, levels=c("True","False", ""))
scores$mis_z = as.numeric(as.character(scores$mis_z))
scores$pLI = as.numeric(as.character(scores$pLI))

scores$gt_score = scores$gt_dup/scores$orthologs_cnt
scores[scores$gt_dup>0, "gt_bool"]="True"
scores[scores$gt_dup==0,"gt_bool"]="False"
scores$gt_bool = factor(scores$gt_bool, levels=c("True","False"))

### Dataset characterisation

## Top 100 
top100= head(scores[order(-scores$prd_score),],n=100)   
scores[,"rank"]="False"
scores[scores$identifier %in% top100$identifier,"rank"]="True"
scores$rank = factor(scores$rank, levels=c("True","False"))

top100_mean=mean(top100$prd_score)
overall_mean = mean(scores$prd_score)

human_lin=scores[scores$human_dup>0,]
human100 = head(human_lin[order(-human_lin$human_frac),],n=100)
scores[,"rank_humlin"]="False"
scores[scores$identifier %in% human_lin$identifier,"rank_humlin"]="True"
scores$rank_humlin = factor(scores$rank_humlin, levels=c("True","False"))
cnt_human_lin = nrow(scores[scores$rank_humlin=="True",])

## Report general statistics/characteristics
print(" ## Dataframe characteristics ")
print(paste0("Total number of OGs: ",nrow(scores)))
print(paste0("Conserved #OGs <=0 ",nrow(scores[scores$prd_score<=0,]),", fraction ",nrow(scores[scores$prd_score<=0,])/nrow(scores)))
print(paste0("Gene dup #OGs >0 ",nrow(scores[scores$gt_dup>0,]),", fraction ",nrow(scores[scores$gt_dup>0,])/nrow(scores)))
print(paste0("Human lineage (full) >1 duplication #OGs ",cnt_human_lin))
print(paste0("Mean PRD score: ",mean(scores$prd_score)))
print(paste0("Number of gene trees: ",nrow(gt_dataframe)))
print(paste0("Mean PRD score, avg gene tree: ",mean(gt_dataframe$prd_score)))
    

## Clan analysis and reporting
print(paste0("## Clan enrichment tests, top ",top_clans_cnt," clans "))
for(i in unique(scores$clan)){
  scores[scores$clan==i,"clan_bool"]="True" 
  scores[scores$clan!=i,"clan_bool"]="False" 
  scores$clan_bool = factor(scores$clan_bool, levels=c("True","False"))
  
  clan_mean = clan_summary[clan_summary$clan==i,c("prd_score")]
  if(clan_mean < overall_mean) {
    w=wilcox.test(prd_score ~ clan_bool, data = scores,alternative="less")
  }else{
    w=wilcox.test(prd_score ~ clan_bool, data = scores,alternative="greater")
  }
  
  fisher.table = matrix(c(nrow(scores[scores$rank=="True"&scores$clan_bool=="True",]),
                          nrow(scores[scores$rank=="True"&scores$clan_bool=="False",]),
                          nrow(scores[scores$rank=="False"&scores$clan_bool=="True",]),
                          nrow(scores[scores$rank=="False"&scores$clan_bool=="False",])
  ), nrow = 2)
  f=fisher.test(fisher.table)
 
  #fill table with values
  clan_summary[clan_summary$clan==i,"fisher.pval"]=f$p.value
  clan_summary[clan_summary$clan==i,"fisher.odss"]=f$estimate
  clan_summary[clan_summary$clan==i,"wilcoxon.pval"]=w$p.value
  
  #report only for top but save for all
  if(i %in% top_clans$clan) {
    print(paste0(i," (",clan_mean,"): testing"))
    
    print(w)
    print(paste0(i,": Fisher enrichment in top 100"))
    print(f)
    
    if(f$p.value<0.05){
      print(fisher.table)
    }
  }
}

write.table(clan_summary, file=clan_summary_file, quote=F, sep="\t", row.names = F, col.names = T)

### Tests on human set, require duplications in human lineage

## Schaper comparison enrichment
print("## Categorial enrichment tests  ")
print("### Schaper ")
print("Schaper comparison enrichment in top 100 human ")
fisher.table = matrix(c(nrow(human100[human100$schaper_pos=="True",]),
                        nrow(human100[human100$schaper_pos=="False",]),
                        nrow(scores[scores$schaper_pos=="True" &! scores$identifier %in% human100$identifier,]),  
                        nrow(scores[scores$schaper_pos=="False" &! scores$identifier %in% human100$identifier,]) 
                        ), nrow = 2)
print(fisher.table)
fisher.test(fisher.table)

print(paste0("Schaper comparison enrichment in  human lin #OGs ",cnt_human_lin))
fisher.table = matrix(c(nrow(human_lin[human_lin$schaper_pos=="True",]),
                        nrow(human_lin[human_lin$schaper_pos=="False",]),
                        nrow(scores[scores$schaper_pos=="True" &! scores$identifier %in% human_lin$identifier,]),  
                        nrow(scores[scores$schaper_pos=="False" &! scores$identifier %in% human_lin$identifier,]) 
), nrow = 2)
print(fisher.table)
fisher.test(fisher.table)
print(paste0("Schaper comparison enrichment in human lin #OGs  ",cnt_human_lin," including unannotated by Schaper as bottom row"))
print(matrix(c(nrow(human_lin[human_lin$schaper_pos=="True",]),
                        nrow(human_lin[human_lin$schaper_pos=="False",]),
                        nrow(human_lin[human_lin$schaper_pos=="",]),
                        nrow(scores[scores$schaper_pos=="True" &! scores$identifier %in% human_lin$identifier,]),  
                        nrow(scores[scores$schaper_pos=="False" &! scores$identifier %in% human_lin$identifier,]),
                        nrow(scores[scores$schaper_pos=="" &! scores$identifier %in% human_lin$identifier,]) 
), nrow = 3))

print("Schaper comparison enrichment False in <0 PRD score")
fisher.table = matrix(c(nrow(scores[scores$schaper_pos=="False" & scores$prd_score<0,]),
                        nrow(scores[scores$schaper_pos=="True" & scores$prd_score<0,]),
                        nrow(scores[scores$schaper_pos=="False" &! scores$prd_score<0,]),  
                        nrow(scores[scores$schaper_pos=="True" &! scores$prd_score<0,]) 
), nrow = 2)
print(fisher.table)
fisher.test(fisher.table)

print("Schaper comparison enrichment False in <0 PRD score, including unannotated by Schaper as bottom row")
fisher.table = matrix(c(nrow(scores[scores$schaper_pos=="False" & scores$prd_score<0,]),
                        nrow(scores[scores$schaper_pos=="True" & scores$prd_score<0,]),
                        nrow(scores[scores$schaper_pos=="" & scores$prd_score<0,]),
                        nrow(scores[scores$schaper_pos=="False" &! scores$prd_score<0,]),  
                        nrow(scores[scores$schaper_pos=="True" &! scores$prd_score<0,]),
                        nrow(scores[scores$schaper_pos=="" &! scores$prd_score<0,]) 
                        
), nrow = 3)
print(fisher.table)
## Selectome enrichment in top 100
print("### Selectome")
print("Selectome enrichment in top 100  PRD score")
fisher.table = matrix(c(nrow(top100[top100$selectome=="True",]),
                        nrow(top100[top100$selectome=="False",]),
                        nrow(scores[scores$selectome=="True" &! scores$identifier %in% top100$identifier,]),  
                        nrow(scores[scores$selectome=="False" &! scores$identifier %in% top100$identifier,]) 
), nrow = 2)
print(fisher.table)
fisher.test(fisher.table)

print("Selectome  comparison enrichment in top 100 human ")
fisher.table = matrix(c(nrow(human100[human100$selectome=="True",]),
                        nrow(human100[human100$selectome=="False",]),
                        nrow(scores[scores$selectome=="True" &! scores$identifier %in% human100$identifier,]),  
                        nrow(scores[scores$selectome=="False" &! scores$identifier %in% human100$identifier,]) 
), nrow = 2)
print(fisher.table)
fisher.test(fisher.table)


print(paste0("Selectome  comparison enrichment in human lin #OGs ",cnt_human_lin))
fisher.table = matrix(c(nrow(human_lin[human_lin$selectome=="True",]),
                        nrow(human_lin[human_lin$selectome=="False",]),
                        nrow(scores[scores$selectome=="True" &! scores$identifier %in% human_lin$identifier,]),  
                        nrow(scores[scores$selectome=="False" &! scores$identifier %in% human_lin$identifier,]) 
), nrow = 2)
print(fisher.table)
fisher.test(fisher.table)

print("Comparing PRD score means in Selectome T/F:")
print(paste0("Selectome True mean: ", mean(scores[scores$selectome=="True",]$prd_score)))
print(paste0("Selectome False mean: ", mean(scores[scores$selectome=="False",]$prd_score)))
wilcox.test(prd_score ~ selectome, data = scores, alternative="less")

## segdup enrichment in top 100
print("### Segmental duplication regions")
print("SegDup enrichment in top 100  PRD score")
fisher.table = matrix(c(nrow(top100[top100$segdup=="True",]),
                        nrow(top100[top100$segdup=="False",]),
                        nrow(scores[scores$segdup=="True" &! scores$identifier %in% top100$identifier,]),  
                        nrow(scores[scores$segdup=="False" &! scores$identifier %in% top100$identifier,]) 
), nrow = 2)
print(fisher.table)
fisher.test(fisher.table)

print("SegDup comparison enrichment in top 100 human ")
fisher.table = matrix(c(nrow(human100[human100$segdup=="True",]),
                        nrow(human100[human100$segdup=="False",]),
                        nrow(scores[scores$segdup=="True" &! scores$identifier %in% human100$identifier,]),  
                        nrow(scores[scores$segdup=="False" &! scores$identifier %in% human100$identifier,]) 
), nrow = 2)
print(fisher.table)
fisher.test(fisher.table)


print(paste0("SegDup comparison enrichment in human lin #OGs ",cnt_human_lin))
fisher.table = matrix(c(nrow(human_lin[human_lin$segdup=="True",]),
                        nrow(human_lin[human_lin$segdup=="False",]),
                        nrow(scores[scores$segdup=="True" &! scores$identifier %in% human_lin$identifier,]),  
                        nrow(scores[scores$segdup=="False" &! scores$identifier %in% human_lin$identifier,]) 
), nrow = 2)
print(fisher.table)
fisher.test(fisher.table)

print("Comparing PRD score means in SegDup T/F:")
print(paste0("SegDup True mean: ", mean(scores[scores$segdup=="True",]$prd_score)))
print(paste0("SegDup False mean: ", mean(scores[scores$segdup=="False",]$prd_score)))
wilcox.test(prd_score ~ segdup, data = scores, alternative="greater")

## chi2 independance 
print("chi2 independance - Selectome and SegDup")
tbl = matrix(c(nrow(scores[scores$segdup=="True" & scores$selectome=="True",]),
         nrow(scores[scores$segdup=="True" & scores$selectome=="False",]),
         nrow(scores[scores$segdup=="False" & scores$selectome=="True",]),  
         nrow(scores[scores$segdup=="False" & scores$selectome=="False",]) 
), nrow = 2)
print(tbl)
print(chisq.test(tbl, correct=F))


print("chi2 independance - SegDup and gt>0")
tbl = matrix(c(nrow(scores[scores$segdup=="True" & scores$gt_dup>0,]),
               nrow(scores[scores$segdup=="True" & scores$gt_dup==0,]),
               nrow(scores[scores$segdup=="False" & scores$gt_dup>0,]),  
               nrow(scores[scores$segdup=="False" & scores$gt_dup==0,]) 
), nrow = 2)
print(tbl)
print(chisq.test(tbl, correct=F))

print("## EXaC comparison missense z-score and pLI score ")
m=matrix(c(mean(scores$mis_z,na.rm=T),
         mean(top100$mis_z,na.rm=T),
         mean(human_lin$mis_z,na.rm=T),
         mean(scores$pLI,na.rm=T),
         mean(top100$pLI,na.rm=T),
         mean(human_lin$pLI,na.rm=T)), nrow = 3)
rownames(m) = c("dataset","top100", "human_lin")
colnames(m) = c("missense","pLI")
print(m)

print("EXac comparision with top 100")
wilcox.test(mis_z ~ rank, data = scores,alternative="less")
wilcox.test(pLI ~ rank, data = scores,alternative="less")

print(paste0("Exac comparision with human lin #OGs  ",cnt_human_lin))
wilcox.test(mis_z ~ rank_humlin, data = scores,alternative="greater")
wilcox.test(pLI ~ rank_humlin, data = scores,alternative="greater")


## Gene duplication and oversplit OGs

print("## Gene duplication and oversplit OGs")

print("Oversplit OGs (non-euteleostomi root t117571) comparison enrichment in top 100 ")
fisher.table = matrix(c(nrow(top100[top100$gt_root!="t117571",]),
                        nrow(top100[top100$gt_root=="t117571",]),
                        nrow(scores[scores$gt_root!="t117571" &! scores$identifier %in% top100$identifier,]),  
                        nrow(scores[scores$gt_root=="t117571" &! scores$identifier %in% top100$identifier,]) 
), nrow = 2)
print(fisher.table)
fisher.test(fisher.table)

print("Mean gene duplication:")
print(paste0("dataset: ", mean(scores$gt_dup,na.rm=T)))
print(paste0("top 100: ", mean(top100$gt_dup,na.rm=T)))
print(paste0("human_lin: ", mean(human_lin$gt_dup,na.rm=T)))

## Gene duplication boolean enrichment in top 100
print("Gene duplication T/F enrichment in top 100")
fisher.table = matrix(c(nrow(top100[top100$gt_bool=="True",]),
                        nrow(top100[top100$gt_bool=="False",]),
                        nrow(scores[scores$gt_bool=="True" &! scores$identifier %in% top100$identifier,]),  
                        nrow(scores[scores$gt_bool=="False" &! scores$identifier %in% top100$identifier,]) 
), nrow = 2)
print(fisher.table)
fisher.test(fisher.table)

print(paste0("Gene duplication T/F enrichment in human lin #OGs ",cnt_human_lin))
fisher.table = matrix(c(nrow(human_lin[human_lin$gt_bool=="True",]),
                        nrow(human_lin[human_lin$gt_bool=="False",]),
                        nrow(scores[scores$gt_bool=="True" &! scores$identifier %in% human_lin$identifier,]),  
                        nrow(scores[scores$gt_bool=="False" &! scores$identifier %in% human_lin$identifier,]) 
), nrow = 2)
print(fisher.table)
fisher.test(fisher.table)


###TODO

### Figures
## Figure 1: PRD score landscape with examples
display = scores[scores$gene_symbol %in% gene_symbol_list, ]
p <- ggplot(data = scores, aes(x=prd_score, ..scaled..))
p <- p + coord_cartesian(clip = 'off')
p <- p + xlab("PRD score") + ylab("Density estimate,scaled") + ggtitle("PRD score landscape")
p <- p + geom_density(adjust = 1.5) +  scale_x_continuous(trans=symlog_trans(thr=1))
p <- p + geom_point(data=display, aes(x=prd_score,y=0.1,color=gene_symbol),show.legend = F) 
p <- p +  geom_segment(data=display, aes(color=gene_symbol,x=prd_score, xend=prd_score, y=0, yend=0.1),show.legend = F)
p <- p + geom_label(data=display,aes(label=gene_symbol,color=gene_symbol,y=0.1),hjust=0, vjust=0, show.legend=F)
fig1 = p
p

# Summarize on gene tree 
gt_display = gt_dataframe[gt_dataframe$genetree_id %in% display$genetree_id, ]
p <- ggplot(data = gt_dataframe, aes(x=prd_score, ..scaled..))
p <- p + coord_cartesian(clip = 'off')
p <- p + xlab("PRD score") + ylab("Density estimate,scaled") + ggtitle("PRD score landscape, gene tree averaged")
p <- p + geom_density(adjust = 1.5) +  scale_x_continuous(trans=symlog_trans(thr=1))
p <- p + geom_point(data=gt_display, aes(x=prd_score,y=0.1,color=gene_symbol),show.legend = F) 
p <- p +  geom_segment(data=gt_display, aes(color=gene_symbol,x=prd_score, xend=prd_score, y=0, yend=0.1),show.legend = F)
p <- p + geom_label(data=gt_display,aes(label=gene_symbol,color=gene_symbol,y=0.1),hjust=0, vjust=0, show.legend=F)
p

## Plots: comparision of top clans

title=paste0("Top ",top_clans_cnt, " clans: PRD score distribution")
p <- ggplot(data=scores[scores$clan %in% top_clans$clan,], aes(x=clan,  y = prd_score,color=clan,fill=clan))
p <- p + coord_cartesian(clip = 'off') + scale_y_continuous(trans=symlog_trans(thr=3))
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1 ), legend.position = "none")
p <- p + xlab("Group") + ylab("PRD score") +  ggtitle(title)
p <- p + geom_violin(alpha=0.5,trim=F)
supfig_clans_full=p
p


title=paste0("Top ",top_clans_cnt, " clans: PRD score distribution (zoom)")
p <- ggplot(data=scores[scores$clan %in% top_clans$clan,], aes(x=clan,  y = prd_score,color=clan,fill=clan))
p <- p + coord_cartesian(clip = 'on',ylim = c(-1,1)) + scale_y_continuous(trans=symlog_trans(thr=3))
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1 ), legend.position = "none")
p <- p + xlab("Group") + ylab("PRD score") +  ggtitle(title)
p <- p + geom_violin(alpha=0.5,trim=T) 
p <- p+stat_summary(fun.y=mean, geom="point", shape=23, size=1,fill="black",color="black")
p <- p+geom_hline(yintercept = 0)
supfig_clans_zoom=p
p


## Figure 2, comparision of 5 top clans including all 
scores_2=scores
scores_2$clan="all"
data_double = rbind(scores,scores_2)
clans = aggregate(identifier ~ clan, data = data_double, FUN=length)
top_clans_cnt=6
data_double$clan = factor(data_double$clan, levels=clans[order(-clans$identifier),c("clan")])
top_clans = head(clans[order(-clans$identifier),],n=top_clans_cnt)

title=paste0("Top 5 clans: PRD score distribution (zoom)")
p <- ggplot(data=data_double[data_double$clan %in% top_clans$clan,], aes(x=clan,  y = prd_score,color=clan,fill=clan))
p <- p + coord_cartesian(clip = 'on',ylim = c(-1,1)) #+ scale_y_continuous(trans=symlog_trans(thr=3))
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1 ), legend.position = "none")
p <- p + xlab("Group") + ylab("PRD score") +  ggtitle(title)
p <- p + geom_violin(alpha=0.5,trim=F)
p <- p+stat_summary(fun.y=mean, geom="point", shape=23, size=1,fill="black",color="black")
p <- p+geom_hline(yintercept = 0)
fig2 = p
p


title=paste0("Top 5 clans: PRD score distribution")
p <- ggplot(data=data_double[data_double$clan %in% top_clans$clan,], aes(x=clan,  y = prd_score,color=clan,fill=clan))
p <- p + coord_cartesian(clip = 'off') + scale_y_continuous(trans=symlog_trans(thr=3))
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1 ), legend.position = "none")
p <- p + xlab("Group") + ylab("PRD score") +  ggtitle(title)
p <- p + geom_violin(alpha=0.5,trim=F)
supfig_clans_all=p
p
my_comparisons <- list(  c("all","motif"), c("all", "Beta_propeller"), c("all","C2H2-zf"), c("all","TPR"), c("all","LRR") )
p <- p+ stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",label="p.signif")
supfig_clans_all_sign=p
p

### Output figures to files
svg(file=Fig1_prd_score_file, width=6,height=3)
fig1
dev.off()

svg(Fig2_topclan_file, width=6,height=4)
fig2
dev.off()

pdf(SupFig_top_clans_file, width=12,height=7)
supfig_clans_full
supfig_clans_zoom
supfig_clans_all
supfig_clans_all_sign
fig2
dev.off()

##Draft figures and wilcox checks


### Zinc finger comparision
scores[scores$clan=="C2H2-zf","zf_bool"]="True" 
scores[scores$clan!="C2H2-zf","zf_bool"]="False" 

##Plot: compare zf c2h2 and not zf
#svg(paste0(path,"ZF_dupscore_scaled.svg"), width=6,height=3)
title="Zinc finger OGs show fast repeat evolution"
p <- ggplot(data = scores, aes(x=prd_score, ..scaled..))
p <- p + coord_cartesian(clip = 'off') + scale_x_continuous(trans=symlog_trans(thr=1))
p <- p + xlab("PRD score") + ylab("Density estimate,scaled") + ggtitle(title)
p <- p + geom_density(data=scores,adjust = 1.5,aes(group=zf_bool, color=zf_bool, fill=zf_bool),alpha=0.1)
p
#dev.off()
#non-scaled
#svg(paste0(path,"ZF_dupscore_nonscaled.svg"), width=6,height=3)
p <- ggplot(data = scores, aes(x=prd_score))
p <- p + coord_cartesian(clip = 'off') + scale_x_continuous(trans=symlog_trans(thr=1))
p <- p + xlab("PRD score") + ylab("Density estimate") + ggtitle(title)
p <- p + geom_density(data=scores,adjust = 1.5,aes(group=zf_bool, color=zf_bool, fill=zf_bool),alpha=0.1)
p
#dev.off()

##Significance
#nonparametric test, because not normally distributed
with(scores,shapiro.test(prd_score)) #shows clearly not normal
wilcox.test(prd_score ~ zf_bool, data = scores,alternative="greater")



p <- ggplot(data = scores, aes(x=prd_score))
p <- p + coord_cartesian(clip = 'off') + scale_x_continuous(trans=symlog_trans(thr=1))
p <- p + xlab("PRD score") + ylab("Density estimate") + ggtitle(title)
p <- p + geom_density(data=scores,adjust = 1.5,aes(group=rank, color=rank, fill=rank),alpha=0.1)
p



### Schaper
title="Comparision to positive set of Schaper et al "
p <- ggplot(data = scores, aes(x=prd_score,))
p <- p + coord_cartesian(clip = 'off') + scale_x_continuous(trans=symlog_trans(thr=1))
p <- p + xlab("PRD score") + ylab("Density estimate,scaled") + ggtitle(title)
p <- p + geom_density(data=scores[scores$schaper_pos!="",],adjust = 1.5,aes(group=schaper_pos, color=schaper_pos, fill=schaper_pos),alpha=0.1)
p

wilcox.test(prd_score ~ schaper_pos, data = scores[scores$schaper_pos!="",],alternative="greater")

##Selectome
title="Comparision to Selectome db of positive selection "
p <- ggplot(data = scores, aes(x=prd_score,..scaled..))
p <- p + coord_cartesian(clip = 'off') + scale_x_continuous(trans=symlog_trans(thr=1))
p <- p + xlab("PRD score") + ylab("Density estimate,scaled") + ggtitle(title)
p <- p + geom_density(data=scores,adjust = 1.5,aes(group=factor(selectome), color=factor(selectome)),alpha=0.1)
p

title="Comparision to Selectome db of positive selection "
p <- ggplot(data = scores[scores$human_dup>0,], aes(x=prd_score,..scaled..))
p <- p + coord_cartesian(clip = 'off') + scale_x_continuous(trans=symlog_trans(thr=1))
p <- p + xlab("PRD score") + ylab("Density estimate,scaled") + ggtitle(title)
p <- p + geom_density(data=scores[scores$human_dup>0,],adjust = 1.5,aes(group=factor(selectome), color=factor(selectome)),alpha=0.1)
p

## Gene duplication
ggplot(data=scores, aes(x=prd_score, y=gt_score))+geom_point()

title="Gene duplication "
p <- ggplot(data = scores, aes(x=prd_score,))
p <- p + coord_cartesian(clip = 'off') + scale_x_continuous(trans=symlog_trans(thr=1))
p <- p + xlab("PRD score") + ylab("Density estimate,scaled") + ggtitle(title)
p <- p + geom_density(data=scores,adjust = 1.5, aes(group=gt_bool, color=gt_bool), alpha=0.1)
p

wilcox.test(prd_score ~ gt_bool, data = scores)




#EXAC

ggplot(data=scores, aes(x=prd_score, y=mis_z))+geom_point()
ggplot(data=scores, aes(x=prd_score, y=pLI))+geom_point()

title="Exac "
p <- ggplot(data = scores, aes(x=mis_z,))
#p <- p + coord_cartesian(clip = 'off') + scale_x_continuous(trans=symlog_trans(thr=1))
p <- p + xlab("ExAc, mis_z") + ylab("Density estimate") + ggtitle(title)
p <- p + geom_density(data=scores,adjust = 1.5, aes(group=rank, color=rank), alpha=0.5)
p

title="Exac "
p <- ggplot(data = scores, aes(x=pLI,))
#p <- p + coord_cartesian(clip = 'off') + scale_x_continuous(trans=symlog_trans(thr=1))
p <- p + xlab("ExAC, pLI") + ylab("Density estimate") + ggtitle(title)
p <- p + geom_density(data=scores,adjust = 1, aes(group=rank, color=rank), alpha=0.5)
p


title="Exac "
p <- ggplot(data = scores, aes(x=mis_z,))
#p <- p + coord_cartesian(clip = 'off') + scale_x_continuous(trans=symlog_trans(thr=1))
p <- p + xlab("ExAc, mis_z") + ylab("Density estimate") + ggtitle(title)
p <- p + geom_density(data=scores,adjust = 1.5, aes(group=rank_humlin, color=rank_humlin), alpha=0.5)
p

title="Exac "
p <- ggplot(data = scores, aes(x=pLI,))
#p <- p + coord_cartesian(clip = 'off') + scale_x_continuous(trans=symlog_trans(thr=1))
p <- p + xlab("ExAC, pLI") + ylab("Density estimate") + ggtitle(title)
p <- p + geom_density(data=scores,adjust = 1, aes(group=rank_humlin, color=rank_humlin), alpha=0.5)
p


wilcox.test(selectome ~ rank_humlin, data = scores,alternative="less")


title="Comparision to segdup db of positive selection "
p <- ggplot(data = scores, aes(x=prd_score,..scaled..))
p <- p + coord_cartesian(clip = 'off') + scale_x_continuous(trans=symlog_trans(thr=1))
p <- p + xlab("PRD score") + ylab("Density estimate,scaled") + ggtitle(title)
p <- p + geom_density(data=scores,adjust = 1.5,aes(group=factor(segdup), color=factor(segdup)),alpha=0.1)
p

wilcox.test(prd_score ~ segdup, data = scores,alternative="greater")

wilcox.test(prd_score ~ segdup, data = scores[scores$clan=="C2H2-zf",],alternative="greater")


title="Comparision to segdup db of positive selection, zf only "
p <- ggplot(data = scores[scores$clan=="C2H2-zf",], aes(x=prd_score,..scaled..))
p <- p + coord_cartesian(clip = 'off') + scale_x_continuous(trans=symlog_trans(thr=1))
p <- p + xlab("PRD score") + ylab("Density estimate,scaled") + ggtitle(title)
p <- p + geom_density(data=scores[scores$clan=="C2H2-zf",],adjust = 1.5,aes(group=factor(segdup), color=factor(segdup)),alpha=0.1)
p


title="Comparision to segdup db of positive selection, motif only "
p <- ggplot(data = scores[scores$clan=="motif",], aes(x=prd_score,..scaled..))
p <- p + coord_cartesian(clip = 'off') + scale_x_continuous(trans=symlog_trans(thr=1))
p <- p + xlab("PRD score") + ylab("Density estimate,scaled") + ggtitle(title)
p <- p + geom_density(data=scores[scores$clan=="motif",],adjust = 1.5,aes(group=factor(segdup), color=factor(segdup)),alpha=0.1)
p


wilcox.test(prd_score ~ segdup, data = scores[scores$clan=="motif",],alternative="greater")

dev.off()
