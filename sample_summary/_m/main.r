suppressMessages({library(SummarizedExperiment)
                  library(tidyverse)
                  library(ggpubr)})

save_ggplots <- function(p, fn, w, h){
    for(ext in c('.pdf', '.svg')){
        ggsave(paste0(fn, ext), plot=p, width=w, height=h)
    }
}

# Load counts and phenotype R variable
load("../../input/counts/_m/caudate_brainseq_phase3_hg38_rseGene_merged_n464.rda")
### Subset and recode
keepIndex = which((rse_gene$Dx %in% c('Control', "Schizo")) & 
                  rse_gene$Race %in% c('CAUC', 'AA'))
rse_gene = rse_gene[, keepIndex]
### Extract phenotypes
pheno_C <- colData(rse_gene) %>% as.data.frame

# Load counts and phenotype R variable
load("../../input/counts/_m/dlpfc_ribozero_brainseq_phase2_hg38_rseGene_merged_n453.rda")
### Subset and recode
keepIndex = which((rse_gene$Dx %in% c('Control', "Schizo")) & 
                  rse_gene$Race %in% c('CAUC', 'AA'))
rse_gene = rse_gene[, keepIndex]
### Extract phenotypes
pheno_D <- colData(rse_gene) %>% as.data.frame

# Load counts and phenotype R variable
load("../../input/counts/_m/hippo_brainseq_phase2_hg38_rseGene_merged_n447.rda")
### Subset and recode
keepIndex = which((rse_gene$Dx %in% c('Control', "Schizo")) & 
                  rse_gene$Race %in% c('CAUC', 'AA'))
rse_gene = rse_gene[, keepIndex]
### Extract phenotypes
pheno_H <- colData(rse_gene) %>% as.data.frame

# Load counts and phenotype R variable
load("../../input/counts/_m/astellas_dg_hg38_rseGene_n263.rda")
### Subset and recode
keepIndex = which((rse_gene$Dx %in% c('Control', "Schizo")) & 
                  rse_gene$Race %in% c('CAUC', 'AA'))
rse_gene = rse_gene[, keepIndex]
### Extract phenotypes
pheno_dg <- colData(rse_gene) %>% as.data.frame

allCols <- intersect(intersect(intersect(colnames(pheno_C), colnames(pheno_D)), 
                               colnames(pheno_H)), 
                     colnames(pheno_dg))
pheno = rbind(pheno_C[, allCols], pheno_D[, allCols], 
              pheno_H[, allCols], pheno_dg[, allCols]) %>% 
    filter(Age > 17) %>% mutate(Race=gsub("CAUC", "EA", Race))

ancestry = data.table::fread("../../input/ancestry_structure/structure.out_ancestry_proportion_raceDemo_compare")
ancestry %>% head(2)

ancestry %>% mutate_if(is.character, as.factor) %>%
    group_by(group) %>% summarize(AA=mean(Afr), EA=mean(Eur))

ancestry %>% inner_join(pheno, by=c("id"="BrNum")) %>%
    filter(Age > 17, Dx == "Control") %>% select(group, Afr, Eur) %>% 
    mutate_if(is.character, as.factor) %>% distinct %>%
    group_by(group) %>% 
    summarize(AA_mean=mean(Afr), AA_sd=sd(Afr), AA_max=max(Afr), AA_min=min(Afr),
              EA_mean=mean(Eur), EA_sd=sd(Eur), EA_max=max(Eur), EA_min=min(Eur))

brp = ancestry %>% inner_join(pheno, by=c("id"="BrNum")) %>%
    filter(Age > 17, Dx == "Control") %>% select(id, Race, Afr, Eur) %>%
    mutate_if(is.character, as.factor) %>% distinct %>%
    pivot_longer(-c("Race", "id"), names_to="Ancestry", values_to="Proportion") %>% 
    mutate_if(is.character, as.factor) %>% group_by(Ancestry) %>% 
    mutate(ID = fct_reorder(id, desc(Proportion))) %>%
    ggbarplot(x="ID", y="Proportion", fill = "Ancestry", color="Ancestry",
              palette="npg", ylab="Admixture", xlab="292 Individuals", 
              ggtheme=theme_pubr(base_size=20), legend="right") +
    geom_hline(yintercept=0.5, linetype="dashed", color="white") +
    geom_hline(yintercept=0.75, linetype="dashed", color="black") +
    geom_hline(yintercept=0.25, linetype="dashed", color="black") +
    font("xy.title", face="bold") + font("legend.title", face="bold") +
    rremove("x.text") + rremove("x.ticks")
save_ggplots(brp, "ancestry_structure_barplot", 12, 5)
brp

bxp = ancestry %>% inner_join(pheno, by=c("id"="BrNum")) %>%
    filter(Age > 17, Dx == "Control") %>% select(id, Race, Afr, Eur) %>%
    mutate_if(is.character, as.factor) %>% distinct %>%
    pivot_longer(-c("Race", "id"), names_to="Ancestry", values_to="Proportion") %>% 
    ggdensity(x="Proportion", color="Race", fill="Race", facet.by="Ancestry", 
              ncol=2, rug=TRUE, add="mean", palette="npg", ylab="Population Density", 
              xlab="Ancestry Proportion", panel.labs.font=list(face='bold'), 
              ggtheme=theme_pubr(base_size=15, border=TRUE)) + 
    font("xy.title", face="bold") + font("legend.title", face="bold")
save_ggplots(bxp, "ancestry_structure_distribution", 10, 5)
bxp

pheno %>% dim

print(paste("There are", unique(pheno$BrNum) %>% length, "unique BrNum."))

pheno %>% select(BrNum, Region) %>% distinct %>%
    mutate_if(is.character, as.factor) %>% 
    group_by(Region) %>% count()

pheno %>% select(BrNum, Race) %>% distinct %>%
    mutate_if(is.character, as.factor) %>% 
    group_by(Race) %>% count()

pheno %>% select(BrNum, Race, Region) %>% distinct %>%
    mutate_if(is.character, as.factor) %>% 
    group_by(Region, Race) %>% count()

pheno %>% select(BrNum, Sex, Region) %>% distinct %>%
    mutate_if(is.character, as.factor) %>% 
    group_by(Region, Sex) %>% count()

pheno %>% group_by(Region) %>% 
  summarise_at(vars(c("Age")), list(mean = mean, sd = sd)) 

pheno %>% group_by(Region, Race) %>% 
  summarise_at(vars(c("Age")), list(mean = mean, sd = sd)) 

pheno %>% filter(RIN != "NA") %>% mutate("RIN"=as.numeric(unlist(RIN))) %>% 
    group_by(Region) %>% summarise_at(vars(c("RIN")), list(mean = mean, sd = sd)) 

pheno %>% filter(RIN != "NA") %>% mutate("RIN"=as.numeric(unlist(RIN))) %>% 
    group_by(Region, Race) %>% summarise_at(vars(c("RIN")), list(mean = mean, sd = sd)) 

plot_pie <- function(tissue){
    pie = pheno %>% mutate_if(is.character, as.factor) %>% group_by(Region, Race) %>%
        count %>% as.data.frame %>% group_by(Region) %>%
        transmute(Race, Percent = round(n/sum(n)*100, 1)) %>%
        mutate(Labels=paste0(Race, " (", Percent, "%)")) %>% filter(Region == tissue) %>%
        ggpie("Percent", label="Labels", fill="Race", color="white", palette="npg", 
              lab.pos="in", lab.font=c(8, "bold", "white"),
              ggtheme=theme_pubr(base_size=20, legend="none"))
    return(pie)
}

## Get and annotate plot
cc_pie = annotate_figure(plot_pie("Caudate"), 
                         top = text_grob("Caudate", face = "bold", size = 26))
gg_pie = annotate_figure(plot_pie("DentateGyrus"), 
                         top = text_grob("Dentate Gyrus", face = "bold", size = 26))
dd_pie = annotate_figure(plot_pie("DLPFC"), 
                         top = text_grob("DLPFC", face = "bold", size = 26))
hh_pie = annotate_figure(plot_pie("HIPPO"), 
                         top = text_grob("Hippocampus", face = "bold", size = 26))
## Arrange figure
figure <- ggarrange(cc_pie, gg_pie, dd_pie, hh_pie, ncol = 2, nrow = 2)
save_ggplots(figure, "ancestry_piecharts", 10, 10)
figure

pheno = pheno %>% filter(Age > 17, Dx == "Control", Race == "AA")
pheno %>% dim

print(paste("There are", unique(pheno$BrNum) %>% length, "unique BrNum."))

pheno %>% select(BrNum, Region) %>% distinct %>%
    mutate_if(is.character, as.factor) %>% 
    group_by(Region) %>% count()

pheno %>% select(BrNum, Race) %>% distinct %>%
    mutate_if(is.character, as.factor) %>% 
    group_by(Race) %>% count()

pheno %>% select(BrNum, Race, Region) %>% distinct %>%
    mutate_if(is.character, as.factor) %>% 
    group_by(Region, Race) %>% count()

pheno %>% select(BrNum, Sex, Region) %>% distinct %>%
    mutate_if(is.character, as.factor) %>% 
    group_by(Region, Sex) %>% count()

pheno %>% group_by(Region) %>% 
  summarise_at(vars(c("Age")), list(mean = mean, sd = sd)) 

pheno %>% group_by(Region, Race) %>% 
  summarise_at(vars(c("Age")), list(mean = mean, sd = sd)) 

pheno %>% filter(RIN != "NA") %>% mutate("RIN"=as.numeric(unlist(RIN))) %>% 
    group_by(Region) %>% summarise_at(vars(c("RIN")), list(mean = mean, sd = sd)) 

pheno %>% filter(RIN != "NA") %>% mutate("RIN"=as.numeric(unlist(RIN))) %>% 
    group_by(Region, Race) %>% summarise_at(vars(c("RIN")), list(mean = mean, sd = sd)) 

Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
