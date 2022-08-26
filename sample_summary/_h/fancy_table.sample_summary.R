## Summarize results and plot correlation with q-value

suppressPackageStartupMessages({
    library(dplyr)
    library(gtsummary)
})

save_table <- function(pp, fn){
    for(ext in c(".pdf", ".tex")){
        gt::gtsave(as_gt(pp), filename=paste0(fn,ext))
    }
}

load_phenotypes <- function(){
    return(data.table::fread("../../input/phenotypes/merged/_m/merged_phenotypes.csv") %>%
           mutate_if(is.character, as.factor))
}
memPHENO <- memoise::memoise(load_phenotypes)

#### MAIN
                                         # Generate pretty tables
fn1 = "sample_breakdown.table"
pp1 = memPHENO() %>% filter(Race == "AA", Age > 17, Dx == "Control") %>%
    select(Region, Sex, Age, RIN) %>%
    mutate(Sex=factor(Sex, labels=c("Female", "Male")),
           Region=factor(Region, labels=c("Caudate", "Dentate Gyrus",
                                          "DLPFC", "Hippocampus"))) %>%
    tbl_summary(by="Region", missing="no", #type=all_continuous() ~ "continuous2",
                statistic = all_continuous() ~ c("{mean} ({sd})")) %>%
    modify_header(all_stat_cols()~"**{level}**<br>N = {n}") %>%
    modify_spanning_header(all_stat_cols()~"**Brain Region**") %>%
    bold_labels() %>% italicize_levels()
save_table(pp1, fn1)
    
fn2 = "sample_breakdown.table.all_samples"
pp2 = memPHENO() %>% filter(Race %in% c("AA", "CAUC"), Age > 17, Dx == "Control") %>%
    select(Region, Race, Sex, Age, RIN) %>%
    mutate(Sex=factor(Sex, labels=c("Female", "Male")),
           Region=factor(Region, labels=c("Caudate", "Dentate Gyrus",
                                          "DLPFC", "Hippocampus")),
           Race=factor(forcats::fct_drop(Race), labels=c("AA", "EA"))) %>%
    tbl_strata(
        strata=Region,
        ~.x %>%
            tbl_summary(by="Race", missing="no",
                        statistic = all_continuous() ~ c("{mean} ({sd})")) %>%
            modify_header(all_stat_cols()~"**{level}**<br>N = {n}") %>%
            modify_spanning_header(all_stat_cols()~"**Brain Region**") %>%
            bold_labels() %>% italicize_levels()
    )
save_table(pp2, fn2)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
