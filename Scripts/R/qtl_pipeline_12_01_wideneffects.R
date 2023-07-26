library(tidyverse)

#Defaults
workflow        <- get0("workflow", ifnotfound="../../Workflows/1")
num_top_qtls    <- get0("num_top_qtls", ifnotfound=2) #Number of top QTLs to show per trait

args = commandArgs(trailingOnly=TRUE)

if(length(args)==0) {
    print("No arguments supplied.")
} else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

effs.tb <- read_csv(file=paste0(workflow,'/traits/effects_collated.csv'), show_col_types=FALSE)
effs.wide.tb <- effs.tb %>%
                    distinct() %>%
                    pivot_wider(names_from="genotype",values_from=c("effect_mean","effect_se"))
write_csv(effs.wide.tb, file=paste0(workflow,'/traits/effects_collated.wide.csv'))
