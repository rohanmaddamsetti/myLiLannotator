## simple-ARG-duplication-analysis.R by Rohan Maddamsetti.
## analyse the distribution of antibiotic resistance genes (ARGs)
## on chromosomes versus plasmids in fully-sequenced genomes and plasmids
## in the NCBI RefSeq database (dated March 26 2021).

## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## CRITICAL TODO: Fix upstream bug where Anthropogenic is Anthropogeni.
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## CRITICAL TODO: Fix upstream bug where Environmental is Environmenta

library(tidyverse)
library(cowplot)
library(data.table)
library(stringi)

################################################################################
## Regular expressions used in this analysis.

## match MGE genes using the following keywords in the "product" annotation
transposon.keywords <- "IS|transpos\\S*|insertion|Tra[A-Z]|Tra[0-9]|tra[A-Z]|conjugate transposon|Transpos\\S*|Tn[0-9]|tranposase|Tnp|Ins|ins"
plasmid.keywords <- "relax\\S*|conjug\\S*|mob\\S*|plasmid|type IV|chromosome partitioning|chromosome segregation|Mob\\S*|Plasmid|Rep|Conjug\\S*"
phage.keywords <- "capsid|phage|Tail|tail|head|tape measure|antiterminatio|Phage|virus|Baseplate|baseplate|coat|entry exclusion"
other.HGT.keywords <- "Integrase|integrase|excision\\S*|exonuclease|recomb|toxin|restrict\\S*|resolv\\S*|topoisomerase|reverse transcrip|intron|antitoxin|toxin|Toxin|Reverse transcriptase|hok|Hok|competence|addiction"

MGE.keywords <- paste(transposon.keywords, plasmid.keywords, phage.keywords, other.HGT.keywords, sep="|")

## antibiotic-specific keywords.
chloramphenicol.keywords <- "chloramphenicol|Chloramphenicol"
tetracycline.keywords <- "tetracycline efflux|Tetracycline efflux|TetA|Tet(A)|tetA|tetracycline-inactivating"
MLS.keywords <- "macrolide|lincosamide|streptogramin"
multidrug.keywords <- "Multidrug resistance|multidrug resistance|antibiotic resistance"
beta.lactam.keywords <- "lactamase|LACTAMASE|beta-lactam|oxacillinase|carbenicillinase|betalactam\\S*"
glycopeptide.keywords <- "glycopeptide resistance|VanZ|vancomycin resistance|VanA|VanY|VanX|VanH|streptothricin N-acetyltransferase"
polypeptide.keywords <- "bacitracin|polymyxin B|phosphoethanolamine transferase|phosphoethanolamine--lipid A transferase"
diaminopyrimidine.keywords <- "trimethoprim|dihydrofolate reductase|dihydropteroate synthase"
sulfonamide.keywords <- "sulfonamide|Sul1|sul1|sulphonamide"
quinolone.keywords <- "quinolone|Quinolone|oxacin|qnr|Qnr"
aminoglycoside.keywords <- "Aminoglycoside|aminoglycoside|streptomycin|Streptomycin|kanamycin|Kanamycin|tobramycin|Tobramycin|gentamicin|Gentamicin|neomycin|Neomycin|16S rRNA (guanine(1405)-N(7))-methyltransferase|23S rRNA (adenine(2058)-N(6))-methyltransferase|spectinomycin 9-O-adenylyltransferase|Spectinomycin 9-O-adenylyltransferase|Rmt"
macrolide.keywords <- "macrolide|ketolide|Azithromycin|azithromycin|Clarithromycin|clarithromycin|Erythromycin|erythromycin|Erm|EmtA"
antimicrobial.keywords <- "QacE|Quaternary ammonium|quaternary ammonium|Quarternary ammonium|quartenary ammonium|fosfomycin|ribosomal protection|rifampin ADP-ribosyl|azole resistance|antimicrob\\S*"

antibiotic.keywords <- paste(chloramphenicol.keywords, tetracycline.keywords, MLS.keywords, multidrug.keywords,
    beta.lactam.keywords, glycopeptide.keywords, polypeptide.keywords, diaminopyrimidine.keywords,
    sulfonamide.keywords, quinolone.keywords, aminoglycoside.keywords, macrolide.keywords, antimicrobial.keywords, sep="|")

antibiotic.or.MGE.keywords <- paste(MGE.keywords,antibiotic.keywords,sep="|")

################################################################################
## Functions

fancy_scientific <- function(x) {
    ## function for plotting better axis labels.
    ## see solution here for nice scientific notation on axes.
    ## https://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales
    ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x)))))
}

## return the first column for several tables.
## shows the number of isolates in each category.
make.isolate.totals.col <- function(gbk.annotation) {
    isolate.totals <- gbk.annotation %>%
        group_by(Annotation) %>%
        summarize(total_isolates = n()) %>%
        arrange(desc(total_isolates))
    return(isolate.totals)
}


read.LLM.gbk.annotation.csv <- function(gbk.annotation.path, ground.truth.gbk.annotation) {
    
    ## only select Annotation_Accession, Organism, Strain, Genus columns from ground_truth_gbk.annotation.
    relevant.ground.truth <- ground.truth.gbk.annotation %>%
        select(Annotation_Accession, Organism, Strain, Genus)
    
    gbk.annotation.path %>%
        read.csv() %>%
        as_tibble() %>%
        ## filter based on ground truth genomes
        inner_join(relevant.ground.truth) %>%
        ## refer to NA annotations as "Unannotated".
        mutate(Annotation = replace_na(Annotation,"Unannotated")) %>%
        ## collapse Annotations into a smaller number of categories as follows:
        ## Marine, Freshwater --> Water
        ## Sediment, Soil, Terrestrial --> Earth
        ## Plants, Agriculture, Animals --> Plants & Animals
        ## Anthropogenic -> Human-impacted
        mutate(Annotation = replace(Annotation, Annotation == "Marine", "Water")) %>%
        mutate(Annotation = replace(Annotation, Annotation == "Freshwater", "Water")) %>%
        mutate(Annotation = replace(Annotation, Annotation == "Sediment", "Earth")) %>%
        mutate(Annotation = replace(Annotation, Annotation == "Soil", "Earth")) %>%
        mutate(Annotation = replace(Annotation, Annotation == "Terrestrial", "Earth")) %>%
        mutate(Annotation = replace(Annotation, Annotation == "Plants", "Plants & Animals")) %>%
        mutate(Annotation = replace(Annotation, Annotation == "Agriculture", "Plants & Animals")) %>%
        mutate(Annotation = replace(Annotation, Annotation == "Animals", "Plants & Animals")) %>%
        ## CRITICAL TODO: FIX THIS UPSTREAM BUG
        ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    mutate(Annotation = replace(Annotation, Annotation == "Anthropogeni", "Human-impacted")) %>%
##        mutate(Annotation = replace(Annotation, Annotation == "Anthropogenic", "Human-impacted")) %>%
        ## And now remove all Unannotated genomes, since these are not analyzed
        ## at all in this first paper.
        filter(Annotation != "Unannotated") %>%
        ## and remove any strains (although none should fall in this category)
        ## that were not annotated by annotate-ecological-category.py.
        filter(Annotation != "blank")
}


##########################################################################
## Functions for Figure 1.

## See Wikipedia reference:
## https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval

## Make Z-distributed confidence intervals for the fraction of isolates with
## duplicated ARGs (panel A),
## the fraction of isolates with single-copy ARGs (panel B),
## the fraction of isolates with duplicated genes (panel C).

## Count data for isolates with duplicated ARGs
## goes into Supplementary Table S1.

calc.isolate.confints <- function(df) {
    df %>%
        ## use the normal approximation for binomial proportion conf.ints
        mutate(se = sqrt(p*(1-p)/total_isolates)) %>%
        ## See Wikipedia reference:
        ## https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
        mutate(Left = p - 1.96*se) %>%
        mutate(Right = p + 1.96*se) %>%
        ## truncate confidence limits to interval [0,1].
        rowwise() %>% mutate(Left = max(0, Left)) %>%
        rowwise() %>% mutate(Right = min(1, Right)) %>%
        ## Sort every table by the total number of isolates.
        arrange(desc(total_isolates))
}


make.TableS1 <- function(gbk.annotation, duplicate.ARGs) {

    ## count the number of isolates with duplicated ARGs in each category.
    ARG.category.counts <- duplicate.ARGs %>%
        ## next two lines is to count isolates rather than genes
        select(Annotation_Accession, Annotation) %>%
        distinct() %>%
        count(Annotation, sort = TRUE) %>%
        rename(isolates_with_duplicated_ARGs = n)
    
    ## join columns to make Table S1.
    TableS1 <- make.isolate.totals.col(gbk.annotation) %>%
        left_join(ARG.category.counts) %>%
        mutate(isolates_with_duplicated_ARGs =
                   replace_na(isolates_with_duplicated_ARGs,0)) %>%
        mutate(p = isolates_with_duplicated_ARGs/total_isolates) %>%
        calc.isolate.confints()
    
    return(TableS1)
}


## generic version of make.TableS1, for examining classes of genes other than
## antibiotic resistance genes.
make.IsolateEnrichmentTable <- function(gbk.annotation, duplicate.proteins, keywords) {
    ## count the number of isolates with duplicated genes of interest in each category.
    category.counts <- duplicate.proteins %>%
        filter(str_detect(.$product, keywords)) %>%
        ## next two lines is to count isolates rather than genes
        select(Annotation_Accession, Annotation) %>%
        distinct() %>%
        count(Annotation, sort = TRUE) %>%
        rename(isolates_with_duplicated_function = n)
    
    ## join columns to make the Table.
    Table <- make.isolate.totals.col(gbk.annotation) %>%
        left_join(category.counts) %>%
        mutate(isolates_with_duplicated_function =
                   replace_na(isolates_with_duplicated_function,0)) %>%
        mutate(p = isolates_with_duplicated_function/total_isolates) %>%
        calc.isolate.confints()
    return(Table)
}


make.IsolateEnrichmentControlTable <- function(gbk.annotation, singleton.proteins, keywords) {
    ## count the number of isolates with singleton genes of interest in each category.
    category.counts <- singleton.proteins %>%
        filter(str_detect(.$product, keywords)) %>%
        ## next two lines is to count isolates rather than genes
        select(Annotation_Accession, Annotation) %>%
        distinct() %>%
        count(Annotation, sort = TRUE) %>%
        rename(isolates_with_singleton_function = n)
    
    ## join columns to make the Table.
    Table <- make.isolate.totals.col(gbk.annotation) %>%
        left_join(category.counts) %>%
        mutate(isolates_with_singleton_function =
                   replace_na(isolates_with_singleton_function,0)) %>%
        mutate(p = isolates_with_singleton_function/total_isolates) %>%
        calc.isolate.confints()
    return(Table)
}


make.confint.figure.panel <- function(Table, order.by.total.isolates, title,
                                      no.category.label = FALSE) {    
    Fig.panel <- Table %>%
        mutate(Annotation = factor(
                   Annotation,
                   levels = rev(order.by.total.isolates))) %>%
        ggplot(aes(y = Annotation, x = p)) +
        geom_point(size=1) +
        ylab("") +
        xlab("Proportion of Isolates") +
        theme_classic() +
        ggtitle(title) +
        ## plot CIs.
        geom_errorbarh(aes(xmin=Left,xmax=Right), height=0.2, linewidth=0.2)
    
    if (no.category.label)
        Fig.panel <- Fig.panel +
            theme(axis.text.y=element_blank())
    
    return(Fig.panel)
}


make.TableS2 <- function(gbk.annotation, singleton.ARGs) {

## count the number of isolates with singleton AR genes in each category.
    ARG.category.counts <- singleton.ARGs %>%
        ## next two lines is to count isolates rather than genes
        select(Annotation_Accession, Annotation) %>%
        distinct() %>%
        group_by(Annotation) %>%
        summarize(isolates_with_singleton_ARGs = n()) %>%
        arrange(desc(isolates_with_singleton_ARGs))
    gc() ## free memory.
    
    ## join columns to make Table S2.
    TableS2 <- make.isolate.totals.col(gbk.annotation) %>%
        left_join(ARG.category.counts) %>%
        mutate(isolates_with_singleton_ARGs =
                   replace_na(isolates_with_singleton_ARGs,0)) %>%
        mutate(p = isolates_with_singleton_ARGs/total_isolates) %>%
        calc.isolate.confints()
    return(TableS2)
}


make.TableS3 <- function(gbk.annotation, duplicate.proteins) {
    ## count the number of isolates with duplicated genes in each category.
    category.counts <- duplicate.proteins %>%
        ## next two lines is to count isolates rather than genes
        select(Annotation_Accession, Annotation) %>%
        distinct() %>%
        group_by(Annotation) %>%
        summarize(isolates_with_duplicated_genes = n()) %>%
        arrange(desc(isolates_with_duplicated_genes))
    
    ## join columns to make Table S3.
    TableS3 <- make.isolate.totals.col(gbk.annotation) %>%
        left_join(category.counts) %>%
        mutate(isolates_with_duplicated_genes =
                   replace_na(isolates_with_duplicated_genes, 0)) %>%
        mutate(p = isolates_with_duplicated_genes/total_isolates) %>%
        calc.isolate.confints()
    return(TableS3)
}


################################################################################
## Set up the key data structures for the analysis.

ground.truth.gbk.annotation <- read.csv(
    "../data/Maddamsetti2024/FileS3-Complete-Genomes-with-Duplicated-ARG-annotation.csv") %>%
    as_tibble()

## This vector is used for ordering axes in figures and tables.
order.by.total.isolates <- make.isolate.totals.col(ground.truth.gbk.annotation)$Annotation

ground.truth.duplicate.proteins <- data.table::fread("../data/Maddamsetti2024/duplicate-proteins.csv",
                                                     drop="sequence") %>%
    ## now merge with gbk annotation.
    inner_join(ground.truth.gbk.annotation, by="Annotation_Accession")
        
ground.truth.duplicate.ARGs <- ground.truth.duplicate.proteins %>%
    filter(str_detect(product, antibiotic.keywords))


## Import llama3.2 ecological annotations
llama3.2.gbk.annotation <- read.LLM.gbk.annotation.csv(
    "../results/llama3.2_latest_gbk-annotation-table.csv",
    ground.truth.gbk.annotation)

llama3.2.duplicate.proteins <- data.table::fread("../data/Maddamsetti2024/duplicate-proteins.csv",
                                                 drop="sequence") %>%
    ## now merge with gbk annotation.
    inner_join(llama3.2.gbk.annotation, by="Annotation_Accession")

llama3.2.duplicate.ARGs <- llama3.2.duplicate.proteins %>%
    filter(str_detect(product, antibiotic.keywords))


## Now look at singleton proteins.
ground.truth.singleton.proteins <- data.table::fread("../data/Maddamsetti2024/all-proteins.csv",
                                                     drop="sequence") %>%
    filter(count == 1) %>%
    inner_join(ground.truth.gbk.annotation, by="Annotation_Accession")

## uses 5.7Gb of memory.

print(object.size(ground.truth.singleton.proteins), units = "auto")

## THIS LINE IS VERY SLOW.
ground.truth.singleton.ARGs <- ground.truth.singleton.proteins[stri_detect_regex(product, antibiotic.keywords)]

## get LifestyleAnnotation, rename to Annotation for existing code to work nicely.
gbk.reannotation <- read.csv("../results/llama3.2_latest_Complete-Genomes-with-lifestyle-annotation.csv") %>%
    select(-hasDuplicatedARGs, -Annotation) %>%
    rename(Annotation = LifestyleAnnotation) %>%
    ## refer to NA annotations as "Unannotated".
    mutate(Annotation = replace_na(Annotation,"Unannotated")) %>%
    ## And now remove all Unannotated genomes, since these are not analyzed
    ## at all in this first paper.
    filter(Annotation != "Unannotated") %>%
    ## CRITICAL TODO: FIX THIS UPSTREAM BUG
    ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    mutate(Annotation = replace(Annotation, Annotation == "Environmenta", "Environmental"))


## and merge.
reannotated.singleton.ARGs <- ground.truth.singleton.ARGs %>%
    select(-Annotation) %>%
    left_join(gbk.reannotation)

reannotated.duplicate.ARGs <- ground.truth.duplicate.ARGs %>%
    select(-Annotation) %>%
    left_join(gbk.reannotation)

reannotated.duplicate.proteins <- ground.truth.duplicate.proteins %>%
    select(-Annotation) %>%
    left_join(gbk.reannotation)


##########################################################################
## Data structures for Figure 1AB.

## Data structure for Figure 1AB:
## normal-approximation confidence intervals for the percentage
## of isolates with duplicated ARGs.
ground.truth.TableS1 <- make.TableS1(ground.truth.gbk.annotation, ground.truth.duplicate.ARGs)

llama3.2.TableS1 <- make.TableS1(llama3.2.gbk.annotation, llama3.2.duplicate.ARGs)

reannotated.TableS1 <- make.TableS1(gbk.reannotation, reannotated.duplicate.ARGs)

######################
## Table S2. Control: does the distribution of ARG singletons
## (i.e. genes that have NOT duplicated) follow the distribution
## of sampled isolates?

## No categories are enriched with ARG singletons,
## as most isolates have a gene that matches an antibiotic keyword.
## Animal-host isolates are depleted (perhaps due to aphid bacteria isolates?)


## This data frame will be used for Figure 1C.
TableS2 <- make.TableS2(ground.truth.gbk.annotation, ground.truth.singleton.ARGs)

## This data frame will be used for Figure 1D.
reannotated.TableS2 <- make.TableS2(gbk.reannotation, reannotated.singleton.ARGs)

#########################################################################
## Table S3. Control: does the number of isolates with duplicate genes
## follow the sampling distribution of isolates?

## Most follow the expected distribution.
## however, isolates from animal-hosts are signficantly depleted
## in duplicate genes: FDR-corrected p = 0.0000314
## while isolates from anthropogenic environments are weakly enriched
## in multi-copy genes: FDR-corrected p = 0.0212.


## Data structure for Figure 1E.
TableS3 <- make.TableS3(ground.truth.gbk.annotation, ground.truth.duplicate.proteins)

## Data structure for Figure 1G
reannotated.TableS3 <- make.TableS3(gbk.reannotation, reannotated.duplicate.proteins)
            
################################################################################
## Save Tables S1, S2, and S3 as Source Data for Fig1ABC.
##write.csv(TableS1, "../results/Source-Data/Fig1A-Source-Data.csv", row.names=FALSE, quote=FALSE)
##write.csv(TableS2, "../results/Source-Data/Fig1B-Source-Data.csv", row.names=FALSE, quote=FALSE)
##write.csv(TableS3, "../results/Source-Data/Fig1C-Source-Data.csv", row.names=FALSE, quote=FALSE)

## make Figure.
## Throughout, add special scales for Figure 1.

Fig1A <- make.confint.figure.panel(llama3.2.TableS1, order.by.total.isolates, "D-ARGs in genomes\nannotated by llama3.2") +
    scale_x_continuous(breaks = c(0, 0.15), limits = c(0,0.16))

Fig1B <- make.confint.figure.panel(ground.truth.TableS1, order.by.total.isolates, "Maddamsetti et al. (2024)\nD-ARGs") +
    scale_x_continuous(breaks = c(0, 0.15), limits = c(0,0.16))

Fig1C <- make.confint.figure.panel(reannotated.TableS1, new.order.by.total.isolates,
                                   "D-ARGs in genomes\nreannotated by llama3.2")

Fig1D <- make.confint.figure.panel(TableS2, order.by.total.isolates,
                                   "Maddamsetti et al. (2024)\nS-ARGs")

## This vector is used for ordering axes in figures and tables.
new.order.by.total.isolates <- make.isolate.totals.col(gbk.reannotation)$Annotation
Fig1E <- make.confint.figure.panel(reannotated.TableS2, new.order.by.total.isolates,
                                   "S-ARGs in genomes\nreannotated by llama3.2")


Fig1F <- make.confint.figure.panel(TableS3, order.by.total.isolates,
                                   "Maddamsetti et al. (2024)\nAll D-genes")

Fig1G <- make.confint.figure.panel(reannotated.TableS3, new.order.by.total.isolates,
                                   "All D-genes in genomes\nreannotated by llama3.2")

Fig1ABCDEFG <- plot_grid(
    Fig1A, NULL, Fig1B, Fig1C, Fig1D, Fig1E, Fig1F, Fig1G,
    labels=c('A', '', 'B', 'C', 'D', 'E', 'F', 'G'), nrow=4)

ggsave("../results/Fig1ABCDEFG.pdf", Fig1ABCDEFG, height=8,width=7)


## Adobe Firefly prompts:
## a muscular cyborg rainbow unicorn with a face stripe like ziggy stardust and a long rainbow mane working hard on a laptop
## a muscular cyborg rainbow llama with a face stripe like ziggy stardust and a long rainbow mane working hard on a laptop
