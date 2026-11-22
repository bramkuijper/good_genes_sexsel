library("tidyverse")
library("patchwork")

#the_data <- read_delim("summary_good_genes.csv",delim=";")
#the_data <- read_delim("summary_goodgenes_nonneg_mutation.csv",delim=";")
the_file <- "summary_good_genes_bias_no_bias.csv"
the_file_no_ext <- tools::file_path_sans_ext(the_file)

the_data <- read_delim(the_file,delim=";")



the_data_l <- pivot_longer(the_data
                           ,cols=c(meant,meanp)
                           ,names_to = "trait"
                           ,values_to = "trait_value")

p1 <- ggplot(data = the_data_l
       ,mapping = aes(x = biasv, y = trait_value)) +
       geom_hline(yintercept=0) +
    geom_point(mapping=aes(colour=trait),alpha=0.5) +
    scale_colour_brewer(palette="Set1") +
    theme_classic() +
    xlab("Probability of biased mutations on viability trait") +
    ylab("Trait value") 

p2 <- ggplot(data = the_data_l
       ,mapping = aes(x = biasv, y = abs(trait_value))) +
    geom_point(mapping=aes(colour=trait)) +
    scale_colour_brewer(palette="Set1") +
    theme_classic() +
    xlab("Probability of biased mutations on viability trait") +
    ylab("Absolute trait value") 

p1 | p2


ggsave(filename=paste0(the_file_no_ext,"_","overview.pdf"),width=12,height=5,device = cairo_pdf)

the_data_l <- pivot_longer(the_data
                           ,cols=c(covtp,covtv,covpv)
                           ,names_to = "trait"
                           ,values_to = "trait_value")

ggplot(data = the_data_l
       ,mapping = aes(x = biasv, y = trait_value)) +
    geom_point(mapping=aes(colour=trait)) +
    scale_colour_brewer(palette="Set1") +
    theme_classic() +
    xlab("Probability of biased mutations on viability trait") +
    ylab("covariance")

ggsave(filename=paste0(the_file_no_ext,"_","overview_covariances.pdf"),width=12,height=5,device = cairo_pdf)
