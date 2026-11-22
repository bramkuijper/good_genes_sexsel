#!/usr/bin/env Rscript

biasv = c(0.55, 0.65, 0.75,0.8,0.85, 0.95)

nrep = 20

maxgen = 50000


# generate a date_time stamp as a character
date_str <- format(Sys.time(), "%d%m%Y_%H%M%S")

# generate names for the output file based on date and time 
# to ensure some sort of uniqueness, to avoiding that output
# files from different runs overwrite eachother
output_file_prefix = paste0("sim_good_genes_",date_str)

counter <- 0

exe = "./good_genes.exe"

batch_file_contents <- ""

a <- 1.0
b <- 0.0025
c <- 0.5

init_t <- 0.0
init_p <- 1.0
init_v <- 5.0

max_mut <- c(0.4) #c(0.05,0.1,0.2,0.3)

max_mut_v <- 0.8

only_positive <- F

for (max_mut_i in max_mut)
{
    for (rep_i in 1:nrep)
    {
        for (biasv_i in biasv)
        {
            counter <- counter + 1
            file_name_i <- paste0(output_file_prefix,"_",counter)

            echo_str <- paste("echo",counter)

            max_mut_t <- max_mut_i
            max_mut_p <- max_mut_i

            command_str <- paste(exe,
                            biasv_i,
                            a,
                            b,
                            c,
                            init_t,
                            init_p,
                            init_v,
                            max_mut_t,
                            max_mut_p,
                            max_mut_v,
                            format(maxgen,scientific=F),
                            as.numeric(only_positive),
                            file_name_i)

            # append to batch file contents
            batch_file_contents <- paste0(batch_file_contents
                    ,"\n"
                    ,echo_str
                    ,"\n"
                    ,command_str)
        }
    }
}

writeLines(text=batch_file_contents)
