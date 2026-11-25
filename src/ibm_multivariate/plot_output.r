#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("tidyverse", warn.conflicts=F))
suppressPackageStartupMessages(library("jsonlite", warn.conflicts=F))
suppressPackageStartupMessages(library("patchwork", warn.conflicts=F))

# from a list of values like x1, x2, x3
# create a reasonable variable name, like x
make.var.name <- function(vars) {
    
    var1 <- vars[[1]]

    return(gsub(pattern="[_0-1]",replacement="",x=var1))
}

find.params <- function(filename) {

    f <- readLines(filename)

    seq <- seq(1,length(f),1)

    for (line_i in seq)
    {
        if (length(grep("^\\d",f[[line_i]])) > 0)
        {
            return(line_i - 1)
        }
    }

    return(NA)
}
# use tidyverse's sym() to transform strings to symbols 
transform.sym  <- function(x) 
{
    if (!is.null(x))
    {
        sym(x)
    }
}

if (!exists("file.name"))
{
    
    # get command line arguments
    args = commandArgs(trailingOnly=TRUE)
    
    # give an error message if you do not provide it with a simulation file name
    if (length(args) < 1)
    {
        print("provide a simulation file name")
        stop()
    }
    
    file.name <- args[[1]]
}

param.line <- find.params(file.name)

data.tibble <- read_delim(file=file.name
        ,delim=";"
        ,skip=param.line-1
        ,col_names=T)

# get the parameters
data.tibble.params <- read_delim(file=file.name
        ,delim=";"
        ,n_max=param.line-3
        ,col_names=c("name","value")
        )

# transpose the tibble with the parameters
params <- data.tibble.params %>% pivot_wider(
        names_from = name
        ,values_from = value)

##### column names for anything involving p #####

names_p <- grep(
        x = names(data.tibble),
        pattern = "^meanp.*",
        value = T)

names_p_str <- paste(names_p, collapse='","')

names_var_p <- grep(x = names(data.tibble),
        pattern = "^varp.*",
        value = T)

names_var_p_str <- paste(names_var_p, collapse='","')
##### column names for anything involving t #####

names_t <- grep(x = names(data.tibble),
        pattern = "^meant.*",
        value = T)

names_t_str <- paste(names_t, collapse='","')

names_var_t <- grep(x = names(data.tibble),
        pattern = "^vart.*",
        value = T)

names_var_t_str <- paste(names_var_t, collapse='","')

##### column names for anything involving x #####

names_x <- grep(x = names(data.tibble),
        pattern = "^meanx.*",
        value = T)

names_x_str <- paste(names_x, collapse='","')

names_var_x <- grep(x = names(data.tibble),
        pattern = "^varx.*",
        value = T)

names_var_x_str <- paste(names_var_x, collapse='","')


##### covariances #####

### cov p_i x v
cov_pv_names <- grep(x = names(data.tibble),
        pattern = "covp\\d+",
        value = T)

cov_pv_names_str <- paste(cov_pv_names, collapse = '","')

### cov t_i x p_i

cov_tp_names <- grep(x = names(data.tibble),
        pattern = "covt\\d+p\\d+",
        value = T)

cov_tp_names_str <- paste(cov_tp_names, collapse = '","')

### cov t_i x v

cov_tv_names <- grep(x = names(data.tibble),
        pattern = "covt\\d+v",
        value = T)

cov_tv_names_str <- paste(cov_tv_names, collapse = '","')


# now make the structure of the plot
jsonstuff <- paste0('[
    {"xvar" : "generation",
        "yvar" : ["',names_p_str,'"]},
    {"xvar" : "generation",
        "yvar" : ["',names_t_str,'"]},
    {"xvar" : "generation",
        "yvar" : ["',names_x_str,'"]},
    {
        "xvar" : "generation",
        "yvar" : "meanv"
    },
    {
        "xvar" : "generation",
        "yvar" : ["',names_var_p_str,'"]
    },
    {
        "xvar" : "generation",
        "yvar" : ["',names_var_t_str,'"]
    },
    {
        "xvar" : "generation",
        "yvar" : ["',names_var_x_str,'"]
    },
    {
        "xvar" : "generation",
        "yvar" : "varv"
    },
    {
        "xvar" : "generation",
        "yvar" : ["nf","nm"]
    },
    {
        "xvar" : "generation",
        "yvar" : ["p_surv_f","p_surv_m"]
    },
    {
        "xvar" : "generation",
        "yvar" : ["nfsurv","nmsurv"]
    },
    {
        "xvar" : "generation",
        "yvar" : ["',cov_tp_names_str,'"]
    },
    {
        "xvar" : "generation",
        "yvar" : ["',cov_tv_names_str,'"]
    },
    {
        "xvar" : "generation",
        "yvar" : ["',cov_pv_names_str,'"]
    }
]
')

plot.structure <- fromJSON(jsonstuff, simplifyVector = F)

plot.structure.l <- length(plot.structure)

# list with all the plots
plot.list <- list(rep(NA,times=plot.structure.l))

plot.list.idx <- 1

for (plot_struct_idx in 1:plot.structure.l)
{
    # get the (potential list of) y variable(s)
    # as this is a list and hence highly structured
    # hence, try to flatten it
    yvar <- unlist(plot.structure[[plot_struct_idx]]$yvar)

    if (length(yvar) > 1)
    {
        yvar_name <- make.var.name(yvar)
        yvar_values <- paste0(yvar_name,"_values")

        sub.data <- pivot_longer(data=data.tibble
                ,cols=yvar
                ,names_to=yvar_name
                ,values_to=yvar_values)

        # get rid of aes_string like this: 
        # https://stackoverflow.com/questions/74414272/how-to-replace-the-deprecated-ggplot2-function-aes-string-accepting-an-arbitrar/74414389#74414389 

        # aes arguments for the ggplot() call
        plot_args <- lapply(X=list(
                        x=plot.structure[[plot_struct_idx]]$xvar
                        ,y=yvar_values),
                FUN=transform.sym)

        # aes arguments for the geom_line() call
        line_args <- lapply(X=list(
                        colour=yvar_name),
                FUN=transform.sym)

        plot.list[[plot.list.idx]] <- ggplot(data=sub.data
                ,mapping=aes(!!!plot_args)) + 
                    geom_line(mapping=aes(!!!line_args))
    } else {
        
        # aes arguments for the ggplot() call
        plot_args <- lapply(X=list(
                        x=plot.structure[[plot_struct_idx]]$xvar
                        ,y=plot.structure[[plot_struct_idx]]$yvar
                        ),
                FUN=transform.sym)

        plot.list[[plot.list.idx]] <- ggplot(data=data.tibble
                ,mapping=aes(!!!plot_args)) + geom_line()
    }

    # add ylim
    if ("ylim" %in% names(plot.structure[[plot_struct_idx]]))
    {
        plot.list[[plot.list.idx]] <- plot.list[[plot.list.idx]] + ylim(
                unlist(
                        plot.structure[[plot.list.idx]]$ylim)
                )
    }
    
    plot.list[[plot.list.idx]] <- plot.list[[plot.list.idx]] + theme_classic()
    
    plot.list.idx <- plot.list.idx + 1
}

sub.data <- pivot_longer(data=data.tibble
        ,cols=all_of(c(names_p,names_x))
        ,names_to=c(".value","trait_idx")
        ,names_pattern = "mean([xp])(\\d+)"
        )

phase_plot <- ggplot(
        data = sub.data,
        mapping = aes(x = p, y = x)) +
    geom_line(mapping = aes(colour = trait_idx))


# now a phase plot to end it all with
plot.list[[plot.list.idx]] <- phase_plot

title <- ""

if (exists("params") && "biasv" %in% names(params))
{
    title <- paste0(
            "bias: ",params["biasv"],
            " a: ",params["a"],
            " b: ",params["b"],
            " c: ",params["c"],
            " theta: ",params["theta"],
            " gamma: ",params["gamma"],
            " weber_k: ",params["weber_k"]
            )
}

wrap_plots(plot.list,ncol=1) + plot_annotation(
        title=title)

file.name.graph <- paste0("graph_",basename(file.name),".pdf")

ggsave(file.name.graph,height= 3 * plot.structure.l)






