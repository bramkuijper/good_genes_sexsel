#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_

#include <string>

class Parameters
{
    public:
        unsigned n{5000}; // population size

        // per-generation, per allele mutation probability of preference alleles
        double mu_p{0.05}; 
        // per-generation, per allele mutation probability of ornament alleles
        double mu_t{0.05};
        // per-generation, per allele mutation probability of viability alleles
        double mu_v{0.05};

        // mutations in p, t and v are drawn from a uniform distribution
        // below are the maximal mutational values p, t and v respectively
        double max_mut_p{0.4};
        double max_mut_t{0.4};
        double max_mut_v{0.8};

        // biased mutation on v; if 0 all mutations are increasing the value of
        // of 0.5, there is no bias and mutations are equally likely to increase vs
        // decrease v; if 1.0 then all mutations decrease the value of v
        double biasv{0.99};

        // efficacy of sexual selection, as odds of mating are given by exp(a * p * x)
        double a{1.0};
        
        // cost of preferences as female survival is given by exp(-b * p * p - |v - v_opt|)
        double b{0.0025};

        // cost of ornaments, as survival is given by exp(-c * x * x - |v - v_opt|)
        double c{0.5};

        // number of males assessed by female
        unsigned choice_sample_size{10};

        // initial values of p, t, v
        double init_t{0.0};
        double init_p{1.0};
        double init_v{5.0};

        // value of v that maximizes survival
        double v_opt{10.0};

        // number of generations until simulation ends
        unsigned max_num_gen{1000};
        // data output to file every numoutgen generations
        unsigned numoutgen{10};

        // file name
        std::string file_name{"sim_good_genes"};
};
#endif
