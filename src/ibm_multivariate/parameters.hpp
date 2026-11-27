#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_

#include <string>

class Parameters
{
    public:
        unsigned n{5000}; // population size

        // per-generation, per allele mutation probability of preference alleles
        std::vector <double> mu_p;
        // per-generation, per allele mutation probability of ornament alleles
        std::vector <double> mu_t;
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

        double weber_k{0.0};
        
        // cost of preferences as female survival is given by exp(-b * p * p - |v - v_opt|)
        double b{0.0025};

        // how costs of multiple preferences should be scaled
        // if theta = 1/gamma then costs of each preference accrue independently
        //
        // if theta -> 1, costs combine in a super additive fashion and females that
        // evaluate males based on multiple preferences face severe fitness costs
        //
        // if theta -> 0, costs are no more than the max cost 
        // if either preference is considered
        double theta{0.5};

        // how quickly survival decays away from the optimum
        double gamma{2};

        // parameter scaling condition-dependence of costs of ornamentation
        // if k = 0, then costs are only dependent on quality in terms of ornament
        // size (i.e., high quality males have higher costs chiefly because they 
        // have larger ornaments). 
        // if k > 0, then higher quality males incur lower costs of ornamentation
        // (i.e., smaller c) than low quality males
        // if k < 0, then higher quality males incur higher costs of ornamentation
        // (i.e., higher c) than low quality males
        double k{1.0};

        // cost of ornaments, as survival is given by exp(-c * x * x - |v - v_opt|)
        double c{0.5};

        // number of males assessed by female
        unsigned choice_sample_size{10};

        // initial values of p, t, v
        std::vector <double> init_t;
        std::vector <double> init_p;
        double init_v{5.0};

        // value of v that maximizes survival
        double v_opt{10.0};

        // number of generations until simulation ends
        unsigned max_num_gen{1000};
        // data output to file every numoutgen generations
        unsigned numoutgen{10};

        // number of traits
        unsigned ntraits{5};

        // file name
        std::string file_name{"sim_good_genes"};

        // whether values of t, p can only be positive (
        // or strictly speaking, nonnegative)
        bool only_positive{false};
};
#endif
