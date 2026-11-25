#include <cassert>
#include <random>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <iostream>
#include "good_genes.hpp"
#include "parameters.hpp"

// the constructor, which initializes 
// and then runs the whole simulation
GoodGenes::GoodGenes(Parameters const &params) :
    par{params}, // initialize the parameter object
    data_file{par.file_name}, // start the output file
    uniform{0.0,1.0}, // initialize a standard uniform distribution
    rd{}, // rd() generates the seed fed to the random number generator
    seed{rd()}, // save the seed
    rng_r{seed}, // initialize the random number generator
    males(par.n/2,Individual(par)), // initialize males
    females(par.n/2,Individual(par)) // initialize females
{

    // finally write parameters to output file
    write_parameters();

    // add headers to the data file
    write_data_headers();

    for (time_step = 0; 
            time_step <= par.max_num_gen; ++time_step)
    {

        // set phenotypes 
        phenotypes();

        // ... then reproduce (which 
        // includes generating a survival dist)
        reproduction();

        // write data to file ever nth timestep
        if (time_step % par.numoutgen == 0)
        {
            write_data();
        }
    }

} // end GoodGenes() constructor

// express phenotypes 
void GoodGenes::phenotypes()
{
    double ti,v;

    for (auto male_iter{males.begin()}; male_iter != males.end(); ++male_iter)
    {
        v = 0.5 * (male_iter->v[0] + male_iter->v[1]);

        for (unsigned trait_idx{0}; trait_idx < par.ntraits; ++trait_idx)
        {
            ti = 0.5 * (male_iter->t[0][trait_idx] + male_iter->t[1][trait_idx]);

            male_iter->x[trait_idx] = ti * std::exp(-std::fabs(par.v_opt - v));
        }
    }
}

void GoodGenes::update_survival_distribution()
{
    // aux variables to store trait values
    double p,v,surv;

    std::vector<double> female_survival{};
    std::vector<double> male_survival{};

    // initialize survival probabilities
    mean_p_survive_f = 0.0;
    mean_p_survive_m = 0.0;

    double sum_preference_costs, sum_trait_costs;

    double theta_inv = 1.0 / par.theta;
    
    for (auto female_iter{females.begin()}; 
            female_iter != females.end(); 
            ++female_iter)
    {
        sum_preference_costs = 0.0;

        for (unsigned trait_idx{0}; 
                trait_idx < par.ntraits;
                ++trait_idx)
        {
            p = std::fabs(0.5 * (female_iter->p[0][trait_idx] + 
                    female_iter->p[1][trait_idx]));

            sum_preference_costs += std::pow(p, theta_inv);
        }

        v = 0.5 * (female_iter->v[0] + female_iter->v[1]);

        surv = std::exp(-par.b * 
                std::pow(sum_preference_costs, par.gamma * par.theta)
                            - std::fabs(par.v_opt - v));

        // update average survival prob
        mean_p_survive_f += surv;

        female_survival.push_back(surv);

    } // end for females

    // param object to update the current fecundity distribution
    std::discrete_distribution<unsigned>::param_type
        female_survival_distribution_param(
                female_survival.begin(),
                female_survival.end());

    // update the female survival distribution
    female_survival_distribution.param(
        female_survival_distribution_param
            );

    for (auto male_iter{males.begin()}; 
            male_iter != males.end(); ++male_iter)
    {
        sum_trait_costs = 0.0;
        
        v = 0.5 * (male_iter->v[0] + male_iter->v[1]);

        for (unsigned trait_idx{0}; trait_idx < par.ntraits;++trait_idx)
        {
            sum_trait_costs += par.c / (1.0 + par.k * v) * 
                male_iter->x[trait_idx] * male_iter->x[trait_idx];
        }
        
        surv = std::exp(-sum_trait_costs - std::fabs(par.v_opt - v)) ;
        
        mean_p_survive_m += surv;

        male_survival.push_back(surv);
    }

    // param object to update the current fecundity distribution
    std::discrete_distribution<unsigned>::param_type
        male_survival_distribution_param(
                male_survival.begin(),
                male_survival.end());

    // update the female survival distribution
    male_survival_distribution.param(
        male_survival_distribution_param
            );

    mean_p_survive_f /= static_cast<unsigned>(females.size());
    mean_p_survive_m /= static_cast<unsigned>(males.size());
} // end survival


void GoodGenes::reproduction()
{
    long unsigned female_idx,male_idx;

    std::vector<Individual> nextgen{};

    // extinction
    if (females.size() < 1 || males.size() < 1)
    {
        write_parameters();
        exit(1);
    }

    for (unsigned newborn_idx{0}; newborn_idx < par.n; ++newborn_idx)
    {
        female_idx = female_survival_distribution(rng_r);

        assert(female_idx >= 0);
        assert(female_idx < females.size());

        male_idx = choose(females[female_idx]);
        
        assert(male_idx >= 0);
        assert(male_idx < males.size());

        Individual Kid(
                females[female_idx],
                males[male_idx],
                rng_r,
                par);

        nextgen.push_back(Kid);
    }

    assert(nextgen.size() == par.n);

    females.clear();
    males.clear();

    for (unsigned newborn_idx{0}; newborn_idx < par.n; ++newborn_idx)
    {
        if (uniform(rng_r) < 0.5)
        {
            males.push_back(nextgen[newborn_idx]);
        } 
        else
        {
            females.push_back(nextgen[newborn_idx]);
        }
    }
} // end reproduction

void GoodGenes::write_parameters()
{
    data_file 
        << std::endl 
        << std::endl
        << "seed;" << seed << ";" << std::endl
        << "type;" << "good_genes_survival" << ";" << std::endl
        << "n;" << par.n << ";" << std::endl
        << "ntraits;" << par.ntraits << ";" << std::endl
        << "mu_p;" << par.mu_p << ";" << std::endl
        << "mu_t;" << par.mu_t << ";" << std::endl
        << "mu_v;" << par.mu_v << ";" << std::endl
        << "max_mut_p;" << par.max_mut_p << ";" << std::endl
        << "max_mut_t;" << par.max_mut_t << ";" << std::endl
        << "max_mut_v;" << par.max_mut_v << ";" << std::endl
        << "biasv;" << par.biasv << ";" << std::endl
        << "a;" << par.a << ";" << std::endl
        << "b;" << par.b << ";" << std::endl
        << "c;" << par.c << ";" << std::endl
        << "theta;" << par.theta << ";" << std::endl
        << "gamma;" << par.gamma << ";" << std::endl
        << "k;" << par.k << ";" << std::endl
        << "weber_k;" << par.weber_k << ";" << std::endl
        << "choice_sample_size;" << par.choice_sample_size << ";" << std::endl
        << "init_t;" << par.init_t << ";" << std::endl
        << "init_p;" << par.init_p << ";" << std::endl
        << "init_v;" << par.init_v << ";" << std::endl
        << "v_opt;" << par.v_opt << ";" << std::endl
        << "only_positive;" << par.only_positive << ";" << std::endl
        << "max_num_gen;" << par.max_num_gen << ";" << std::endl
        << "numoutgen;" << par.numoutgen << ";" << std::endl;
}

void GoodGenes::write_data()
{

    std::vector <double> meanp(par.ntraits, 0.0);
    std::vector <double> ssp(par.ntraits, 0.0);
    std::vector <double> meant(par.ntraits, 0.0);
    std::vector <double> sst(par.ntraits, 0.0);
    std::vector <double> meanx(par.ntraits, 0.0);
    std::vector <double> ssx(par.ntraits, 0.0);
    double meanv{0.0};
    double ssv{0.0};
    
    //  covariances
    std::vector <double> stv(par.ntraits, 0.0);
    std::vector <double> stp(par.ntraits, 0.0);
    std::vector <double> spv(par.ntraits, 0.0);

    // keep track of population sizes
    unsigned long nf{females.size()};
    unsigned long nm{males.size()};

    // aux variables to store trait values
    double p,t,v,x;

    double nf_plus_nm = static_cast<double>(nf + nm);
    
    data_file << time_step << ";";

    for (auto female_iter{females.begin()};
            female_iter != females.end();
            ++female_iter)
    {
        for (unsigned trait_idx{0}; 
                trait_idx < par.ntraits; ++trait_idx)
        {
            p = 0.5 * (female_iter->p[0][trait_idx] + female_iter->p[1][trait_idx]);
            meanp[trait_idx] += p;
            ssp[trait_idx] += p*p;
        
            t = 0.5 * (female_iter->t[0][trait_idx] + female_iter->t[1][trait_idx]);
            meant[trait_idx] += t;
            sst[trait_idx] += t*t;

            stp[trait_idx] += t * p;
            spv[trait_idx] += p * v;
            stv[trait_idx] += t * v;
        }

        v = 0.5 * (female_iter->v[0] + female_iter->v[1]);
        meanv += v;
        ssv += v*v;
    }
    
    for (auto male_iter{males.begin()};
            male_iter != males.end();
            ++male_iter)
    {
        for (unsigned trait_idx{0}; 
                trait_idx < par.ntraits; ++trait_idx)
        {
            p = 0.5 * (male_iter->p[0][trait_idx] + male_iter->p[1][trait_idx]);
            meanp[trait_idx] += p;
            ssp[trait_idx] += p*p;
        
            t = 0.5 * (male_iter->t[0][trait_idx] + male_iter->t[1][trait_idx]);
            meant[trait_idx] += t;
            sst[trait_idx] += t*t;
            
            x = male_iter->x[trait_idx];
            meanx[trait_idx] += x;
            ssx[trait_idx] += x*x;

            stp[trait_idx] += t * p;
            spv[trait_idx] += p * v;
            stv[trait_idx] += t * v;
        }

        v = 0.5 * (male_iter->v[0] + male_iter->v[1]);
        meanv += v;
        ssv += v*v;
    }

    for (unsigned trait_idx{0}; 
            trait_idx < par.ntraits; ++trait_idx)
    {
        meanp[trait_idx] /= nf_plus_nm;
        meant[trait_idx] /= nf_plus_nm;
        meanx[trait_idx] /= static_cast<double>(nm);
    
        double varp = ssp[trait_idx] / nf_plus_nm - 
            meanp[trait_idx] * meanp[trait_idx];

        double vart = sst[trait_idx] / nf_plus_nm - 
            meant[trait_idx] * meant[trait_idx];

        double varx = ssx[trait_idx] / static_cast<double>(nm) - 
            meanx[trait_idx] * meanx[trait_idx];
    
        double covtp = stp[trait_idx] / nf_plus_nm - 
            meant[trait_idx] * meanp[trait_idx];
        
        double covtv = stv[trait_idx] / nf_plus_nm - 
            meant[trait_idx] * meanv;

        double covpv = spv[trait_idx] / nf_plus_nm - 
            meanp[trait_idx] * meanv;

        data_file << meanp[trait_idx] << ";"
            << meant[trait_idx] << ";"
            << meanx[trait_idx] << ";"
            << varp << ";"
            << vart << ";"
            << varx << ";"
            << covtp << ";"
            << covtv << ";"
            << covpv << ";";
    }

    meanv /= nf_plus_nm;
    double varv = ssv / nf_plus_nm - meanv * meanv;
    

    data_file << meanv << ";"
        << varv << ";"
        << mean_p_survive_f << ";"
        << mean_p_survive_m << ";"
        << nfsurv << ";"
        << nmsurv << ";"
        << nf << ";"
        << nm << ";" << std::endl;

} // write_data()

void GoodGenes::write_data_headers()
{
    data_file << "generation;";

    for (unsigned trait_idx{0}; trait_idx < par.ntraits; ++trait_idx)
    {
        data_file << "meanp" << (trait_idx + 1) << ";"
            << "meant" << (trait_idx + 1) << ";"
            << "meanx" << (trait_idx + 1) << ";"
            << "varp" << (trait_idx + 1) << ";"
            << "vart" << (trait_idx + 1) << ";"
            << "varx" << (trait_idx + 1) << ";"
            << "covt" << (trait_idx + 1) << "p" << (trait_idx + 1) << ";"
            << "covt" << (trait_idx + 1) << "v;" 
            << "covp" << (trait_idx + 1) << "v;";
    }

    data_file << "meanv;varv;p_surv_f;p_surv_m;nfsurv;nmsurv;nf;nm;" << std::endl;
}

// choose surviving male according to its ornament 
unsigned GoodGenes::choose(Individual const &female)
{
    std::vector <double> male_fitness{};
    std::vector <unsigned> male_idxs{};

    unsigned sampled_male_idx;
    
    double sum_fitness;

    double x;

    std::vector p(par.ntraits, 0.0);

    for (unsigned trait_idx{0}; 
            trait_idx < par.ntraits; ++trait_idx)
    {
        p[trait_idx] = 0.5 *(female.p[0][trait_idx] + 
                female.p[1][trait_idx]);
    }

    for (unsigned inspected_male_idx{0}; 
            inspected_male_idx < par.choice_sample_size; 
            ++inspected_male_idx)
    {
        // sample male according to his survival value
        sampled_male_idx = 
            male_survival_distribution(rng_r);

        sum_fitness = 0.0;

        for (unsigned trait_idx{0}; trait_idx < par.ntraits; ++trait_idx)
        {
            x = males[sampled_male_idx].x[trait_idx];

            sum_fitness += par.a * p[trait_idx] * (x + 
                    x * par.weber_k * normal(rng_r));

        }

        sum_fitness = std::exp(sum_fitness);

        assert(std::isnormal(sum_fitness));

        male_idxs.push_back(sampled_male_idx);
        male_fitness.push_back(sum_fitness);
    }

    // now make distribution of the fitnesses to choose from
    std::discrete_distribution<unsigned> choose_male(male_fitness.begin(),
            male_fitness.end());

    unsigned the_chosen_male = male_idxs[choose_male(rng_r)];

    return(the_chosen_male);
} // choose()



