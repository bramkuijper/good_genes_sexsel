#include <cassert>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <iostream>
#include "good_genes.hpp"
#include "parameters.hpp"

GoodGenes::GoodGenes(Parameters const &params) :
    par{params},
    data_file{par.file_name},
    uniform{0.0,1.0},
    rd{},
    seed{rd()},
    rng_r{seed},
    males(par.n/2,Individual(par)),
    females(par.n/2,Individual(par))
{
    write_data_headers();

    for (time_step = 0; 
            time_step <= par.max_num_gen; ++time_step)
    {
        phenotypes();
        reproduction();

        if (time_step % par.numoutgen == 0)
        {
            std::cout << time_step << std::endl;
            write_data();
        }
    }

    write_parameters();

} // end GoodGenes() constructor

void GoodGenes::reproduction()
{
    update_survival_distribution();

    long unsigned female_idx,male_idx;

    std::vector<Individual> nextgen{};

    for (unsigned newborn_idx{0}; newborn_idx < par.n; ++newborn_idx)
    {
        female_idx = survival_distribution(rng_r);

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
        } else
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
        << "n;" << par.n << ";" << std::endl
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
        << "choice_sample_size;" << par.choice_sample_size << ";" << std::endl
        << "init_t;" << par.init_t << ";" << std::endl
        << "init_p;" << par.init_p << ";" << std::endl
        << "init_v;" << par.init_v << ";" << std::endl
        << "v_opt;" << par.v_opt << ";" << std::endl
        << "max_num_gen;" << par.max_num_gen << ";" << std::endl
        << "numoutgen;" << par.numoutgen << ";" << std::endl;
}

void GoodGenes::write_data()
{
    double meanp{0.0};
    double ssp{0.0};
    double meant{0.0};
    double sst{0.0};
    double meanv{0.0};
    double ssv{0.0};
    double meanx{0.0};
    double ssx{0.0};


    // keep track of population sizes
    unsigned long nf{females.size()};
    unsigned long nm{males.size()};

    // aux variables to store trait values
    double p,t,v,x;

    for (auto female_iter{females.begin()};
            female_iter != females.end();
            ++female_iter)
    {
        p = 0.5 * (female_iter->p[0] + female_iter->p[1]);
        meanp += p;
        ssp += p*p;
        
        t = 0.5 * (female_iter->t[0] + female_iter->t[1]);
        meant += t;
        sst += t*t;

        v = 0.5 * (female_iter->v[0] + female_iter->v[1]);
        meanv += v;
        ssv += v*v;
    }
    
    for (auto male_iter{males.begin()};
            male_iter != males.end();
            ++male_iter)
    {
        p = 0.5 * (male_iter->p[0] + male_iter->p[1]);
        meanp += p;
        ssp += p*p;
        
        t = 0.5 * (male_iter->t[0] + male_iter->t[1]);
        meant += t;
        sst += t*t;

        v = 0.5 * (male_iter->v[0] + male_iter->v[1]);
        meanv += v;
        ssv += v*v;
        
        x = male_iter->x;
        meanx += x;
        ssx += x*x;
    }


    meanp /= (nf + nm);
    meant /= (nf + nm);
    meanv /= (nf + nm);
    meanx /= nm;
    
    double varp = ssp / (nf + nm) - meanp * meanp;
    double vart = sst / (nf + nm) - meant * meant;
    double varv = ssv / (nf + nm) - meanv * meanv;
    double varx = ssx / nm - meanx * meanx;

    data_file << time_step << ";"
        << meanp << ";"
        << meant << ";"
        << meanv << ";"
        << meanx << ";"
        << varp << ";"
        << vart << ";"
        << varv << ";"
        << varx << ";"
        << nf << ";"
        << nm << ";" << std::endl;
} // write_data()

void GoodGenes::write_data_headers()
{
    data_file << "generation;meanp;meant;meanv;meanx;varp;vart;varv;varx;nf;nm;" << std::endl;
}
void GoodGenes::phenotypes()
{
    double t,v;
    for (auto male_iter{males.begin()};
            male_iter != males.end();
            ++male_iter)
    {
        t = 0.5 * (male_iter->t[0] + male_iter->t[1]);
        v = 0.5 * (male_iter->v[0] + male_iter->v[1]);

        male_iter->x = t * std::exp(-std::fabs(par.v_opt - v));
    }
} // end phenotypes()

// choose surviving male according to its ornament 
unsigned GoodGenes::choose(Individual const &female)
{
    std::vector <double> male_fitness{};
    std::vector <unsigned> male_idxs{};

    // distribution to sample males from
    std::uniform_int_distribution<unsigned> 
        male_sampler(0, males.size() - 1);

    unsigned sampled_male_idx;
    
    double fitness;

    double p = 0.5 * (female.p[0] + female.p[1]);

    double x,v,survival_odds;

    for (unsigned inspected_male_idx{0}; 
            inspected_male_idx < par.choice_sample_size; 
            ++inspected_male_idx)
    {
        sampled_male_idx = male_sampler(rng_r);

        x = males[sampled_male_idx].x;
        v = 0.5 * (males[sampled_male_idx].v[0] + males[sampled_male_idx].v[1]);

        survival_odds = std::exp(-par.c * x * x - std::fabs(par.v_opt - v));

        fitness = survival_odds * std::exp(par.a * p * x);
        
        male_idxs.push_back(sampled_male_idx);
        male_fitness.push_back(fitness);

    }

    // now make distribution of the fitnesses to choose from
    std::discrete_distribution<unsigned> choose_male(male_fitness.begin(),
            male_fitness.end());

    unsigned the_chosen_male = male_idxs[choose_male(rng_r)];

    return(the_chosen_male);
} // choose()



void GoodGenes::update_survival_distribution()
{
    std::vector<double> female_survival{}; 

    double surv;

    double p, v;

    for (auto female_iter{females.begin()}; 
            female_iter != females.end();
            ++female_iter
            )
    {
        p = 0.5 * (female_iter->p[0] + female_iter->p[1]);
        v = 0.5 * (female_iter->v[0] + female_iter->v[1]);

        surv = std::exp(-par.b * p * p - std::fabs(par.v_opt - v)) ;

        female_survival.push_back(surv);
    }

    assert(female_survival.size() == females.size());

    // param object to override the current fecundity distribution
    // that is associated to the patch
    std::discrete_distribution<unsigned>::param_type 
        survival_distribution_param(
                female_survival.begin()
                ,female_survival.end());

    // update the survival distribution
    survival_distribution.param(
            survival_distribution_param);
}//end update_survival_distributions
