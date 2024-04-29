#include <cassert>
#include <stdexcept>
#include <vector>
#include <cmath>
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
    for (time_step = 0; 
            time_step < par.max_num_gen; ++time_step)
    {
        phenotypes();
        reproduction();
    }
} // end GoodGenes() constructor

void GoodGenes::reproduction()
{
    update_survival_distribution();

    int female_idx,male_idx;

    std::vector<Individual> nextgen;

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

    data_file << time_step << ";";

    unsigned nf{females.size()};
    unsigned nm{males.size()};

    for (auto female_iter{females.begin()};
            female_iter != females.end();
            ++female_iter)
    {
        p = female_iter->p[0] + female_iter->p[1];
        meanp += p;
        ssp += p*p;
        
        t = female_iter->t[0] + female_iter->t[1];
        meant += t;
        sst += t*t;

        v = female_iter->v[0] + female_iter->v[1];
        meanv += v;
        ssv += v*v;
    }
    
    for (auto male_iter{males.begin()};
            male_iter != males.end();
            ++male_iter)
    {
        p = male_iter->p[0] + male_iter->p[1];
        meanp += p;
        ssp += p*p;
        
        t = male_iter->t[0] + male_iter->t[1];
        meant += t;
        sst += t*t;

        v = male_iter->v[0] + male_iter->v[1];
        meanv += v;
        ssv += v*v;
        
        x = male_iter->x[0] + male_iter->x[1];
        meanx += x;
        ssx += x*x;
    }


    meanp /= (nf + nm)
    meant /= (nf + nm)
    meanv /= (nf + nm)
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
        << nm << ";";
} // write_data()

void write_data_headers()
{

}
void GoodGenes::phenotypes()
{
    double t,v;
    for (auto male_iter{males.begin()};
            male_iter != males.end();
            ++male_iter)
    {
        t = male_iter->t[0] + male_iter->t[1];
        v = male_iter->v[0] + male_iter->v[1];

        male_iter->x = t * std::exp(-abs(par.v_opt - v));
    }
} // end phenotypes()

// choose surviving male according to its ornament 
unsigned GoodGenes::choose(Individual const &female)
{
    std::vector <double> male_fitness;
    std::vector <unsigned> male_idxs;

    // distribution to sample males from
    std::uniform_int_distribution<unsigned> 
        male_sampler(0, males.size() - 1);

    unsigned sample_idx;
    
    double fitness;

    double p = female.p[0] + female.p[1];

    double x,v,survival_odds;
    for (unsigned sample_idx{0}; 
            sample_idx < par.choice_sample_size; 
            ++sample_idx)
    {
        sample_idx = male_sampler(rng_r);

        x = males[sample_idx].x;
        v = males[sample_idx].v[0] + males[sample_idx].v[1];

        survival_odds = std::exp(-par.c * x * x - abs(par.v_opt - v));

        fitness = survival_odds * std::exp(par.a * p * x);
        
        male_idxs.push_back(sample_idx);
        male_fitness.push_back(sample_idx);
    }

    // now make distribution of the fitnesses to choose from
    std::discrete_distribution<unsigned> choose_male(male_fitness.begin(),
            male_fitness.end());

    return(choose_male(rng_r));
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
        p = female_iter->p[0] + female_iter->p[1];
        v = female_iter->v[0] + female_iter->v[1];

        surv = std::exp(-par.b * p*p - abs(par.v_opt - v));

        female_survival.push_back(surv);
    }

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
