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
    for (unsigned time_step{0}; 
            time_step < par.max_num_gen; ++time_step)
    {
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

        male_idx = choose(female_idx);
        
        assert(male_idx >= 0);
        assert(male_idx < males.size());

        Individual Kid(
                females[female_idx],
                males[male_idx],
                rng_r,
                params);

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


void GoodGenes::choose()
{
    std::vector male_fitness;

    unsigned male_idx;

    for (unsigned sample_idx{0}; 
            sample_idx < par.choice_sample_size; 
            ++sample_idx)
    {
        sample_idx = 
    }
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
