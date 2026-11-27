#include <cmath>
#include <random>
#include "individual.hpp"
#include "parameters.hpp"

// standard constructor
Individual::Individual(Parameters const &params) :
    t(2, std::vector<double> (params.ntraits, 0.0)), // initialize t at 0
    p(2, std::vector<double> (params.ntraits, 0.0)), // initialize p at 0
    v{params.init_v, params.init_v},
    x(params.ntraits, 0.0) // initialize x at 0
{
    // initialize vectors for t and p
    for (unsigned trait_idx{0}; trait_idx < 2; ++trait_idx)
    {
        for (unsigned allele_idx{0}; allele_idx < 2; ++allele_idx)
        {
            p[allele_idx][trait_idx] = params.init_p[trait_idx];
            t[allele_idx][trait_idx] = params.init_t[trait_idx];
        }
    }
} // end Individual::Individual

Individual::Individual(Individual const &other) :
    t{other.t[0], other.t[1]},
    p{other.p[0], other.p[1]},
    v{other.v[0], other.v[1]},
    x{other.x}
{}



Individual::Individual(
                Individual const &mother,
                Individual const &father,
                std::mt19937 &rng_r,
                Parameters const &params) :
    t(2, std::vector<double> (params.ntraits, 0.0)),
    p(2, std::vector<double> (params.ntraits, 0.0)),
    x(params.ntraits, 0.0)

{
    std::bernoulli_distribution segregator{0.5};
    std::uniform_real_distribution<double> uniform{0.0,1.0};

    v[0] = mother.v[segregator(rng_r)];
    v[1] = father.v[segregator(rng_r)];

    // aux variable for multational increments
    double incr;

    for (unsigned trait_idx{0}; trait_idx < params.ntraits; ++trait_idx)
    {
        p[0][trait_idx] = mother.p[segregator(rng_r)][trait_idx];
        p[1][trait_idx] = father.p[segregator(rng_r)][trait_idx];
        
        t[0][trait_idx] = mother.t[segregator(rng_r)][trait_idx];
        t[1][trait_idx] = father.t[segregator(rng_r)][trait_idx];

        for (unsigned allele_idx{0}; allele_idx < 2; ++allele_idx)
        {
            if (uniform(rng_r) < params.mu_p[trait_idx])
            {
                incr = uniform(rng_r) * params.max_mut_p;

                p[allele_idx][trait_idx] = p[allele_idx][trait_idx] + 
                    (uniform(rng_r) < 0.5 ? -incr : incr);

                if (params.only_positive && p[allele_idx][trait_idx] < 0.0)
                {
                    p[allele_idx][trait_idx] = 0.0;
                }
            }

            if (uniform(rng_r) < params.mu_t[trait_idx])
            {
                incr = uniform(rng_r) * params.max_mut_t;
                
                t[allele_idx][trait_idx] = t[allele_idx][trait_idx] +
                     (uniform(rng_r) < 0.5 ? -incr : incr);
                
                if (params.only_positive && t[allele_idx][trait_idx] < 0.0)
                {
                    t[allele_idx][trait_idx] = 0.0;
                }
            } //end if mu_t
        } // end for allele_idx
    } // end for trait_idx 
        
    for (unsigned allele_idx{0}; allele_idx < 2; ++allele_idx)
    {
        if (uniform(rng_r) < params.mu_v)
        {
            incr = uniform(rng_r) * params.max_mut_v;
            v[allele_idx] = v[allele_idx] + (uniform(rng_r) < params.biasv ? -incr : incr);
        }
    }
} // end Individual()

void Individual::operator=(Individual const &other)
{
    for (unsigned allele_idx{0}; allele_idx < 2; ++allele_idx)
    {
        p[allele_idx] = other.p[allele_idx];
        t[allele_idx] = other.t[allele_idx];
        v[allele_idx] = other.v[allele_idx];
    }

    x = other.x;
}//end operator=
