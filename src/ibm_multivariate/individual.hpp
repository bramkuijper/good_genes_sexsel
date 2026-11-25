#ifndef _INDIVIDUAL_HPP_
#define _INDIVIDUAL_HPP_

#include <random>
#include "parameters.hpp"

class Individual
{
    public:
        std::vector < std::vector < double > > t{};             
        std::vector < std::vector < double > > p{};             
        double v[2]{0.0,0.0};

        // realized ornament
        std::vector < double > x{}; 

        // standard constructor
        Individual(Parameters const &param);

        // copy ctor
        Individual(Individual const &other);

        // birth constructor
        Individual(Individual const &mother,
                Individual const&father,
                std::mt19937 &rng_r,
                Parameters const &params);

        void operator=(Individual const &other);
};

#endif 
