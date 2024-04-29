#ifndef _GOOD_GENES_HPP_
#define _GOOD_GENES_HPP_

#include <vector>
#include <random>

#include <fstream>
#include "individual.hpp"

class GoodGenes
{
    private:
        Parameters par;

        std::ofstream data_file;
        std::uniform_real_distribution<double> uniform;

        std::random_device rd;
        unsigned seed;
        std::mt19937 rng_r;

        unsigned time_step{0};


        std::discrete_distribution<unsigned> survival_distribution{1,1};

        std::vector<Individual> males;
        std::vector<Individual> females;

        void update_survival_distribution();
        void reproduction();
        void phenotypes();
        unsigned choose(Individual const &female);

        void write_data();

    public:
        GoodGenes(Parameters const &parameters);



};


#endif
