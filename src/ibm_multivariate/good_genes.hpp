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
        std::normal_distribution <double> normal{0.0, 1.0};

        std::discrete_distribution<unsigned> female_survival_distribution{};
        std::discrete_distribution<unsigned> male_survival_distribution{};

        std::random_device rd;
        unsigned seed;
        std::mt19937 rng_r;

        unsigned time_step{0};

        double mean_p_survive_f{0.0};
        double mean_p_survive_m{0.0};

        std::vector<Individual> males;
        std::vector<Individual> females;

        unsigned nfsurv{0};
        unsigned nmsurv{0};

        void reproduction();
        void update_survival_distribution();
        void phenotypes();
        unsigned choose(Individual const &female);

        void write_data();
        void write_data_headers();
        void write_parameters();

    public:
        GoodGenes(Parameters const &parameters);



};


#endif
