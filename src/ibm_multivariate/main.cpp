#include <string>
#include "good_genes.hpp"
#include "parameters.hpp"

int main(int argc, char **argv)
{
    Parameters parameters;

    unsigned arg_ctr{1};

    parameters.biasv = std::stod(argv[arg_ctr++]);
    parameters.a = std::stod(argv[arg_ctr++]);
    parameters.b = std::stod(argv[arg_ctr++]);
    parameters.c = std::stod(argv[arg_ctr++]);
    parameters.theta = std::stod(argv[arg_ctr++]);
    parameters.gamma = std::stod(argv[arg_ctr++]);
    parameters.k = std::stod(argv[arg_ctr++]);
    parameters.weber_k = std::stod(argv[arg_ctr++]);
    parameters.init_v = std::stod(argv[arg_ctr++]);
    parameters.max_mut_t = std::stod(argv[arg_ctr++]);
    parameters.max_mut_p = std::stod(argv[arg_ctr++]);
    parameters.max_mut_v = std::stod(argv[arg_ctr++]);
    parameters.max_num_gen = static_cast<unsigned>(std::stoi(argv[arg_ctr++]));
    parameters.only_positive = std::stoi(argv[arg_ctr++]);
    parameters.ntraits = static_cast<unsigned>(std::stoi(argv[arg_ctr++]));
    
    // initialize vectors with mutation rates
    parameters.mu_p = std::vector<double> (parameters.ntraits, 0.0);
    parameters.mu_t = std::vector<double> (parameters.ntraits, 0.0);
    parameters.init_p = std::vector<double> (parameters.ntraits, 0.0);
    parameters.init_t = std::vector<double> (parameters.ntraits, 0.0);

    
    for (unsigned trait_idx{0}; 
            trait_idx < parameters.ntraits;
            ++trait_idx)
    {
        parameters.mu_p[trait_idx] = std::stod(argv[arg_ctr++]);
    }
        
    for (unsigned trait_idx{0}; 
            trait_idx < parameters.ntraits;
            ++trait_idx)
    {
        parameters.mu_t[trait_idx] = std::stod(argv[arg_ctr++]);
    }

    for (unsigned trait_idx{0}; 
            trait_idx < parameters.ntraits;
            ++trait_idx)
    {
        parameters.init_p[trait_idx] = std::stod(argv[arg_ctr++]);
    }

    for (unsigned trait_idx{0}; 
            trait_idx < parameters.ntraits;
            ++trait_idx)
    {
        parameters.init_t[trait_idx] = std::stod(argv[arg_ctr++]);
    }


    parameters.file_name = argv[arg_ctr++];

    GoodGenes gg(parameters);

} // end main
