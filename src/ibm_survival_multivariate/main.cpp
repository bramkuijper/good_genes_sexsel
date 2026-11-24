#include <string>
#include "good_genes.hpp"
#include "parameters.hpp"

int main(int argc, char **argv)
{
    Parameters parameters;

    parameters.biasv = std::stod(argv[1]);
    parameters.a = std::stod(argv[2]);
    parameters.b = std::stod(argv[3]);
    parameters.c = std::stod(argv[4]);
    parameters.theta = std::stod(argv[5]);
    parameters.gamma = std::stod(argv[6]);
    parameters.k = std::stod(argv[7]);
    parameters.weber_k = std::stod(argv[8]);
    parameters.init_t = std::stod(argv[9]);
    parameters.init_p = std::stod(argv[10]);
    parameters.init_v = std::stod(argv[11]);
    parameters.max_mut_t = std::stod(argv[12]);
    parameters.max_mut_p = std::stod(argv[13]);
    parameters.max_mut_v = std::stod(argv[14]);
    parameters.max_num_gen = static_cast<unsigned>(std::stoi(argv[15]));
    parameters.only_positive = std::stoi(argv[16]);
    parameters.ntraits = static_cast<unsigned>(std::stoi(argv[17]));
    parameters.file_name = argv[18];

    GoodGenes gg(parameters);

} // end main
