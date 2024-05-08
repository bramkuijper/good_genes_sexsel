# Good-genes sexual selection: individual-based simulations

## Versions

See the different folders in the `src` directory for the different versions of the simulation code.

 * `ibm`: the individual-based simulations as used in Fawcett et al 2011 PNAS https://doi.org/10.1073/pnas.1105721108 (without sex allocation)

 * `ibm_survival`: same as above, but in which we incorporate an explicit survival stage rather than having survival contribute to relative fitness, in line with older C++ code of Fisher, [here](https://github.com/bramkuijper/fisher_sexsel/blob/fe805f7c89d8ea99eb77adac71dbbc8ac291f939/src/ibm/fisherian.cpp)
 Adding such realism increases the strength of survival selection on ornaments and preferences, which reduces absolute exaggeration of `p` and `t`


## Downloading the code
On mac or Linux with `git`, `make` and `g++` or `clang` installed, download the code:

```
git clone https://github.com/bramkuijper/good_genes_sexsel.git 
```

## Running the simulation

### Compiling the code
First, compile the `C++` code:

```
cd good_genes_sexsel/src/ibm
make
```

There should now be an executable named `good_genes.exe` in this folder. However, we cannot just run this as a standalone programme, as running the executable requires a bunch of command-line parameters.


### Generate a shell file to run the simulation
The programme `good_genes.exe` requires a few command line parameters (you can add more by modifying the code in `main.cpp` and subsequently recompiling). 

To generate all permutations of these command line parameters, I have created a standalone `R` script `vary_parameters.r`. This script generates all permutations of parameters and then creates a `bash` shell command to run simulations for each permutation. The script `vary_parameters.r` writes these commands to the bash file `batch_file.sh`:

```
./vary_parameters.r > batch_file.sh
```
It might be insightful to inspect the contents of `vary_parameters.r` and `batch_file.sh` in an editor, so that you know what is going on. You can edit `vary_parameters.r` to change the values of the parameter combinations.

If you look at `batch_file.sh`, you will see that every second line contains a command to run a single simulation (for a single permutation of the parameters). The other line contains a print statement with the number of the current simulation.

### Actually running the simulation
Once you have a `batch_file.sh` with instructions to run the code, you can simply type the following on the command line:

```
bash batch_file.sh
```

after which each command in the file `batch_file.sh` will be run sequentially.

### Modifying the code to add different parameters etc
If you would like to add more permutations of different parameters that are present in `parameters.cpp`, you will need to change both `vary_parameters.r` (to add your parameter of choice to the permuation generating code, which is a simple `for` loop) as well as `main.cpp` to have the simulation overwrite the defaults via the contents of the corresponding `argv` array.

Don't forget to re-run `make` once you have accommodated for additional command line parameters in `main.cpp`.



