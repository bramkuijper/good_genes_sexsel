# Good-genes sexual selection: individual-based simulations

## Versions

See the different folders in the `src` directory for the different versions of the simulation code.

 * `ibm`: the individual-based simulations as used in Fawcett et al 2011 PNAS https://doi.org/10.1073/pnas.1105721108 (without sex allocation)

 * `ibm_survival`: same as above, but in which we incorporate an explicit survival stage rather than having survival contribute to relative fitness, in line with older C++ code of Fisher, [here](https://github.com/bramkuijper/fisher_sexsel/blob/fe805f7c89d8ea99eb77adac71dbbc8ac291f939/src/ibm/fisherian.cpp)
 Adding such realism increases the strength of survival selection on ornaments and preferences, which reduces absolute exaggeration of `p` and `t`
