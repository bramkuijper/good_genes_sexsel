# Individual-based simulations

The heart of the code is available in the file `good_genes.cpp`. In the constructor function `GoodGenes::GoodGenes()` (constructor means: a function that sets up the whole simulation), you will see a life cycle (reflected by a `for` loop). Each time step, the simulation expresses phenotypes, using the function `phenotypes()` and then we have mate choice and reproduction in the function `reproduction()`.

See the top-level `README` file about how to compile the function.


