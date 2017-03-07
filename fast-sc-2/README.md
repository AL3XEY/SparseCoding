#README

##Disclaimer

The following code is a modification of [fast-sc-2](http://ai.stanford.edu/~hllee/softwares/nips06-sparsecoding.htm), a code provided by Honglak Lee, Alexis Battle, Rajat Raina and Andrew Y. Ng from Stanford University. The file `README` is their original instruction file.

Files added :

* `reconstruction.m`
* `is_octave.m`
* `main.m`
* `getdata_imagearray_some.m`
* `getdata_imagearray_all.m`

Files updates :

* `l2ls_learn_basis_dual.m`
* `l1ls_featuresign.m`
* `display_network_nonsquare2.m`
* `demo_fast_sc.m`
* `sc2/getObjective2.m`

##Usage

To run the learning algorithm, use the function `sparse_coding(...)`. There are some parameters to this function, so you can take a look at the example `demo_fast_sc(...)` function to have an idea about the values to give the parameters. You can also run the script file `main.m`, that will simply call the `demo_fast_sc(...)` function.

By default, the algorithm will run on the dataset in the file `data/IMAGES_RAW.mat`.

To recover an image, use the `reconstruction(...)` function. Results of the dictionnary learning and the reconstruction will be saved under the folder `results/`.

If you want to test the reconstruction without learning a new dictionnary, you can use the ones under the folders `results/<number>/`, where a README file will explain how the dictionnary was learn for each.

This software is still work in progress but should work under Matlab R2013a and Octave 4.2.0 and 4.0.2
