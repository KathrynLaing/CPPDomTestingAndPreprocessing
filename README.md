# CPPDomTestingAndPreprocessing

We start by randomly generating a set of CP-nets and outcomes in the exact same way as in the DQ-Pruning repository. These are the dominance queries we wish to test. They are stored as .RData objects, as before.

Then, we use the CPNConvert function to take these RObject CP-nets and convert them into a form that we can pass to our Rcpp functions (within R)
These functions return a list in R, which we can then manipulate as we like including storing the results.
The necessary set up to call the Rcpp functions in R is:...

