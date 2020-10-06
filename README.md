# C++ Dominance Testing and Preprocessing Functions

This repository gives the code for the Chapter 2 and 3 experiments from my thesis:

*Laing, K. (2020). Conditional Preference Networks: Efficient Dominance Testing and Learning. PhD Thesis. School of Mathematics, University of Leeds, UK.*

In particular, C++ versions of the dominance testing functions (in R) given in my DQ-Pruning repository as well as CP-net preprocessing functions.

All of these functions are given in the file `DQFunctions.cpp`. 

These functions use the Rcpp package. Rcpp is a tool for integrating C++ and R, which means that these functions can be called directly from R if one uses the following preamble in the R script:
```
library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-std=c++0x")
sourceCpp("DQFunctions.cpp")
```
As we are now using C++ rather than R, we had to develop a new encoding structure for CP-nets. We discuss this and how to move between this and the R encoding from the DQ-Pruning repository in the following section. We then discuss the updated (C++) dominance testing functions in the following section. In the final section, we discuss the prepocessing functions.

## C++ CP-net Encoding
The R encoding of CP-nets we used can be found in the DQ-Pruning repository. In both of our Chapter 2 and 3 experiments, we use the CP-net and dominance query generators given in the DQ-Pruning repository and then pass these to our C++ function. However, these R objects do not translate directly to C++ inputs, so we convert CP-nets into a form we can utilise in C++ functions.

In C++ we input CP-nets to a function as 4 inputs `(A, N, ENTRIES, BREAK)`.
- `A` is a matrix giving the CP-net structure adjacency matrix. This is the first element of our old (R) CP-net encoding
- `N` is an n length vector (n - number of variables) giving the domain size of each variable
- `ENTRIES` and `BREAKS` are vectors which specify our CPTs as follows
DQ-Pruning repository
`ENTRIES`: Suppose variable X has 3 parents. If we let all variable domains be {1,2,...,|Domain|}, then we can view parental assignments as length 3vectors e.g (1,2,2) (given a parental order). In all of our functions variables are enumerated 1 to n, so we always have a variable (and parent) order. In general, we can view Pa(X) assignments as vectors and so we can put them into lexicographic order e.g. (1,1,1), (1,1,2), (1,2,1),(1,2,2),...., (2,2,2) if it is 3 binary variables. The entries for CPT(X) are the CPT(X) rows given in lexicographic order of the parent assignments. That is, we have the CPT row corresponding to (1,1,1) followed by the CPT row corresponding to (1,1,2) etc. The CPT row corresponding to (1,1,1) is a vector of |Dom(X)| values, where the ith entry is the preference position of X=i under parent assignment (1,1,1). Note that this is the entries of the CPT[X] in the R encoding corresponding to the (1,1,1) parent assignment. Thus, to obtain the new CPT(X) encoding we are essentially flattening the R CPT array into a vector, ordering rows lexicographically by parent assignment. The CPT(X) representation is then a |Dom(X)|\*(# parent assignments) length vector. `ENTRIES` is the CPT(1) vector followed by the CPT(2) vector etc.

`BREAKS`: This gives us the starting indices of each CPT in `ENTRIES`. That is, the ith entry of `BREAKS` is the index (k) such that CPT(i) starts at the kth entry of `ENTRIES`. If c_i is the length of the CPT(i) vector, then `BREAKS` is (1, 1+c_1, 1+c_1+c_2, ..., 1+c_1+...+c_n). The last element of `BREAKS` is redundant (no CPT starts at this index). However, it is useful in practice as the length of `ENTRIES` is c_1+c_2+...+c_n. `BREAKS` is useful as it allows us to extract CPT(i) from `ENTRIES` efficiently

The file `CPNConvert.R` contains the function `CPNConvert`. This R function takes an R structured CP-net at transforms it into the C++ encoding (outputting the elements in a list) so we can pass it to our C++ functions as follows. Suppose `CPN` is an R encoded CP-net, generated by our DQ-Pruning repository generator for example. Let `Function` be one of our C++ functions that requires input `(A, N, ENTRIES, BREAK)`. Then we apply `Function` to `N` in R as follows:
```
M=CPNConvert(CPN)
Function(M[[1]], M[[2]], M[[3]], M[[4]])
```

## Dominance Testing Functions:
In `DQFunctions.cpp` We give the following functions, which are the C++ translation of the dominance testing functions given in the DQ-Pruning repository (see here for further details). These are the functions we use in the Chapter 2 experiments of my thesis.
```
RankDQRankPriority(A, N, ENTRIES, BREAK, o1, o2)
RankDQRankDiffPriority(A, N, ENTRIES, BREAK, o1, o2)
PenaltyDQ(A, N, ENTRIES, BREAK, o1, o2)
SFDQ(A, N, ENTRIES, BREAK, o1, o2)
RSDQRankPriority(A, N, ENTRIES, BREAK, o1, o2)
RSDQRankDiffPriority(A, N, ENTRIES, BREAK, o1, o2)
RPDQRankPriority(A, N, ENTRIES, BREAK, o1, o2)
RPDQRankDiffPriority(A, N, ENTRIES, BREAK, o1, o2)
RPDQPenaltyPriority(A, N, ENTRIES, BREAK, o1, o2)
PSDQ(A, N, ENTRIES, BREAK, o1, o2)
RPSDQRankDiffPriority(A, N, ENTRIES, BREAK, o1, o2)
RPSDQRankPriority(A, N, ENTRIES, BREAK, o1, o2)
RPSDQPenaltyPriority(A, N, ENTRIES, BREAK, o1, o2)
```
Each function takes a CP-net and two outcomes, `o1` and `o2`. This is asking the dominance query "Is `o1` preferred to `o2`?". Note that outcomes are formated as vectors, the same as in DQ-Pruning (R outcomes can be passed straight to C++ functions without conversion).

The functions return the answer of the query (true/false) and the number of outcomes considered in the process of answering the query.

Functions starting with:
- `Rank` use rank pruning only
- `Penalty` use penalty pruning only
- `SF` use suffix fixing pruning only
- `RS` use rank pruning and suffix fixing
- `RP` use rank pruning and penalty pruning
- `PS` use penalty pruning and suffix fixing
- `RPS` use all three pruning methods

The priority method is either specified in the function title (Rank, Rank + Diff, Penalty) except in the cases where only 1 prioritisation is possible (see DQ-Pruning or Thesis Chapter 2 for details).

For each dominance testing function, we also have a corresponding timed version. For example, `RankDQRankPriority` has the corresponding function `RankDQRankPriorityTimed`. This function takes the same input and simply applies `RankDQRankPriority` and records the time it takes to answer the query. This function then returns the answer of the query (true/false), the number of outcomes considered, and the start and end time of the function (allowing us to calculate time taken). These are the functions we used in our Chapter 2 experiments in order to evaluate time elapsed results as well.

The following functions in `DQFunctions.cpp` are called within the dominance testing functions
```
CPTRow(A, N, ENTRIES, BREAK, x, PA)
Rank(A, N, ENTRIES, BREAK, o)
Penalty(A, N, ENTRIES, BREAK, o)
RankDifferences(A, N, ENTRIES, BREAK)
```

`CPTRow` takes a CP-net, a variable `x` (index between 1 and #variables), and `PA` a parental assignment to the parents of x (or 0 is x has no parents). This returns the row of CPT(X) that corresponds to this parental assignment (or just CPT(X) if x has no parents). See discussion of `ENTRIES` for what we mean by CPT row.

`Rank` and `Penalty` both take a CP-net and outcome, `o`, as inputs. They return the rank of `o` and the penalty of `o` respectively. 

`RankDifferences` takes a CP-net as input. It returns an n length vector (n=#variables). The ith element of this vector is the least rank improvement of variable i, L(i) (see thesis, section 2.4.1 for details).

Note that in these functions we enumerate our variables from 1, that is we have variables {1,2,...,n} and similarly, the domain of variable X is {1,2,...,|Dom(X)|}. In some cases we may use enumeration from 0 ({0,1,...,n-1}) for calculation convenience within the structure.

As we discussed above, parental assignments can be viewed as vectors and ordered lexicographically. This gives an enumeration of parent assignments e.g (1,1,1) is 1, (1,1,2) is 2 and so on. We will often move between parental assignments and their lexicographic enumeration within these functions (and the later preprocessing functions). We may also use the equivalent enumeration (1,1,1) is 0, (1,1,2) is 1 and so on. In certain cases, we may consider a restricted set of parental assignments, these can be put in lexicographic order and then enumerated sequentially and we may use this enumeration also (which will assign  values to assinments differently than the previous enumerations). These processes of enumeration can also be applied to outcomes viewed as vectors.

## Preprocessing Functions
The `DQFunctions.cpp` gives the following functions which we use to apply and test preprocessing in our Chapter 3 experiments:
```
UvRemove(A, N, ENTRIES, BREAK, o1, o2)
```
This function takes a CP-net and two outcomes - i.e. the dominance query "Is `o1` preferred to `o2`?". These outcomes, o1 and o2, must be distinct.

The function removes all variables from the CP-net that are unimportant to the query in the manner we discuss in Chapter 3. Any variable that has a parent (or several) removed by this process has all remaining parents tested for degeneracy. Any degenerate parents are removed (that is, the parent-child edge is removed from the structure). The function then returns the reduced CP-net and dominance query
```
ConnectedComponents(A)
```
This function takes an adjacency matrix, A, which we assume to be acyclic. 

The function obtains the connected components of A. It returns the number of connected components, k. If k>1, then it also returns a  kxn matrix (n=#variables), C. The ith row of C is a vector of 0s and 1s. the jth entry is 1 if and only if variable j is in the ith connected component. Note that there is no specific ordering of the connected components.
```
UVRSAndRPSDQRankPriority(A, N, ENTRIES, BREAK, o1, o2)
```
  //checks Are applied and if found false, DQFALSE is returned with 0 outcomes considered and start
  //and end times
  //Otherwise, once reduced it is separated into connected components. 
  //Each subDQ (if >1) is tested at formation with numerical checks.
  //If any found false DQFALSE is returned with 0 outcomes considered and start and end times
  //If all pass this is the point at with DonePreProc Time is recorded.
  //Each SubDQ is answered in order of increasing size. The outcomes considered each time are added on
  //If any is false, proccess immediately stops and DQFALSE is returned with
  //outcomes considered and time steps START, DONEPP, END
  //If all are true then proccess returns DQTRUE with total outcomes considered and 
  //time steps START, DONEPP, END
  //Further, if the function gets past PP then the reduced #outcomes is output
  
 uses numerical checks and logs outcomes traversed and times and #outcomes as explained in thesis
 using increasing size to order subCPN
 
 
 ForwardPruning(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2)
  //This function applies forward pruning. It then removes all size 1 variables and normalises the remaining domains so
  //that they are back to being enumerated 1,2,3,..
  //OUTPUTS: - If any domain goes to 0 (automatic false) then outputs the variable of interest and the Domains vector 
  //(pruned up to that point)
  // -OTHERWISEit outputs the pruned domains matrix, the reduced (pruned, normalised, degeneracies rmvd) CPN and assoc reduced DQ
  
ForwardPruning(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2)
Mention WHY we have to normalise and get rid of 1 value variables
output formats...

FPAndRPSDQRankPriority(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2)
applies FP, answers queries
outputs formated akin to ForwardPruning
applying FP with numerical checks, as in thesis. Outcomes considered etc. defined as in thesis

CombAndRPSDQRankPriority(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2)
diff to algorithm - thesis remark - add comment in function?
outputs same form as above
uses numerical checks too 

NumericalCheck(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2)
checks the three initial condition
false - any 1 of them holds
true - none hold (need to answer DQ, MIGHT be true)

ImpVar(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2)
input DQ
return a 0/1 vector of important variables
if o1=o2, returns a vector of 2s

Degenerate(IntegerVector NRed, IntegerVector ENTRIES)
//NRed - vector giving the domain sizes of Xs parents in order, then the domain size of X
  // ENTRIES - the CPT(X) section of the CP-net entries vector
  //Returns -1 if CPT(X) is non degenerate
  // Returns i (integer between 0 and |Pa| -1) giving the index of a degenerate parent (if CPT is degenerate)
  //Note that if we return that i is a degenerate parent, this does not mean that there aren't more degenerate parents
  
how we determine degen (theoretically), tho some easy cases
Parental asst enumeration - use partials to dtermine total

not consistent re variable enumeration tho always have to convert back to 0 indexing

RemovePa(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, int PA, int CH)
PA-CH is a degenerate edge we remove
PA, CH must be indexed from 0
Again using parent lex encodings to identify the right rows to extract more easily and cyclicng through is easier
have to upd entries and breaks at the end

