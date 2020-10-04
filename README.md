# CPPDomTestingAndPreprocessing

We start by randomly generating a set of CP-nets and outcomes in the exact same way as in the DQ-Pruning repository. These are the dominance queries we wish to test. They are stored as .RData objects, as before.

Then, we use the CPNConvert function to take these RObject CP-nets and convert them into a form that we can pass to our Rcpp functions (within R)
These functions return a list in R, which we can then manipulate as we like including storing the results.
The necessary set up to call the Rcpp functions in R is:...

CPTRow(A, N, ENTRIES, BREAK, x, PA)
Rank(A, N, ENTRIES, BREAK, o)
Penalty(A, N, ENTRIES, BREAK, o)
RankDifferences(A, N, ENTRIES, BREAK)

enumeration from 1 domains, variables

using Rcpp objects

lexicographic parent assts
parental lex enumerations actually the same.... No its weird. Made to start at 1.. make sure other repo starts at 0 in example, weird conversion  

CP-net form
outcomes - vectors

accurancy of ranks and encodings here is good enough for d<=5 and n<=20

RankDQRankPriority(A, N, ENTRIES, BREAK, o1, o2)
RankDQRankDiffPriority(A, N, ENTRIES, BREAK, o1, o2)
PenaltyDQ(A, N, ENTRIES, BREAK, o1, o2)
SFDQ(A, N, ENTRIES, BREAK, o1, o2)
RSDQRankPriority(A, N, ENTRIES, BREAK, o1, o2)
RSDQRankDiffPriority(A, N, ENTRIES, BREAK, o1, o2)
RPDQRankPriority(A, N, ENTRIES, BREAK, o1, o2)
RPDQRankDiffPriority
RPDQPenaltyPriority
PSDQ
RPSDQRankDiffPriority
RPSDQRankPriority
RPSDQPenaltyPriority

return answer and OT

NumericalCheck(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2)
checks the three initial condition
false - any 1 of them holds
true - none hold (need to answer DQ, MIGHT be true)

ImpVar(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2)
input DQ
return a 0/1 vector of important variables
if o1=o2, returns a vector of 2s

update ref-> paper!
