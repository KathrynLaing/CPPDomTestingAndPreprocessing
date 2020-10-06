#include <Rcpp.h>
#include <Rcpp/Benchmark/Timer.h>
#include <chrono>
#include <vector>
#include <limits>
using namespace Rcpp;

// [[Rcpp::export]]
/* Take CP- net, variable X and parent assignment as inputs, output preference row corresponding to parental assignment */
/* No Parents - parent asst input should just be the number 0*/
/*variables and variable assignments (including CPT entries) indexed from 1*/
IntegerVector CPTRow(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, int x, IntegerVector PA){
  //extract CPT(X) from entries
  int x0=BREAK[(x-1)];
  int x1=BREAK[(x)];
  IntegerVector CPT((x1-x0));
  for(int i=0; i<CPT.size();i++){
    CPT[i]=ENTRIES[(x0+i-1)];
  }
  //CPT is the CPT of variable x - the entries only, flattened out into a vector
  //nPa - number of parents of X
  int nPa=PA.size();
  if(PA[0]==0){
    // If X has no parents
    //CPT only has one parent assignment and only one row  so return this row
    return CPT;
  }
  else{
    //PaDom is the domain sizes of the parents of X in order
    IntegerVector PaDom(nPa);
    int j=0;
    //Note A.nrow is the number of variables
    for(int i=0;i<A.nrow();i++){
      if(A(i,(x-1))!=0){
        PaDom[j]=N[i];
        j+=1;
      }
    }
    //nPaAsst is the number associated with the input parent assignment when parent assignments are enumerated lexicographically
    int nPaAsst;
    if(PA.size()==1){
      //If X has only 1 parent
      nPaAsst=PA[0];
    }
    else{
      //If X has multiple parents
      nPaAsst=1;
      //iter is the lexicographic multiplier vector that we use to change between parental assignments and their lexicographic enumeration
      IntegerVector iter(nPa);
      for(int i=0;i<nPa;i++){
        int prod=1;
        for(int j=i+1;j<nPa;j++){
          prod=prod*PaDom[j];
        }
        iter[i]=prod;
      }
      for(int i=0;i<nPa;i++){
        if(PA[i]!=1){
          int a=PA[i]-1;
          nPaAsst+=(a*iter[i]);
        }
      }
    }
    nPaAsst-=1;
    //dom - domain size of X
    int dom=N[(x-1)];
    IntegerVector cptRow(dom);
    //Extract the nPaAsst-th row of CPT(X). This will be the preference row corresponding to the input parent assignment
    for(int i=0;i<dom;i++){
      cptRow[i]=CPT[(nPaAsst*dom)+i];
    }
    return cptRow;
  }
  /*We have returned the preference row of CPT(X) corresp to given parent asst */
}

// [[Rcpp::export]]
double Rank(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o){
  // Input a CPN and an outcome o
  //Output r(o)
  //Start by setting rank to 0
  double rank=0.0;
  //AF - vector of ancestral factors
  //Note A.nrow is the number of variables
  std::vector<double> AF(A.nrow());
  //AncSizes - vector of #ancestors for each variable
  IntegerVector AncSizes(A.nrow());
  //For each variable: 
  for(int i=0;i<A.nrow();i++){
    //calculate the ancestor set of variable i (Anc)
    //Start by adding parents and then repeatedly add the parents of every variable in the set until there is no change
    IntegerVector Anc(A.nrow(),0);
    int AncCheck=0;
    int AncSum=0;
    for(int j=0;j<A.nrow();j++){
      if(A(j,i)==1){
        Anc[j]=1;
        AncSum+=1;
      }
    }
    while(AncCheck!=AncSum){
      AncCheck=AncSum;
      for(int j=0;j<A.nrow();j++){
        if(Anc[j]==1){
          for(int k=0;k<A.nrow();k++){
            if(A(k,j)==1){
              Anc[k]=1;
            }
          }
        }
      }
      AncSum=0;
      for(int j=0;j<A.nrow();j++){
        AncSum+=Anc[j];
      }
    }
    //Anc is now a 0/1 vector giving the ancestors of variable i
    //AncSum - number of ancestors (non-zero entries in Anc) of variable i
    AncSizes[i]=AncSum;
    //Calculate the ancestral factor of variable i (AFi) by dividing 1 by the domain size of every ancestor of i
    double AFi=1.0;
    for(int j=0;j<A.nrow();j++){
      if(Anc[j]==1){
        AFi=AFi/N[j];
      }
    }
    AF[i]=AFi;
  }
  //DP- vector giving the number of descendent paths for each variable
  IntegerVector DP(A.nrow(),-1);
  int counter=0;
  int index;
  //until all DP values are calculated:
  while(counter<A.nrow()){
    //identify a variable (index) with maximal ancestor set which has not yet had descendent paths calculated
    //by definition all of the children of index have already had DP value calculated.
    index=-1;
    int ancsize=-1;
    for(int i=0;i<A.nrow();i++){
      if(DP[i]==(-1)){
        if(AncSizes[i]>ancsize){
          index=i;
          ancsize=AncSizes[i];
        }
      }
    }
    //calculate DP of index by summing (DP+1) over all children of index
    int dp=0;
    for(int i=0;i<A.nrow();i++){
      if(A(index,i)==1){
        dp+=(DP[i]+1);
      }
    }
    DP[index]=dp;
    //increase counter of #DP values calculated
    counter+=1;
  }
  //For each variable, calculate the weight contributed to rank
  for(int i=0;i<A.nrow();i++){
    //First, we must determine the preference position of the value taken by variable i in o
    //Parents - 0/1 vector of the parents of i
    //counter - number of parents of i
    IntegerVector Parents(A.nrow(),0);
    counter=0;
    for(int j=0;j<A.nrow();j++){
      if(A(j,i)==1){
        Parents[j]=1;
        counter+=1;
      }
    }
    IntegerVector pref;
    if(counter==0){
      //if theres no parents:
      //pref - the CPT of i (only one row as no parents)
      pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), 0);
    }
    else{
      // if i does have parents
      //PA values assigned to parents of i in o
      IntegerVector PA(counter);
      counter=0;
      for(int j=0;j<A.nrow();j++){
        if(Parents[j]==1){
          PA[counter]=o[j];
          counter+=1;
        }
      }
      //pref - The row of CPT(i) corresponding to the parents (of i) assigment given in o (PA)
      pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), PA);
    }
    //pp - position of preference term in rank formula
    //This is calculable with |Dom(i)|  (N[i]) because we can extract the preference position of i in o from pref (as this is the preference rule over i corresponding to the parental assignment in o)
    double pp=((N[i]-pref[o[i]-1]+1));
    pp=pp/(N[i]);
    //we can now calculate the weight contributed to the rank by variable i as follows (using the rank formula)
    rank+=AF[i]*(DP[i]+1)*pp;
  }
  //rank is now r(o) as we started at 0 and summed the weights contributed by all variables.
  return rank;
}

// [[Rcpp::export]]
long long int Penalty(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o){
  // Input a CPN and an outcome o
  //Output penalty(o)
  //Start by setting penalty to 0
  long long int penalty=0;
  //DecSizes - vector of #descendants for each variable
  //Note A.nrow is the number of variables
  IntegerVector DecSizes(A.nrow());
  IntegerVector Dec(A.nrow(),0);
  //For each variable: 
  for(int i=0;i<A.nrow();i++){
    //calculate the descendent set of variable i (Dec)
    for(int j=0;j<A.nrow();j++){
      Dec[j]=0;
    }
    int DecCheck=0;
    int DecSum=0;
    //Start by adding children of i to the set Dec
    for(int j=0;j<A.nrow();j++){
      if(A(i,j)==1){
        Dec[j]=1;
        DecSum+=1;
      }
    }
    //Then repeatedly add the children of every variable in the Dec until there is no change
    while(DecCheck!=DecSum){
      DecCheck=DecSum;
      for(int j=0;j<A.nrow();j++){
        if(Dec[j]==1){
          for(int k=0;k<A.nrow();k++){
            if(A(j,k)==1){
              Dec[k]=1;
            }
          }
        }
      }
      DecSum=0;
      for(int j=0;j<A.nrow();j++){
        DecSum+=Dec[j];
      }
    }
    //Dec is now a 0/1 vector showingthe descendants of i
    //DecSum - number of descendents of i (non-zero entries in Dec)
    DecSizes[i]=DecSum;
  }
  //ImpWeights- vector giving the importance weights for each variable
  std::vector<long long int> ImpWeights(A.nrow(),0);
  int counter=A.nrow();
  //for variables with no children, set importance weights to 1
  for(int i=0;i<A.nrow();i++){
    if(DecSizes[i]==0){
      ImpWeights[i]=1;
      counter-=1;
    }
  }
  //until all importance weights have been calculated:
  while(counter>0){
    //find a variable with minimum descendants that has not yet had its weight calculated (index)
    int LeastDec=A.nrow();
    int index;
    for(int i=0;i<A.nrow();i++){
      if(ImpWeights[i]==0){
        if(DecSizes[i]<LeastDec){
          LeastDec=DecSizes[i];
          index=i;
        }
      }
    }
    //calculate the importance weight of index according to the weight formula (using the importance weights of its children which must already be calculated)
    long long int impweight=0;
    for(int i=0;i<A.nrow();i++){
      if(A(index,i)==1){
        //if i is a child of index:
        impweight+=ImpWeights[i]*(N[i]-1);
      }
    }
    ImpWeights[index]=(impweight+1);
    //reduce the counter as another importance weight has been calculated
    counter-=1;
  }
  //For each variable, calculate the weight contributed to penalty
  for(int i=0;i<A.nrow();i++){
    //First, we must determine the preference position of the value taken by variable i in o
    //Parents 0/1 vector of the parents of i
    //counter - number of parents of i
    IntegerVector Parents(A.nrow(),0);
    counter=0;
    for(int j=0;j<A.nrow();j++){
      if(A(j,i)==1){
        Parents[j]=1;
        counter+=1;
      }
    }
    IntegerVector pref;
    if(counter==0){
      //if theres no parents:
      //pref - the CPT of i (only one row as no parents)
      pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), 0);
    }
    else{
      // if i does have parents
      //PA - values assigned to parents of i in o
      IntegerVector PA(counter);
      counter=0;
      for(int j=0;j<A.nrow();j++){
        if(Parents[j]==1){
          PA[counter]=o[j];
          counter+=1;
        }
      }
      //pref - The row of CPT(i) corresponding to the parent assigment given in o (to the parents of i)
      pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), PA);
    }
    //calclulate the weight added to penalty by variable i
    //This is importance weight x degree of penalty of assignment to variable i in o
    //The degree of penalty can be extracted from pref as this is the preference rule over i corresponding to the parental assignment in o
    penalty+=ImpWeights[i]*(pref[o[i]-1]-1);
  }
  //penalty is now penalty(o) as we started at 0 and summed the weights contributed by all variables.
  return penalty;
}

// [[Rcpp::export]]
std::vector<double> RankDifferences(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK){
  //Input - CP-net
  //Output - vector of the least rank improvements
  //We will store the least rank improvement in Diff
  std::vector<double> Diff(A.nrow());
  //AF - will be a vector of ancestral factors
  std::vector<double> AF(A.nrow());
  //AncSizes - will store the number of ancestors of each variable
  IntegerVector AncSizes(A.nrow());
  //For each variable, calculate the ancestral factor:
  for(int i=0;i<A.nrow();i++){
    //First calculate the ancestor set (Anc) of variable i as in the rank function as we did in the rank function
    IntegerVector Anc(A.nrow(),0);
    int AncCheck=0;
    int AncSum=0;
    for(int j=0;j<A.nrow();j++){
      if(A(j,i)==1){
        Anc[j]=1;
        AncSum+=1;
      }
    }
    while(AncCheck!=AncSum){
      AncCheck=AncSum;
      for(int j=0;j<A.nrow();j++){
        if(Anc[j]==1){
          for(int k=0;k<A.nrow();k++){
            if(A(k,j)==1){
              Anc[k]=1;
            }
          }
        }
      }
      AncSum=0;
      for(int j=0;j<A.nrow();j++){
        AncSum+=Anc[j];
      }
    }
    //AncSum is now the #ancestors of variable i
    AncSizes[i]=AncSum;
    //Calculate ancestral factor of i (AFi) as we did in the rank function:
    double AFi=1.0;
    for(int j=0;j<A.nrow();j++){
      if(Anc[j]==1){
        AFi=AFi/N[j];
      }
    }
    AF[i]=AFi;
  }
  //AF now contains ancestral factors and AncSizes contains #ancestors for each variable
  //DP - number of descendent paths for each variable, calculated in the same way as in the rank function
  IntegerVector DP(A.nrow(),-1);
  int counter=0;
  int index;
  while(counter<A.nrow()){
    index=-1;
    int ancsize=-1;
    for(int i=0;i<A.nrow();i++){
      if(DP[i]==(-1)){
        if(AncSizes[i]>ancsize){
          index=i;
          ancsize=AncSizes[i];
        }
      }
    }
    int dp=0;
    for(int i=0;i<A.nrow();i++){
      if(A(index,i)==1){
        dp+=(DP[i]+1);
      }
    }
    DP[index]=dp;
    counter+=1;
  }
  //For each variable, calculate the least rank improvement term
  for(int i=0;i<A.nrow();i++){
    //difference - least rank improvement of variable i
    double difference=AF[i]*(DP[i]+1);
    difference=difference/(N[i]);
    double ChildTerms=0.0;
    //for each child of i, calculate the deduction term
    for(int j=0;j<A.nrow();j++){
      if(A(i,j)==1){
        double term=AF[j]*(DP[j]+1)*(N[j]-1);
        term=term/(N[j]);
        ChildTerms+=term;
      }
    }
    difference-=ChildTerms;
    //difference is now the least rank improvement of variable i, calculated according to the given formula
    Diff[i]=difference;
  }
  //Diff - vector of least rank improvement terms for each variable.
  return Diff;
}

// [[Rcpp::export]]
List RankDQRankPriority(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2) {
  //Function takes a dominance query (o1>o2?) and answers it using rank + diff pruning with rank prioritisation
  //Omegais number of outcomes
  double Omega=1.0;
  //Multiply all domain sizes together
  //Note A.nrow is the number of variables
  for(int i=0;i<A.nrow();i++){
    Omega*=N[i];
  }
  int n=A.nrow();
  //diff- Hamming distance between o1 and o2 as vectors
  int diff=0;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      diff+=1;
    }
  }
  if(diff==0){
    //If o1=o2 then return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //rank1=r(o1)
  //rank2=r(o2)
  //RDiffs - vector of least rank improvement terms
  double rank1=Rank(A, N, ENTRIES, BREAK, o1);
  double rank2=Rank(A, N, ENTRIES, BREAK, o2);
  std::vector<double> RDiffs=RankDifferences(A, N, ENTRIES, BREAK);
  //comp - r(o1)-r(o2)-L_D(o1,o2) [least rank difference term]
  double comp=rank1-rank2;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      //for each i where o1 and o2 differ:
      comp-=RDiffs[i];
    }
  }
  double prec=-0.5/Omega;
  //if comp<0 our initial rank check fails.
  //However, as we are comparing the value of doubles and comp=0 does not fail the check, we use the prec term to account for possible rounding/precision errors.
  if(comp<prec){
    //If initial rank check fails then return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //LexMult - lexicographic multipliers used to transform outcomes (as vectors) to/from their lexicographic enumeration
  std::vector<long long int> LexMult(A.nrow());
  LexMult[(n-1)]=1;
  for(int i=2;i<=n;i++){
    LexMult[(n-i)]=LexMult[(n-i+1)]*N[(n-i+1)];
  }
  //SearchTree - list of outcomes in the search tree (lexicographic enumeration) - currently just o2
  //SearchTreeRanks - For each outcome in SearchTree, this gives rank value of outcomes that havent been considered and -1 for those that have
  std::vector<long long int> SearchTree(1);
  std::vector<double> SearchTreeRanks(1);
  //lex - lexicographic enumeration of outcome o2
  long long int lex=0;
  for(int i=0;i<n;i++){
    if(i!=(n-1)){
      lex+=(o2[i]-1)*LexMult[i];
    }
    else{
      lex+=(o2[i])*LexMult[i];
    }
  }
  //Add o2 (and its rank) to search tree
  SearchTree[0]=lex;
  SearchTreeRanks[0]=rank2;
  //unconsidered - number of outcomes in the tree not yet considered
  long long int unconsidered=1;
  //While there are unconsidered outcomes in the tree (i.e. while it is possible to search further):
  while(unconsidered>0){
    //identify the outcome in the tree that is not yet considered and has maximal rank (call this outcome 3)
    //index - index of this outcome in SearchTree
    //MaxRank is the maximum rank of the unconsidered outcomes
    double MaxRank=-1;
    long long int index;
    for(long long int i=0;i<SearchTree.size();i++){
      if(SearchTreeRanks[i]>MaxRank){
        MaxRank=SearchTreeRanks[i];
        index= i;
      }
    }
    //rank3 - rank of outcome 3
    double rank3=MaxRank;
    //Change the SearchTreeRanks value to -1 as this outcome is now `considered'
    //take 1 from #unconsidered
    SearchTreeRanks[index]=-1.0;
    unconsidered -= 1;
    // lex3 is the lexicographic enumeration of  outcome3
    long long int lex3=SearchTree[index];
    //Convert lex3 into the vector form of outcome 3
    IntegerVector o3(n);
    for(int i=0;i<n;i++){
      if((lex3 % LexMult[i])==0){
        o3[i]=lex3/LexMult[i];
        lex3-=(o3[i]-1)*LexMult[i];
      }
      else{
        o3[i]=(int)lex3/LexMult[i];
        o3[i]+=1;
        lex3-=(o3[i]-1)*LexMult[i];
      }
    }
    for(int i=0;i<n;i++){
      //for each variable i
      //Parents is  0/1 vecor of the parents of variable i
      //NPa is the number of parents
      IntegerVector Parents(n,0);
      int NPa=0;
      for(int j=0;j<n;j++){
        if(A(j,i)==1){
          Parents[j]=1;
          NPa+=1;
        }
      }
      //use CPTRow to extract the preference rule from CPT(i) corresponding to parent assignment in o3 as done in previous functions
      IntegerVector pref;
      IntegerVector PaAsst(NPa);
      if(NPa==0){
        pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), 0);
      }
      else{
        int counter=0;
        for(int j=0;j<n;j++){
          if(Parents[j]==1){
            PaAsst[counter]=o3[j];
            counter+=1;
          }
        }
        pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), PaAsst);
      }
      //pp - position of preference of variable i in o3 according to preference rule pref
      int pp=pref[o3[i]-1];
      if(pp>1){
        //if i can be improved (i.e. o3[i] is not the optimal assignment for i):
        //flips - vector of values that variable i could be flipped to and improve preference (according to preference rule pref)
        //obtain flips by finding which values have a lower preference position than pp in pref
        IntegerVector flips((pp-1));
        int counter=0;
        for(int j=0;j<N[i];j++){
          if(pref[j]<(pp)){
            flips[counter]=(j+1);
            counter+=1;
          }
        }
        //for each `better value' for i in flips
        for(int j=0;j<(pp-1);j++){
          //o4 is outcome 3 with variable i changed to take one of the better values in flips
          IntegerVector o4(n);
          for(int k=0;k<n;k++){
            o4[k]=o3[k];
          }
          o4[i]=flips[j];
          //diffs - hamming distance between o4 and o1
          int diffs=0;
          for(int k=0;k<n;k++){
            if(o1[k]!=o4[k]){
              diffs+=1;
            }
          }
          if(diffs==0){
            //if o4=o1 then we have reached o1 via an improving flip (o4) of outcome 3
            //return true and the number of outcomes in the search tree
            long long int Count=SearchTree.size();
            List R=List::create(true,Count);
            return R;
          }
          //rank4 - rank of o4
          double rank4=Rank(A, N, ENTRIES, BREAK, o4);
          //comp = r(o1)-r(o4)-L_D(o1,o4) [least rand difference between o1 and o4]
          double comp=rank1-rank4;
          for(int k=0;k<n;k++){
            if(o1[k]!=o4[k]){
              comp-=RDiffs[k];
            }
          }
          //If comp<0, then we prune o4 from the tree according to the rank pruning condition
          if(!(comp<prec)){
            //If o4 is not pruned by rank pruning:
            //lex4 - the lexicographic enumeration of o4
            long long int lex4=0;
            for(int k=0;k<n;k++){
              if(k!=(n-1)){
                lex4+=(o4[k]-1)*LexMult[k];
              }
              else{
                lex4+=(o4[k])*LexMult[k];
              }
            }
            //repeat - true if and only if o4 is already present in the tree
            bool repeat=false;
            for(long long int k=0;k<SearchTree.size();k++){
              if(SearchTree[k]==lex4){
                repeat=true;
              }
            }
            if(!repeat){
              //if o4 is not already in the tree:
              //add o4 to the search Tree and put its rank into SearchTreeRanks
              SearchTree.push_back (lex4);
              SearchTreeRanks.push_back (rank4);
              //add one to the number of unconsidered outcomes to account for the addition of o4 to the tree
              unconsidered +=1;
            }
          }
        }
      }
      //We have now assessed each improving i flip of outcome 3
      //If it a flip was not pruned by ranks and not already present in the tree, then this improving flip was added to the tree
    }
    //We have now added to the tree all improving flips of outcome 3 that were not pruned by rank pruning or already present
  }
  //if there are no more leaves in the tree to consider (and o1 was not found), then we cannot search further and the query must be false
  //return that DQ is false and the size of the tree
  long long int Count=SearchTree.size();
  List R=List::create(false,Count);
  return R;
}

// [[Rcpp::export]]
List RankDQRankDiffPriority(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2) {
  //Function takes a dominance query (o1>o2?) and answers it using rank + diff pruning with rank + diff prioritisation
  //Omegais number of outcomes
  double Omega=1.0;
  //Multiply all domain sizes together
  //Note A.nrow is the number of variables
  for(int i=0;i<A.nrow();i++){
    Omega*=N[i];
  }
  int n=A.nrow();
  //diff- Hamming distance between o1 and o2 as vectors
  int diff=0;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      diff+=1;
    }
  }
  if(diff==0){
    //If o1=o2 then return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //rank1=r(o1)
  //rank2=r(o2)
  //RDiffs - vector of least improvement terms
  double rank1=Rank(A, N, ENTRIES, BREAK, o1);
  double rank2=Rank(A, N, ENTRIES, BREAK, o2);
  double rankdiff=0.0;
  std::vector<double> RDiffs=RankDifferences(A, N, ENTRIES, BREAK);
  //comp - r(o1)-r(o2)-L_D(o1,o2) [least rank difference term]
  double comp=rank1-rank2;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      //for each i where o1 and o2 differ:
      comp-=RDiffs[i];
      rankdiff+=RDiffs[i];
    }
  }
  double prec=-0.5/Omega;
  //if comp<0 our initial rank check fails.
  //However, as we are comparing the value of doubles and comp=0 does not fail the check, we use the prec term to account for possible rounding/precision errors.
  if(comp<prec){
    //If initial rank check fails then return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //LexMult - lexicographic multipliers used to transform outcomes (as vectors) to/from their lexicographic enumeration
  std::vector<long long int> LexMult(A.nrow());
  LexMult[(n-1)]=1;
  for(int i=2;i<=n;i++){
    LexMult[(n-i)]=LexMult[(n-i+1)]*N[(n-i+1)];
  }
  //SearchTree - list of outcomes in the search tree (lexicographic enumeration) - currently just o2
  //SearchTreeRanks - For each outcome in SearchTree, this gives r(o)+L_D(o,o1) value of outcomes, o, that havent been considered and -1 for those that have
  std::vector<long long int> SearchTree(1);
  std::vector<double> SearchTreeRanks(1);
  //lex - lexicographic enumeration of outcome o2
  long long int lex=0;
  for(int i=0;i<n;i++){
    if(i!=(n-1)){
      lex+=(o2[i]-1)*LexMult[i];
    }
    else{
      lex+=(o2[i])*LexMult[i];
    }
  }
  //Add o2 (and its rank+L_D term) to search tree
  SearchTree[0]=lex;
  SearchTreeRanks[0]=(rank2+rankdiff);
  //unconsidered - number of outcomes in the tree not yet considered
  long long int unconsidered=1;
  //While there are unconsidered outcomes in the tree (i.e. while it is possible to search further):
  while(unconsidered>0){
    //identify the outcome in the tree that is not yet considered and has maximal rank+L_D term (call this outcome 3)
    //index - index of this outcome in SearchTree
    //MaxRank is the maximum rank+L_D of the unconsidered outcomes
    double MaxRank=-1;
    long long int index;
    for(long long int i=0;i<SearchTree.size();i++){
      if(SearchTreeRanks[i]>MaxRank){
        MaxRank=SearchTreeRanks[i];
        index= i;
      }
    }
    //rank3 - rank+L_D of outcome 3
    double rank3=MaxRank;
    //Change the SearchTreeRanks value to -1 as outcome 3 is now `considered'
    //take 1 from #unconsidered
    SearchTreeRanks[index]=-1.0;
    unconsidered -= 1;
    // lex3 is the lexicographic enumeration of  outcome3
    long long int lex3=SearchTree[index];
    //Convert lex3 into the vector form of outcome 3, o3
    IntegerVector o3(n);
    for(int i=0;i<n;i++){
      if((lex3 % LexMult[i])==0){
        o3[i]=lex3/LexMult[i];
        lex3-=(o3[i]-1)*LexMult[i];
      }
      else{
        o3[i]=(int)lex3/LexMult[i];
        o3[i]+=1;
        lex3-=(o3[i]-1)*LexMult[i];
      }
    }
    for(int i=0;i<n;i++){
      //for each variable i
      //Parents is  0/1 vecor of the parents of variable i
      //NPa is the number of parents
      IntegerVector Parents(n,0);
      int NPa=0;
      for(int j=0;j<n;j++){
        if(A(j,i)==1){
          Parents[j]=1;
          NPa+=1;
        }
      }
      //use CPTRow to extract the preference rule from CPT(i) corresponding to parental assignment (for i) in o3, as we have done in previous functions
      IntegerVector pref;
      IntegerVector PaAsst(NPa);
      if(NPa==0){
        pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), 0);
      }
      else{
        int counter=0;
        for(int j=0;j<n;j++){
          if(Parents[j]==1){
            PaAsst[counter]=o3[j];
            counter+=1;
          }
        }
        pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), PaAsst);
      }
      //pp - position of preference of variable i in o3 according to preference rule pref
      int pp=pref[o3[i]-1];
      if(pp>1){
        //if i can be improved (i.e. o3[i] is not the optimal assignment for i):
        //flips - vector of values that variable i could be flipped to and improve preference (according to preference rule pref)
        //obtain flips by finding which values have a lower preference position than pp in pref
        IntegerVector flips((pp-1));
        int counter=0;
        for(int j=0;j<N[i];j++){
          if(pref[j]<(pp)){
            flips[counter]=(j+1);
            counter+=1;
          }
        }
        //for each `better value' for i in flips
        for(int j=0;j<(pp-1);j++){
          //o4 is outcome 3 with variable i changed to take one of the better values in flips
          IntegerVector o4(n);
          for(int k=0;k<n;k++){
            o4[k]=o3[k];
          }
          o4[i]=flips[j];
          //diffs - hamming distance between o4 and o1
          int diffs=0;
          for(int k=0;k<n;k++){
            if(o1[k]!=o4[k]){
              diffs+=1;
            }
          }
          if(diffs==0){
            //if o4=o1 then we have reached o1 via an improving flip (o4) of outcome 3
                //return true and the number of outcomes in the search tree
            long long int Count=SearchTree.size();
            List R=List::create(true,Count);
            return R;
          }
          //rank4 - rank of o4
          double rank4=Rank(A, N, ENTRIES, BREAK, o4);
          //comp = r(o1)-r(o4)-L_D(o1,o4) [least rand difference between o1 and o4]
          double comp=rank1-rank4;
          double rankdiff=0.0;
          for(int k=0;k<n;k++){
            if(o1[k]!=o4[k]){
              comp-=RDiffs[k];
              rankdiff+=RDiffs[k];
            }
          }
          //If comp<0, then we prune o4 from the tree according to the rank pruning condition
          if(!(comp<prec)){
            //If o4 is not pruned by rank pruning:
            //lex4 - the lexicographic enumeration of o4
            long long int lex4=0;
            for(int k=0;k<n;k++){
              if(k!=(n-1)){
                lex4+=(o4[k]-1)*LexMult[k];
              }
              else{
                lex4+=(o4[k])*LexMult[k];
              }
            }
            //repeat - true if and only if o4 is already present in the tree
            bool repeat=false;
            for(long long int k=0;k<SearchTree.size();k++){
              if(SearchTree[k]==lex4){
                repeat=true;
              }
            }
            if(!repeat){
              //if o4 is not already in the tree:
              //add o4 to the search Tree and put its rank+L_D term into SearchTreeRanks
              SearchTree.push_back (lex4);
              SearchTreeRanks.push_back ((rank4+rankdiff));
              //add one to the number of unconsidered outcomes to account for the addition of o4 to the tree
              unconsidered +=1;
            }
          }
        }
      }
      //We have now assessed each improving i flip of outcome 3
      //If it was not pruned by ranks and not already present in the tree, then this improving flip was added to the tree
    }
    //We have now added to the tree all improving flips of outcome 3 that were not pruned by rank pruning or already present
  }
  //if there are no more leaves in the tree to consider (and o1 was not found), then we cannot search further and the query must be false
  //return that DQ is false and the size of the tree
  long long int Count=SearchTree.size();
  List R=List::create(false,Count);
  return R;
}

// [[Rcpp::export]]
List PenaltyDQ(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2) {
  //Function takes a dominance query and answers it using penalty pruning
  //n - number of variables
  int n=A.nrow();
  //diff- Hamming distance between o1 and o2 as vectors
  int diff=0;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      diff+=1;
    }
  }
  if(diff==0){
    //If o1=o2 then return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //pen1=penalty(o1)
  //pen2=penalty(o2)
  long long int pen1=Penalty(A, N, ENTRIES, BREAK, o1);
  long long int pen2=Penalty(A, N, ENTRIES, BREAK, o2);
  //f - evaluation function
  long long int f=pen2-pen1-diff;
  if(f<0){
    //if f<0 initial penalty check fails
    //return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //LexMult - lexicographic multipliers used to transform outcomes (as vectors) to/from their lexicographic enumeration
  std::vector<long long int> LexMult(A.nrow());
  LexMult[(n-1)]=1;
  for(int i=2;i<=n;i++){
    LexMult[(n-i)]=LexMult[(n-i+1)]*N[(n-i+1)];
  }
  //SearchTree - list of outcomes in the search tree (lexicographic enumeration) - currently just o2
  //SearchTreePenalties - For each outcome in SearchTree, this gives evaluation function (f) value of outcomes that havent been considered and -1 for those that have
  std::vector<long long int> SearchTree(1);
  std::vector<long long int> SearchTreePenalties(1);
  //lex - lexicographic enumeration of outcome o2
  long long int lex=0;
  for(int i=0;i<n;i++){
    if(i!=(n-1)){
      lex+=(o2[i]-1)*LexMult[i];
    }
    else{
      lex+=(o2[i])*LexMult[i];
    }
  }
  //Add o2 (and its f value) to search tree
  SearchTree[0]=lex;
  SearchTreePenalties[0]=f;
  //unconsidered - number of outcomes in the tree not yet considered
  long long int unconsidered=1;
  //While there are unconsidered outcomes in the tree (i.e. while it is possible to search further):
  while(unconsidered>0){
    //identify the outcome in the tree that is not yet considered and has minimal f value (call this outcome 3)
    //index - index of this outcome in SearchTree
    //Minf is the minimum f value of the unconsidered outcomes
    long long int Minf;
    long long int index;
    for(long long int i=0;i<SearchTree.size();i++){
      if(SearchTreePenalties[i]!=(-1)){
        Minf=SearchTreePenalties[i];
        index=i;
        break;
      }
    }
    for(long long int i=0;i<SearchTree.size();i++){
      if(SearchTreePenalties[i]!=(-1)){
        if(SearchTreePenalties[i]<Minf){
          Minf=SearchTreePenalties[i];
          index= i;
        }
      }
    }
    //f3 - f value for outcome 3
    long long int f3=Minf;
    //Change the SearchTreePenalties value to -1 as this outcome is now `considered'
    //take 1 from #unconsidered
    SearchTreePenalties[index]=-1;
    unconsidered -= 1;
    // lex3 is the lexicographic enumeration of  outcome3
    long long int lex3=SearchTree[index];
    //Convert lex3 into the vector form of outcome 3, o3
    IntegerVector o3(n);
    for(int i=0;i<n;i++){
      if((lex3 % LexMult[i])==0){
        o3[i]=lex3/LexMult[i];
        lex3-=(o3[i]-1)*LexMult[i];
      }
      else{
        o3[i]=(int)lex3/LexMult[i];
        o3[i]+=1;
        lex3-=(o3[i]-1)*LexMult[i];
      }
    }
    for(int i=0;i<n;i++){
      //for each variable i
      //Parents is  0/1 vecor of the parents of variable i
      //NPa is the number of parents
      IntegerVector Parents(n,0);
      int NPa=0;
      for(int j=0;j<n;j++){
        if(A(j,i)==1){
          Parents[j]=1;
          NPa+=1;
        }
      }
      //use CPTRow to extract the preference rule from CPT(i) corresponding to parental assignment (for i) in o3, as we have done in previous functions
      IntegerVector pref;
      IntegerVector PaAsst(NPa);
      if(NPa==0){
        pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), 0);
      }
      else{
        int counter=0;
        for(int j=0;j<n;j++){
          if(Parents[j]==1){
            PaAsst[counter]=o3[j];
            counter+=1;
          }
        }
        pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), PaAsst);
      }
      //pp - position of preference of variable i in o3 according to preference rule pref
      int pp=pref[o3[i]-1];
      if(pp>1){
        //if i can be improved (i.e. o3[i] is not the optimal assignment for i):
        //flips - vector of values that variable i could be flipped to and improve preference (according to preference rule pref)
        //obtain flips by finding which values have a lower preference position than pp in pref
        IntegerVector flips((pp-1));
        int counter=0;
        for(int j=0;j<N[i];j++){
          if(pref[j]<(pp)){
            flips[counter]=(j+1);
            counter+=1;
          }
        }
        //for each `better value' for i in flips
        for(int j=0;j<(pp-1);j++){
          //o4 is outcome 3 with variable i changed to take one of the better values in flips
          IntegerVector o4(n);
          for(int k=0;k<n;k++){
            o4[k]=o3[k];
          }
          o4[i]=flips[j];
          //diffs - hamming distance between o4 and o1
          int diffs=0;
          for(int k=0;k<n;k++){
            if(o1[k]!=o4[k]){
              diffs+=1;
            }
          }
          if(diffs==0){
            //if o4=o1 then we have reached o1 via an improving flip (o4) of outcome 3
            //return true and the number of outcomes in the search tree
            long long int Count=SearchTree.size();
            List R=List::create(true,Count);
            return R;
          }
          //pen4 - penalty of o4
          long long int pen4=Penalty(A, N, ENTRIES, BREAK, o4);
          //f - evaluation function value for o4
          long long int f=pen4-pen1-diffs;
          //If f<0, then we prune o4 from the tree according to the penalty pruning condition
          if(f>=0){
            //If o4 is not pruned by penalty pruning:
            //lex4 - the lexicographic enumeration of o4
            long long int lex4=0;
            for(int k=0;k<n;k++){
              if(k!=(n-1)){
                lex4+=(o4[k]-1)*LexMult[k];
              }
              else{
                lex4+=(o4[k])*LexMult[k];
              }
            }
            //repeat - true if and only if o4 is already present in the tree
            bool repeat=false;
            for(long long int k=0;k<SearchTree.size();k++){
              if(SearchTree[k]==lex4){
                repeat=true;
              }
            }
            if(!repeat){
              //if o4 is not already in the tree
              //add o4 to the search Tree and put its f value into SearchTreePenalties
              SearchTree.push_back (lex4);
              SearchTreePenalties.push_back (f);
              //add one to the number of unconsidered outcomes to account for the addition of o4 to the tree
              unconsidered +=1;
            }
          }
        }
      }
      //We have now assessed each improving i flip of outcome 3 
      //If it was not pruned by penalty pruning and not already present in the tree, then this improving flip was added to the tree
    }
    //We have now added to the tree all improving flips of outcome 3 that were not pruned by penalty pruning or already present
  }
  //if there are no more leaves in the tree to consider (and and o1 was not found), then we cannot search further and the query must be false
  //return that DQ is false and the size of the tree
  long long int Count=SearchTree.size();
  List R=List::create(false,Count);
  return R;
}

// [[Rcpp::export]]
List SFDQ(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2) {
  //Function takes a dominance query and answers it using Suffix Fixing pruning
  //n - number of variables
  int n=A.nrow();
  //diff- Hamming distance between o1 and o2 as vectors
  int diff=0;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      diff+=1;
    }
  }
  if(diff==0){
    //If o1=o2 then return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //AncSizes - vector of #ancestors for each variable
  IntegerVector AncSizes(A.nrow());
  for(int i=0;i<A.nrow();i++){
    //For variable i, calculate the ancestor set (Anc) in the same way as in rank function
    IntegerVector Anc(A.nrow(),0);
    int AncCheck=0;
    int AncSum=0;
    for(int j=0;j<A.nrow();j++){
      if(A(j,i)==1){
        Anc[j]=1;
        AncSum+=1;
      }
    }
    while(AncCheck!=AncSum){
      AncCheck=AncSum;
      for(int j=0;j<A.nrow();j++){
        if(Anc[j]==1){
          for(int k=0;k<A.nrow();k++){
            if(A(k,j)==1){
              Anc[k]=1;
            }
          }
        }
      }
      AncSum=0;
      for(int j=0;j<A.nrow();j++){
        AncSum+=Anc[j];
      }
    }
    //add #ancestors of i to AncSizes vector
    AncSizes[i]=AncSum;
  }
  // TO - vector giving variables in topological order
  // calculated by putting variables in increasing order according to AncSizes
  std::vector<int> TO(A.nrow());
  int count=n;
  while(count>0){
    int MinAnc=n;
    int index;
    for(int i=0;i<A.nrow();i++){
      if(AncSizes[i]<MinAnc){
        MinAnc=AncSizes[i];
        index=i;
      }
    }
    TO[(n-count)]=index;
    AncSizes[index]=n;
    count-=1;
  }
  //Note - AncSizes is no longer accurate due to modification in above calculation
  //LexMult - lexicographic multipliers used to transform outcomes (as vectors) to/from their lexicographic enumeration
  std::vector<long long int> LexMult(A.nrow());
  LexMult[(n-1)]=1;
  for(int i=2;i<=n;i++){
    LexMult[(n-i)]=LexMult[(n-i+1)]*N[(n-i+1)];
  }
  ///SearchTree - list of outcomes in the search tree (lexicographic enumeration) - currently just o2
  std::vector<long long int> SearchTree(1);
  //Position will be the next outcome index in the tree to be considered
  //We will simply consider the outcomes in SearchTree in order as they must be in increasing order of depth in the tree by construction
  long long int Position=0;
  //lex - lexicographic enumeration of outcome o2
  long long int lex=0;
  for(int i=0;i<n;i++){
    if(i!=(n-1)){
      lex+=(o2[i]-1)*LexMult[i];
    }
    else{
      lex+=(o2[i])*LexMult[i];
    }
  }
  //Add o2  to search tree
  SearchTree[0]=lex;
  //unconsidered - number of outcomes in the tree not yet considered
  long long int unconsidered=1;
  //While there are unconsidered outcomes in the tree (i.e. while it is possible to search further):
  while(unconsidered>0){
    // we next consider the outcome located at the Position'th entry of SearchTree, call this outcome 3
    //decrease the number of unconsidered outcomes by to account for this
    unconsidered -= 1;
    //lex3 - lexicographic enumeration of outcome 3
    long long int lex3=SearchTree[Position];
    //increase the Position varaible for the next iteration
    Position+=1;
    //Convert lex3 into the vector form of outcome 3, o3
    IntegerVector o3(n);
    for(int i=0;i<n;i++){
      if((lex3 % LexMult[i])==0){
        o3[i]=lex3/LexMult[i];
        lex3-=(o3[i]-1)*LexMult[i];
      }
      else{
        o3[i]=(int)lex3/LexMult[i];
        o3[i]+=1;
        lex3-=(o3[i]-1)*LexMult[i];
      }
    }
    //SuffN is the position of the first element of the matching suffix between o1 and o3 (according to TO ordering)
    //if no matching suffix then SuffN=n
    int SuffN=n;
    for(int i=0;i<n;i++){
      if(o1[TO[(n-i-1)]]==o3[TO[(n-1-i)]]){
        SuffN=n-i-1;
      }
      else{
        break;
      }
    }
    for(int i=0;i<n;i++){
      //for each variable i
      //TOi is the position of i in the TO
      int TOi;
      for(int j=0;j<n;j++){
        if(TO[j]==i){
          TOi=j;
        }
      }
      if(TOi<SuffN){
        //if i is NOT in the matching suffix:
        //Parents is  0/1 vecor of the parents of variable i
        //NPa is the number of parents
        IntegerVector Parents(n,0);
        int NPa=0;
        for(int j=0;j<n;j++){
          if(A(j,i)==1){
            Parents[j]=1;
            NPa+=1;
          }
        }
        //use CPTRow to extract the preference rule from CPT(i) corresponding to parental assignment (for i) in o3, as we have done in previous functions
        IntegerVector pref;
        IntegerVector PaAsst(NPa);
        if(NPa==0){
          pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), 0);
        }
        else{
          int counter=0;
          for(int j=0;j<n;j++){
            if(Parents[j]==1){
              PaAsst[counter]=o3[j];
              counter+=1;
            }
          }
          pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), PaAsst);
        }
        //pp - position of preference of variable i in o3 according to preference rule pref
        int pp=pref[o3[i]-1];
        if(pp>1){
          //if i can be improved (i.e. o3[i] is not the optimal assignment for i):
          //flips - vector of values that variable i could be flipped to and improve preference (according to preference rule pref)
          //obtain flips by finding which values have a lower preference position than pp in pref
          IntegerVector flips((pp-1));
          int counter=0;
          for(int j=0;j<N[i];j++){
            if(pref[j]<(pp)){
              flips[counter]=(j+1);
              counter+=1;
            }
          }
          //for each `better value' for i in flips
          for(int j=0;j<(pp-1);j++){
            //o4 is outcome 3 with variable i changed to take one of the better values in flips
            IntegerVector o4(n);
            for(int k=0;k<n;k++){
              o4[k]=o3[k];
            }
            o4[i]=flips[j];
            //diffs - hamming distance between o4 and o1
            int diffs=0;
            for(int k=0;k<n;k++){
              if(o1[k]!=o4[k]){
                diffs+=1;
              }
            }
            if(diffs==0){
              //if o4=o1 then we have reached o1 via an improving flip (o4) of outcome 3
                //return true and the number of outcomes in the search tree
              long long int Count=SearchTree.size();
              List R=List::create(true,Count);
              return R;
            }
            //lex4 - the lexicographic enumeration of o4
            long long int lex4=0;
            for(int k=0;k<n;k++){
              if(k!=(n-1)){
                lex4+=(o4[k]-1)*LexMult[k];
              }
              else{
                lex4+=(o4[k])*LexMult[k];
              }
            }
            //repeat - true if and only if o4 is already present in the tree
            bool repeat=false;
            for(long long int k=0;k<SearchTree.size();k++){
              if(SearchTree[k]==lex4){
                repeat=true;
              }
            }
            if(!repeat){
              //if o4 is not already in the tree:
              //add o4 to the search Tree
              SearchTree.push_back (lex4);
              //add one to the number of unconsidered outcomes to account for the addition of o4 to the tree
              unconsidered +=1;
            }
          }
        }
        //We have now assessed each improving i flip of outcome 3 for i not in a matching suffix
        //If it was not already present in the tree, then this improving flip was added to the tree
      }
    }
    //We have now added to the tree all improving flips of outcome 3 that were not pruned by suffix fixing or already present
  }
  //if there are no more leaves in the tree to consider (and and o1 was not found), then we cannot search further and the query must be false
  //return that DQ is false and the size of the tree
  long long int Count=SearchTree.size();
  List R=List::create(false,Count);
  return R;
}

// [[Rcpp::export]]
List RSDQRankPriority(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2) {
  //Function takes a dominance query (o1>o2?) and answers it using rank + diff pruning with rank prioritisation
  //Omegais number of outcomes
  double Omega=1.0;
  //Multiply all domain sizes together
  //Note A.nrow is the number of variables
  for(int i=0;i<A.nrow();i++){
    Omega*=N[i];
  }
  int n=A.nrow();
  //diff- Hamming distance between o1 and o2 as vectors
  int diff=0;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      diff+=1;
    }
  }
  if(diff==0){
    //If o1=o2 then return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //rank1=r(o1)
  //rank2=r(o2)
  //RDiffs - vector of least improvement terms
  double rank1=Rank(A, N, ENTRIES, BREAK, o1);
  double rank2=Rank(A, N, ENTRIES, BREAK, o2);
  std::vector<double> RDiffs=RankDifferences(A, N, ENTRIES, BREAK);
  //comp - r(o1)-r(o2)-L_D(o1,o2) [least rank difference term]
  double comp=rank1-rank2;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      //for each i where o1 and o2 differ:
      comp-=RDiffs[i];
    }
  }
  double prec=-0.5/Omega;
  //if comp<0 our initial rank check fails.
  //However, as we are comparing the value of doubles and comp=0 does not fail the check, we use the prec term to account for possible rounding/precision errors.
  if(comp<prec){
    //If initial rank check fails then return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //LexMult - lexicographic multipliers used to transform outcomes (as vectors) to/from their lexicographic enumeration
  std::vector<long long int> LexMult(A.nrow());
  LexMult[(n-1)]=1;
  for(int i=2;i<=n;i++){
    LexMult[(n-i)]=LexMult[(n-i+1)]*N[(n-i+1)];
  }
  //AncSizes - vector of #ancestors for each variable
  IntegerVector AncSizes(A.nrow());
  for(int i=0;i<A.nrow();i++){
    //For variable i, calculate the ancestor set (Anc) in the same way as in rank function
    IntegerVector Anc(A.nrow(),0);
    int AncCheck=0;
    int AncSum=0;
    for(int j=0;j<A.nrow();j++){
      if(A(j,i)==1){
        Anc[j]=1;
        AncSum+=1;
      }
    }
    while(AncCheck!=AncSum){
      AncCheck=AncSum;
      for(int j=0;j<A.nrow();j++){
        if(Anc[j]==1){
          for(int k=0;k<A.nrow();k++){
            if(A(k,j)==1){
              Anc[k]=1;
            }
          }
        }
      }
      AncSum=0;
      for(int j=0;j<A.nrow();j++){
        AncSum+=Anc[j];
      }
    }
    //add #ancestors of i to AncSizes vector
    AncSizes[i]=AncSum;
  }
  // TO - vector giving variables in topological order
  // calculated by putting variables in increasing order according to AncSizes
  std::vector<int> TO(A.nrow());
  int count=n;
  while(count>0){
    int MinAnc=n;
    int index;
    for(int i=0;i<A.nrow();i++){
      if(AncSizes[i]<MinAnc){
        MinAnc=AncSizes[i];
        index=i;
      }
    }
    TO[(n-count)]=index;
    AncSizes[index]=n;
    count-=1;
  }
  //Note - AncSizes is no longer accurate due to modification in above calculation
  //SearchTree - list of outcomes in the search tree (lexicographic enumeration) - currently just o2
  //SearchTreeRanks - For each outcome in SearchTree, this gives rank value of outcomes that havent been considered and -1 for those that have
  std::vector<long long int> SearchTree(1);
  std::vector<double> SearchTreeRanks(1);
  //lex - lexicographic enumeration of outcome o2
  long long int lex=0;
  for(int i=0;i<n;i++){
    if(i!=(n-1)){
      lex+=(o2[i]-1)*LexMult[i];
    }
    else{
      lex+=(o2[i])*LexMult[i];
    }
  }
  //Add o2 (and its rank) to search tree
  SearchTree[0]=lex;
  SearchTreeRanks[0]=rank2;
  //unconsidered - number of outcomes in the tree not yet considered
  long long int unconsidered=1;
  //While there are unconsidered outcomes in the tree (i.e. while it is possible to search further):
  while(unconsidered>0){
    //identify the outcome in the tree that is not yet considered and has maximal rank (call this outcome 3)
    //index - index of this outcome in SearchTree
    //MaxRank is the maximum rank of the unconsidered outcomes
    double MaxRank=-1;
    long long int index;
    for(long long int i=0;i<SearchTree.size();i++){
      if(SearchTreeRanks[i]>MaxRank){
        MaxRank=SearchTreeRanks[i];
        index= i;
      }
    }
    //rank3 - rank of outcome 3
    double rank3=MaxRank;
    //Change the SearchTreeRanks value to -1 as this outcome is now `considered'
    //take 1 from #unconsidered
    SearchTreeRanks[index]=-1.0;
    unconsidered -= 1;
    // lex3 is the lexicographic enumeration of  outcome3
    long long int lex3=SearchTree[index];
    //Convert lex3 into the vector form of outcome 3, o3
    IntegerVector o3(n);
    for(int i=0;i<n;i++){
      if((lex3 % LexMult[i])==0){
        o3[i]=lex3/LexMult[i];
        lex3-=(o3[i]-1)*LexMult[i];
      }
      else{
        o3[i]=(int)lex3/LexMult[i];
        o3[i]+=1;
        lex3-=(o3[i]-1)*LexMult[i];
      }
    }
    //SuffN is the position of the first element of the matching suffix between o1 and o3 (according to TO ordering)
    //if no matching suffix then SuffN=n
    int SuffN=n;
    for(int i=0;i<n;i++){
      if(o1[TO[(n-i-1)]]==o3[TO[(n-1-i)]]){
        SuffN=n-i-1;
      }
      else{
        break;
      }
    }
    for(int i=0;i<n;i++){
      //for each variable i
      //TOi is the position of i in the TO
      int TOi;
      for(int j=0;j<n;j++){
        if(TO[j]==i){
          TOi=j;
        }
      }
      if(TOi<SuffN){
        //if i is NOT in the matching suffix:
        //Parents is  0/1 vecor of the parents of variable i
        //NPa is the number of parents
        IntegerVector Parents(n,0);
        int NPa=0;
        for(int j=0;j<n;j++){
          if(A(j,i)==1){
            Parents[j]=1;
            NPa+=1;
          }
        }
        //use CPTRow to extract the preference rule from CPT(i) corresponding to parental assignment (for i) in o3, as we have done in previous functions
        IntegerVector pref;
        IntegerVector PaAsst(NPa);
        if(NPa==0){
          pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), 0);
        }
        else{
          int counter=0;
          for(int j=0;j<n;j++){
            if(Parents[j]==1){
              PaAsst[counter]=o3[j];
              counter+=1;
            }
          }
          pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), PaAsst);
        }
        //pp - position of preference of variable i in o3 according to preference rule pref
        int pp=pref[o3[i]-1];
        if(pp>1){
          //if i can be improved (i.e. o3[i] is not the optimal assignment for i):
          //flips - vector of values that variable i could be flipped to and improve preference (according to preference rule pref)
          //obtain flips by finding which values have a lower preference position than pp in pref
          IntegerVector flips((pp-1));
          int counter=0;
          for(int j=0;j<N[i];j++){
            if(pref[j]<(pp)){
              flips[counter]=(j+1);
              counter+=1;
            }
          }
          //for each `better value' for i in flips
          for(int j=0;j<(pp-1);j++){
            //o4 is outcome 3 with variable i changed to take one of the better values in flips
            IntegerVector o4(n);
            for(int k=0;k<n;k++){
              o4[k]=o3[k];
            }
            o4[i]=flips[j];
            //diffs - hamming distance between o4 and o1
            int diffs=0;
            for(int k=0;k<n;k++){
              if(o1[k]!=o4[k]){
                diffs+=1;
              }
            }
            if(diffs==0){
              //if o4=o1 then we have reached o1 via an improving flip (o4) of outcome 3
              //return true and the number of outcomes in the search tree
              long long int Count=SearchTree.size();
              List R=List::create(true,Count);
              return R;
            }
            //rank4 - rank of o4
            double rank4=Rank(A, N, ENTRIES, BREAK, o4);
            //comp = r(o1)-r(o4)-L_D(o1,o4) [least rand difference between o1 and o4]
            double comp=rank1-rank4;
            for(int k=0;k<n;k++){
              if(o1[k]!=o4[k]){
                comp-=RDiffs[k];
              }
            }
            //If comp<0, then we prune o4 from the tree according to the rank pruning condition
            if(!(comp<prec)){
              //If o4 is not pruned by rank pruning:
              //lex4 - the lexicographic enumeration of o4
              long long int lex4=0;
              for(int k=0;k<n;k++){
                if(k!=(n-1)){
                  lex4+=(o4[k]-1)*LexMult[k];
                }
                else{
                  lex4+=(o4[k])*LexMult[k];
                }
              }
              //repeat - true if and only if o4 is already present in the tree
              bool repeat=false;
              for(long long int k=0;k<SearchTree.size();k++){
                if(SearchTree[k]==lex4){
                  repeat=true;
                }
              }
              if(!repeat){
                //if o4 is not already in the tree
                //add o4 to the search Tree and put its rank into SearchTreeRanks
                SearchTree.push_back (lex4);
                SearchTreeRanks.push_back (rank4);
                //add one to the number of unconsidered outcomes to account for the addition of o4 to the tree
                unconsidered +=1;
              }
            }
          }
        }
        //We have now assessed each improving i flip of outcome 3 for i not in a matching suffix
        //If it was not pruned by ranks and not already present in the tree, then this improving flip was added to the tree
      }
    }
    //We have now added to the tree all improving flips of outcome 3 that were not pruned by suffix fixing, rank pruning, or already present
  }
  //if there are no more leaves in the tree to consider (and and o1 was not found), then we cannot search further and the query must be false
  //return that DQ is false and the size of the tree
  long long int Count=SearchTree.size();
  List R=List::create(false,Count);
  return R;
}

// [[Rcpp::export]]
List RSDQRankDiffPriority(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2) {
  //Function takes a dominance query (o1>o2?) and answers it using rank pruning + SF with rank + diff prioritisation
  //Omegais number of outcomes
  double Omega=1.0;
  //Multiply all domain sizes together
  //Note A.nrow is the number of variables
  for(int i=0;i<A.nrow();i++){
    Omega*=N[i];
  }
  int n=A.nrow();
  //diff- Hamming distance between o1 and o2 as vectors
  int diff=0;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      diff+=1;
    }
  }
  if(diff==0){
    //If o1=o2 then return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //rank1=r(o1)
  //rank2=r(o2)
  //RDiffs - vector of least improvement terms
  double rank1=Rank(A, N, ENTRIES, BREAK, o1);
  double rank2=Rank(A, N, ENTRIES, BREAK, o2);
  double rankdiff=0.0;
  std::vector<double> RDiffs=RankDifferences(A, N, ENTRIES, BREAK);
  //comp - r(o1)-r(o2)-L_D(o1,o2) [least rank difference term]
  double comp=rank1-rank2;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      //for each i where o1 and o2 differ:
      comp-=RDiffs[i];
      rankdiff+=RDiffs[i];
    }
  }
  double prec=-0.5/Omega;
  //if comp<0 our initial rank check fails.
  //However, as we are comparing the value of doubles and comp=0 does not fail the check, we use the prec term to account for possible rounding/precision errors.
  if(comp<prec){
    //If initial rank check fails then return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //LexMult - lexicographic multipliers used to transform outcomes (as vectors) to/from their lexicographic enumeration
  std::vector<long long int> LexMult(A.nrow());
  LexMult[(n-1)]=1;
  for(int i=2;i<=n;i++){
    LexMult[(n-i)]=LexMult[(n-i+1)]*N[(n-i+1)];
  }
  //AncSizes - vector of #ancestors for each variable
  IntegerVector AncSizes(A.nrow());
  for(int i=0;i<A.nrow();i++){
    //For variable i, calculate the ancestor set (Anc) in the same way as in rank function
    IntegerVector Anc(A.nrow(),0);
    int AncCheck=0;
    int AncSum=0;
    for(int j=0;j<A.nrow();j++){
      if(A(j,i)==1){
        Anc[j]=1;
        AncSum+=1;
      }
    }
    while(AncCheck!=AncSum){
      AncCheck=AncSum;
      for(int j=0;j<A.nrow();j++){
        if(Anc[j]==1){
          for(int k=0;k<A.nrow();k++){
            if(A(k,j)==1){
              Anc[k]=1;
            }
          }
        }
      }
      AncSum=0;
      for(int j=0;j<A.nrow();j++){
        AncSum+=Anc[j];
      }
    }
    ///add #ancestors of i to AncSizes vector
    AncSizes[i]=AncSum;
  }
  / TO - vector giving variables in topological order
  // calculated by putting variables in increasing order according to AncSizes
  std::vector<int> TO(A.nrow());
  int count=n;
  while(count>0){
    int MinAnc=n;
    int index;
    for(int i=0;i<A.nrow();i++){
      if(AncSizes[i]<MinAnc){
        MinAnc=AncSizes[i];
        index=i;
      }
    }
    TO[(n-count)]=index;
    AncSizes[index]=n;
    count-=1;
  }
  //Note - AncSizes is no longer accurate due to modification in above calculation
  //SearchTree - list of outcomes in the search tree (lexicographic enumeration) - currently just o2
  //SearchTreeRanks - For each outcome in SearchTree, this gives r(o)+L_D(o,o1) [rank + least rank difference] value of outcomes, o, that havent been considered and -1 for those that have
  std::vector<long long int> SearchTree(1);
  std::vector<double> SearchTreeRanks(1);
  //lex - lexicographic enumeration of outcome o2
  long long int lex=0;
  for(int i=0;i<n;i++){
    if(i!=(n-1)){
      lex+=(o2[i]-1)*LexMult[i];
    }
    else{
      lex+=(o2[i])*LexMult[i];
    }
  }
  //Add o2 (and its rank+L_D term) to search tree
  SearchTree[0]=lex;
  SearchTreeRanks[0]=(rank2+rankdiff);
  //unconsidered - number of outcomes in the tree not yet considered
  long long int unconsidered=1;
  //While there are unconsidered outcomes in the tree (i.e. while it is possible to search further):
  while(unconsidered>0){
    //identify the outcome in the tree that is not yet considered and has maximal rank+L_D (call this outcome 3)
    //index - index of this outcome in SearchTree
    //MaxRank is the maximum rank+L_D of the unconsidered outcomes
    double MaxRank=-1;
    long long int index;
    for(long long int i=0;i<SearchTree.size();i++){
      if(SearchTreeRanks[i]>MaxRank){
        MaxRank=SearchTreeRanks[i];
        index= i;
      }
    }
    //rank3 - rank+L_D of outcome 3
    double rank3=MaxRank;
    //Change the SearchTreeRanks value to -1 as this outcome is now `considered'
    //take 1 from #unconsidered
    SearchTreeRanks[index]=-1.0;
    unconsidered -= 1;
    // lex3 is the lexicographic enumeration of  outcome3
    long long int lex3=SearchTree[index];
    //Convert lex3 into the vector form of outcome 3, o3
    IntegerVector o3(n);
    for(int i=0;i<n;i++){
      if((lex3 % LexMult[i])==0){
        o3[i]=lex3/LexMult[i];
        lex3-=(o3[i]-1)*LexMult[i];
      }
      else{
        o3[i]=(int)lex3/LexMult[i];
        o3[i]+=1;
        lex3-=(o3[i]-1)*LexMult[i];
      }
    }
    //SuffN is the position of the first element of the matching suffix between o1 and o3 (according to TO ordering)
    //if no matching suffix then SuffN=n
    int SuffN=n;
    for(int i=0;i<n;i++){
      if(o1[TO[(n-i-1)]]==o3[TO[(n-1-i)]]){
        SuffN=n-i-1;
      }
      else{
        break;
      }
    }
    for(int i=0;i<n;i++){
      //for each variable i
      //TOi is the position of i in the TO
      int TOi;
      for(int j=0;j<n;j++){
        if(TO[j]==i){
          TOi=j;
        }
      }
      if(TOi<SuffN){
        //if i is NOT in the matching suffix:
        //Parents is  0/1 vecor of the parents of variable i
        //NPa is the number of parents
        IntegerVector Parents(n,0);
        int NPa=0;
        for(int j=0;j<n;j++){
          if(A(j,i)==1){
            Parents[j]=1;
            NPa+=1;
          }
        }
        //use CPTRow to extract the preference rule from CPT(i) corresponding to parental assignment (for i) in o3, as we have done in previous functions
        IntegerVector pref;
        IntegerVector PaAsst(NPa);
        if(NPa==0){
          pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), 0);
        }
        else{
          int counter=0;
          for(int j=0;j<n;j++){
            if(Parents[j]==1){
              PaAsst[counter]=o3[j];
              counter+=1;
            }
          }
          pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), PaAsst);
        }
        //pp - position of preference of variable i in o3 according to preference rule pref
        int pp=pref[o3[i]-1];
        if(pp>1){
          //if i can be improved (i.e. o3[i] is not the optimal assignment for i):
          //flips - vector of values that variable i could be flipped to and improve preference (according to preference rule pref)
          //obtain flips by finding which values have a lower preference position than pp in pref
          IntegerVector flips((pp-1));
          int counter=0;
          for(int j=0;j<N[i];j++){
            if(pref[j]<(pp)){
              flips[counter]=(j+1);
              counter+=1;
            }
          }
          //for each `better value' for i in flips
          for(int j=0;j<(pp-1);j++){
            //o4 is outcome 3 with variable i changed to take one of the better values in flips
            IntegerVector o4(n);
            for(int k=0;k<n;k++){
              o4[k]=o3[k];
            }
            o4[i]=flips[j];
            //diffs - hamming distance between o4 and o1
            int diffs=0;
            for(int k=0;k<n;k++){
              if(o1[k]!=o4[k]){
                diffs+=1;
              }
            }
            if(diffs==0){
              //if o4=o1 then we have reached o1 via an improving flip (o4) of outcome 3
              //return true and the number of outcomes in the search tree
              long long int Count=SearchTree.size();
              List R=List::create(true,Count);
              return R;
            }
            //rank4 - rank of o4
            double rank4=Rank(A, N, ENTRIES, BREAK, o4);
            //comp = r(o1)-r(o4)-L_D(o1,o4) [least rand difference between o1 and o4]
            double comp=rank1-rank4;
            double rankdiff=0.0;
            for(int k=0;k<n;k++){
              if(o1[k]!=o4[k]){
                comp-=RDiffs[k];
                rankdiff+=RDiffs[k];
              }
            }
            //If comp<0, then we prune o4 from the tree according to the rank pruning condition
            if(!(comp<prec)){
              //If o4 is not pruned by rank pruning:
              //lex4 - the lexicographic enumeration of o4
              long long int lex4=0;
              for(int k=0;k<n;k++){
                if(k!=(n-1)){
                  lex4+=(o4[k]-1)*LexMult[k];
                }
                else{
                  lex4+=(o4[k])*LexMult[k];
                }
              }
              //repeat - true if and only if o4 is already present in the tree
              bool repeat=false;
              for(long long int k=0;k<SearchTree.size();k++){
                if(SearchTree[k]==lex4){
                  repeat=true;
                }
              }
              if(!repeat){
                //if o4 is not already in the tree
                //add o4 to the search Tree and put its rank+L_D into SearchTreeRanks
                SearchTree.push_back (lex4);
                SearchTreeRanks.push_back ((rank4+rankdiff));
                //add one to the number of unconsidered outcomes to account for the addition of o4 to the tree
                unconsidered +=1;
              }
            }
          }
        }
        //We have now assessed each improving i flip of outcome 3 for i not in a matching suffix
        //If it was not pruned by ranks and not already present in the tree, then this improving flip was added to the tree
      }
    }
    //We have now added to the tree all improving flips of outcome 3 that were not pruned by suffix fixing, rank pruning, or already present
  }
  //if there are no more leaves in the tree to consider (and and o1 was not found), then we cannot search further and the query must be false
  //return that DQ is false and the size of the tree
  long long int Count=SearchTree.size();
  List R=List::create(false,Count);
  return R;
}

// [[Rcpp::export]]
List RPDQRankPriority(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2) {
  //Function takes a dominance query and answers it using rank + penalty pruning with rank prioritisation
  //Omegais number of outcomes
  double Omega=1.0;
  //Multiply all domain sizes together
  //Note A.nrow is the number of variables
  for(int i=0;i<A.nrow();i++){
    Omega*=N[i];
  }
  int n=A.nrow();
  //diff- Hamming distance between o1 and o2 as vectors
  int diff=0;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      diff+=1;
    }
  }
  if(diff==0){
    //If o1=o2 then return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //pen1=penalty(o1)
  //pen2=penalty(o2)
  long long int pen1=Penalty(A, N, ENTRIES, BREAK, o1);
  long long int pen2=Penalty(A, N, ENTRIES, BREAK, o2);
  //f - evaluation function
  long long int f=pen2-pen1-diff;
  if(f<0){
    //if f<0 initial penalty check fails
    //return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //rank1=r(o1)
  //rank2=r(o2)
  //RDiffs - vector of least improvement terms
  double rank1=Rank(A, N, ENTRIES, BREAK, o1);
  double rank2=Rank(A, N, ENTRIES, BREAK, o2);
  std::vector<double> RDiffs=RankDifferences(A, N, ENTRIES, BREAK);
  //comp - r(o1)-r(o2)-L_D(o1,o2) [least rank difference term]
  double comp=rank1-rank2;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      //for each i where o1 and o2 differ:
      comp-=RDiffs[i];
    }
  }
  double prec=-0.5/Omega;
  //if comp<0 our initial rank check fails.
  //However, as we are comparing the value of doubles and comp=0 does not fail the check, we use the prec term to account for possible rounding/precision errors.
  if(comp<prec){
    //If initial rank check fails then return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //LexMult - lexicographic multipliers used to transform outcomes (as vectors) to/from their lexicographic enumeration
  std::vector<long long int> LexMult(A.nrow());
  LexMult[(n-1)]=1;
  for(int i=2;i<=n;i++){
    LexMult[(n-i)]=LexMult[(n-i+1)]*N[(n-i+1)];
  }
  //SearchTree - list of outcomes in the search tree (lexicographic enumeration) - currently just o2
  //SearchTreeRanks - For each outcome in SearchTree, this gives rank value of outcomes that havent been considered and -1 for those that have
  std::vector<long long int> SearchTree(1);
  std::vector<double> SearchTreeRanks(1);
  //lex - lexicographic enumeration of outcome o2
  long long int lex=0;
  for(int i=0;i<n;i++){
    if(i!=(n-1)){
      lex+=(o2[i]-1)*LexMult[i];
    }
    else{
      lex+=(o2[i])*LexMult[i];
    }
  }
  //Add o2 (and its rank) to search tree
  SearchTree[0]=lex;
  SearchTreeRanks[0]=rank2;
  //unconsidered - number of outcomes in the tree not yet considered
  long long int unconsidered=1;
  //While there are unconsidered outcomes in the tree (i.e. while it is possible to search further):
  while(unconsidered>0){
    //identify the outcome in the tree that is not yet considered and has maximal rank (call this outcome 3)
    //index - index of this outcome in SearchTree
    //MaxRank is the maximum rank of the unconsidered outcomes
    double MaxRank=-1;
    long long int index;
    for(long long int i=0;i<SearchTree.size();i++){
      if(SearchTreeRanks[i]>MaxRank){
        MaxRank=SearchTreeRanks[i];
        index= i;
      }
    }
    //rank3 - rank of outcome 3
    double rank3=MaxRank;
    //Change the SearchTreeRanks value to -1 as this outcome is now `considered'
    //take 1 from #unconsidered
    SearchTreeRanks[index]=-1.0;
    unconsidered -= 1;
    // lex3 is the lexicographic enumeration of  outcome3
    long long int lex3=SearchTree[index];
    //Convert lex3 into the vector form of outcome 3, o3
    IntegerVector o3(n);
    for(int i=0;i<n;i++){
      if((lex3 % LexMult[i])==0){
        o3[i]=lex3/LexMult[i];
        lex3-=(o3[i]-1)*LexMult[i];
      }
      else{
        o3[i]=(int)lex3/LexMult[i];
        o3[i]+=1;
        lex3-=(o3[i]-1)*LexMult[i];
      }
    }
    for(int i=0;i<n;i++){
      //for each variable i
      //Parents is  0/1 vecor of the parents of variable i
      //NPa is the number of parents
      IntegerVector Parents(n,0);
      int NPa=0;
      for(int j=0;j<n;j++){
        if(A(j,i)==1){
          Parents[j]=1;
          NPa+=1;
        }
      }
      //use CPTRow to extract the preference rule from CPT(i) corresponding to parental assignment (for i) in o3, as we have done in previous functions
      IntegerVector pref;
      IntegerVector PaAsst(NPa);
      if(NPa==0){
        pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), 0);
      }
      else{
        int counter=0;
        for(int j=0;j<n;j++){
          if(Parents[j]==1){
            PaAsst[counter]=o3[j];
            counter+=1;
          }
        }
        pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), PaAsst);
      }
      //pp - position of preference of variable i in o3 according to preference rule pref
      int pp=pref[o3[i]-1];
      if(pp>1){
        //if i can be improved (i.e. o3[i] is not the optimal assignment for i):
        //flips - vector of values that variable i could be flipped to and improve preference (according to preference rule pref)
        //obtain flips by finding which values have a lower preference position than pp in pref
        IntegerVector flips((pp-1));
        int counter=0;
        for(int j=0;j<N[i];j++){
          if(pref[j]<(pp)){
            flips[counter]=(j+1);
            counter+=1;
          }
        }
        //for each `better value' for i in flips
        for(int j=0;j<(pp-1);j++){
          //o4 is outcome 3 with variable i changed to take one of the better values in flips
          IntegerVector o4(n);
          for(int k=0;k<n;k++){
            o4[k]=o3[k];
          }
          o4[i]=flips[j];
          //diffs - hamming distance between o4 and o1
          int diffs=0;
          for(int k=0;k<n;k++){
            if(o1[k]!=o4[k]){
              diffs+=1;
            }
          }
          if(diffs==0){
            //if o4=o1 then we have reached o1 via an improving flip (o4) of outcome 3
            //return true and the number of outcomes in the search tree
            long long int Count=SearchTree.size();
            List R=List::create(true,Count);
            return R;
          }
          //pen4 - penalty of o4
          long long int pen4=Penalty(A, N, ENTRIES, BREAK, o4);
          //f - evaluation function value for o4
          long long int f=pen4-pen1-diffs;
          //If f<0, then we prune o4 from the tree according to the penalty pruning condition
          if(f>=0){
            //If o4 is not pruned by penalty pruning:
            //rank4 - rank of o4
            double rank4=Rank(A, N, ENTRIES, BREAK, o4);
            //comp = r(o1)-r(o4)-L_D(o1,o4) [least rand difference between o1 and o4]
            double comp=rank1-rank4;
            for(int k=0;k<n;k++){
              if(o1[k]!=o4[k]){
                comp-=RDiffs[k];
              }
            }
            //If comp<0, then we prune o4 from the tree according to the rank pruning condition
            if(!(comp<prec)){
              //If o4 is not pruned by rank pruning:
              //lex4 - the lexicographic enumeration of o4
              long long int lex4=0;
              for(int k=0;k<n;k++){
                if(k!=(n-1)){
                  lex4+=(o4[k]-1)*LexMult[k];
                }
                else{
                  lex4+=(o4[k])*LexMult[k];
                }
              }
              //repeat - true if and only if o4 is already present in the tree
              bool repeat=false;
              for(long long int k=0;k<SearchTree.size();k++){
                if(SearchTree[k]==lex4){
                  repeat=true;
                }
              }
              if(!repeat){
                //if o4 is not already in the tree:
                //add o4 to the search Tree and put its rank into SearchTreeRanks
                SearchTree.push_back (lex4);
                SearchTreeRanks.push_back (rank4);
                //add one to the number of unconsidered outcomes to account for the addition of o4 to the tree
                unconsidered +=1;
              }
            }
          }
        }
      }
      //We have now assessed each improving i flip of outcome 3
      //If it was not pruned by rank or penalty pruning and not already present in the tree, then this improving flip was added to the tree
    }
    //We have now added to the tree all improving flips of outcome 3 that were not pruned by penalty pruning, rank pruning, or already present
  }
  //if there are no more leaves in the tree to consider (and and o1 was not found), then we cannot search further and the query must be false
  //return that DQ is false and the size of the tree
  long long int Count=SearchTree.size();
  List R=List::create(false,Count);
  return R;
}

// [[Rcpp::export]]
List RPDQRankDiffPriority(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2) {
  //Function takes a dominance query and answers it using rank + penalty pruning with rank + diff prioritisation
  //Omegais number of outcomes
  double Omega=1.0;
  //Multiply all domain sizes together
  //Note A.nrow is the number of variables
  for(int i=0;i<A.nrow();i++){
    Omega*=N[i];
  }
  int n=A.nrow();
  //diff- Hamming distance between o1 and o2 as vectors
  int diff=0;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      diff+=1;
    }
  }
  if(diff==0){
    //If o1=o2 then return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //pen1=penalty(o1)
  //pen2=penalty(o2)
  long long int pen1=Penalty(A, N, ENTRIES, BREAK, o1);
  long long int pen2=Penalty(A, N, ENTRIES, BREAK, o2);
  //f - evaluation function
  long long int f=pen2-pen1-diff;
  if(f<0){
    //if f<0 initial penalty check fails
    //return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //rank1=r(o1)
  //rank2=r(o2)
  //RDiffs - vector of least improvement terms
  double rank1=Rank(A, N, ENTRIES, BREAK, o1);
  double rank2=Rank(A, N, ENTRIES, BREAK, o2);
  double rankdiff=0.0;
  std::vector<double> RDiffs=RankDifferences(A, N, ENTRIES, BREAK);
  //comp - r(o1)-r(o2)-L_D(o1,o2) [least rank difference term]
  double comp=rank1-rank2;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      //for each i where o1 and o2 differ:
      comp-=RDiffs[i];
      rankdiff+=RDiffs[i];
    }
  }
  double prec=-0.5/Omega;
  //if comp<0 our initial rank check fails.
  //However, as we are comparing the value of doubles and comp=0 does not fail the check, we use the prec term to account for possible rounding/precision errors.
  if(comp<prec){
    //If initial rank check fails then return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //LexMult - lexicographic multipliers used to transform outcomes (as vectors) to/from their lexicographic enumeration
  std::vector<long long int> LexMult(A.nrow());
  LexMult[(n-1)]=1;
  for(int i=2;i<=n;i++){
    LexMult[(n-i)]=LexMult[(n-i+1)]*N[(n-i+1)];
  }
  //SearchTree - list of outcomes in the search tree (lexicographic enumeration) - currently just o2
  //SearchTreeRanks - For each outcome in SearchTree, this gives r(o)+L_D(o,o1) [rank + Least rank difference] value of outcomes, o, that havent been considered and -1 for those that have
  std::vector<long long int> SearchTree(1);
  std::vector<double> SearchTreeRanks(1);
  //lex - lexicographic enumeration of outcome o2
  long long int lex=0;
  for(int i=0;i<n;i++){
    if(i!=(n-1)){
      lex+=(o2[i]-1)*LexMult[i];
    }
    else{
      lex+=(o2[i])*LexMult[i];
    }
  }
  //Add o2 (and its rank+L_D term) to search tree
  SearchTree[0]=lex;
  SearchTreeRanks[0]=(rank2+rankdiff);
  //unconsidered - number of outcomes in the tree not yet considered
  long long int unconsidered=1;
  //While there are unconsidered outcomes in the tree (i.e. while it is possible to search further):
  while(unconsidered>0){
    //identify the outcome in the tree that is not yet considered and has maximal rank+L_D (call this outcome 3)
    //index - index of this outcome in SearchTree
    //MaxRank is the maximum rank+L_D of the unconsidered outcomes
    double MaxRank=-1;
    long long int index;
    for(long long int i=0;i<SearchTree.size();i++){
      if(SearchTreeRanks[i]>MaxRank){
        MaxRank=SearchTreeRanks[i];
        index= i;
      }
    }
    //rank3 - r(o3)+L_D(o1,o3) 
    double rank3=MaxRank;
    //Change the SearchTreeRanks value to -1 as this outcome is now `considered'
    //take 1 from #unconsidered
    SearchTreeRanks[index]=-1.0;
    unconsidered -= 1;
    // lex3 is the lexicographic enumeration of  outcome3
    long long int lex3=SearchTree[index];
    //Convert lex3 into the vector form of outcome 3, o3
    IntegerVector o3(n);
    for(int i=0;i<n;i++){
      if((lex3 % LexMult[i])==0){
        o3[i]=lex3/LexMult[i];
        lex3-=(o3[i]-1)*LexMult[i];
      }
      else{
        o3[i]=(int)lex3/LexMult[i];
        o3[i]+=1;
        lex3-=(o3[i]-1)*LexMult[i];
      }
    }
    for(int i=0;i<n;i++){
      //for each variable i
      //Parents is  0/1 vecor of the parents of variable i
      //NPa is the number of parents
      IntegerVector Parents(n,0);
      int NPa=0;
      for(int j=0;j<n;j++){
        if(A(j,i)==1){
          Parents[j]=1;
          NPa+=1;
        }
      }
      //use CPTRow to extract the preference rule from CPT(i) corresponding to parental assignment (for i) in o3, as we have done in previous functions
      IntegerVector pref;
      IntegerVector PaAsst(NPa);
      if(NPa==0){
        pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), 0);
      }
      else{
        int counter=0;
        for(int j=0;j<n;j++){
          if(Parents[j]==1){
            PaAsst[counter]=o3[j];
            counter+=1;
          }
        }
        pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), PaAsst);
      }
      //pp - position of preference of variable i in o3 according to preference rule pref
      int pp=pref[o3[i]-1];
      if(pp>1){
        //if i can be improved (i.e. o3[i] is not the optimal assignment for i):
        //flips - vector of values that variable i could be flipped to and improve preference (according to preference rule pref)
        //obtain flips by finding which values have a lower preference position than pp in pref
        IntegerVector flips((pp-1));
        int counter=0;
        for(int j=0;j<N[i];j++){
          if(pref[j]<(pp)){
            flips[counter]=(j+1);
            counter+=1;
          }
        }
        //for each `better value' for i in flips
        for(int j=0;j<(pp-1);j++){
          //o4 is outcome 3 with variable i changed to take one of the better values in flips
          IntegerVector o4(n);
          for(int k=0;k<n;k++){
            o4[k]=o3[k];
          }
          o4[i]=flips[j];
          //diffs - hamming distance between o4 and o1
          int diffs=0;
          for(int k=0;k<n;k++){
            if(o1[k]!=o4[k]){
              diffs+=1;
            }
          }
          if(diffs==0){
            //if o4=o1 then we have reached o1 via an improving flip (o4) of outcome 3
            //return true and the number of outcomes in the search tree
            long long int Count=SearchTree.size();
            List R=List::create(true,Count);
            return R;
          }
          //pen4 - penalty of o4
          long long int pen4=Penalty(A, N, ENTRIES, BREAK, o4);
          //f - evaluation function value for o4
          long long int f=pen4-pen1-diffs;
          //If f<0, then we prune o4 from the tree according to the penalty pruning condition
          if(f>=0){
            //If o4 is not pruned by penalty pruning:
            //rank4 - rank of o4
            double rank4=Rank(A, N, ENTRIES, BREAK, o4);
            //comp = r(o1)-r(o4)-L_D(o1,o4) [least rand difference between o1 and o4]
            double comp=rank1-rank4;
            double rankdiff=0.0;
            for(int k=0;k<n;k++){
              if(o1[k]!=o4[k]){
                comp-=RDiffs[k];
                rankdiff+=RDiffs[k];
              }
            }
            //If comp<0, then we prune o4 from the tree according to the rank pruning condition
            if(!(comp<prec)){
              //If o4 is not pruned by rank pruning:
              //lex4 - the lexicographic enumeration of o4
              long long int lex4=0;
              for(int k=0;k<n;k++){
                if(k!=(n-1)){
                  lex4+=(o4[k]-1)*LexMult[k];
                }
                else{
                  lex4+=(o4[k])*LexMult[k];
                }
              }
              //repeat - true if and only if o4 is already present in the tree
              bool repeat=false;
              for(long long int k=0;k<SearchTree.size();k++){
                if(SearchTree[k]==lex4){
                  repeat=true;
                }
              }
              if(!repeat){
                //if o4 is not already in the tree
                //add o4 to the search Tree and put its rank+L_D into SearchTreeRanks
                SearchTree.push_back (lex4);
                SearchTreeRanks.push_back ((rank4+rankdiff));
                //add one to the number of unconsidered outcomes to account for the addition of o4 to the tree
                unconsidered +=1;
              }
            }
          }
        }
      }
      //We have now assessed each improving i flip of outcome 3 
      //If it was not pruned by ranks or penalty pruning and not already present in the tree, then this improving flip was added to the tree
    }
    //We have now added to the tree all improving flips of outcome 3 that were not pruned by penalty or rank pruning, or already present
  }
  //if there are no more leaves in the tree to consider (and and o1 was not found), then we cannot search further and the query must be false
  //return that DQ is false and the size of the tree
  long long int Count=SearchTree.size();
  List R=List::create(false,Count);
  return R;
}

// [[Rcpp::export]]
List RPDQPenaltyPriority(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2) {
  //Function takes a dominance query and answers it using rank + penalty pruning with penalty priority
  //n - number of variables
  int n=A.nrow();
  //Omega is number of outcomes
  double Omega=1.0;
  //Multiply all domain sizes together
  //Note A.nrow is the number of variables
  for(int i=0;i<A.nrow();i++){
    Omega*=N[i];
  }
  //diff- Hamming distance between o1 and o2 as vectors
  int diff=0;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      diff+=1;
    }
  }
  if(diff==0){
    //If o1=o2 then return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //pen1=penalty(o1)
  //pen2=penalty(o2)
  long long int pen1=Penalty(A, N, ENTRIES, BREAK, o1);
  long long int pen2=Penalty(A, N, ENTRIES, BREAK, o2);
  //f - evaluation function
  long long int f=pen2-pen1-diff;
  if(f<0){
    //if f<0 initial penalty check fails
    //return false and 0 outcomes considere
    List R=List::create(false,0);
    return R;
  }
  //rank1=r(o1)
  //rank2=r(o2)
  //RDiffs - vector of least improvement terms
  double rank1=Rank(A, N, ENTRIES, BREAK, o1);
  double rank2=Rank(A, N, ENTRIES, BREAK, o2);
  std::vector<double> RDiffs=RankDifferences(A, N, ENTRIES, BREAK);
  //comp - r(o1)-r(o2)-L_D(o1,o2) [least rank difference term]
  double comp=rank1-rank2;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      //for each i where o1 and o2 differ:
      comp-=RDiffs[i];
    }
  }
  double prec=-0.5/Omega;
  //if comp<0 our initial rank check fails.
  //However, as we are comparing the value of doubles and comp=0 does not fail the check, we use the prec term to account for possible rounding/precision errors.
  if(comp<prec){
    //If initial rank check fails then return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //LexMult - lexicographic multipliers used to transform outcomes (as vectors) to/from their lexicographic enumeration
  std::vector<long long int> LexMult(A.nrow());
  LexMult[(n-1)]=1;
  for(int i=2;i<=n;i++){
    LexMult[(n-i)]=LexMult[(n-i+1)]*N[(n-i+1)];
  }
  //SearchTree - list of outcomes in the search tree (lexicographic enumeration) - currently just o2
  //SearchTreePenalties - For each outcome in SearchTree, this gives evaluation function (f) value of outcomes that havent been considered and -1 for those that have
  std::vector<long long int> SearchTree(1);
  std::vector<long long int> SearchTreePenalties(1);
  //lex - lexicographic enumeration of outcome o2
  long long int lex=0;
  for(int i=0;i<n;i++){
    if(i!=(n-1)){
      lex+=(o2[i]-1)*LexMult[i];
    }
    else{
      lex+=(o2[i])*LexMult[i];
    }
  }
  //Add o2 (and its f value) to search tree
  SearchTree[0]=lex;
  SearchTreePenalties[0]=f;
  //unconsidered - number of outcomes in the tree not yet considered
  long long int unconsidered=1;
  //While there are unconsidered outcomes in the tree (i.e. while it is possible to search further):
  while(unconsidered>0){
    //identify the outcome in the tree that is not yet considered and has minimal f value (call this outcome 3)
    //index - index of this outcome in SearchTree
    //Minf is the minimum f value of the unconsidered outcomes
    long long int Minf;
    long long int index;
    for(long long int i=0;i<SearchTree.size();i++){
      if(SearchTreePenalties[i]!=(-1)){
        Minf=SearchTreePenalties[i];
        index=i;
        break;
      }
    }
    for(long long int i=0;i<SearchTree.size();i++){
      if(SearchTreePenalties[i]!=(-1)){
        if(SearchTreePenalties[i]<Minf){
          Minf=SearchTreePenalties[i];
          index= i;
        }
      }
    }
    //f3 - f value for outcome 3
    long long int f3=Minf;
    //Change the SearchTreePenalties value to -1 as this outcome is now `considered'
    //take 1 from #unconsidered
    SearchTreePenalties[index]=-1;
    unconsidered -= 1;
    // lex3 is the lexicographic enumeration of  outcome3
    long long int lex3=SearchTree[index];
    //Convert lex3 into the vector form of outcome 3, o3
    IntegerVector o3(n);
    for(int i=0;i<n;i++){
      if((lex3 % LexMult[i])==0){
        o3[i]=lex3/LexMult[i];
        lex3-=(o3[i]-1)*LexMult[i];
      }
      else{
        o3[i]=(int)lex3/LexMult[i];
        o3[i]+=1;
        lex3-=(o3[i]-1)*LexMult[i];
      }
    }
    for(int i=0;i<n;i++){
      //for each variable i
      //Parents is  0/1 vecor of the parents of variable i
      //NPa is the number of parents
      IntegerVector Parents(n,0);
      int NPa=0;
      for(int j=0;j<n;j++){
        if(A(j,i)==1){
          Parents[j]=1;
          NPa+=1;
        }
      }
      //use CPTRow to extract the preference rule from CPT(i) corresponding to parental assignment (for i) in o3, as we have done in previous functions
      IntegerVector pref;
      IntegerVector PaAsst(NPa);
      if(NPa==0){
        pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), 0);
      }
      else{
        int counter=0;
        for(int j=0;j<n;j++){
          if(Parents[j]==1){
            PaAsst[counter]=o3[j];
            counter+=1;
          }
        }
        pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), PaAsst);
      }
      //pp - position of preference of variable i in o3 according to preference rule pref
      int pp=pref[o3[i]-1];
      if(pp>1){
        //if i can be improved (i.e. o3[i] is not the optimal assignment for i):
        //flips - vector of values that variable i could be flipped to and improve preference (according to preference rule pref)
        //obtain flips by finding which values have a lower preference position than pp in pref
        IntegerVector flips((pp-1));
        int counter=0;
        for(int j=0;j<N[i];j++){
          if(pref[j]<(pp)){
            flips[counter]=(j+1);
            counter+=1;
          }
        }
        //for each `better value' for i in flips
        for(int j=0;j<(pp-1);j++){
          //o4 is outcome 3 with variable i changed to take one of the better values in flips
          IntegerVector o4(n);
          for(int k=0;k<n;k++){
            o4[k]=o3[k];
          }
          o4[i]=flips[j];
          //diffs - hamming distance between o4 and o1
          int diffs=0;
          for(int k=0;k<n;k++){
            if(o1[k]!=o4[k]){
              diffs+=1;
            }
          }
          if(diffs==0){
            //if o4=o1 then we have reached o1 via an improving flip (o4) of outcome 3
            //return true and the number of outcomes in the search tree
            long long int Count=SearchTree.size();
            List R=List::create(true,Count);
            return R;
          }
          //pen4 - penalty of o4
          long long int pen4=Penalty(A, N, ENTRIES, BREAK, o4);
          //f - evaluation function value for o4
          long long int f=pen4-pen1-diffs;
          //If f<0, then we prune o4 from the tree according to the penalty pruning condition
          if(f>=0){
            //If o4 is not pruned by penalty pruning:
            //rank4 - rank of o4
            double rank4=Rank(A, N, ENTRIES, BREAK, o4);
            //comp = r(o1)-r(o4)-L_D(o1,o4) [least rand difference between o1 and o4]
            double comp=rank1-rank4;
            for(int k=0;k<n;k++){
              if(o1[k]!=o4[k]){
                comp-=RDiffs[k];
              }
            }
            //If comp<0, then we prune o4 from the tree according to the rank pruning condition
            if(!(comp<prec)){
              //If o4 is not pruned by rank pruning:
              //lex4 - the lexicographic enumeration of o4
              long long int lex4=0;
              for(int k=0;k<n;k++){
                if(k!=(n-1)){
                  lex4+=(o4[k]-1)*LexMult[k];
                }
                else{
                  lex4+=(o4[k])*LexMult[k];
                }
              }
              //repeat - true if and only if o4 is already present in the tree
              bool repeat=false;
              for(long long int k=0;k<SearchTree.size();k++){
                if(SearchTree[k]==lex4){
                  repeat=true;
                }
              }
              if(!repeat){
                //if o4 is not already in the tree
              //add o4 to the search Tree and put its f value into SearchTreePenalties
                SearchTree.push_back (lex4);
                SearchTreePenalties.push_back (f);
                //add one to the number of unconsidered outcomes to account for the addition of o4 to the tree
                unconsidered +=1;
              }
            }
          }
        }
      }
      //We have now assessed each improving i flip of outcome 3 
      //If it was not pruned by rank or penalty pruning and not already present in the tree, then this improving flip was added to the tree
    }
    //We have now added to the tree all improving flips of outcome 3 that were not pruned by rank or penalty pruning or already present
  }
  //if there are no more leaves in the tree to consider (and and o1 was not found), then we cannot search further and the query must be false
  //return that DQ is false and the size of the tree
  long long int Count=SearchTree.size();
  List R=List::create(false,Count);
  return R;
}

// [[Rcpp::export]]
List PSDQ(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2) {
  //Function takes a dominance query (o1>o2?) and answers it using penalty pruning + SF
  //n - number of variables
  int n=A.nrow();
  //diff- Hamming distance between o1 and o2 as vectors
  int diff=0;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      diff+=1;
    }
  }
  if(diff==0){
    //If o1=o2 then return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //pen1=penalty(o1)
  //pen2=penalty(o2)
  long long int pen1=Penalty(A, N, ENTRIES, BREAK, o1);
  long long int pen2=Penalty(A, N, ENTRIES, BREAK, o2);
  //f - evaluation function
  long long int f=pen2-pen1-diff;
  if(f<0){
    //if f<0 initial penalty check fails
    //return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //AncSizes - vector of #ancestors for each variable
  IntegerVector AncSizes(A.nrow());
  for(int i=0;i<A.nrow();i++){
    //For variable i, calculate the ancestor set (Anc) in the same way as in rank function
    IntegerVector Anc(A.nrow(),0);
    int AncCheck=0;
    int AncSum=0;
    for(int j=0;j<A.nrow();j++){
      if(A(j,i)==1){
        Anc[j]=1;
        AncSum+=1;
      }
    }
    while(AncCheck!=AncSum){
      AncCheck=AncSum;
      for(int j=0;j<A.nrow();j++){
        if(Anc[j]==1){
          for(int k=0;k<A.nrow();k++){
            if(A(k,j)==1){
              Anc[k]=1;
            }
          }
        }
      }
      AncSum=0;
      for(int j=0;j<A.nrow();j++){
        AncSum+=Anc[j];
      }
    }
    //add #ancestors of i to AncSizes vector
    AncSizes[i]=AncSum;
  }
  // TO - vector giving variables in topological order
  // calculated by putting variables in increasing order according to AncSizes
  std::vector<int> TO(A.nrow());
  int count=n;
  while(count>0){
    int MinAnc=n;
    int index;
    for(int i=0;i<A.nrow();i++){
      if(AncSizes[i]<MinAnc){
        MinAnc=AncSizes[i];
        index=i;
      }
    }
    TO[(n-count)]=index;
    AncSizes[index]=n;
    count-=1;
  }
  //Note - AncSizes is no longer accurate due to modification in above calculation
  //LexMult - lexicographic multipliers used to transform outcomes (as vectors) to/from their lexicographic enumeration
  std::vector<long long int> LexMult(A.nrow());
  LexMult[(n-1)]=1;
  for(int i=2;i<=n;i++){
    LexMult[(n-i)]=LexMult[(n-i+1)]*N[(n-i+1)];
  }
  //SearchTree - list of outcomes in the search tree (lexicographic enumeration) - currently just o2
  //SearchTreePenalties - For each outcome in SearchTree, this gives evaluation function (f) value of outcomes that havent been considered and -1 for those that have
  std::vector<long long int> SearchTree(1);
  std::vector<long long int> SearchTreePenalties(1);
  //lex - lexicographic enumeration of outcome o2
  long long int lex=0;
  for(int i=0;i<n;i++){
    if(i!=(n-1)){
      lex+=(o2[i]-1)*LexMult[i];
    }
    else{
      lex+=(o2[i])*LexMult[i];
    }
  }
  //Add o2 (and its f value) to search tree
  SearchTree[0]=lex;
  SearchTreePenalties[0]=f;
  //unconsidered - number of outcomes in the tree not yet considered
  long long int unconsidered=1;
  //While there are unconsidered outcomes in the tree (i.e. while it is possible to search further):
  while(unconsidered>0){
    //identify the outcome in the tree that is not yet considered and has minimal f value (call this outcome 3)
    //index - index of this outcome in SearchTree
    //Minf is the minimum f value of the unconsidered outcomes
    long long int Minf;
    long long int index;
    for(long long int i=0;i<SearchTree.size();i++){
      if(SearchTreePenalties[i]!=(-1)){
        Minf=SearchTreePenalties[i];
        index=i;
        break;
      }
    }
    for(long long int i=0;i<SearchTree.size();i++){
      if(SearchTreePenalties[i]!=(-1)){
        if(SearchTreePenalties[i]<Minf){
          Minf=SearchTreePenalties[i];
          index= i;
        }
      }
    }
    //f3 - f value for outcome 3
    long long int f3=Minf;
    //Change the SearchTreePenalties value to -1 as this outcome is now `considered'
    //take 1 from #unconsidered
    SearchTreePenalties[index]=-1;
    unconsidered -= 1;
    // lex3 is the lexicographic enumeration of  outcome3
    long long int lex3=SearchTree[index];
    //Convert lex3 into the vector form of outcome 3, o3
    IntegerVector o3(n);
    for(int i=0;i<n;i++){
      if((lex3 % LexMult[i])==0){
        o3[i]=lex3/LexMult[i];
        lex3-=(o3[i]-1)*LexMult[i];
      }
      else{
        o3[i]=(int)lex3/LexMult[i];
        o3[i]+=1;
        lex3-=(o3[i]-1)*LexMult[i];
      }
    }
    //SuffN is the position of the first element of the matching suffix between o1 and o3 (according to TO ordering)
    //if no matching suffix then SuffN=n
    int SuffN=n;
    for(int i=0;i<n;i++){
      if(o1[TO[(n-i-1)]]==o3[TO[(n-1-i)]]){
        SuffN=n-i-1;
      }
      else{
        break;
      }
    }
    for(int i=0;i<n;i++){
      //for each variable i
      //TOi is the position of i in the TO
      int TOi;
      for(int j=0;j<n;j++){
        if(TO[j]==i){
          TOi=j;
        }
      }
      if(TOi<SuffN){
        //if i is NOT in the matching suffix:
        //Parents is  0/1 vecor of the parents of variable i
        //NPa is the number of parents
        IntegerVector Parents(n,0);
        int NPa=0;
        for(int j=0;j<n;j++){
          if(A(j,i)==1){
            Parents[j]=1;
            NPa+=1;
          }
        }
        //use CPTRow to extract the preference rule from CPT(i) corresponding to parental assignment (for i) in o3, as we have done in previous functions
        IntegerVector pref;
        IntegerVector PaAsst(NPa);
        if(NPa==0){
          pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), 0);
        }
        else{
          int counter=0;
          for(int j=0;j<n;j++){
            if(Parents[j]==1){
              PaAsst[counter]=o3[j];
              counter+=1;
            }
          }
          pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), PaAsst);
        }
        //pp - position of preference of variable i in o3 according to preference rule pref
        int pp=pref[o3[i]-1];
        if(pp>1){
          //if i can be improved (i.e. o3[i] is not the optimal assignment for i):
          //flips - vector of values that variable i could be flipped to and improve preference (according to preference rule pref)
          //obtain flips by finding which values have a lower preference position than pp in pref
          IntegerVector flips((pp-1));
          int counter=0;
          for(int j=0;j<N[i];j++){
            if(pref[j]<(pp)){
              flips[counter]=(j+1);
              counter+=1;
            }
          }
          //for each `better value' for i in flips
          for(int j=0;j<(pp-1);j++){
            //o4 is outcome 3 with variable i changed to take one of the better values in flips
            IntegerVector o4(n);
            for(int k=0;k<n;k++){
              o4[k]=o3[k];
            }
            o4[i]=flips[j];
            //diffs - hamming distance between o4 and o1
            int diffs=0;
            for(int k=0;k<n;k++){
              if(o1[k]!=o4[k]){
                diffs+=1;
              }
            }
            if(diffs==0){
              //if o4=o1 then we have reached o1 via an improving flip (o4) of outcome 3
            //return true and the number of outcomes in the search tree
              long long int Count=SearchTree.size();
              List R=List::create(true,Count);
              return R;
            }
            //pen4 - penalty of o4
            long long int pen4=Penalty(A, N, ENTRIES, BREAK, o4);
            //f - evaluation function value for o4
            long long int f=pen4-pen1-diffs;
            //If f<0, then we prune o4 from the tree according to the penalty pruning condition
            if(f>=0){
              //If o4 is not pruned by penalty pruning:
            //lex4 - the lexicographic enumeration of o4
              long long int lex4=0;
              for(int k=0;k<n;k++){
                if(k!=(n-1)){
                  lex4+=(o4[k]-1)*LexMult[k];
                }
                else{
                  lex4+=(o4[k])*LexMult[k];
                }
              }
              //repeat - true if and only if o4 is already present in the tree
              bool repeat=false;
              for(long long int k=0;k<SearchTree.size();k++){
                if(SearchTree[k]==lex4){
                  repeat=true;
                }
              }
              if(!repeat){
                //if o4 is not already in the tree
                //add o4 to the search Tree and put its f value into SearchTreePenalties
                SearchTree.push_back (lex4);
                SearchTreePenalties.push_back (f);
                //add one to the number of unconsidered outcomes to account for the addition of o4 to the tree
                unconsidered +=1;
              }
            }
          }
        }
        //We have now assessed each improving i flip of outcome 3 for i not in a matching suffix
        //If it was not pruned by penalties and not already present in the tree, then this improving flip was added to the tree
      }
    }
    //We have now added to the tree all improving flips of outcome 3 that were not pruned by suffix fixing, penalty pruning, or already present
  }
  //if there are no more leaves in the tree to consider (and and o1 was not found), then we cannot search further and the query must be false
  //return that DQ is false and the size of the tree
  long long int Count=SearchTree.size();
  List R=List::create(false,Count);
  return R;
}

// [[Rcpp::export]]
List RPSDQRankDiffPriority(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2) {
  //Function takes a dominance query (o1>o2?) and answers it using all pruning methods with rank + diff prioritisation
  //Omegais number of outcomes
  double Omega=1.0;
  //Multiply all domain sizes together
  //Note A.row is the number of variables
  for(int i=0;i<A.nrow();i++){
    Omega*=N[i];
  }
  int n=A.nrow();
  //diff- Hamming distance between o1 and o2 as vectors
  int diff=0;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      diff+=1;
    }
  }
  if(diff==0){
    //If o1=o2 then return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //pen1=penalty(o1)
  //pen2=penalty(o2)
  long long int pen1=Penalty(A, N, ENTRIES, BREAK, o1);
  long long int pen2=Penalty(A, N, ENTRIES, BREAK, o2);
  //f - evaluation function
  long long int f=pen2-pen1-diff;
  if(f<0){
    //if f<0 initial penalty check fails
    //return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //rank1=r(o1)
  //rank2=r(o2)
  //RDiffs - vector of least improvement terms
  double rank1=Rank(A, N, ENTRIES, BREAK, o1);
  double rank2=Rank(A, N, ENTRIES, BREAK, o2);
  double rankdiff=0.0;
  std::vector<double> RDiffs=RankDifferences(A, N, ENTRIES, BREAK);
  //comp - r(o1)-r(o2)-L_D(o1,o2) [least rank difference term]
  double comp=rank1-rank2;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      //for each i where o1 and o2 differ:
      comp-=RDiffs[i];
      rankdiff+=RDiffs[i];
    }
  }
  double prec=-0.5/Omega;
  //if comp<0 our initial rank check fails.
  //However, as we are comparing the value of doubles and comp=0 does not fail the check, we use the prec term to account for possible rounding/precision errors.
  if(comp<prec){
    //If initial rank check fails then return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //AncSizes - vector of #ancestors for each variable
  IntegerVector AncSizes(A.nrow());
  for(int i=0;i<A.nrow();i++){
    //For variable i, calculate the ancestor set (Anc) in the same way as in rank function
    IntegerVector Anc(A.nrow(),0);
    int AncCheck=0;
    int AncSum=0;
    for(int j=0;j<A.nrow();j++){
      if(A(j,i)==1){
        Anc[j]=1;
        AncSum+=1;
      }
    }
    while(AncCheck!=AncSum){
      AncCheck=AncSum;
      for(int j=0;j<A.nrow();j++){
        if(Anc[j]==1){
          for(int k=0;k<A.nrow();k++){
            if(A(k,j)==1){
              Anc[k]=1;
            }
          }
        }
      }
      AncSum=0;
      for(int j=0;j<A.nrow();j++){
        AncSum+=Anc[j];
      }
    }
    //add #ancestors of i to AncSizes vector
    AncSizes[i]=AncSum;
  }
  // TO - vector giving variables in topological order
  // calculated by putting variables in increasing order according to AncSizes
  std::vector<int> TO(A.nrow());
  int count=n;
  while(count>0){
    int MinAnc=n;
    int index;
    for(int i=0;i<A.nrow();i++){
      if(AncSizes[i]<MinAnc){
        MinAnc=AncSizes[i];
        index=i;
      }
    }
    TO[(n-count)]=index;
    AncSizes[index]=n;
    count-=1;
  }
  //Note - AncSizes is no longer accurate due to modification in above calculation
  //LexMult - lexicographic multipliers used to transform outcomes (as vectors) to/from their lexicographic enumeration
  std::vector<long long int> LexMult(A.nrow());
  LexMult[(n-1)]=1;
  for(int i=2;i<=n;i++){
    LexMult[(n-i)]=LexMult[(n-i+1)]*N[(n-i+1)];
  }
  //SearchTree - list of outcomes in the search tree (lexicographic enumeration) - currently just o2
  //SearchTreeRanks - For each outcome in SearchTree, this gives r(o)+L_D(o,o1) [rank+least rank difference] value of outcomes, o, that havent been considered and -1 for those that have
  std::vector<long long int> SearchTree(1);
  std::vector<double> SearchTreeRanks(1);
  //lex - lexicographic enumeration of outcome o2
  long long int lex=0;
  for(int i=0;i<n;i++){
    if(i!=(n-1)){
      lex+=(o2[i]-1)*LexMult[i];
    }
    else{
      lex+=(o2[i])*LexMult[i];
    }
  }
  //Add o2 (and its rank+L_D term) to search tree
  SearchTree[0]=lex;
  SearchTreeRanks[0]=(rank2+rankdiff);
  //unconsidered - number of outcomes in the tree not yet considered
  long long int unconsidered=1;
  //While there are unconsidered outcomes in the tree (i.e. while it is possible to search further):
  while(unconsidered>0){
    //identify the outcome in the tree that is not yet considered and has maximal rank+L_D term (call this outcome 3)
    //index - index of this outcome in SearchTree
    //MaxRank is the maximum rank+L_D of the unconsidered outcomes
    double MaxRank=-1;
    long long int index;
    for(long long int i=0;i<SearchTree.size();i++){
      if(SearchTreeRanks[i]>MaxRank){
        MaxRank=SearchTreeRanks[i];
        index= i;
      }
    }
    //rank3 - rank+L_D of outcome 3
    double rank3=MaxRank;
    //Change the SearchTreeRanks value to -1 as this outcome is now `considered'
    //take 1 from #unconsidered
    SearchTreeRanks[index]=-1.0;
    unconsidered -= 1;
    // lex3 is the lexicographic enumeration of  outcome3
    long long int lex3=SearchTree[index];
    //Convert lex3 into the vector form of outcome 3, o3
    IntegerVector o3(n);
    for(int i=0;i<n;i++){
      if((lex3 % LexMult[i])==0){
        o3[i]=lex3/LexMult[i];
        lex3-=(o3[i]-1)*LexMult[i];
      }
      else{
        o3[i]=(int)lex3/LexMult[i];
        o3[i]+=1;
        lex3-=(o3[i]-1)*LexMult[i];
      }
    }
    //SuffN is the position of the first element of the matching suffix between o1 and o3 (according to TO ordering)
    //if no matching suffix then SuffN=n
    int SuffN=n;
    for(int i=0;i<n;i++){
      if(o1[TO[(n-i-1)]]==o3[TO[(n-1-i)]]){
        SuffN=n-i-1;
      }
      else{
        break;
      }
    }
    for(int i=0;i<n;i++){
      //for each variable i
      //TOi is the position of i in the TO
      int TOi;
      for(int j=0;j<n;j++){
        if(TO[j]==i){
          TOi=j;
        }
      }
      if(TOi<SuffN){
        //if i is NOT in the matching suffix:
        //Parents is  0/1 vecor of the parents of variable i
        //NPa is the number of parents
        IntegerVector Parents(n,0);
        int NPa=0;
        for(int j=0;j<n;j++){
          if(A(j,i)==1){
            Parents[j]=1;
            NPa+=1;
          }
        }
        //use CPTRow to extract the preference rule from CPT(i) corresponding to parental assignment (for i) in o3, as we have done in previous functions
        IntegerVector pref;
        IntegerVector PaAsst(NPa);
        if(NPa==0){
          pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), 0);
        }
        else{
          int counter=0;
          for(int j=0;j<n;j++){
            if(Parents[j]==1){
              PaAsst[counter]=o3[j];
              counter+=1;
            }
          }
          pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), PaAsst);
        }
        //pp - position of preference of variable i in o3 according to preference rule pref
        int pp=pref[o3[i]-1];
        if(pp>1){
          //if i can be improved (i.e. o3[i] is not the optimal assignment for i):
          //flips - vector of values that variable i could be flipped to and improve preference (according to preference rule pref)
          //obtain flips by finding which values have a lower preference position than pp in pref
          IntegerVector flips((pp-1));
          int counter=0;
          for(int j=0;j<N[i];j++){
            if(pref[j]<(pp)){
              flips[counter]=(j+1);
              counter+=1;
            }
          }
          //for each `better value' for i in flips
          for(int j=0;j<(pp-1);j++){
            //o4 is outcome 3 with variable i changed to take one of the better values in flips
            IntegerVector o4(n);
            for(int k=0;k<n;k++){
              o4[k]=o3[k];
            }
            o4[i]=flips[j];
            //diffs - hamming distance between o4 and o1
            int diffs=0;
            for(int k=0;k<n;k++){
              if(o1[k]!=o4[k]){
                diffs+=1;
              }
            }
            if(diffs==0){
              //if o4=o1 then we have reached o1 via an improving flip (o4) of outcome 3
                //return true and the number of outcomes in the search tree
              long long int Count=SearchTree.size();
              List R=List::create(true,Count);
              return R;
            }
            //pen4 - penalty of o4
            long long int pen4=Penalty(A, N, ENTRIES, BREAK, o4);
            //f - evaluation function value for o4
            long long int f=pen4-pen1-diffs;
            //If f<0, then we prune o4 from the tree according to the penalty pruning condition
            if(f>=0){
              //If o4 is not pruned by penalty pruning:
              //rank4 - rank of o4
              double rank4=Rank(A, N, ENTRIES, BREAK, o4);
              //comp = r(o1)-r(o4)-L_D(o1,o4) [least rand difference between o1 and o4]
              double comp=rank1-rank4;
              double rankdiff=0.0;
              for(int k=0;k<n;k++){
                if(o1[k]!=o4[k]){
                  comp-=RDiffs[k];
                  rankdiff+=RDiffs[k];
                }
              }
              //If comp<0, then we prune o4 from the tree according to the rank pruning condition
              if(!(comp<prec)){
                //If o4 is not pruned by rank pruning:
                //lex4 - the lexicographic enumeration of o4
                long long int lex4=0;
                for(int k=0;k<n;k++){
                  if(k!=(n-1)){
                    lex4+=(o4[k]-1)*LexMult[k];
                  }
                  else{
                    lex4+=(o4[k])*LexMult[k];
                  }
                }
                //repeat - true if and only if o4 is already present in the tree
                bool repeat=false;
                for(long long int k=0;k<SearchTree.size();k++){
                  if(SearchTree[k]==lex4){
                    repeat=true;
                  }
                }
                if(!repeat){
                  //if o4 is not already in the tree
                  //add o4 to the search Tree and put its rank+L_D into SearchTreeRanks
                  SearchTree.push_back (lex4);
                  SearchTreeRanks.push_back ((rank4+rankdiff));
                  //add one to the number of unconsidered outcomes to account for the addition of o4 to the tree
                  unconsidered +=1;
                }
              }
            }
          }
        }
        //We have now assessed each improving i flip of outcome 3 for i not in a matching suffix
        //If it was not pruned by ranks or penalites and not already present in the tree, then this improving flip was added to the tree
      }
    }
    //We have now added to the tree all improving flips of outcome 3 that were not pruned by suffix fixing, penalty pruning, or rank pruning, or already present
  }
  //if there are no more leaves in the tree to consider (and and o1 was not found), then we cannot search further and the query must be false
  //return that DQ is false and the size of the tree
  long long int Count=SearchTree.size();
  List R=List::create(false,Count);
  return R;
}

// [[Rcpp::export]]
List RPSDQRankPriority(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2) {
  //Function takes a dominance query (o1>o2?) and answers it using all3 pruning methods with rank prioritisation
  //Omegais number of outcomes
  double Omega=1.0;
  //Multiply all domain sizes together
  //Note A.nrow is the number of variables
  for(int i=0;i<A.nrow();i++){
    Omega*=N[i];
  }
  int n=A.nrow();
  //diff- Hamming distance between o1 and o2 as vectors
  int diff=0;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      diff+=1;
    }
  }
  if(diff==0){
    //If o1=o2 then return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //pen1=penalty(o1)
  //pen2=penalty(o2)
  long long int pen1=Penalty(A, N, ENTRIES, BREAK, o1);
  long long int pen2=Penalty(A, N, ENTRIES, BREAK, o2);
  //f - evaluation function
  long long int f=pen2-pen1-diff;
  if(f<0){
    //if f<0 initial penalty check fails
    //return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //rank1=r(o1)
  //rank2=r(o2)
  //RDiffs - vector of least improvement terms
  double rank1=Rank(A, N, ENTRIES, BREAK, o1);
  double rank2=Rank(A, N, ENTRIES, BREAK, o2);
  std::vector<double> RDiffs=RankDifferences(A, N, ENTRIES, BREAK);
  //comp - r(o1)-r(o2)-L_D(o1,o2) [least rank difference term]
  double comp=rank1-rank2;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      //for each i where o1 and o2 differ:
      comp-=RDiffs[i];
    }
  }
  double prec=-0.5/Omega;
  //if comp<0 our initial rank check fails.
  //However, as we are comparing the value of doubles and comp=0 does not fail the check, we use the prec term to account for possible rounding/precision errors.
  if(comp<prec){
    //If initial rank check fails then return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //AncSizes - vector of #ancestors for each variable
  IntegerVector AncSizes(A.nrow());
  for(int i=0;i<A.nrow();i++){
    //For variable i, calculate the ancestor set (Anc) in the same way as in rank function
    IntegerVector Anc(A.nrow(),0);
    int AncCheck=0;
    int AncSum=0;
    for(int j=0;j<A.nrow();j++){
      if(A(j,i)==1){
        Anc[j]=1;
        AncSum+=1;
      }
    }
    while(AncCheck!=AncSum){
      AncCheck=AncSum;
      for(int j=0;j<A.nrow();j++){
        if(Anc[j]==1){
          for(int k=0;k<A.nrow();k++){
            if(A(k,j)==1){
              Anc[k]=1;
            }
          }
        }
      }
      AncSum=0;
      for(int j=0;j<A.nrow();j++){
        AncSum+=Anc[j];
      }
    }
    //add #ancestors of i to AncSizes vector
    AncSizes[i]=AncSum;
  }
  // TO - vector giving variables in topological order
  // calculated by putting variables in increasing order according to AncSizes
  std::vector<int> TO(A.nrow());
  int count=n;
  while(count>0){
    int MinAnc=n;
    int index;
    for(int i=0;i<A.nrow();i++){
      if(AncSizes[i]<MinAnc){
        MinAnc=AncSizes[i];
        index=i;
      }
    }
    TO[(n-count)]=index;
    AncSizes[index]=n;
    count-=1;
  }
  //Note - AncSizes is no longer accurate due to modification in above calculation
  //LexMult - lexicographic multipliers used to transform outcomes (as vectors) to/from their lexicographic enumeration
  std::vector<long long int> LexMult(A.nrow());
  LexMult[(n-1)]=1;
  for(int i=2;i<=n;i++){
    LexMult[(n-i)]=LexMult[(n-i+1)]*N[(n-i+1)];
  }
  //SearchTree - list of outcomes in the search tree (lexicographic enumeration) - currently just o2
  //SearchTreeRanks - For each outcome in SearchTree, this gives rank value of outcomes that havent been considered and -1 for those that have
  std::vector<long long int> SearchTree(1);
  std::vector<double> SearchTreeRanks(1);
  //lex - lexicographic enumeration of outcome o2
  long long int lex=0;
  for(int i=0;i<n;i++){
    if(i!=(n-1)){
      lex+=(o2[i]-1)*LexMult[i];
    }
    else{
      lex+=(o2[i])*LexMult[i];
    }
  }
  //Add o2 (and its rank) to search tree
  SearchTree[0]=lex;
  SearchTreeRanks[0]=rank2;
  //unconsidered - number of outcomes in the tree not yet considered
  long long int unconsidered=1;
  //While there are unconsidered outcomes in the tree (i.e. while it is possible to search further):
  while(unconsidered>0){
    //identify the outcome in the tree that is not yet considered and has maximal rank (call this outcome 3)
    //index - index of this outcome in SearchTree
    //MaxRank is the maximum rank of the unconsidered outcomes
    double MaxRank=-1;
    long long int index;
    for(long long int i=0;i<SearchTree.size();i++){
      if(SearchTreeRanks[i]>MaxRank){
        MaxRank=SearchTreeRanks[i];
        index= i;
      }
    }
    //rank3 - rank of outcome 3
    double rank3=MaxRank;
    //Change the SearchTreeRanks value to -1 as this outcome is now `considered'
    //take 1 from #unconsidered
    SearchTreeRanks[index]=-1.0;
    unconsidered -= 1;
    // lex3 is the lexicographic enumeration of  outcome3
    long long int lex3=SearchTree[index];
    //Convert lex3 into the vector form of outcome 3, o3
    IntegerVector o3(n);
    for(int i=0;i<n;i++){
      if((lex3 % LexMult[i])==0){
        o3[i]=lex3/LexMult[i];
        lex3-=(o3[i]-1)*LexMult[i];
      }
      else{
        o3[i]=(int)lex3/LexMult[i];
        o3[i]+=1;
        lex3-=(o3[i]-1)*LexMult[i];
      }
    }
    //SuffN is the position of the first element of the matching suffix between o1 and o3 (according to TO ordering)
    //if no matching suffix then SuffN=n
    int SuffN=n;
    for(int i=0;i<n;i++){
      if(o1[TO[(n-i-1)]]==o3[TO[(n-1-i)]]){
        SuffN=n-i-1;
      }
      else{
        break;
      }
    }
    for(int i=0;i<n;i++){
      //for each variable i
      //TOi is the position of i in the TO
      int TOi;
      for(int j=0;j<n;j++){
        if(TO[j]==i){
          TOi=j;
        }
      }
      if(TOi<SuffN){
        //if i is NOT in the matching suffix:
        //Parents is  0/1 vecor of the parents of variable i
        //NPa is the number of parents
        IntegerVector Parents(n,0);
        int NPa=0;
        for(int j=0;j<n;j++){
          if(A(j,i)==1){
            Parents[j]=1;
            NPa+=1;
          }
        }
        //use CPTRow to extract the preference rule from CPT(i) corresponding to parental assignment (for i) in o3, as we have done in previous functions
        IntegerVector pref;
        IntegerVector PaAsst(NPa);
        if(NPa==0){
          pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), 0);
        }
        else{
          int counter=0;
          for(int j=0;j<n;j++){
            if(Parents[j]==1){
              PaAsst[counter]=o3[j];
              counter+=1;
            }
          }
          pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), PaAsst);
        }
        //pp - position of preference of variable i in o3 according to preference rule pref
        int pp=pref[o3[i]-1];
        if(pp>1){
          //if i can be improved (i.e. o3[i] is not the optimal assignment for i):
          //flips - vector of values that variable i could be flipped to and improve preference (according to preference rule pref)
          //obtain flips by finding which values have a lower preference position than pp in pref
          IntegerVector flips((pp-1));
          int counter=0;
          for(int j=0;j<N[i];j++){
            if(pref[j]<(pp)){
              flips[counter]=(j+1);
              counter+=1;
            }
          }
          //for each `better value' for i in flips
          for(int j=0;j<(pp-1);j++){
            //o4 is outcome 3 with variable i changed to take one of the better values in flips
            IntegerVector o4(n);
            for(int k=0;k<n;k++){
              o4[k]=o3[k];
            }
            o4[i]=flips[j];
            //diffs - hamming distance between o4 and o1
            int diffs=0;
            for(int k=0;k<n;k++){
              if(o1[k]!=o4[k]){
                diffs+=1;
              }
            }
            if(diffs==0){
              //if o4=o1 then we have reached o1 via an improving flip (o4) of outcome 3
                //return true and the number of outcomes in the search tree
              long long int Count=SearchTree.size();
              List R=List::create(true,Count);
              return R;
            }
            //pen4 - penalty of o4
            long long int pen4=Penalty(A, N, ENTRIES, BREAK, o4);
            //f - evaluation function value for o4
            long long int f=pen4-pen1-diffs;
            //If f<0, then we prune o4 from the tree according to the penalty pruning condition
            if(f>=0){
              //If o4 is not pruned by penalty pruning:
              //rank4 - rank of o4
              double rank4=Rank(A, N, ENTRIES, BREAK, o4);
              //comp = r(o1)-r(o4)-L_D(o1,o4) [least rand difference between o1 and o4]
              double comp=rank1-rank4;
              for(int k=0;k<n;k++){
                if(o1[k]!=o4[k]){
                  comp-=RDiffs[k];
                }
              }
              //If comp<0, then we prune o4 from the tree according to the rank pruning condition
              if(!(comp<prec)){
                //If o4 is not pruned by rank pruning:
                //lex4 - the lexicographic enumeration of o4
                long long int lex4=0;
                for(int k=0;k<n;k++){
                  if(k!=(n-1)){
                    lex4+=(o4[k]-1)*LexMult[k];
                  }
                  else{
                    lex4+=(o4[k])*LexMult[k];
                  }
                }
                //repeat - true if and only if o4 is already present in the tree
                bool repeat=false;
                for(long long int k=0;k<SearchTree.size();k++){
                  if(SearchTree[k]==lex4){
                    repeat=true;
                  }
                }
                if(!repeat){
                  //if o4 is not already in the tree:
                  //add o4 to the search Tree and put its rank into SearchTreeRanks
                  SearchTree.push_back (lex4);
                  SearchTreeRanks.push_back (rank4);
                  //add one to the number of unconsidered outcomes to account for the addition of o4 to the tree
                  unconsidered +=1;
                }
              }
            }
          }
        }
        //We have now assessed each improving i flip of outcome 3 for i not in a matching suffix
        //If it was not pruned by ranks or penalties and not already present in the tree, then this improving flip was added to the tree
      }
    }
    //We have now added to the tree all improving flips of outcome 3 that were not pruned by suffix fixing, rank pruning, or penalty pruning or already present
  }
  //if there are no more leaves in the tree to consider (and and o1 was not found), then we cannot search further and the query must be false
  //return that DQ is false and the size of the tree
  long long int Count=SearchTree.size();
  List R=List::create(false,Count);
  return R;
}

// [[Rcpp::export]]
List RPSDQPenaltyPriority(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2) {
  //Function takes a dominance query and answers it using all3 pruning methods with penalty priority
  //n - number of variables
  int n=A.nrow();
  //Omegais number of outcomes
  double Omega=1.0;
  //Multiply all domain sizes together
  //Note A.nrow is the number of variables
  for(int i=0;i<A.nrow();i++){
    Omega*=N[i];
  }
  //diff- Hamming distance between o1 and o2 as vectors
  int diff=0;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      diff+=1;
    }
  }
  if(diff==0){
    //If o1=o2 then return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //pen1=penalty(o1)
  //pen2=penalty(o2)
  long long int pen1=Penalty(A, N, ENTRIES, BREAK, o1);
  long long int pen2=Penalty(A, N, ENTRIES, BREAK, o2);
  //f - evaluation function
  long long int f=pen2-pen1-diff;
  if(f<0){
    //if f<0 initial penalty check fails
    //return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //rank1=r(o1)
  //rank2=r(o2)
  //RDiffs - vector of least improvement terms
  double rank1=Rank(A, N, ENTRIES, BREAK, o1);
  double rank2=Rank(A, N, ENTRIES, BREAK, o2);
  std::vector<double> RDiffs=RankDifferences(A, N, ENTRIES, BREAK);
  //comp - r(o1)-r(o2)-L_D(o1,o2) [least rank difference term]
  double comp=rank1-rank2;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      //for each i where o1 and o2 differ:
      comp-=RDiffs[i];
    }
  }
  double prec=-0.5/Omega;
  //if comp<0 our initial rank check fails.
  //However, as we are comparing the value of doubles and comp=0 does not fail the check, we use the prec term to account for possible rounding/precision errors.
  if(comp<prec){
    //If initial rank check fails then return false and 0 outcomes considered
    List R=List::create(false,0);
    return R;
  }
  //AncSizes - vector of #ancestors for each variable
  IntegerVector AncSizes(A.nrow());
  for(int i=0;i<A.nrow();i++){
    //For variable i, calculate the ancestor set (Anc) in the same way as in rank function
    IntegerVector Anc(A.nrow(),0);
    int AncCheck=0;
    int AncSum=0;
    for(int j=0;j<A.nrow();j++){
      if(A(j,i)==1){
        Anc[j]=1;
        AncSum+=1;
      }
    }
    while(AncCheck!=AncSum){
      AncCheck=AncSum;
      for(int j=0;j<A.nrow();j++){
        if(Anc[j]==1){
          for(int k=0;k<A.nrow();k++){
            if(A(k,j)==1){
              Anc[k]=1;
            }
          }
        }
      }
      AncSum=0;
      for(int j=0;j<A.nrow();j++){
        AncSum+=Anc[j];
      }
    }
    //add #ancestors of i to AncSizes vector
    AncSizes[i]=AncSum;
  }
  // TO - vector giving variables in topological order
  // calculated by putting variables in increasing order according to AncSizes
  std::vector<int> TO(A.nrow());
  int count=n;
  while(count>0){
    int MinAnc=n;
    int index;
    for(int i=0;i<A.nrow();i++){
      if(AncSizes[i]<MinAnc){
        MinAnc=AncSizes[i];
        index=i;
      }
    }
    TO[(n-count)]=index;
    AncSizes[index]=n;
    count-=1;
  }
  //Note - AncSizes is no longer accurate due to modification in above calculation
  //LexMult - lexicographic multipliers used to transform outcomes (as vectors) to/from their lexicographic enumeration
  std::vector<long long int> LexMult(A.nrow());
  LexMult[(n-1)]=1;
  for(int i=2;i<=n;i++){
    LexMult[(n-i)]=LexMult[(n-i+1)]*N[(n-i+1)];
  }
  //SearchTree - list of outcomes in the search tree (lexicographic enumeration) - currently just o2
  //SearchTreePenalties - For each outcome in SearchTree, this gives evaluation function (f) value of outcomes that havent been considered and -1 for those that have
  std::vector<long long int> SearchTree(1);
  std::vector<long long int> SearchTreePenalties(1);
  //lex - lexicographic enumeration of outcome o2
  long long int lex=0;
  for(int i=0;i<n;i++){
    if(i!=(n-1)){
      lex+=(o2[i]-1)*LexMult[i];
    }
    else{
      lex+=(o2[i])*LexMult[i];
    }
  }
  //Add o2 (and its f value) to search tree
  SearchTree[0]=lex;
  SearchTreePenalties[0]=f;
  //unconsidered - number of outcomes in the tree not yet considered
  long long int unconsidered=1;
  //While there are unconsidered outcomes in the tree (i.e. while it is possible to search further):
  while(unconsidered>0){
    //identify the outcome in the tree that is not yet considered and has minimal f value (call this outcome 3)
    //index - index of this outcome in SearchTree
    //Minf is the minimum f value of the unconsidered outcomes
    long long int Minf;
    long long int index;
    for(long long int i=0;i<SearchTree.size();i++){
      if(SearchTreePenalties[i]!=(-1)){
        Minf=SearchTreePenalties[i];
        index=i;
        break;
      }
    }
    for(long long int i=0;i<SearchTree.size();i++){
      if(SearchTreePenalties[i]!=(-1)){
        if(SearchTreePenalties[i]<Minf){
          Minf=SearchTreePenalties[i];
          index= i;
        }
      }
    }
    //f3 - f value for outcome 3
    long long int f3=Minf;
    //Change the SearchTreePenalties value to -1 as this outcome is now `considered'
    //take 1 from #unconsidered
    SearchTreePenalties[index]=-1;
    unconsidered -= 1;
    // lex3 is the lexicographic enumeration of  outcome3
    long long int lex3=SearchTree[index];
    //Convert lex3 into the vector form of outcome 3, o3
    IntegerVector o3(n);
    for(int i=0;i<n;i++){
      if((lex3 % LexMult[i])==0){
        o3[i]=lex3/LexMult[i];
        lex3-=(o3[i]-1)*LexMult[i];
      }
      else{
        o3[i]=(int)lex3/LexMult[i];
        o3[i]+=1;
        lex3-=(o3[i]-1)*LexMult[i];
      }
    }
    //SuffN is the position of the first element of the matching suffix between o1 and o3 (according to TO ordering)
    //if no matching suffix then SuffN=n
    int SuffN=n;
    for(int i=0;i<n;i++){
      if(o1[TO[(n-i-1)]]==o3[TO[(n-1-i)]]){
        SuffN=n-i-1;
      }
      else{
        break;
      }
    }
    for(int i=0;i<n;i++){
      //for each variable i
      //TOi is the position of i in the TO
      int TOi;
      for(int j=0;j<n;j++){
        if(TO[j]==i){
          TOi=j;
        }
      }
      if(TOi<SuffN){
        //if i is NOT in the matching suffix:
        //Parents is  0/1 vecor of the parents of variable i
        //NPa is the number of parents
        IntegerVector Parents(n,0);
        int NPa=0;
        for(int j=0;j<n;j++){
          if(A(j,i)==1){
            Parents[j]=1;
            NPa+=1;
          }
        }
        //use CPTRow to extract the preference rule from CPT(i) corresponding to parental assignment (for i) in o3, as we have done in previous functions
        IntegerVector pref;
        IntegerVector PaAsst(NPa);
        if(NPa==0){
          pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), 0);
        }
        else{
          int counter=0;
          for(int j=0;j<n;j++){
            if(Parents[j]==1){
              PaAsst[counter]=o3[j];
              counter+=1;
            }
          }
          pref=CPTRow(A, N, ENTRIES, BREAK, (i+1), PaAsst);
        }
        //pp - position of preference of variable i in o3 according to preference rule pref
        int pp=pref[o3[i]-1];
        if(pp>1){
          //if i can be improved (i.e. o3[i] is not the optimal assignment for i):
          //flips - vector of values that variable i could be flipped to and improve preference (according to preference rule pref)
          //obtain flips by finding which values have a lower preference position than pp in pref
          IntegerVector flips((pp-1));
          int counter=0;
          for(int j=0;j<N[i];j++){
            if(pref[j]<(pp)){
              flips[counter]=(j+1);
              counter+=1;
            }
          }
          //for each `better value' for i in flips
          for(int j=0;j<(pp-1);j++){
            //o4 is outcome 3 with variable i changed to take one of the better values in flips
            IntegerVector o4(n);
            for(int k=0;k<n;k++){
              o4[k]=o3[k];
            }
            o4[i]=flips[j];
            //diffs - hamming distance between o4 and o1
            int diffs=0;
            for(int k=0;k<n;k++){
              if(o1[k]!=o4[k]){
                diffs+=1;
              }
            }
            if(diffs==0){
              //if o4=o1 then we have reached o1 via an improving flip (o4) of outcome 3
              //return true and the number of outcomes in the search tree
              long long int Count=SearchTree.size();
              List R=List::create(true,Count);
              return R;
            }
            //pen4 - penalty of o4
            long long int pen4=Penalty(A, N, ENTRIES, BREAK, o4);
            //f - evaluation function value for o4
            long long int f=pen4-pen1-diffs;
            //If f<0, then we prune o4 from the tree according to the penalty pruning condition
            if(f>=0){
              //If o4 is not pruned by penalty pruning:
              //rank4 - rank of o4
              double rank4=Rank(A, N, ENTRIES, BREAK, o4);
              //comp = r(o1)-r(o4)-L_D(o1,o4) [least rand difference between o1 and o4]
              double comp=rank1-rank4;
              for(int k=0;k<n;k++){
                if(o1[k]!=o4[k]){
                  comp-=RDiffs[k];
                }
              }
              //If comp<0, then we prune o4 from the tree according to the rank pruning condition
              if(!(comp<prec)){
                //If o4 is not pruned by rank pruning:
                //lex4 - the lexicographic enumeration of o4
                long long int lex4=0;
                for(int k=0;k<n;k++){
                  if(k!=(n-1)){
                    lex4+=(o4[k]-1)*LexMult[k];
                  }
                  else{
                    lex4+=(o4[k])*LexMult[k];
                  }
                }
                //repeat - true if and only if o4 is already present in the tree
                bool repeat=false;
                for(long long int k=0;k<SearchTree.size();k++){
                  if(SearchTree[k]==lex4){
                    repeat=true;
                  }
                }
                if(!repeat){
                  //if o4 is not already in the tree
                  //add o4 to the search Tree and put its f value into SearchTreePenalties
                  SearchTree.push_back (lex4);
                  SearchTreePenalties.push_back (f);
                  //add one to the number of unconsidered outcomes to account for the addition of o4 to the tree
                  unconsidered +=1;
                }
              }
            }
          }
        }
        //We have now assessed each improving i flip of outcome 3 for i not in a matching suffix
        //If it was not pruned by ranks or penalties and not already present in the tree, then this improving flip was added to the tree
      }
    }
    //We have now added to the tree all improving flips of outcome 3 that were not pruned by suffix fixing, rank pruning, or penalty pruning or already present
  }
  //if there are no more leaves in the tree to consider (and and o1 was not found), then we cannot search further and the query must be false
  //return that DQ is false and the size of the tree
  long long int Count=SearchTree.size();
  List R=List::create(false,Count);
  return R;
}

// [[Rcpp::export]]
bool NumericalCheck(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2){
  //Function that takes a CP-net and Dominance query (is o1>o2?)
  // This function checks whether o1=o2  and whether the initial rank and penalty conditions hold
  // If any of these conditions hold false is returned - the DQ is FALSE
  // Otherwise true is returned - meaning the DQ MIGHT be true - needs answering
  //n - number of variables 
  int n=A.nrow();
  //Omega is number of outcomes
  double Omega=1.0;
  //Multiply all domain sizes together
  //Note A.nrow is the number of variables
  for(int i=0;i<A.nrow();i++){
    Omega*=N[i];
  }
  //diff- Hamming distance between o1 and o2 as vectors
  int diff=0;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      diff+=1;
    }
  }
  if(diff==0){
    //If o1=o2 then return false
    return false;
  }
  //pen1=penalty(o1)
  //pen2=penalty(o2)
  long long int pen1=Penalty(A, N, ENTRIES, BREAK, o1);
  long long int pen2=Penalty(A, N, ENTRIES, BREAK, o2);
  ////f - evaluation function
  long long int f=pen2-pen1-diff;
  if(f<0){
    //if f<0 initial penalty check fails
    //return false
    return false;
  }
  //rank1=r(o1)
  //rank2=r(o2)
  //RDiffs - vector of least improvement terms
  double rank1=Rank(A, N, ENTRIES, BREAK, o1);
  double rank2=Rank(A, N, ENTRIES, BREAK, o2);
  std::vector<double> RDiffs=RankDifferences(A, N, ENTRIES, BREAK);
  //comp - r(o1)-r(o2)-L_D(o1,o2) [least rank difference term]
  double comp=rank1-rank2;
  for(int i=0;i<n;i++){
    if(o1[i]!=o2[i]){
      //for each i where o1 and o2 differ:
      comp-=RDiffs[i];
    }
  }
  double prec=-0.5/Omega;
  //if comp<0 our initial rank check fails.
  //However, as we are comparing the value of doubles and comp=0 does not fail the check, we use the prec term to account for possible rounding/precision errors.
  if(comp<prec){
    //If initial rank check fails then return false
    return false;
  }
  //if none of the three conditions hold, return true
  return true;
}

// [[Rcpp::export]]
IntegerVector ImpVar(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2){
  //Function that takes a CP-net and Dominance query (is o1>o2?)
  //returns a 0/1 vector indicating the important variables to the query
  //If o1=o2, return a vector with all 2 entries
  //A.nrow - number of variables
  //D - 0/1 vector indicating where o1 and o2 are the same/differ
  IntegerVector D(A.nrow(),0);
  //HD - hamming distance between o1 and o2
  int HD=0;
  for(int i=0;i<A.nrow();i++){
    if(o1[i]!=o2[i]){
      D[i]=1;
      HD+=1;
    }
  }
  if(HD==0){
    //If o1=o2, return a vector of 2s
    IntegerVector DUM(A.nrow(),2);
    return DUM;
  }
  else{
    //if o1 and o2 are distinct:
    //ANC - vector indicating the variables that are ancestors to those variables in D (if the ith entry is non zero, variable i is an ancestor)
    IntegerVector ANC(A.nrow(),0);
    int ancSum=0;
    int ancCheck=0;
    //DES - vector indicating the variables that are descendants of those variables in D (if the ith entry is non zero, variable i is a descendant)
    IntegerVector DES(A.nrow(),0);
    int desSum=0;
    int desCheck=0;
    for(int i=0; i<A.nrow();i++){
      if(D[i]!=0){
        //for each variable, i, in D:
        //Add the parents of i to ANC
        //Add the children of i to DES
        for(int j=0; j<A.nrow();j++){
          ANC[j]+=A(j,i);
          DES[j]+=A(i,j);
        }
      }
    }
    //ancCheck - sum of indices with non-zero entries in ANC
    //desCheck - sum of indices with non-zero entries in DES
    for(int i=0; i<A.nrow();i++){
      if(ANC[i]!=0){
        ancCheck+=(i+1);
      }
      if(DES[i]!=0){
        desCheck+=(i+1);
      }
    }
    //until there is no change to which variables are in ANC
    while(ancSum!=ancCheck){
      ancSum=ancCheck;
      ancCheck=0;
      for(int i=0; i<A.nrow();i++){
        if(ANC[i]!=0){
          //for every variable, i, in ANC:
          //Add the parents of i to ANC
          for(int j=0; j<A.nrow();j++){
            ANC[j]+=A(j,i);
          }
        }
      }
      //calculate ancCheck
      for(int i=0; i<A.nrow();i++){
        if(ANC[i]!=0){
          ancCheck+=(i+1);
        }
      }
    }
    //until there is no change to which variables are in DEC
    while(desSum!=desCheck){
      desSum=desCheck;
      desCheck=0;
      for(int i=0; i<A.nrow();i++){
        if(DES[i]!=0){
          //for every variable, i, in DES:
          //Add the children of i to DES
          for(int j=0; j<A.nrow();j++){
            DES[j]+=A(i,j);
          }
        }
      }
      //calculate desCheck
      for(int i=0; i<A.nrow();i++){
        if(DES[i]!=0){
          desCheck+=(i+1);
        }
      }
    }
    // ANC and DES now have non-zero entries for exactly the ancestors/descendants of the variables in D
    //ImpVar - 0/1 vector indicating which variables are important to our query
    IntegerVector ImpVar(A.nrow(),0);
    for(int i=0; i<A.nrow();i++){
      //For each variable, i
      //If i is in D, then add it to important variables set
      if(D[i]!=0){
        ImpVar[i]=1;
      }
      else{
        //If i is not in D but it is both an ancestor and descendant of the variables in D, add it to important variables set
        if((ANC[i]!=0)&&(DES[i]!=0)){
          ImpVar[i]=1;
        }
      }
    }
    //By definition, ImpVar has a 1 entry for important variables to our query and 0 everywhere else
    return ImpVar;
  }
}

// [[Rcpp::export]]
int Degenerate(IntegerVector NRed, IntegerVector ENTRIES){
  //NRed - vector giving the domain sizes of Xs parents in order, then the domain size of X
  // ENTRIES - the CPT(X) section of the CP-net entries vector
  //Returns -1 if CPT(X) is non degenerate
  // Returns i (integer between 0 and |Pa| -1) giving the index of a degenerate parent (if CPT is degenerate)
  //Note that if we return that i is a degenerate parent, this does not mean that there aren't more degenerate parents
  //NPa - number of parents of X
  int NPa=(NRed.size()-1);
  //If X has no parents, then its CPT must be non degenerate, return -1
  if(NPa==0){
    return -1;
  }
  else{
    //LexMult - vector of lexicographic multipliers for obtaining the lexicographic enumeration of a parental assignment
    IntegerVector LexMult(NPa);
    int mult=1;
    for(int i=1; i<=NPa; i++){
      LexMult[(NPa-i)]=mult;
      if(i<NPa){
        mult*=NRed[(NPa-i)];
      }
    }
    if(NPa==1){
      //If there is one parent:
      //NRed[1] - size of Dom(X)
      IntegerVector pref(NRed[1],0);
      //pref - first |Dom(X)| elements of ENTRIES
      // this is the preference row of the CPT under parent=1
      for(int i=0; i<NRed[1]; i++){
        pref[i]=ENTRIES[i];
      }
      //for parent=2,3,4,...,|Dom(parent)|, similarly extract the CPT rows from ENTRIES
      for(int i=2; i<=NRed[0];i++){
        for(int j=0; j<NRed[1]; j++){
          //If the jth entry of the parent=i row is different from the jth entry of pref
        if(pref[j]!=ENTRIES[((i-1)*NRed[1]+j)]){
          R//Then the CPT contains at least two distinct rows and so the single parent is valid
          //Thus, the CP-net is non-degenerate, return -1
          return -1;
        }
      }
      }
      //If all of the parent=2,3,4,...,|Dom(parent)| rows were identical to pref, then all CPT rows are the same
      //Thus the only parent is degenerate, return its index (0)
      return 0;
    }
    else{//If there is >1 parent
      //TPa - product of all parent domain sizes - this is the total number of possible parent assignments
      int TPa=LexMult[0]*NRed[0];
      for(int i=0; i<NPa; i++){
        //for each parent, i, check degeneracy of CP-net with respect to i
        //ValidPA-indicates whether i is a valid (non-degenerate) parent
        bool ValidPa=false;
        //alpha - #value assignments of parents that come after i
        int alpha=LexMult[i];
        //PaDom - size of domain of parent i
        int PaDom=NRed[i];
        //beta - #value assignments of parents that come before i
        int beta=TPa/(alpha*PaDom);
        //loop over all possible assignments of values to the set of parents that come before i and the set of parents that come after
        //We cycle through these assignments in lexicographic order
        for(int j=0;j<beta;j++){
          for(int k=0;k<alpha;k++){
            //assume the parents that come before i are in the jth assignment, and those that come after are in their kth assignment
            //pref - the CPT row corresponding to this assignment and parent i =1
            IntegerVector pref(NRed[NPa]);
            for(int l=0;l<NRed[NPa];l++){
              pref[l]=ENTRIES[(j*(TPa/beta)+k)*NRed[NPa]+l];
            }
            //loop through the cases of parent i =2,3,4,..
            for(int l=1;l<NRed[i];l++){
              //Given parent i= (l+1)
              //extract the relevant preference row from CPT
              for(int m=0;m<NRed[NPa];m++){
                //If this row differs from pref on the mth entry:
                if(pref[m]!=ENTRIES[(j*(TPa/beta)+l*alpha+k)*NRed[NPa]+m]){
                  //Then there are 2 distinct preference rows where the parent assignment differs only on parent i
                  //This makes parent i valid, and we can move onto the next i so we break out of all loops until we get back to the i loop
                  ValidPa=true;
                  break;
                }
              }
              if(ValidPa){
                break;
              }
            }
            if(ValidPa){
                break;
              }
          }
          if(ValidPa){
                break;
              }
        }
        //If we cycled through all assignments to the set of parents that come before i and the set of parents that come after and didn't find i valid for any of them
        if(!ValidPa){
          //Then, given any assignment to the remaining parent, all parent i assignments give the same preference rule
          //By definition, parent i is degenerate, thus we return this index
          return i;
        }
        //If we found parent i to be valid, then we simply move onto the next parent
      }
      //If we exit the loop, then every parent must have been found valid
      //Thus, the CPT is non-degenerate, return -1
      return -1;
    }
  }
}

// [[Rcpp::export]]
List RemovePa(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, int PA, int CH){
  //PA and CH  - parent anc child in our cp-net - indexed from 0
  //We take the CP-net, remove this (degenerate) parent-child relation and return the simplified CP-net
  //StartCPTLength - number of entries in the CPT of the child in the input CP-net
  //N[CH] - domain size of child
  int StartCPTLength=N[CH];
  //NPa - number of parents of child (must be at least 1)
  int NPa=0;
  for(int i=0; i<A.nrow(); i++){
    if(A(i,CH)!=0){
      //For each parent of CH, multiply CPT length by parent domain size
      StartCPTLength*=N[i];
      NPa+=1;
    }
  }
  //EndCPTLength Length of CPT(CH) after we remove the parent
  int EndCPTLength=StartCPTLength/N[PA];
  //StartCPT - CPT of CH in input CP-net, extracted from entries vector
  IntegerVector StartCPT(StartCPTLength);
  for(int i=0;i<StartCPTLength;i++){
    StartCPT[i]=ENTRIES[BREAK[CH]-1+i];
  }
  //EndCPT will be the CPT of CH after the parent is removed
  IntegerVector EndCPT(EndCPTLength);
  if(NPa==1){
    //If PA is the only parent then all rows are the same as it is degenerate
    //Thus, the reduced CP-net can be obtained by taking any row, so we take the first (given by the first |Dom(CH)| entries of start CPT)
    for(int i=0;i<N[CH];i++){
      EndCPT[i]=StartCPT[i];
    }
  }
  else{
    //If CH has multiple parents:
    //PaDoms -vector of domain sizes of the parents of CH
    IntegerVector PaDoms(NPa);
    int counter=0;
    //PaNum - the position of PA within the set of parents of CH (ordered)
    int PANum=0;
    for(int i=0; i<A.nrow(); i++){
      if(A(i,CH)!=0){
        PaDoms[counter]=N[i];
        if(i==PA){
          PANum=counter;
        }
        counter+=1;
      }
    }
    //LexMult - vector of multipliers used to obtain the lexicographic enumeration of a parental assignment (to parents of CH)
    IntegerVector LexMult(NPa);
    int mult=1;
    for(int i=1; i<=NPa; i++){
      LexMult[(NPa-i)]=mult;
      if(i<NPa){
        mult*=PaDoms[(NPa-i)];
      }
    }
    //TPA - total number of possible assignments to the parents of CH
    int TPa=LexMult[0]*PaDoms[0];
    //alpha - total number of assignments to the parents that come after PA
    int alpha=LexMult[PANum];
    //beta - total number of assignments to the parents that come before PA
    int beta=TPa/(alpha*N[PA]);
    //Cycle through the possible assignments to the parents before PA and the assinments to the parents after PA (in both cases, moving through the assignments in lexicographic order)
    for(int i=0;i<beta;i++){
      for(int j=0;j<alpha;j++){
        //Given this assignments to the parents of CH other than PA, all assignments to PA will give the same preference order as PA is a degenerate parent
        //extract the CPT row from StartCPT that corresponds to this parent assignment with PA=1
        //enter this row in the appropriate location in EndCPT
        //(this is the row corresponding to the parents before PA taking assignment #i and those after taking assignment #j - as PA is not a parent here, this is a full parental assignment)
        for(int k=0;k<N[CH];k++){
          EndCPT[(i*alpha+j)*N[CH]+k]=StartCPT[(i*alpha*N[PA]+j)*N[CH]+k];
        }
      }
    }
    //This process fully specifies the new CPT EndCPT
  }
  //EndCPT is now the correct CPT for CH, with degen parent PA removed
  //Remove the parent-child edge from the structure
  A(PA,CH)=0;
  //replace StartCPT with EndCPT in the entries vector - this requires changing the length of the vector - New Entries
  IntegerVector NewEntries((ENTRIES.size()+EndCPTLength-StartCPTLength));
  //entries before CPT(CH) copied over
  for(int i=0; i<(BREAK[CH]-1);i++){
    NewEntries[i]=ENTRIES[i];
  }
  //Enter EndCPT
  for(int i=0;i<(EndCPTLength);i++){
    NewEntries[(BREAK[CH]-1)+i]=EndCPT[i];
  }
  //Copy the entries that came after CPT(CH)
  for(int i=0;i<(BREAK[A.nrow()]-BREAK[(CH+1)]);i++){
    NewEntries[(BREAK[CH]-1+EndCPTLength)+i]=ENTRIES[(BREAK[(CH+1)]-1)+i];
  }
  //Recalculate the breaks vector - NewBreaks
  IntegerVector NewBreak(BREAK.size());
  for(int i=0; i<BREAK.size();i++){
    //up to CH, all CPT break values are the same
    if(i<=CH){
      NewBreak[i]=BREAK[i];
    }
    else{
      //As CPT(CH) halved in size, all break points from this point must be reduced by this decrease amount
      NewBreak[i]=BREAK[i]+EndCPTLength-StartCPTLength;
    }
  }
  //L is the new CP-net, with the PA-CH relation removed from the structure and CPTs
  List L=List::create(A,NewEntries,NewBreak);
  return L;
}


// [[Rcpp::export]]
List UvRemove(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2){
  //Input dominance query
  //Function removes the unimportant variables of this query and then removes any degenerate edges
  //We assume that o1=/=o2
  //IV - 0/1 vector indicating the important variables
  IntegerVector IV(A.nrow());
  IV=ImpVar(A,N,ENTRIES,BREAK,o1,o2);
  // V - number of important variables
  int V=0;
  for(int i=0;i<A.nrow();i++){
    V+=IV[i];
  }
  //A2 - matrix made up of the important variable rows of A
  IntegerMatrix A2(V,A.nrow());
  //IVNum - the iindices of the important variables
  IntegerVector IvNum(V,0);
  int n=0;
  int m=0;
  for(int i=0;i<A.nrow();i++){
    if(IV[i]==1){
      IvNum[n]=i;
      n+=1;
      //For each important variable, i, add the ith row of A to A2
      for(int j=0;j<A.nrow();j++){
        A2(m,j)=A(i,j);
      }
      m+=1;
    }
  }
  // RedO1 - outcome 1 restricted to important variables
  // RedO2 - outcome 2 restricted to important variables
  int count=0;
  IntegerVector RedO1(IvNum.size());
  IntegerVector RedO2(IvNum.size());
  for(int i=0;i<A.nrow();i++){
    if(IV[i]==1){
      RedO1[count]=o1[i];
      RedO2[count]=o2[i];
      count+=1;
    }
  }
  //A3 - the important variable rows and columns of A only (i.e. A restricted to important variables)
  IntegerMatrix A3(V,V);
  for(int i=0;i<V;i++){
    //For each important variable i
    //next (ith) row of A3 is the ith A2 row, giving only the entries from important variable columns
    for(int j=0;j<V;j++){
      A3(i,j)=A2(i,IvNum[j]);
    }
  }
  //A3 is the new structure matrix
  //CptLen - vector of the lengths of the CPTs of the important variables (can extract from breaks vector)
  //NewCptLen is the sum of these lengths
  IntegerVector CptLen(V);
  int NewCptLen=0;
  //Break2 vector of break points for a CP-net with only the important variables
  //Calculable given the lengths of each important variable CPT
  IntegerVector Break2((V+1));
  Break2[0]=1;
  for(int i=0;i<V;i++){
    int var=IvNum[i];
    int len=(BREAK[(var+1)]-BREAK[var]);
    NewCptLen+=len;
    CptLen[i]=len;
    Break2[(i+1)]=Break2[i]+len;
  }
  //CPT2 - List of the CPTs of important variables in order  - extractable from the original Entries vector
  IntegerVector Cpt2(NewCptLen);
  int counter=0;
  for(int i=0; i<V;i++){
    //for the ith important variable, IvNum[i]:
    //Extract the CPT of IvNum[i] from Entries and add it to Cpt2
    for(int j=0; j<CptLen[i];j++){
      Cpt2[counter]=ENTRIES[(BREAK[IvNum[i]]-1+j)];
      counter+=1;
    }
  }
  //Cpt2 and Break 2 are now the CPT entries for the CP-net with unimportant variables removed
  //A3 is the corresponding structure
  // It remains to reduce each CPT of variables that have lost a parent in the unimportant variable removal (thenn we will have a well defined reduced CP-net)
  //Cpt3 will be the CPT entries after CPTs are reduced (plus some 0s at the end)
  IntegerVector Cpt3(NewCptLen,0);
  //Cpt3Count - the next element of cpt3 to enter
  int Cpt3Count=0;
  //For each important variable
  for(int i=0; i<V; i++){
    //Let X be the ith important variable
    //nPa - #parents of X in original structure
    int nPa=0;
    for(int j=0; j<A.nrow(); j++){
      nPa+=A(j,IvNum[i]);
    }
    if(nPa==0){
      //If X had no parents originally, it cannot have lost any in the variable removal so its CPT does not need reducing
      //Enter the CPT of X as it is in Cpt2
      for(int j=0; j<CptLen[i]; j++){
        Cpt3[(Cpt3Count+j)]=Cpt2[(Break2[i]-1+j)];
      }
      //update Cpt3Count to the next unfilled element
      Cpt3Count+=CptLen[i];
    }
    else{
      //RemovedPa is a 0/1 vector. If an original parent of X is important, it has a 1 entry, otherwise it has 0 entry
      IntegerVector RemovedPa(nPa,0);
      //PaDomSize - domain sizes of the origina parent set of X
      IntegerVector PaDomSize(nPa,0);
      int counter=0;
      for(int j=0; j<A.nrow(); j++){
        if(A(j,IvNum[i])==1){
          //For each parent of X in A, j
          //Add domain size to PaDomSize
          PaDomSize[counter]=N[j];
          //Check if j is an important variable (b is result)
          bool b=false;
          for(int k=0;k<V;k++){
            if(IvNum[k]==j){
              b=true;
            }
          }
          if(b){
            //if j is important, enter a 1 into the RemovedPa vector, otherwise leave it as 0
            RemovedPa[counter]=1;
          }
          counter+=1;
        }
      }
      //nPaAsst - number of possible assignments to the remaining (important) parents of X - obtained by multiplying their domain sizes
      int nPaAsst=1;
      for(int j=0; j<nPa; j++){
        if(RemovedPa[j]==1){
          nPaAsst*=PaDomSize[j];
        }
      }
      if(nPaAsst==1){
        //If all parents of X were removed:
        //SetPa - the values assigned to the parents of X in both o1 and o2 (must be the same in both as unimportant)
        IntegerVector SetPa(nPa);
        int counter=0;
        for(int j=0; j<A.nrow(); j++){
          if(A(j,IvNum[i])==1){
            //For each parent of X, j:
            //Add o1 assignment to j to SetPa vector
            SetPa[counter]=o1[j];
            counter+=1;
          }
        }
        //Cpt - the (only) row of the original CPT(X) corresponding to the SetPa parental assignment
        IntegerVector Cpt(N[IvNum[i]]);
        Cpt=CPTRow(A, N, ENTRIES, BREAK, (IvNum[i]+1), SetPa);
        //Cpt is the reduced CPT for X by our preprocessing procedure
        //Enter Cpt into Cpt3 as the new CPT(X)
        for(int j=0; j<N[IvNum[i]]; j++){
          Cpt3[(Cpt3Count+j)]=Cpt[j];
        }
        //update Cpt3Count
        Cpt3Count+=N[IvNum[i]];
      }
      //CptLen[i]/N[IvNum[i]] - CPT(X) length/size of Dom(X) - this is the number of parental assignments in original CPT
      if(nPaAsst==(CptLen[i]/N[IvNum[i]])){
        //If no parents of X have been removed, then no reduction is required
        //Enter the CPT(X) into Cpt3 exactly as it is in Cpt2:
        for(int j=0; j<CptLen[i]; j++){
          Cpt3[(Cpt3Count+j)]=Cpt2[(Break2[i]-1+j)];
        }
        //Update Cpt3Count
        Cpt3Count+=CptLen[i];
      }
      if((nPaAsst>1)&&(nPaAsst<(CptLen[i]/N[IvNum[i]]))){
        //If X lost some, but not all parents:
        //PaAsst will list every possible assignment to the original parents of X (viewed as a vector) where the removed variables take the values they have in o1 (and o2)
        //assingments listed end to end in a vector
        IntegerVector PaAsst((nPaAsst*nPa));
        //First, go through and assign the removed parents the values they take in o1 in every parental assignment vector (within PaAsst)
        for(int j=0;j<(nPaAsst*nPa);j+=nPa){
          int counter=0;
          for(int k=0;k<A.nrow();k++){
            if(A(k,IvNum[i])==1){
              if(RemovedPa[counter]==0){
                PaAsst[(j+counter)]=o1[k];
              }
              counter+=1;
            }
          }
        }
        //nPaLeft - number of parents of X in reduced structure (must be positive)
        int nPaLeft=0;
        for(int j=0;j<nPa;j++){
          if(RemovedPa[j]==1){
            nPaLeft+=1;
          }
        }
        //PaLeftDomSize - vector of domain sizes of the remaining parents of X
        IntegerVector PaLeftDomSize(nPaLeft);
        int counter=0;
        for(int j=0;j<nPa;j++){
          if(RemovedPa[j]==1){
            PaLeftDomSize[counter]=PaDomSize[j];
            counter+=1;
          }
        }
        //Pa - gives the indices of the parents of X in original structure
        IntegerVector Pa(nPa);
        counter=0;
        for(int j=0;j<A.nrow();j++){
          if(A(j,IvNum[i])==1){
            Pa[counter]=j;
            counter+=1;
          }
        }
        //For each original parent of X, say Y
        for(int j=0;j<nPa;j++){
          //If this parent is still present after unimportant variable removal
          if(RemovedPa[j]==1){
            //beta is the lexicographic multiplier for Y (viewed as part of a parent assignment of X in the reduced structure)
            int beta=1;
            for(int k=(j+1);k<nPa;k++){
              if(RemovedPa[k]==1){
                beta*=PaDomSize[k];
              }
            }
            //Use the lexicographic enumeration of assignments in PaAsst to assign the Y values to the parental assignments in PaAsst
            for(int x=0;x<N[Pa[j]];x++){
              for(int y=x;y<(nPaAsst/beta);y+=N[Pa[j]]){
                if(beta==1){
                  PaAsst[(y*nPa + j)]=(x+1);
                }
                else{
                  for(int z=0;z<beta;z++){
                    PaAsst[(y*beta*nPa + z*nPa + j)]=(x+1);
                  }
                }
              }
            }
          }
        }
        //PaAsst is now a complete list of the assignments to parents of X such that removed variables take their o1 values
        //These are listed in the lexicographic order of the assignments to the remaining variables of X (as the others are fixed in each assignment)
        //For each assignment in PaAsst
        for(int j=0;j<(nPaAsst*nPa);j+=nPa){
          //PaSet is the next parent assignment
          IntegerVector PaSet(nPa);
          for(int k=0;k<nPa;k++){
            PaSet[k]=PaAsst[(j+k)];
          }
          //CptRow is the CPT(X) row corresponding to this preference assignment in the original CPT(X)
          IntegerVector CptRow(N[IvNum[i]]);
          CptRow=CPTRow(A, N, ENTRIES, BREAK, (IvNum[i]+1), PaSet);
          //Add CptRow to Cpt3
          for(int k=0;k<N[IvNum[i]];k++){
            Cpt3[(Cpt3Count +k)]=CptRow[k];
          }
          //Update Cpt3Count
          Cpt3Count+=N[IvNum[i]];
          //As we are adding the CPT(X) rows where unimportant parent variables take their o1 values, and we are adding them in the lexicographic order of the remaining variables
          //By our preprocessing procedure, this is correctly adding the reduced CPT(X) to Cpt3
        }
      }
    }
  }
  //Cpt3 is now a list of the reduced CPTs of the reduced structure (possibly with some 0s at the end)
  //RedA- the new, reduced structure
  IntegerMatrix RedA=A3;
  //RedN - the Domain sizes of the remaining variables
  IntegerVector RedN(IvNum.size());
  for(int i=0; i<RedA.nrow();i++){
    RedN[i]=N[IvNum[i]];
  }
  //RedBreaks - the breaks vector of the reduced structure (calculable from structure and domain sizes)
  //These breaks will match up with the lengths of the reduced CPTs we obtained in making Cpt3
  IntegerVector RedBreaks((RedA.nrow()+1));
  RedBreaks[0]=1;
  for(int i=0;i<RedA.nrow();i++){
    int CptSize=RedN[i];
    for(int j=0;j<RedA.nrow();j++){
      if(RedA(j,i)!=0){
        CptSize*=RedN[j];
      }
    }
    RedBreaks[(i+1)]=RedBreaks[i]+CptSize;
  }
  //RedEntries - the list of reduced CPts (from Cpt3) with no 0s at the end.
  //We can use breaks value to tell us how long this vector should be and simply extract these first n entries from Cpt3
  IntegerVector RedEntries((RedBreaks[(RedA.nrow())]-1));
  for(int i=0;i<(RedBreaks[(RedA.nrow())]-1);i++){
    RedEntries[i]=Cpt3[i];
  }
  //The Red... items make up the well defined CP-net obtained by removing the unimportant variables (and reducing CPTs accordingly) from the imput CP-net
  //UV - a vector giving the indices of the unimportant variables, which we removed
  IntegerVector UV(A.nrow()-RedA.nrow());
  counter=0;
  for(int i=0;i<A.nrow();i++){
    if(IV[i]==0){
      UV[counter]=i;
      counter+=1;
    }
  }
  //LostParents - for each variable in the reduced structure this entry is the #remaining parents if the variable lost some but not all parents in the reduction
  //otherwise the entry is 0
  //The non-zero entries are the variables for which we need to check degeneracy of remaining parents
  IntegerVector LostParents(RedA.nrow(),0);
  for(int i=0;i<RedA.nrow();i++){
    //for each remaining variable
    //remainingParents - #parents of i in the new structure, redA
    int remainingParents=0;
    for(int k=0;k<RedA.nrow();k++){
      remainingParents+=RedA(k,i);
    }
    if(remainingParents>0){
      //if the variable i has parents in the reduced structure
      //x - index of the variable in original structure, i
      int x=IvNum[i];
      //For each removed variable
      for(int j=0; j<UV.size();j++){
        if(A(UV[j],x)==1){
          //If one of the removed variables was a parent of i in A:
          //Add the number of parents of i in redA to LostParents
          LostParents[i]=remainingParents;
          break;
        }
      }
    }
  }
  
  //We now check for degeneracy of remaining parents where appropriate (LostParents non-zero variables)
  //If we find a degenerate parent, we remove it
  //For each variable in the reduced structure:
  for(int i=0;i<RedA.nrow();i++){
    //If we need to check the degeneracy of the parents of i:
    if(LostParents[i]>0){
      //RemainingParents - #parents of i in redA
      int RemainingParents=LostParents[i];
      //NewPa - indices pf  parents of i in redA
      //NewPaDom - domain sizes of the parents, with domain of variable i at the end
      IntegerVector NewPa(LostParents[i]);
      IntegerVector NewPaDom(LostParents[i]+1);
      counter=0;
      for(int j=0;j<RedA.nrow();j++){
        if(RedA(j,i)==1){
          NewPa[counter]=j;
          NewPaDom[counter]=RedN[j];
          counter+=1;
        }
      }
      NewPaDom[LostParents[i]]=RedN[i];
      //Extract the CPT of i (in reduced CP-net) from RedEntries - this is CPTi
      IntegerVector CPTi(RedBreaks[(i+1)]-RedBreaks[i]);
      for(int j=0;j<(RedBreaks[(i+1)]-RedBreaks[i]);j++){
        CPTi[j]=RedEntries[RedBreaks[i]-1+j];
      }
      //Evaluate the degeneracy of CPTi
      int D=Degenerate(NewPaDom, CPTi);
      //If we find the CPT contains degenerate parents, repeat:
      while(D>=0){
        //D is an invalid parent of i
        //L - the CP-net obtained from RedA by removing degenerate edge D-> i
        List L=RemovePa(RedA, RedN, RedEntries, RedBreaks, NewPa[D], i);
        //Update the Red... items to be the CP-net given by L
        IntegerMatrix RedA=L[0];
        RedEntries=L[1];
        RedBreaks=L[2];
        //i now has one less parent
        RemainingParents-=1;
        //Obtain the new CPT of i from RedEntries (which has been updated) - this is NewCPTi
        IntegerVector NewCPTi(RedBreaks[(i+1)]-RedBreaks[i]);
        for(int j=0;j<(RedBreaks[(i+1)]-RedBreaks[i]);j++){
          NewCPTi[j]=RedEntries[RedBreaks[i]-1+j];
        }
        //RedNewPa - indices of the parents of i in the updated RedA
        //RedNewPaDom - domain sizes of these parents, with the domain of i at the end
        IntegerVector RedNewPa(RemainingParents);
        IntegerVector RedNewPaDom(RemainingParents+1);
        counter=0;
        for(int j=0;j<RedA.nrow();j++){
          if(RedA(j,i)==1){
            RedNewPa[counter]=j;
            RedNewPaDom[counter]=RedN[j];
            counter+=1;
          }
        }
        RedNewPaDom[RemainingParents]=RedN[i];
        //pass these vectors to NewPaDom and NewPa so the next iteration of the while loop works correctly
        NewPaDom=RedNewPaDom;
        NewPa=RedNewPa;
        //Assess the degeneracy of the new CPTi
        D=Degenerate(NewPaDom, NewCPTi);
        //If the CPT is still degenerate, we repeat the while loop, removing the newly identified degenerate parent of i from the structure
      }
      //the CPT of i must now be degenerate, having all degenerate parents been removed
    }
  }
  //We have now removed all degenerate parents (of all variables) from the reduced structure
  //The resulting CP-net is summarised by T
  List T=List::create(RedA,RedN,RedEntries,RedBreaks,RedO1,RedO2);
  //Return the CP-net with all unimportant variables (and subsequent degenerate edges) removed
  return T;
}

// [[Rcpp::export]]
List ConnectedComponents(IntegerMatrix A){
  //input adjacency matrix A (assumed acyclic)
  //returns the number of connected components
  //If there are k>1 connected components, also outputs a matrix with k rows. The ith row is a 0/1 vector indicating the variables in the ith connected component
  //n - number of nodes/variables in the graph
  int n=A.nrow();
  if(n==1){
    //If there is only 1 variable, it is a signle connected component, return 1
    int dummy=1;
    List L=List::create(dummy);
    return L;
  }
  //Opt - 0/1 vector indicating the variables which have 0 parents in A
  IntegerVector Opt(n,0);
  for(int i=0; i<n; i++){
    //for each variabl i:
    //parents - number of parents of i in A
    int parents=0;
    for(int j=0; j<n; j++){
      parents+=A(j,i);
    }
    if(parents==0){
      //if i has no parents, put a 1 entry in Opt
      Opt[i]=1;
    }
  }
  //NOpt - number of variables with no parents
  int NOpt=0;
  for(int i=0; i<n; i++){
    NOpt+=Opt[i];
  }
  if(NOpt==1){
    //If there is only one parent with 0 parents, the whole graph must be connected (due to acyclicity)
    //Thus, only 1 connected component, return 1
    int dummy=1;
    List L=List::create(dummy);
    return L;
  }
  // If there is more than one variable without parents:
  //DescendantSets - matrix whose ith row gives a 0/1 vector indicating the descendants of the ith variable, x, with 0 parents and a 1 entry x
  IntegerMatrix DescendantSets(NOpt,n);
  //OptVar - indices of the variables with 0 parents
  IntegerVector OptVar(NOpt);
  int counter=0;
  for(int i=0; i<n; i++){
    if(Opt[i]==1){
      //if variable i has 0 parents:
      //enter index i in OptVar
      OptVar[counter]=i;
      // put a 1 in the next row of DescendantSets for each child of i and for i itself
      for(int j=0;j<n;j++){
        if(i==j){
          DescendantSets(counter,j)=1;
        }
        else{
          DescendantSets(counter,j)=A(i,j);
        }
      }
      counter+=1;
    }
  }
  for(int i=0;i<NOpt;i++){
    // For each variable with 0 parents:
    //Add the children of the variabls in row i of DescendantSets to the row until no change happens - this will now be the descendent set
    bool Same=false;
    while(!Same){
      //Start - row i of DescendantSets
      IntegerVector Start(n);
      for(int j=0;j<n;j++){
        Start[j]=DescendantSets(i,j);
      }
      //For each variable in Start (with entry 1), add any children of that variable to Start
      for(int j=0;j<n;j++){
        if(Start[j]==1){
          for(int k=0;k<n;k++){
            if(DescendantSets(i,k)==0){
              DescendantSets(i,k)+=A(j,k);
            }
          }
        }
      }
      //End - row i of DescendantSets
      IntegerVector End(n);
      for(int j=0;j<n;j++){
        End[j]=DescendantSets(i,j);
      }
      //Diff- hamming distance between end and start
      int Diff=0;
      for(int j=0;j<n;j++){
        if(Start[j]!=End[j]){
          Diff+=1;
        }
      }
      //If we have made no changes, row i is now the descendant set (plus x) of the ith variable with no parents, x
      if(Diff==0){
        Same=true;
      }
    }
  }
  //We repeatedly assess whether any 2 rows of DescendantSets have a common 1 entry. 
  //If so, merge such a pair of rows by replacing them with a row that has a 1 in every position that had a 1 in either of the two rows
  //This continues until no two rows have a common 1 entry
  //The idea here is that, say x and y have 0 parents. If (x and its descendents) and (y and its descendents) have an overlap, then all of these variables are in 1 connected component
  //We thus group them together and continue looking for overlaps. Once the sets have no more overlaps, they must all be in distinct connected components
  bool Crossover=true;
  int counter2=0;
  while(Crossover){
    //Detect - are there two overlapping rows
    bool Detect=false;
    //OVerlap - 2 such overlapping rows
    IntegerVector Overlap(2);
    //For each variable:
    for(int i=0;i<n;i++){
      //OptAnc - number of rows in DescendantSets that have a 1 in position i
      int OptAnc=0;
      for(int j=0;j<NOpt;j++){
        OptAnc+=DescendantSets(j,i);
      }
      if(OptAnc>1){
        //If more than row of DescendantSets contains a 1 at position i
        //Overlap - (the indices of) two such rows
        counter=0;
        for(int j=0;j<NOpt;j++){
          if(counter<2){
            if(DescendantSets(j,i)==1){
              Overlap[counter]=j;
              counter+=1;
            }
          }
        }
        Detect=true;
        break;
      }
    }
    if(Detect){
      //If we found some case of two rows of DescendantSets having variables (those with entry 1) in common
      if(DescendantSets.nrow()==2){
        //If there are only 2 rows in DescendantSets, they must have a value in common and so the whole structure is connected
        //return 1
        int dummy=1;
        List L=List::create(dummy);
        return L;
      }
      //We now merge the two rows that have an element in common
      //Reduce the number of DescendantSets rows by 1
      NOpt-=1;
      //NewDescendantSets - matrix with the same rows as DescendantSets except rows a and b (where Overlap=(a,b))
      //These two rows are replaced by 1 in NewDescendantSets - this row has an entry 1 if either row a or b had a 1
      IntegerMatrix NewDescendantSets(NOpt,n);
      //Add rows that came before row a
      for(int i=0; i<Overlap[0];i++){
        for(int j=0; j<n;j++){
          NewDescendantSets(i,j)=DescendantSets(i,j);
        }
      }
      //Add the merged row a and b row
      for(int i=0; i<n;i++){
        if((DescendantSets((Overlap[0]),i)+DescendantSets((Overlap[1]),i))==0){
          NewDescendantSets(Overlap[0],i)=0;
        }
        else{
          NewDescendantSets(Overlap[0],i)=1;
        }
      }
      //Add the rows between rows a and b
      for(int i=(Overlap[0]+1); i<Overlap[1];i++){
        for(int j=0; j<n;j++){
          NewDescendantSets(i,j)=DescendantSets(i,j);
        }
      }
      //Ass the rows that came after row b
      for(int i=(Overlap[1]+1);i<(NOpt+1);i++){
        for(int j=0; j<n;j++){
          NewDescendantSets((i-1),j)=DescendantSets(i,j);
        }
      }
      //replace DescendantSets with NewDescendantSets
      DescendantSets=NewDescendantSets;
    }
    else{
      //If no two rows in DescendantSets have a common 1 entry:
      Crossover=false;
    }
    counter2+=1;
  }
  //If the function has not yet exited, then DescendantSets now has >=2 rows, meaning we have multiple connected components
  //The rows of DescendantSets consititue 0/1  vectors of the connected components
  //return the number of connected components (#rows) and the matrix
  int rows=DescendantSets.nrow();
  List L=List::create(rows,DescendantSets);
  return L;
}

// [[Rcpp::export]]
List UVRSAndRPSDQRankPriority(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2){
  //This function takes a CP-net and Dominance query (o1>o2?)
  //Apply UVRS to the CP-net with Numerical Checks, as described in thesis experiment section
  //The resulting reduced query is then answered using rank pruning, penalty pruning and suffix fixing with rank priority
  //If there are multiple reduced queries, they are answered in order on increasing size (number of variables)
  //return the answer of the dominance query and the number of outcomes traversed and the following timestamps
  //start and end of function if the query is answered via preprocessing
  //start, end of preprocessing, end of function otherwise
  //#Outcomes in UVRS reduced CP-net if preprocessing does not answer the query
  //Mark the time of commencing the function
  Timer time;
  time.step("Start");
  //Apply numerical checks to original input query
  bool Check=NumericalCheck(A, N, ENTRIES, BREAK, o1, o2);
  if(!Check){
    // If Numerical Checks are failed then query is false
    //return false, 0 outcomes considered and start/end times
    //log end time 
    time.step("End");
    NumericVector TimeJumps(time);
    List R=List::create(false,0,TimeJumps);
    return R;
  }
  //Reduced - CP-net obtained by removing unimportant variables
  List Reduced=UvRemove(A,N,ENTRIES,BREAK,o1,o2);
  //Apply numerical checks to reduced CP-net
  Check=NumericalCheck(Reduced[0],Reduced[1],Reduced[2],Reduced[3],Reduced[4],Reduced[5]);
  if(!Check){
    // If Numerical Checks after initial reduction are failed then query is false
    //return false, 0 outcomes considered and start/end times
    //log end time 
    time.step("End");
    NumericVector TimeJumps(time);
    List R=List::create(false,0,TimeJumps);
    return R;
  }
  //Continue to remove unimportant variables until an iteration does not alter the CP-net structure
  //Same - are the structures before and after removal the same?
  bool Same=false;
  while(!Same){
    //NextIter - CP-net obtained by removing unimportant variables from reduced structure
    List NextIter=UvRemove(Reduced[0],Reduced[1],Reduced[2],Reduced[3],Reduced[4],Reduced[5]);
    //Original - structure before removal
    IntegerMatrix Original=Reduced[0];
    //New - structure after
    IntegerMatrix New=NextIter[0];
    if(Original.nrow()==New.nrow()){
      //If structure unchanged by removal:
      Same=true;
    }
    else{
      //If structure changed by removal:
      //replace Reduced by the further reduced CP-net, NextIter
      Reduced=NextIter;
      //Evaluate numerical checks for the newly reduced CP-net
      Check=NumericalCheck(Reduced[0],Reduced[1],Reduced[2],Reduced[3],Reduced[4],Reduced[5]);
      if(!Check){
        // If Numerical Checks after succesive reduction fail then original query was false
        //return false, 0 outcomes considered and start/end times
        //log end time 
        time.step("End");
        NumericVector TimeJumps(time);
        List R=List::create(false,0,TimeJumps);
        return R;
      }
    }
  }
  //Reduced is the now the CP-net and query after iterative Uimplortant Variable (UV) removal
  //(and no further reduction can be done this way as last attempt did not affect structure)
  //Conn - connected components in the reduced CP-net structure
  List Conn=ConnectedComponents(Reduced[0]);
  //SubCPN - the reduced CP-net split into its connected component
  List SubCPN;
  //nAfter - number of outcomes after we separate the CP-net into connected components
  long long int nAfter=1;
  if(Conn[0]==1){
    //if there is only 1 connected component, we don't need to split up the CP-net
    //SubCPN is a list of 1 - the reduced CPN and DQ
    SubCPN=List::create(Reduced);
    //SmallN - domain size vector for Reduced
    IntegerVector SmallN=Reduced[1];
    //nAfter - #outcomes in Reduced
    for(int i=0;i<SmallN.size();i++){
      nAfter*=SmallN[i];
    }
  }
  else{
    //For each connected component, restrict the reduced CP-net and query to this component and add this sub-CP-net to our list SubCPN
    for(int k=0;k<Conn[0];k++){
      // For each connected component
      //B - structure of the reduced CP-net
      IntegerMatrix B=Reduced[0];
      //CCM - matrix of connected component entries
      IntegerMatrix CCM=Conn[1];
      //CC - 0/1 vector of kth connected component (row k of CCM)
      IntegerVector CC(B.nrow(),0);
      //NVar - number of variables in the kth connected component
      int NVar=0;
      for(int i=0;i<B.nrow();i++){
        CC[i]+=CCM(k,i);
        NVar+=CCM(k,i);
      }
      //SubA - B reduced to variables in kth connected component (Adjacency matrix of the kth connected component)
      IntegerMatrix SubA(NVar,NVar);
      int counter1=0;
      int counter2=0;
      for(int i=0;i<B.nrow();i++){
        if(CC[i]==1){
          //For each row of B corresponding to a variable in kth connected component
          for(int j=0;j<B.nrow();j++){
            if(CC[j]==1){
              //For each column of B corresponding to a variable in kth connected component
              //Add this entry to SubA
              SubA(counter1,counter2)=B(i,j);
              counter2+=1;
            }
          }
          counter2=0;
          counter1+=1;
        }
      }
      //SubN - domain sizes vector for the kth connected component variables
      IntegerVector SubN(NVar);
      //Subo1 - outcome 1 restricted to the kth connected component variables
      IntegerVector Subo1(NVar);
      //Subo2 - outcome 2 restricted to the kth connected component variable
      IntegerVector Subo2(NVar);
      //Obtain these by taking the relevant element of Reduced and restricting to the kth connected component variables
      int counter=0;
      for(int i=0;i<B.nrow();i++){
        if(CC[i]==1){
          IntegerVector RedN=Reduced[1];
          SubN[counter]=RedN[i];
          IntegerVector Redo1=Reduced[4];
          Subo1[counter]=Redo1[i];
          IntegerVector Redo2=Reduced[5];
          Subo2[counter]=Redo2[i];
          counter+=1;
        }
      }
      //Afterni- number of variable assingments to the variables in the kth connected components (multiply ther domain sizes together)
      long long int Afterni=1;
      for(int i=0;i<SubN.size();i++){
        Afterni*=SubN[i];
      }
      //Add Afterni to nAfter as nAfter is the sum of the #outcomes in each connected component in the case of >1 connected component
      nAfter+=Afterni;
      //Calculate the breaks vector for the CP-net restricted to the kth connected component variables
      IntegerVector SubBreaks((NVar+1));
      SubBreaks[0]=1;
      counter=1;
      //Use breaks vector from Reduced to calculate the CPT lengths of the kth connected component variables
      IntegerVector RedBreaks=Reduced[3];
      IntegerVector RedEntries=Reduced[2];
      for(int i=0;i<B.nrow();i++){
        if(CC[i]==1){
          //For each variable in kth component, add the CPT length to obtain the next breaks value
          SubBreaks[counter]=SubBreaks[(counter-1)]+RedBreaks[(i+1)]-RedBreaks[(i)];
          counter+=1;
        }
      }
      //SubEntries - entries vector for the CP-net reduced to the kth connected component variables
      //Obtained by copying the CPTs for the kth connected component variables over from the entries vector for Reduced
      IntegerVector SubEntries(SubBreaks[NVar]-1);
      counter=0;
      for(int i=0;i<B.nrow();i++){
        if(CC[i]==1){
          //For each of the kth connected component variables:
          //Extract CPT length from Reduced breaks vector - CPTSize
          int CPTSize=RedBreaks[(i+1)]-RedBreaks[(i)];
          //Copy the CPT entries from the Reduced entries vector into SubEntries
          for(int j=0;j<CPTSize;j++){
            SubEntries[counter]=RedEntries[(RedBreaks[(i)]-1+j)];
            counter+=1;
          }
        }
      }
      //Sub... constitute the CP-net and dominance query in reduced restricted to the kth connected component
      //Apply numerical checks to this Sub... dominance query
      bool Check=NumericalCheck(SubA,SubN,SubEntries,SubBreaks,Subo1,Subo2);
      if(!Check){
        //If the checks fail on this sub-CP-net then the original query is also false
        //return false, 0 outcomes considered and start/end times
        //log end time 
        time.step("End");
        NumericVector TimeJumps(time);
        List R=List::create(false,0,TimeJumps);
        return R;
      }
      //SubRed - the Sub CP-net (and query) we obtained by restructing to the kth connected component
      List SubRed=List::create(SubA,SubN,SubEntries,SubBreaks,Subo1,Subo2);
      //Add SubRed to out list of sub CP-nets
      SubCPN.push_back (SubRed);
    }
  }
  //SubCPN - list of reduced dominance queries we must now answer
  //log the end time of preprocessing (UVRS is now complete)
  time.step("DonePreProc");
  //Order - indeces 0 to (#connected components -1) in order of increasing size (number of variables) of the associated connected component
  IntegerVector Order(SubCPN.size());
  IntegerVector Record(SubCPN.size(),1);
  int counter=0;
  while(counter<SubCPN.size()){
    //Find the smallest component not yet listed and add index to Order
    int MinSize=(A.nrow())+1;
    int index;
    for(int i=0;i<SubCPN.size();i++){
      //For each connected component:
      List DQi=SubCPN[i];
      //Ai - #variables in this component/CP-net
      IntegerMatrix Ai=DQi[0];
      //Record[i]==1 if and only if this connected component not yet added to Order
      if((Ai.nrow()<MinSize)&(Record[i]==1)){
        MinSize=Ai.nrow();
        index=i;
      }
    }
    //Add index of minimal remaining component to Order
    Order[counter]=index;
    //Change corresponding Record entry to 0 so that we know it has been added
    Record[index]=0;
    counter+=1;
  }
  //Order is the order in which the reduced queries should be answered - smallest first
  //OutCons - counts the number of outcomes we consider in answering these reduced queries
  long long int OutCons=0;
  for(int i=0;i<SubCPN.size();i++){
    //for each reduced CP-net and query:
    //DQi - ith reduced query (in order of size)
    List DQi=SubCPN[Order[i]];
    //Answer the dominance query
    List Answeri=RPSDQRankPriority(DQi[0],DQi[1],DQi[2],DQi[3],DQi[4],DQi[5]);
    //Add #outcomes considered to the count
    OutCons+=Answeri[1];
    if(!Answeri[0]){
      //If the result of the query is false, then the original (input) query is also false
      //return false, #outcomes considered, start/end preprocessing/end times, and the number of outcomes in the UVRS reduced CP-net
      //log end time 
      time.step("End");
      NumericVector TimeJumps(time);
      List R=List::create(false,OutCons,TimeJumps,nAfter);
      return R;
    }
  }
  //If all reduced querys were answered and found true, the original query must also be true
  //return true, #outcomes considered, start/end preprocessing/end times, and the number of outcomes in the UVRS reduced CP-net
  //log end time 
  time.step("End");
  NumericVector TimeJumps(time);
  List R=List::create(true,OutCons,TimeJumps,nAfter);
  return R;
}

// [[Rcpp::export]]
List ForwardPruning(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2){
  //Input CP-net and dominance query (o1>o2?)
  //Assume that o1=/=o2
  //Function applies forward pruning to this query as described in my thesis
  //If forward pruning answers the query by reducing a variable domain to 0 - output variable of interest and pruned domains matrix
  //Otherwise, the pruned domains matrix and reduced CP-net and query are output
  //n-number of variables
  int n=A.nrow();
  //NAnce - #ancestors for each variable. Calculated in the same way as rank function does
  IntegerVector NAnc(n);
  for(int i=0;i<n;i++){
    int anc=0;
    int check=0;
    IntegerVector Anc(n,0);
    for(int j=0;j<n;j++){
      Anc[j]=A(j,i);
      check+=A(j,i);
    }
    while(anc!=check){
      anc=check;
      for(int j=0;j<n;j++){
        if(Anc[j]==1){
          for(int k=0;k<n;k++){
            if(Anc[k]==0){
              Anc[k]+=A(k,j);
            }
          }
        }
      }
      check=0;
      for(int j=0;j<n;j++){
        check+=Anc[j];
      }
    }
    NAnc[i]=anc;
  }
  //TopOrder - topological order of the variables (according to CP-net structure)
  //Obtained by puting the variables in increasing order of the #ancestors of a variable
  IntegerVector TopOrder(n);
  int counter=0;
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      if(NAnc[j]==i){
        TopOrder[counter]=j;
        counter+=1;
      }
    }
  }
  //M - maximum domain size of variables in CP-net
  int M=0;
  for(int i=0;i<n;i++){
    if(N[i]>M){
      M=N[i];
    }
  }
  //Domains - matrix where ith row gives the possible domain values for variable i (indexed from 1 i.e 1,2 for binary vectors). Any remaining elements of a row are filled with 0s
  IntegerMatrix Domains(n,M);
  for(int i=0;i<n;i++){
    for(int j=0;j<N[i];j++){
      Domains(i,j)=(j+1);
    }
  }
  //We now prune the domain of each variable (in topological order) according to forward pruning
  for(int i=0;i<n;i++){
    //x - ith variable in topological order
    int x=TopOrder[i];
    //DTG - domain transition graph for x (as an adjacency matrix)
    IntegerMatrix DTG(N[x],N[x]);
    //Pa - 0/1 vector showing the parents of x
    //PA - number of parents
    IntegerVector Pa(n,0);
    int PA=0;
    for(int j=0;j<n;j++){
      Pa[j]=A(j,x);
      PA+=A(j,x);
    }

    if(PA==0){
      //If x has no parents
      //pref - the single row that constitutes CPT(x)
      IntegerVector pref=CPTRow(A, N, ENTRIES, BREAK, (x+1), 0);
      //For j from 1 to |Dom(x)|-1, add an edge to the DTG from the domain value in preference position j+a, to that in position j (in pref)
      for(int j=1;j<(N[x]);j++){
        int target;
        int source;
        for(int k=0;k<N[x];k++){
          if(pref[k]==j){
            target=k;
          }
          if(pref[k]==(j+1)){
            source=k;
          }
        }
        DTG(source,target)=1;
      }
    }
    if(PA>0){
      //If x has parents
      //NPa - 0/1 vector iving the parent variables of x
      IntegerVector NPa(PA);
      counter=0;
      for(int j=0;j<n;j++){
        if(Pa[j]==1){
          NPa[counter]=j;
          counter+=1;
        }
      }
      //PrunedPaDom - current domain sizes of the parents (note that the parent domains have been pruned already)
      IntegerVector PrunedPaDom(PA);
      for(int j=0;j<PA;j++){
        int pa=NPa[j];
        int dom=0;
        //size of domain of pa is the number of non-zero entries in the corresponding Domain Matrix row 
        for(int k=0;k<M;k++){
          if(Domains(pa,k)!=0){
            dom+=1;
          }
        }
        PrunedPaDom[j]=dom;
      }
      //Once variables are pruned we have domains like {1,3,4} with values missing. For ease we will still consider these to be assignments 1,2, and 3 and then, where necessary translate these to be the actual domain values
      //We call these normalised and actual assignments
      //NormAsst - first (lexicographically) normalised assignment to the parents of X
      IntegerVector NormAsst(PA,1);
      //ActualAsst - find the actual assignment corresponding to NormAsst
      IntegerVector ActualAsst(PA);
      for(int j=0;j<PA;j++){
        //For each parent
        counter=1;
        //cycle through the possible domain value until you find the mth no-pruned value (if NormAsst assigns variable j the value m)
        //this domain value is the actual assignment for variable j
        for(int k=0;k<M;k++){
          if(Domains(NPa[j],k)!=0){
            if(counter==NormAsst[j]){
              ActualAsst[j]=Domains(NPa[j],k);
              break;
            }
            else{
              counter+=1;
            }
          }
        }
      }
      //Cycle through all NormAsst (normalised assignments to the parents of x) lexicographically [that is, all non-pruned parental assignments normalised]
      //For each assignment, we will add the associated edges to DTG
      bool Complete=false;
      while(!Complete){
        //Diff - hamming distance between NormAsst and PrunedPaDom
        int Diff=0;
        for(int j=0;j<PA;j++){
          if(NormAsst[j]!=PrunedPaDom[j]){
            Diff+=1;
          }
        }
        if(Diff==0){
          //If NormAsst = PrunedPaDom, then NormAsst is the largest/last parental assignment so we exit the while loop after this iteration
          Complete=true;
        }
        //pref - CPT(x) row corresponding to the current assignment to parents (ActualAsst)
        IntegerVector pref=CPTRow(A, N, ENTRIES, BREAK, (x+1), ActualAsst);
        //For j from 1 to |Dom(x)|-1, add an edge to the DTG from the domain value in preference position j+a, to that in position j (in pref)
        for(int j=1;j<(N[x]);j++){
          int target;
          int source;
          for(int k=0;k<N[x];k++){
            if(pref[k]==j){
              target=k;
            }
            if(pref[k]==(j+1)){
              source=k;
            }
          }
          DTG(source,target)=1;
        }
        if(!Complete){
          //If we are performing another iteration, move NormAsst to the next (lexicographically) assignment
          int MaxDiff=0;
          for(int j=0;j<PA;j++){
            if(NormAsst[j]!=PrunedPaDom[j]){
              if(j>MaxDiff){
                MaxDiff=j;
              }
            }
          }
          NormAsst[MaxDiff]+=1;
          if(MaxDiff<(PA-1)){
            for(int j=(MaxDiff +1);j<PA;j++){
              NormAsst[j]=1;
            }
          }
          //Find the actual assignment corresponding to NormAsst (ActualAsst)
          for(int j=0;j<PA;j++){
          counter=1;
          for(int k=0;k<M;k++){
            if(Domains(NPa[j],k)!=0){
              if(counter==NormAsst[j]){
                ActualAsst[j]=Domains(NPa[j],k);
                break;
              }
              else{
                counter+=1;
              }
            }
          }
        }
        }
      }
    }
    //We have now built the domain transition graph
    //We now prune the domain of x
    //Target - index of domain value taken by x in o1
    //Source - index of domain value taken by x in o2
    int Target=o1[x]-1;
    int Source=o2[x]-1;
    //Forward - 0/1 vector of the domain values reachable in DTS from source (this trivially includes source)
    //Backward - 0/1 vector of the domain values that can reach target in DTS (this trivially includes target)
    IntegerVector Forward(N[x],0);
    IntegerVector Backward(N[x],0);
    Forward[Source]=1;
    Backward[Target]=1;
    int FTot=1;
    //Add all children of Source in DTG to Forward
    for(int j=0;j<N[x];j++){
      if(DTG(Source,j)==1){
        Forward[j]=1;
      }
    }
    int FCheck=0;
    for(int j=0;j<N[x];j++){
      FCheck+=Forward[j];
    }
    //Add the children of the values in Forward to Forward repeatedly until this makes no difference to the vector
    while(FTot!=FCheck){
      FTot=FCheck;
      for(int j=0;j<N[x];j++){
        if(Forward[j]==1){
          for(int k=0;k<N[x];k++){
            if(DTG(j,k)==1){
              Forward[k]=1;
            }
          }
        }
      }
      FCheck=0;
      for(int j=0;j<N[x];j++){
        FCheck+=Forward[j];
      }
    }
    if(Forward[Target]==0){
      //If Target not reachable from Source then we would prune the whole domain of x in forward pruning and thus the original dominance query is false
      //Return the variable we pruned completely and the matrix of domains as they are currently (mid way through pruning)
      //Note - output format indicates the query to be false
      List L=List::create(x,Domains);
      return L;
    }
    int BTot=1;
    //Add all parents of Target in DTG to Backward
    for(int j=0;j<N[x];j++){
      if(DTG(j,Target)==1){
        Backward[j]=1;
      }
    }
    int BCheck=0;
    for(int j=0;j<N[x];j++){
      BCheck+=Backward[j];
    }
    //Add the parents of the values in Backward to Backward repeatedly until this makes no difference to the vector
    while(BTot!=BCheck){
      BTot=BCheck;
      for(int j=0;j<N[x];j++){
        if(Backward[j]==1){
          for(int k=0;k<N[x];k++){
            if(DTG(k,j)==1){
              Backward[k]=1;
            }
          }
        }
      }
      BCheck=0;
      for(int j=0;j<N[x];j++){
        BCheck+=Backward[j];
      }
    }
    //Intersection - 0/1 vector indicating the x values that lie in both Forward and Backward
    //These are the values we DONT prune
    IntegerVector Intersection(N[x],0);
    for(int j=0;j<N[x];j++){
      if((Forward[j]+Backward[j])==2){
        Intersection[j]=1;
      }
    }
    //Prune the domain values of x that are not in Intersection
    for(int j=0;j<N[x];j++){
      if(Intersection[j]==0){
        Domains(x,j)=0;
      }
    }
  }
  //Each variable has now had its domain pruned, all domains are still of size >0
  //the pruned domains are given in the matrix Domains
  //We now need to reduce the CP-net according to these value removals
  //FPN - list of pruned domain sizes
  IntegerVector FPN(n);
  for(int i=0;i<n;i++){
    int dom=0;
    for(int j=0;j<M;j++){
      if(Domains(i,j)!=0){
        dom+=1;
      }
    }
    FPN[i]=dom;
  }
  //NVar - number of variables that have >1 value in their pruned domains (variables with only 1 value will be removed)
  int NVar=0;
  for(int i=0;i<n;i++){
    if(FPN[i]>1){
      NVar+=1;
    }
  }
  //FPVar - indices  of variables that have >1 value in their pruned domains
  IntegerVector FPVar(NVar);
  counter=0;
  for(int i=0;i<n;i++){
    if(FPN[i]>1){
      FPVar[counter]=i;
      counter+=1;
    }
  }
  //FPA - CP-net adjacency matrix reduced to the variables in FPVar
  IntegerMatrix FPA(NVar,NVar);
  for(int i=0;i<NVar;i++){
    for(int j=0;j<NVar;j++){
      FPA(i,j)=A(FPVar[i],FPVar[j]);
    }
  }
  //FPo1 - o1 reduced to the variables in FPVar and normalised (as previously)
  //FPo2 - o2 reduced to the variables in FPVar and normalised (as previously)
  IntegerVector FPo1(NVar);
  IntegerVector FPo2(NVar);
  for(int i=0;i<NVar;i++){
    //For ith remaining variables (FPVar), v
    //x1- value taken by v in o1
    //x2- value taken by v in o2
    int x1=o1[FPVar[i]];
    int x2=o2[FPVar[i]];
    counter=1;
    //Now go through the domain of v and if x1 is the kth non-zero entry, we enter k into FPo1, similarly for FPo2
    for(int j=0;j<M;j++){
      if(Domains(FPVar[i],j)>0){
        if(Domains(FPVar[i],j)==x1){
          FPo1[i]=counter;
        }
        if(Domains(FPVar[i],j)==x2){
          FPo2[i]=counter;
        }
        counter+=1;
      }
    }
  }
  //To finish reducing the CP-net (with respect to pruned romains and removing variables with 1 value) we need to obtain the reduced entries and breaks vectors
  //FPBreaks will be the reduced CP-net breaks vector
  IntegerVector FPBreaks((NVar+1));
  //Calculate as usual, using the reduced structure and pruned domain sizes (CPT size is the multipilcation of its domain size and all its parents domain sizes)
  FPBreaks[0]=1;
  for(int i=0;i<NVar;i++){
    int CPTSize=FPN[FPVar[i]];
    for(int j=0;j<NVar;j++){
      if(FPA(j,i)==1){
        CPTSize*=FPN[FPVar[j]];
      }
    }
    FPBreaks[(i+1)]=CPTSize+FPBreaks[i];
  }
  //FPNRed - domain sizes for the remaining variables
  IntegerVector FPNRed(NVar);
  for(int i=0;i<NVar;i++){
    FPNRed[i]=FPN[FPVar[i]];
  }
  //FPBreaks will be the reduced CP-net entries vector
  IntegerVector FPEntries((FPBreaks[NVar]-1));
  //LostParents- #parents in reduced structure for each variable in reduced structure that lost a parent (0 otherwise)
  IntegerVector LostParents(NVar,0);
  //FPParents - #parents of each variable in reduced structure
  IntegerVector FPParents(NVar,0);
  for(int i=0;i<NVar;i++){
    //For each remaining variable, X
    //OldPa - #parents X had in the original structure
    //NwPa - #parents X has in the reduced structure
    int OldPa=0;
    int NewPa=0;
    for(int j=0;j<n;j++){
      OldPa+=A(j,FPVar[i]);
    }
    for(int j=0;j<NVar;j++){
      NewPa+=FPA(j,i);
    }
    //Add values to FPParents and LostParents
    FPParents[i]=NewPa;
    if(OldPa>NewPa){
      LostParents[i]=NewPa;
    }
    if(FPEntries.size()==ENTRIES.size()){
      //if new entries vector has same size, nothing is removed and so we can simply copy over the entries vector
      FPEntries=ENTRIES;
    }
    else{
      if(NewPa==0){
        //If x has no parents in the reduced structure
        //PaAsst - Parental assignment of x (in original CP-net) in o1 (and o2)
        IntegerVector PaAsst(OldPa);
        counter=0;
        for(int j=0;j<n;j++){
          if(A(j,FPVar[i])==1){
            PaAsst[counter]=o1[j];
            counter+=1;
          }
        }
        //pref- CPT(X) row corresponding to the parental assignment in o1
        //if all parents of x are removed by reducing domain to 1, then they are all fixed at their o1 values
        //Thus, pref is preference row for CPT(X) in the reduced structure (no parents means the CPT is a single row)
        IntegerVector pref=CPTRow(A, N, ENTRIES, BREAK, (FPVar[i]+1), PaAsst);
        //RedPref reduces this preference order to the remaining (unpruned) domain values and normalises it so that preference positions are values from 1 to #unpruned values
        IntegerVector RedPref(FPNRed[i]);
        //Reduce pref to the unpruned value entries
        counter=0;
        for(int j=0;j<M;j++){
          if(Domains(FPVar[i],j)>0){
            RedPref[counter]=pref[j];
            counter+=1;
          }
        }
        //normalise RefPref
        counter=1;
        for(int j=1;j<=N[FPVar[i]];j++){
          for(int k=0;k<FPNRed[i];k++){
            if(RedPref[k]==j){
              RedPref[k]=counter;
              counter+=1;
            }
          }
        }
        //RedPref is now the CPT(X) well defined for our reduced CP-net
        //Add RedPref to our entries vector
        for(int j=0;j<FPNRed[i];j++){
          FPEntries[FPBreaks[i]-1+j]=RedPref[j];
        }
      }
      else{
        // if X has lost some but not all parents:
        //PaDom - reduced parent (of x) domain sizes for original structure
        //PaVar - parents of x in original structure
        IntegerVector PaDom(OldPa);
        IntegerVector PaVar(OldPa);
        counter=0;
        for(int j=0;j<n;j++){
          if(A(j,FPVar[i])==1){
            PaDom[counter]=FPN[j];
            PaVar[counter]=j;
            counter+=1;
          }
        }
        //We again cycle through all possible assignments to Pa(X) in original structure (only taking pruned values) in lexicographic order
        //LexPaAsst - normalised assignment (this is the one we cycle through)
        //ActualPaAsst -actual assignment
        int CPTcounter=0;
        //LexPaAsst first possible assignment
        IntegerVector LexPaAsst(OldPa,1);
        //determine the actual corresponding assignment
        IntegerVector ActualPaAsst(OldPa);
        for(int j=0;j<OldPa;j++){
          int pa=PaVar[j];
          counter=1;
          for(int k=0;k<M;k++){
            if(Domains(pa,k)>0){
              if(counter==LexPaAsst[j]){
                ActualPaAsst[j]=Domains(pa,k);
                break;
              }
              counter+=1;
            }
          }
        }
        bool Comp=false;
        while(!Comp){
          //Diff - Hamming distance between LexPaAsst and PaDom
          int Diff=0;
          for(int j=0;j<OldPa;j++){
            if(LexPaAsst[j]!=PaDom[j]){
              Diff+=1;
            }
          }
          if(Diff==0){
            //if LexPaAsst = PaDom, then LexPaAsst is the last/largest parental assignment and this is the last iteration of the loop
            Comp=true;
          }
          //extract the (original) CPT(X) row corresponding to ActualPaAsst - pref
          IntegerVector pref=CPTRow(A, N, ENTRIES, BREAK, (FPVar[i]+1), ActualPaAsst);
          //RedPref - pref restricted to unpruned domain values and normalised
          IntegerVector RedPref(FPNRed[i]);
          counter=0;
          for(int j=0;j<M;j++){
            if(Domains(FPVar[i],j)>0){
              RedPref[counter]=pref[j];
              counter+=1;
            }
          }
          counter=1;
          for(int j=1;j<=N[FPVar[i]];j++){
            for(int k=0;k<FPNRed[i];k++){
              if(RedPref[k]==j){
                RedPref[k]=counter;
                counter+=1;
              }
            }
          }
          //Add RedPref to our new entries vector
          for(int j=0;j<FPNRed[i];j++){
            FPEntries[FPBreaks[i]-1+CPTcounter]=RedPref[j];
            CPTcounter+=1;
          }
          //This process adds the preference rule for each remaining (unpruned) parent assignment of X, changed to be well defined for our new CP-net
          //this gives the forward pruning reduced CPT(X) once all parental assignments completed
          if(!Comp){
            //If we are performing another iteration, change LexPaAsst to the (lexicographically) next (not pruned) parent assignment
            int D=0;
            for(int j=0; j<OldPa;j++){
              if(LexPaAsst[j]!=PaDom[j]){
                if(j>D){
                  D=j;
                }
              }
            }
            LexPaAsst[D]+=1;
            if(D<(OldPa-1)){
              for(int j=(D+1);j<OldPa;j++){
                LexPaAsst[j]=1;
              }
            }
            //Determine the actual assignment that LexPaAsst corresponds to - this is now ActualPaAsst
            for(int j=0;j<OldPa;j++){
              int pa=PaVar[j];
              counter=1;
              for(int k=0;k<M;k++){
                if(Domains(pa,k)>0){
                  if(counter==LexPaAsst[j]){
                    ActualPaAsst[j]=Domains(pa,k);
                    break;
                  }
                  counter+=1;
                }
              }
            }
          }
        }
      }
    }
  }
  //We have now added the forward pruning reduced CPTs of each variable to FPEntries, making it the entries vector for the forward pruning reduced CP-net
  //Thus, the forward pruning reduced CP-net is not fully defined and in standard CP-net form
  //Finally, we check that no remaining parent-child edges have become degenerate in this process
  //ParentCheck is a 0/1 vector indicating which of the remaining variables we need to check the parent degeneracy for
  //i.e. which parent relations could plausibly have become degenerate
  IntegerVector ParentCheck(NVar,0);
  for(int i=0;i<NVar;i++){
    //For each remaining variable, x
    if(FPParents[i]>0){
       if(LostParents[i]>0){
        //If X lost a parent in forward pruning, but not all parents, then we must check parent validity
        ParentCheck[i]=1;
       }
       else{
        if(FPN[FPVar[i]]!=N[FPVar[i]]){
          //If X lost domain values, then we must check parent validity
          ParentCheck[i]=1;
        }
        else{
          for(int j=0;j<NVar;j++){
            if(FPA(j,i)==1){
              //for any of the remaining parents of x, y
              if(FPN[FPVar[j]]!=N[FPVar[j]]){
              //If a y lost domain values, then we must check x parent validity
              ParentCheck[i]=1;
             }
            }
          }
        }
      }
    }
  }
  for(int i=0;i<FPA.nrow();i++){
    if(ParentCheck[i]>0){
      //For each remainig variable, x, for which we need to check parent validity:
      //RemainingParents - # parents of x in reduced structure
      int RemainingParents=FPParents[i];
      //NewPa - indices of parents of x in reduced structure
      //NewPaDom - (pruned) domain sizes of the parents of x, with domain size of x at end
      IntegerVector NewPa(FPParents[i]);
      IntegerVector NewPaDom(FPParents[i]+1);
      counter=0;
      for(int j=0;j<NVar;j++){
        if(FPA(j,i)==1){
          NewPa[counter]=j;
          NewPaDom[counter]=FPNRed[j];
          counter+=1;
        }
      }
      NewPaDom[FPParents[i]]=FPNRed[i];
      //CPTi - CPT(X) in reduced CP-net, extracted from reduced CP-net entries vector
      IntegerVector CPTi(FPBreaks[(i+1)]-FPBreaks[i]);
      for(int j=0;j<(FPBreaks[(i+1)]-FPBreaks[i]);j++){
        CPTi[j]=FPEntries[FPBreaks[i]-1+j];
      }
      //Check whether this CPT(X) is degenerate
      int D=Degenerate(NewPaDom, CPTi);
      //If it is degenerate, we obtain an invalid (degenerate) parent
      //We repeatedly remove the invalid parent-child edge and then assess again if the new CPT(X) is degenerate until all invalid parents are removed and so the CPT is non-degenerate
      while(D>=0){
        //if CPT(X) is degenerate:
        //D is an invalid parent
        //L - CP-net obtained by removing the edge from degenerate parent to x in the FP... CP-net
        List L=RemovePa(FPA, FPNRed, FPEntries, FPBreaks, NewPa[D], i);
        //Update FP... CP-net to be this CP-net, L
        IntegerMatrix FPA=L[0];
        FPEntries=L[1];
        FPBreaks=L[2];
        // Number of parents of X in FP structure has reduced by 1
        RemainingParents-=1;
        //NewCPTi - CPT(X) in this newly reduced (edge removed) CP-net, extracted from entries vector
        IntegerVector NewCPTi(FPBreaks[(i+1)]-FPBreaks[i]);
        for(int j=0;j<(FPBreaks[(i+1)]-FPBreaks[i]);j++){
          NewCPTi[j]=FPEntries[FPBreaks[i]-1+j];
        }
        //RedNewPa - indices of new parent set of X after edge removal
        //RedNewPaDom - (pruned) domain sizes of these parents, followed by domain size of x
        IntegerVector RedNewPa(RemainingParents);
        IntegerVector RedNewPaDom(RemainingParents+1);
        counter=0;
        for(int j=0;j<FPA.nrow();j++){
          if(FPA(j,i)==1){
            RedNewPa[counter]=j;
            RedNewPaDom[counter]=FPNRed[j];
            counter+=1;
          }
        }
        RedNewPaDom[RemainingParents]=FPNRed[i];
        //Replace NewPaDom and NewPa with RedNewPaDom and RedNewPa so the next iteration works correctly
        NewPaDom=RedNewPaDom;
        NewPa=RedNewPa;
        //Assess whether the new CPT(X) (after edge removal) is still degenerate
        D=Degenerate(NewPaDom, NewCPTi);
      }
    }
  }
  //We have now removed all degenerate parent-child relations that could have been caused by our forward pruning reduction
  //The result is a CP-net that is reduced by forward pruning, with no degenerate edges, and in the standard (normalised) form of a CP-net - this is L
  List L=List::create(Domains,FPA,FPNRed,FPEntries,FPBreaks,FPo1,FPo2);
  //return this CP-net
  return L;
}

// [[Rcpp::export]]
List FPAndRPSDQRankPriority(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2){
  //Takes a CP-net and a dominance query (is o1>o2?)
  //Applies forward pruning and answers the reduced query using rank pruning, penalty pruning and suffix fixing with rank priority
  //Returns the answer of the query, the number of outcomes considered and the following timestamps
  //start and end of function if the query is answered via forward pruning
  //start, end of preprocessing, and end of function otherwise
  //Also output #Outcomes in forward pruning reduced CP-net if preprocessing does not answer the query
  //Log start time of the function
  Timer time;
  time.step("Start");
  //Evaluate numerical checks on the query
  bool Check=NumericalCheck(A, N, ENTRIES, BREAK, o1, o2);
  if(!Check){
    // If Numerical Checks fail then the query is false
    //  return false, 0 outcomes comsidered and start/end times
    //Log the end time
    time.step("End");
    NumericVector TimeJumps(time);
    List R=List::create(false,0,TimeJumps);
    return R;
  }
  //Apply forward pruning to our CP-net and query - result given by PrunedQuery
  List PrunedQuery=ForwardPruning(A, N, ENTRIES, BREAK, o1, o2);
  if(PrunedQuery.size()==2){
    //If forward pruning answers the query (finds it false):
    //  return false, 0 outcomes comsidered and start/end times
    //log the end time
    time.step("End");
    NumericVector TimeJumps(time);
    List R=List::create(false,0,TimeJumps);
    return R;
  }
  //If forward pruning doesnt answer the query then PrunedQuery is the reduced CP-net and query obtained by applying forward pruning
  //Apply numerical checks to reduced query
  Check=NumericalCheck(PrunedQuery[1], PrunedQuery[2], PrunedQuery[3], PrunedQuery[4], PrunedQuery[5], PrunedQuery[6]);
  if(!Check){
    // If Numerical Checks fail then the original input query is false
    //  return false, 0 outcomes comsidered and start/end times
    //Log the end time
    time.step("End");
    NumericVector TimeJumps(time);
    List R=List::create(false,0,TimeJumps);
    return R;
  }
  //Log the time we finish applying preprocessing (forward pruning)
  time.step("DonePreProc");
  //NRed - vector of domain sizes of the forward pruning reduced CP-net
  IntegerVector NRed=PrunedQuery[2];
  //nAfter - number of outcomes in the forward pruning reduced CP-net
  long long int nAfter=1;
  for(int i=0;i<NRed.size();i++){
    nAfter*=NRed[i];
  }
  //Answer the dominance query obtained by reducing CP-net via forward pruning
  List Answer=RPSDQRankPriority(PrunedQuery[1], PrunedQuery[2], PrunedQuery[3], PrunedQuery[4], PrunedQuery[5], PrunedQuery[6]);
  //log the end time of the function
  time.step("End");
  //Return the answer of the dominance query (which is also the answer of the input query), number of outcomes considered in dominance testing,
  //the start/end of forward pruning/end timestamps, and the #outcomes in forward pruning reduced CP-net
  NumericVector TimeJumps(time);
  if(!Answer[0]){
    List R=List::create(false,Answer[1],TimeJumps,nAfter);
    return R;
  }
  List R=List::create(true,Answer[1],TimeJumps,nAfter);
  return R;
}

// [[Rcpp::export]]
List CombAndRPSDQRankPriority(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2){
  //input - CP-net and dominance query
  //Function applies the combination of UVRS and Forward pruning (described in the thesis, see page 122 remark about the distinction between this code and the algorithm presented) to the query
  //If this preprocessing results in reduced queries, we answer these in size order using rank pruning, penalty pruning and suffix fixing with rank priority
  //Returns the answer of the original query, the number of outcomes considered and the following timestamps
  //start and end of function if the query is answered via preprocessing
  //start, end of preprocessing, and end of function otherwise
  //Also output #Outcomes in preprocessed (reduced) CP-net if preprocessing does not answer the query
  //Log start time of the function
  Timer time;
  time.step("Start");
  //Apply numerical checks to original input query
  bool Check=NumericalCheck(A, N, ENTRIES, BREAK, o1, o2);
  if(!Check){
    // If Numerical Checks are failed then query is false
    //return false, 0 outcomes considered and start/end times
    //log end time 
    time.step("End");
    NumericVector TimeJumps(time);
    List R=List::create(false,0,TimeJumps);
    return R;
  }
  //Reduced - CP-net obtained by removing unimportant variables
  List Reduced=UvRemove(A,N,ENTRIES,BREAK,o1,o2);
  //Apply numerical checks to reduced CP-net
  Check=NumericalCheck(Reduced[0],Reduced[1],Reduced[2],Reduced[3],Reduced[4],Reduced[5]);
  if(!Check){
    // If Numerical Checks after initial reduction are failed then query is false
    //return false, 0 outcomes considered and start/end times
    //log end time 
    time.step("End");
    NumericVector TimeJumps(time);
    List R=List::create(false,0,TimeJumps);
    return R;
  }
  //Continue to remove unimportant variables until an iteration does not alter the CP-net structure
  //Same - are the structures before and after removal the same?
  bool Same=false;
  while(!Same){
    //NextIter - CP-net obtained by removing unimportant variables from reduced structure
    List NextIter=UvRemove(Reduced[0],Reduced[1],Reduced[2],Reduced[3],Reduced[4],Reduced[5]);
    //Original - structure before removal
    IntegerMatrix Original=Reduced[0];
    //New - structure after
    IntegerMatrix New=NextIter[0];
    if(Original.nrow()==New.nrow()){
      //If structure unchanged by removal:
      Same=true;
    }
    else{
      //If structure changed by removal:
      //replace Reduced by the further reduced CP-net, NextIter
      Reduced=NextIter;
      //Evaluate numerical checks for the newly reduced CP-net
      Check=NumericalCheck(Reduced[0],Reduced[1],Reduced[2],Reduced[3],Reduced[4],Reduced[5]);
      if(!Check){
        // If Numerical Checks after succesive reduction fail then original query was false
        //return false, 0 outcomes considered and start/end times
        //log end time 
        time.step("End");
        NumericVector TimeJumps(time);
        List R=List::create(false,0,TimeJumps);
        return R;
      }
    }
  }
  //Reduced is the now the CP-net and query after iterative Uimplortant Variable (UV) removal
  //(and no further reduction can be done this way as last attempt did not affect structure)
  //Conn - connected components in the reduced CP-net structure
  List Conn=ConnectedComponents(Reduced[0]);
  //SubCPN - the reduced CP-net split into its connected component
  List SubCPN;
  //CurrentStruc - list of the adjacency matrices of the (sub-)CP-nets given in SubCPN
  List CurrentStruc;
  if(Conn[0]==1){
    //if there is only 1 connected component, we don't need to split up the CP-net
    //SubCPN is a list of 1 - the reduced CPN and DQ
    SubCPN=List::create(Reduced);
    //CurrentStruc in this case is simply a list of 1 - the adjacency matrix of Reduced CP-net
    CurrentStruc=List::create(Reduced[0]);
    //SmallN - domain size vector for Reduced
    IntegerVector SmallN=Reduced[1];
  }
  else{
    //For each connected component, restrict the reduced CP-net and query to this component and add this sub-CP-net to our list SubCPN
    for(int k=0;k<Conn[0];k++){
      // For each connected component
      //B - structure of the reduced CP-net
      IntegerMatrix B=Reduced[0];
      //CCM - matrix of connected component entries
      IntegerMatrix CCM=Conn[1];
      //CC - 0/1 vector of kth connected component (row k of CCM)
      IntegerVector CC(B.nrow(),0);
      //NVar - number of variables in the kth connected component
      int NVar=0;
      for(int i=0;i<B.nrow();i++){
        CC[i]+=CCM(k,i);
        NVar+=CCM(k,i);
      }
      //SubA - B reduced to variables in kth connected component (Adjacency matrix of the kth connected component)
      IntegerMatrix SubA(NVar,NVar);
      int counter1=0;
      int counter2=0;
      for(int i=0;i<B.nrow();i++){
        if(CC[i]==1){
          for(int j=0;j<B.nrow();j++){
            if(CC[j]==1){
              SubA(counter1,counter2)=B(i,j);
              counter2+=1;
            }
          }
          counter2=0;
          counter1+=1;
        }
      }
      //Add SubA to CurrentStruc - this is the adjacency matrix of the kth connected component and thus of the kth (sub-)CP-net in SubCPN
      CurrentStruc.push_back (SubA);
      //SubN - domain sizes vector for the kth connected component variables
      IntegerVector SubN(NVar);
      //Subo1 - outcome 1 restricted to the kth connected component variables
      IntegerVector Subo1(NVar);
      //Subo2 - outcome 2 restricted to the kth connected component variable
      IntegerVector Subo2(NVar);
      //Obtain these by taking the relevant element of Reduced and restricting to the kth connected component variables
      int counter=0;
      for(int i=0;i<B.nrow();i++){
        if(CC[i]==1){
          IntegerVector RedN=Reduced[1];
          SubN[counter]=RedN[i];
          IntegerVector Redo1=Reduced[4];
          Subo1[counter]=Redo1[i];
          IntegerVector Redo2=Reduced[5];
          Subo2[counter]=Redo2[i];
          counter+=1;
        }
      }
      //Calculate the breaks vector for the CP-net restricted to the kth connected component variables
      IntegerVector SubBreaks((NVar+1));
      SubBreaks[0]=1;
      counter=1;
      //Use breaks vector from Reduced to calculate the CPT lengths of the kth connected component variables
      IntegerVector RedBreaks=Reduced[3];
      IntegerVector RedEntries=Reduced[2];
      for(int i=0;i<B.nrow();i++){
        if(CC[i]==1){
          //For each variable in kth component, add the CPT length to obtain the next breaks value
          SubBreaks[counter]=SubBreaks[(counter-1)]+RedBreaks[(i+1)]-RedBreaks[(i)];
          counter+=1;
        }
      }
      //SubEntries - entries vector for the CP-net reduced to the kth connected component variables
      //Obtained by copying the CPTs for the kth connected component variables over from the entries vector for Reduced
      IntegerVector SubEntries(SubBreaks[NVar]-1);
      counter=0;
      for(int i=0;i<B.nrow();i++){
        if(CC[i]==1){
          //For each of the kth connected component variables, x:
          //Extract CPT length from Reduced breaks vector - CPTSize
          int CPTSize=RedBreaks[(i+1)]-RedBreaks[(i)];
          //Copy the CPT(x) entries from the Reduced entries vector into SubEntries
          for(int j=0;j<CPTSize;j++){
            SubEntries[counter]=RedEntries[(RedBreaks[(i)]-1+j)];
            counter+=1;
          }
        }
      }
      //Sub... constitute the CP-net and dominance query in Reduced restricted to the kth connected component
      //Apply numerical checks to this Sub... dominance query
      bool Check=NumericalCheck(SubA,SubN,SubEntries,SubBreaks,Subo1,Subo2);
      if(!Check){
        //If the checks fail on this sub-CP-net then the original query is also false
        //return false, 0 outcomes considered and start/end times
        //log end time 
        time.step("End");
        NumericVector TimeJumps(time);
        List R=List::create(false,0,TimeJumps);
        return R;
      }
      //SubRed - the Sub CP-net (and query) we obtained by restricting to the kth connected component
      List SubRed=List::create(SubA,SubN,SubEntries,SubBreaks,Subo1,Subo2);
      //Add SubRed to out list of sub-CP-nets
      SubCPN.push_back (SubRed);
    }
  }
  //SubCPN - list of reduced CP-net and dominance queries produced by applying UVRS to input query
  //FinishedCPNs will be the list of reduced CP-net and dominance queries produced by applying UVRS+Forward pruning combination to input query
  //When a query cannot be further reduced by UVRS or forward pruning, we add it to FinishedCPNs
  List FinishedCPNs=List::create();
  //STOP - boolean that tells us when we should stop applying preprocessing (each iteration applies forward pruning and then UVRS to the remaining queries - first iteration starts with SubCPN queries)
  //Essentially, STOP becomes true when there are no queries left that can be potentially further reduced (they have all been added to FinishedCPNs)
  bool STOP=false;
  while(!STOP){
    //NewSubCPN - The sub-CP-nets (and queries) in SubCPN after being reduced by forward pruning
    //NewStruc - structures of the Sub-CP-nets in NewSubCPN
    List NewSubCPN=List::create();
    List NewStruc=List::create();
    //For each Sub-CP-net (and associated dominance query) in SubCPN, apply forward pruning to the query
    //If forward pruning finds the query false, the original input query is false (function will terminate)
    //Otherwise, record the reduced CP-net (and query) in NewSubCPN and the structure of this reduced CP-net in NewStruc
    for(int cc=0;cc<SubCPN.size();cc++){
      //For each Sub-CP-net and associated dominance query in SubCPN, DQi
      List DQi=SubCPN[cc];
      //Apply Forward pruning to DQi
      List FPred=ForwardPruning(DQi[0],DQi[1],DQi[2],DQi[3],DQi[4],DQi[5]);
      if(FPred.size()==2){
      //If forward pruning answers the query (finds it false) then the original input query is false:
      //  return false, 0 outcomes comsidered and start/end times
      //log the end time
      time.step("End");
      NumericVector TimeJumps(time);
      List R=List::create(false,0,TimeJumps);
      return R;
      }
      //Evaluate numerical checks for forward pruning reduced query
      bool Check=NumericalCheck(FPred[1],FPred[2],FPred[3],FPred[4],FPred[5],FPred[6]);
      if(!Check){
        // If Numerical Checks fail then the original input query is false
        //  return false, 0 outcomes comsidered and start/end times
        //Log the end time
        time.step("End");
        NumericVector TimeJumps(time);
        List R=List::create(false,0,TimeJumps);
        return R;
      }
      //RedDQi - reduced CP-net and dominance query produced by applying forward pruning
      List RedDQi=List::create(FPred[1],FPred[2],FPred[3],FPred[4],FPred[5],FPred[6]);
      //Add RedDQi to NewSubCPN
      NewSubCPN.push_back (RedDQi);
      //Add structure of RedDQi CP-net to NewStruc
      NewStruc.push_back (RedDQi[0]);
    }
    //Diffs - 0/1 vector indicating, for each sub-CP-net, whether the above forward pruning application altered the structure (1 for yes, 0 for no)
    IntegerVector Diffs(NewStruc.size());
    for(int i=0;i<NewStruc.size();i++){
      //For each of the Sub-CP-nets (after the above forward pruning reduction):
      //StrucDiff - did applying forward pruning change structure of ith sub-CP-net?
      bool StrucDiff=false;
      //New - structure of ith sub-CP-net after forward pruning reduction
      IntegerMatrix New=NewStruc[i];
      //Old -structure prior to forward pruning reduction
      IntegerMatrix Old=CurrentStruc[i];
      if(New.nrow()!=Old.nrow()){
        //If new and old have different dimensions, the structure was changed
        StrucDiff=true;
      }
      else{
        //If new and old have same dimension, check equality of each matrix entry
        for(int j=0;j<New.nrow();j++){
          for(int k=0;k<New.nrow();k++){
            if(New(j,k)!=Old(j,k)){
              //If any entries are different, the structure was changed and we can stop checking (break out of this loop)
              StrucDiff=true;
              break;
            }
          }
          if(StrucDiff){
            break;
          }
        }
      }
      //If structure changed, ith entry of Diffs is 1, otherwise its 0
      if(StrucDiff){
        Diffs[i]=1;
      }
      else{
        Diffs[i]=0;
      }
    }
    //Note only the SubCP-nets with Diff=1 can be further reduced by UVRS or forward pruning
    //SubCPN - list of the sub-CP-nets (and queries), after being reduced by forward pruning, for which forward pruning changed the structure
    SubCPN=List::create();
    //CurrentStruc - structures of the CP-nets in SubCPN
    CurrentStruc=List::create();
    //totaldiffs - number of sub-CP-nets for which forward pruning changed the structure (i.e. #CP-nets in SubCPN)
    int totaldiffs=0;
    //We add the sub-cp-nets where forward pruning did not affect the structure (i.e those not added to subCPN) to FinishedCPNs (as they cannot be further reduced)
    for(int i=0;i<NewStruc.size();i++){
      //For each of the Sub-CP-nets
      if(Diffs[i]==1){
        //If forward pruning changed the structure:
        //Add the sub-CP-net reduced by forward pruning (NewSubDQi) to SubCPN
        SubCPN.push_back (NewSubCPN[i]);
        List NewSubDQi=NewSubCPN[i];
        //Add the structure of NewSubDQi to CurrentStruc
        IntegerMatrix NewA=NewSubDQi[0];
        CurrentStruc.push_back (NewA);
        //count this Sub CP-net in our  totaldiffs count
        totaldiffs+=1;
      }
      if(Diffs[i]==0){
        //For those Sub-CP-nets for which forward pruning did not affect the structure:
        //Add the sub-CP-net reduced by forward pruning to FinishedCPNs (as it cannot be further reduced)
        FinishedCPNs.push_back (NewSubCPN[i]);
      }
    }
    if(totaldiffs==0){
      //If all sub-CP-nets added to FinishedCPNs then there are no more CP-nets left to be reduced further and so our preprocessing terminates
      STOP=true;
    }
    if(!STOP){
      //If there are Sub CP-nets in SubCPN (that is, there are still sub CP-nets we can reduce further):
      //Apply UVRS to each Sub-CP-net
      //If this makes any change to the Sub CP-net structure, we put the resulting reduced CP-net(s) in NewSubCPN
      NewSubCPN=List::create();
      //NewStruc - structures of the CP-nets in NewSubCPN
      NewStruc=List::create();
      //changes - #SubCP-nets such that UVRS makes a structural change and so reduced CP-nets are added to NewSubCPN
      int changes=0;
      //If UVRS does not alter the structure of a sub-cp-net then this query cannot be reduced any further by UVRS  or forward pruning
      //Thus, we add it to FinishedCPNs (note that in this case, UVRS has no affect at all so the sub CP-net and UVRS reduction of it are the same CP-net)
      for(int i=0;i<CurrentStruc.size();i++){
        //For each of the remaining sub CP-nets, DQi:
        List DQi=SubCPN[i];
        //Apply UVRS to DQi exactly as we did to the original input query
        //Remove unimportant variables from DQi 
        List Reducedi=UvRemove(DQi[0],DQi[1],DQi[2],DQi[3],DQi[4],DQi[5]);
        //Evaluate numerical checks for reduced CP-net
        Check=NumericalCheck(Reducedi[0],Reducedi[1],Reducedi[2],Reducedi[3],Reducedi[4],Reducedi[5]);
        if(!Check){
          time.step("End");
          NumericVector TimeJumps(time);
          List R=List::create(false,0,TimeJumps);
          return R;
        }
        //Repeatedly remove unimportant variables (and evaluate numerical checks) until an iteration does not change the structure
        bool Same=false;
        while(!Same){
          List NextIteri=UvRemove(Reducedi[0],Reducedi[1],Reducedi[2],Reducedi[3],Reducedi[4],Reducedi[5]);
          IntegerMatrix OriginalA=Reducedi[0];
          IntegerMatrix NewA=NextIteri[0];
          if(OriginalA.nrow()==NewA.nrow()){
            Same=true;
          }
          else{
            Reducedi=NextIteri;
            Check=NumericalCheck(Reducedi[0],Reducedi[1],Reducedi[2],Reducedi[3],Reducedi[4],Reducedi[5]);
            if(!Check){
              time.step("End");
              NumericVector TimeJumps(time);
              List R=List::create(false,0,TimeJumps);
              return R;
            }
          }
        }
        //Reducedi - DQi after iterative unimportant variable removal
        //Conn - connected components of Reducedi structure
        List Conn=ConnectedComponents(Reducedi[0]);
        if(Conn[0]==1){
          //If Reducedi has only 1 connected component:
          //Newi - structure of Reducedi
          IntegerMatrix Newi=Reducedi[0];
          //Oldi - structure of DQi
          IntegerMatrix Oldi=DQi[0];
          if(Newi.nrow()==Oldi.nrow()){
            //If iterative variable removal removed no variables (dimensions are the same), then no change was made - DQi=Reducedi
            //Add Sub-CP-net (and query) to FinishedCPNs as no further reduction by UVRS or forward pruning is possible
            FinishedCPNs.push_back (Reducedi);
          }
          else{
            //IF iterative variable removal did have an affect:
            //Add a count to changes
            changes+=1;
            //if there is only 1 connected component, we don't need to split up the CP-net into connected components
            //SubCPN is a list of 1 - the reduced CPN and DQ
            NewSubCPN.push_back (Reducedi);
            //Add structure to NewStruc
            NewStruc.push_back (Reducedi[0]);
          }
        }
        else{
          //If Reducedi has multiple connected components:
          //UVRS will have an effect on DQi so add a count to changes
          changes+=1;
          //Split Reduced i up into its connected component sub-CP-nets
          //Check numerical checks for each connected component CP-net - if they fail, original input query false and function terminates
          //If not  - Add these sub-cp-nets to NewSubCPN and their structures to NewStruc
          //These steps are perfomed identically to the corresponding section of our original application of UVRS to input query:
          for(int k=0;k<Conn[0];k++){
            IntegerMatrix B=Reducedi[0];
            IntegerMatrix CCM=Conn[1];
            IntegerVector CC(B.nrow(),0);
            int NVar=0;
            for(int i=0;i<B.nrow();i++){
              CC[i]+=CCM(k,i);
              NVar+=CCM(k,i);
            }
            IntegerMatrix SubA(NVar,NVar);
            int counter1=0;
            int counter2=0;
            for(int i=0;i<B.nrow();i++){
              if(CC[i]==1){
                for(int j=0;j<B.nrow();j++){
                  if(CC[j]==1){
                    SubA(counter1,counter2)=B(i,j);
                    counter2+=1;
                  }
                }
                counter2=0;
                counter1+=1;
              }
            }
            NewStruc.push_back (SubA);
            IntegerVector SubN(NVar);
            IntegerVector Subo1(NVar);
            IntegerVector Subo2(NVar);
            int counter=0;
            for(int i=0;i<B.nrow();i++){
              if(CC[i]==1){
                IntegerVector RedN=Reducedi[1];
                SubN[counter]=RedN[i];
                IntegerVector Redo1=Reducedi[4];
                Subo1[counter]=Redo1[i];
                IntegerVector Redo2=Reducedi[5];
                Subo2[counter]=Redo2[i];
                counter+=1;
              }
            }
            IntegerVector SubBreaks((NVar+1));
            SubBreaks[0]=1;
            counter=1;
            IntegerVector RedBreaks=Reducedi[3];
            IntegerVector RedEntries=Reducedi[2];
            for(int i=0;i<B.nrow();i++){
              if(CC[i]==1){
                SubBreaks[counter]=SubBreaks[(counter-1)]+RedBreaks[(i+1)]-RedBreaks[(i)];
                counter+=1;
              }
            }
            IntegerVector SubEntries(SubBreaks[NVar]-1);
            counter=0;
            for(int i=0;i<B.nrow();i++){
              if(CC[i]==1){
                int CPTSize=RedBreaks[(i+1)]-RedBreaks[(i)];
                for(int j=0;j<CPTSize;j++){
                  SubEntries[counter]=RedEntries[(RedBreaks[(i)]-1+j)];
                  counter+=1;
                }
              }
            }
            bool Check=NumericalCheck(SubA,SubN,SubEntries,SubBreaks,Subo1,Subo2);
            if(!Check){
              time.step("End");
              NumericVector TimeJumps(time);
              List R=List::create(false,0,TimeJumps);
              return R;
            }
            List SubRed=List::create(SubA,SubN,SubEntries,SubBreaks,Subo1,Subo2);
            NewSubCPN.push_back (SubRed);
          }
        }
      }
      //NewSubCPN is now the list of the CP-nets resulting form UVRS reducing the sub-cp-nets  (in the cases where UVRS had an effect and so the reduced CP-nets can be potentially further reduced)
      //Replace SubCPN and CurrentStruc so that the next iteration of the while loop performs correctly
      SubCPN=NewSubCPN;
      CurrentStruc=NewStruc;
      if(changes==0){
        //If changes ==0 then UVRS had no affect on any of the sub CP-nets and, thus, no queries can be reduced further by UVRS or forward pruning and so we do not attempt further preprocessing
        STOP=true;
      }
      //If UVRS has an effect and so SubCPN contains CP-nets (and queries) that may be reduced further, we apply this loop again, attempting to reduce them further by applying forward pruning and UVRS in turn
    }
  }
  //Once queries have been reduced as far as possible (resulting in FinishedCPNs), we cannot preprocess further
  //Log the time we finishe preprocessing
  time.step("DonePreProc");
  //Order - indeces 0 to (#reduced queries in FinishedCPNs) in order of increasing size (number of variables) of the associated CP-net in FinishedCPNs
  IntegerVector Order(FinishedCPNs.size());
  IntegerVector Record(FinishedCPNs.size(),1);
  int counter=0;
  while(counter<FinishedCPNs.size()){
    //Find the smallest component not yet listed and add index to Order
    int MinSize=(A.nrow())+1;
    int index;
    for(int i=0;i<FinishedCPNs.size();i++){
      //For each sub CP-net:
      List DQi=FinishedCPNs[i];
      //Ai - #variables in this CP-net
      IntegerMatrix Ai=DQi[0];
      //Record[i]==1 if and only if this sub CP-net not yet added to Order
      if((Ai.nrow()<MinSize)&(Record[i]==1)){
        MinSize=Ai.nrow();
        index=i;
      }
    }
    //Add index of minimal remaining sub CP-net to Order
    Order[counter]=index;
    //Change corresponding Record entry to 0 so that we know it has been added
    Record[index]=0;
    counter+=1;
  }
  //Order is the order in which the reduced queries should be answered - smallest first
  //nAfter - sum of #outcomes in each CP-net in FinishedCPNs
  long long int nAfter=0;
  for(int i=0;i<FinishedCPNs.size();i++){
    //For each CP-net in FinishedCPNs, DQi:
    List DQi=FinishedCPNs[i];
    //Ni - domain size vector of DQi
    IntegerVector Ni=DQi[1];
    //ni - #outcomes in DQi (product of Ni)
    long long int ni=1;
    for(int j=0;j<Ni.size();j++){
      ni*=Ni[j];
    }
    //Add ni to nAfter sum
    nAfter+=ni;
  }
  //OutCons - counts the number of outcomes we consider in answering the reduced queries
  long long int OutCons=0;
  for(int i=0;i<FinishedCPNs.size();i++){
    //for each reduced CP-net and query:
    //DQi - ith reduced query (in order of size)
    List DQi=FinishedCPNs[Order[i]];
    //Answer the dominance query
    List Answeri=RPSDQRankPriority(DQi[0],DQi[1],DQi[2],DQi[3],DQi[4],DQi[5]);
    //Add #outcomes considered to the count
    OutCons+=Answeri[1];
    if(!Answeri[0]){
      //If the result of the query is false, then the original (input) query is also false
      //return false, #outcomes considered, start/end preprocessing/end times, and the total number of outcomes in the reduced CP-nets
      //log end time 
      time.step("End");
      NumericVector TimeJumps(time);
      List R=List::create(false,OutCons,TimeJumps,nAfter);
      return R;
    }
  }
  //If all reduced querys were answered and found true, the original query must also be true
  //return true, #outcomes considered, start/end preprocessing/end times, and the total number of outcomes in the reduced CP-nets
  //log end time 
  time.step("End");
  NumericVector TimeJumps(time);
  List R=List::create(true,OutCons,TimeJumps,nAfter);
  return R;
}

// [[Rcpp::export]]
List RPSDQRankPriorityTimed(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2){
  //Function takes a dominance query and times how long (the function above) takes to 
  //answer it with RPS and rank priority
  //This outputs the DQ result, outcomes considered, and time taken
  Timer time;
  time.step("Start");
  List Result;
  Result=RPSDQRankPriority( A,  N,  ENTRIES,  BREAK,  o1,  o2);
  time.step("End");
  if(!Result[0]){
    NumericVector TimeJumps(time);
    List R=List::create(false,Result[1],TimeJumps);
    return R;
  }
  NumericVector TimeJumps(time);
  List R=List::create(true,Result[1],TimeJumps);
  return R;
}

// [[Rcpp::export]]
List RPSDQPenaltyPriorityTimed(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2){
  //Function takes a dominance query and times how long (the function above) takes to 
  //answer it with RPS and penalty priority
  //This outputs the DQ result, outcomes considered, and time taken
  Timer time;
  time.step("Start");
  List Result;
  Result=RPSDQPenaltyPriority( A,  N,  ENTRIES,  BREAK,  o1,  o2);
  time.step("End");
  if(!Result[0]){
    NumericVector TimeJumps(time);
    List R=List::create(false,Result[1],TimeJumps);
    return R;
  }
  NumericVector TimeJumps(time);
  List R=List::create(true,Result[1],TimeJumps);
  return R;
}

// [[Rcpp::export]]
List RPSDQRankDiffPriorityTimed(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2){
  //Function takes a dominance query and times how long (the function above) takes to 
  //answer it with RPS and rank+diff priority
  //This outputs the DQ result, outcomes considered, and time taken
  Timer time;
  time.step("Start");
  List Result;
  Result=RPSDQRankDiffPriority( A,  N,  ENTRIES,  BREAK,  o1,  o2);
  time.step("End");
  if(!Result[0]){
    NumericVector TimeJumps(time);
    List R=List::create(false,Result[1],TimeJumps);
    return R;
  }
  NumericVector TimeJumps(time);
  List R=List::create(true,Result[1],TimeJumps);
  return R;
}

// [[Rcpp::export]]
List PSDQTimed(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2){
  //Function takes a dominance query and times how long (the function above) takes to 
  //answer it with PS and penalty priority
  //This outputs the DQ result, outcomes considered, and time taken
  Timer time;
  time.step("Start");
  List Result;
  Result=PSDQ( A,  N,  ENTRIES,  BREAK,  o1,  o2);
  time.step("End");
  if(!Result[0]){
    NumericVector TimeJumps(time);
    List R=List::create(false,Result[1],TimeJumps);
    return R;
  }
  NumericVector TimeJumps(time);
  List R=List::create(true,Result[1],TimeJumps);
  return R;
}

// [[Rcpp::export]]
List RPDQPenaltyPriorityTimed(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2){
  //Function takes a dominance query and times how long (the function above) takes to 
  //answer it with RP and penalty priority
  //This outputs the DQ result, outcomes considered, and time taken
  Timer time;
  time.step("Start");
  List Result;
  Result=RPDQPenaltyPriority( A,  N,  ENTRIES,  BREAK,  o1,  o2);
  time.step("End");
  if(!Result[0]){
    NumericVector TimeJumps(time);
    List R=List::create(false,Result[1],TimeJumps);
    return R;
  }
  NumericVector TimeJumps(time);
  List R=List::create(true,Result[1],TimeJumps);
  return R;
}

// [[Rcpp::export]]
List RPDQRankDiffPriorityTimed(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2){
  //Function takes a dominance query and times how long (the function above) takes to 
  //answer it with RP and rank+diff priority
  //This outputs the DQ result, outcomes considered, and time taken
  Timer time;
  time.step("Start");
  List Result;
  Result=RPDQRankDiffPriority( A,  N,  ENTRIES,  BREAK,  o1,  o2);
  time.step("End");
  if(!Result[0]){
    NumericVector TimeJumps(time);
    List R=List::create(false,Result[1],TimeJumps);
    return R;
  }
  NumericVector TimeJumps(time);
  List R=List::create(true,Result[1],TimeJumps);
  return R;
}

// [[Rcpp::export]]
List RPDQRankPriorityTimed(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2){
  //Function takes a dominance query and times how long (the function above) takes to 
  //answer it with RP and rank priority
  //This outputs the DQ result, outcomes considered, and time taken
  Timer time;
  time.step("Start");
  List Result;
  Result=RPDQRankPriority( A,  N,  ENTRIES,  BREAK,  o1,  o2);
  time.step("End");
  if(!Result[0]){
    NumericVector TimeJumps(time);
    List R=List::create(false,Result[1],TimeJumps);
    return R;
  }
  NumericVector TimeJumps(time);
  List R=List::create(true,Result[1],TimeJumps);
  return R;
}

// [[Rcpp::export]]
List RSDQRankDiffPriorityTimed(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2){
  //Function takes a dominance query and times how long (the function above) takes to 
  //answer it with RS and rank+diff priority
  //This outputs the DQ result, outcomes considered, and time taken
  Timer time;
  time.step("Start");
  List Result;
  Result=RSDQRankDiffPriority( A,  N,  ENTRIES,  BREAK,  o1,  o2);
  time.step("End");
  if(!Result[0]){
    NumericVector TimeJumps(time);
    List R=List::create(false,Result[1],TimeJumps);
    return R;
  }
  NumericVector TimeJumps(time);
  List R=List::create(true,Result[1],TimeJumps);
  return R;
}

// [[Rcpp::export]]
List RSDQRankPriorityTimed(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2){
  //Function takes a dominance query and times how long (the function above) takes to 
  //answer it with RS and rank priority
  //This outputs the DQ result, outcomes considered, and time taken
  Timer time;
  time.step("Start");
  List Result;
  Result=RSDQRankPriority( A,  N,  ENTRIES,  BREAK,  o1,  o2);
  time.step("End");
  if(!Result[0]){
    NumericVector TimeJumps(time);
    List R=List::create(false,Result[1],TimeJumps);
    return R;
  }
  NumericVector TimeJumps(time);
  List R=List::create(true,Result[1],TimeJumps);
  return R;
}

// [[Rcpp::export]]
List SFDQTimed(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2){
  //Function takes a dominance query and times how long (the function above) takes to 
  //answer it with SF
  //This outputs the DQ result, outcomes considered, and time taken
  Timer time;
  time.step("Start");
  List Result;
  Result=SFDQ( A,  N,  ENTRIES,  BREAK,  o1,  o2);
  time.step("End");
  if(!Result[0]){
    NumericVector TimeJumps(time);
    List R=List::create(false,Result[1],TimeJumps);
    return R;
  }
  NumericVector TimeJumps(time);
  List R=List::create(true,Result[1],TimeJumps);
  return R;
}

// [[Rcpp::export]]
List PenaltyDQTimed(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2){
  //Function takes a dominance query and times how long (the function above) takes to 
  //answer it with Penalty
  //This outputs the DQ result, outcomes considered, and time taken
  Timer time;
  time.step("Start");
  List Result;
  Result=PenaltyDQ( A,  N,  ENTRIES,  BREAK,  o1,  o2);
  time.step("End");
  if(!Result[0]){
    NumericVector TimeJumps(time);
    List R=List::create(false,Result[1],TimeJumps);
    return R;
  }
  NumericVector TimeJumps(time);
  List R=List::create(true,Result[1],TimeJumps);
  return R;
}

// [[Rcpp::export]]
List RankDQRankDiffPriorityTimed(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2){
  //Function takes a dominance query and times how long (the function above) takes to 
  //answer it with Rank and Rank+diff priority
  //This outputs the DQ result, outcomes considered, and time taken
  Timer time;
  time.step("Start");
  List Result;
  Result=RankDQRankDiffPriority( A,  N,  ENTRIES,  BREAK,  o1,  o2);
  time.step("End");
  if(!Result[0]){
    NumericVector TimeJumps(time);
    List R=List::create(false,Result[1],TimeJumps);
    return R;
  }
  NumericVector TimeJumps(time);
  List R=List::create(true,Result[1],TimeJumps);
  return R;
}

// [[Rcpp::export]]
List RankDQRankPriorityTimed(IntegerMatrix A, IntegerVector N, IntegerVector ENTRIES, IntegerVector BREAK, IntegerVector o1, IntegerVector o2){
  //Function takes a dominance query and times how long (the function above) takes to 
  //answer it with Rank and Rank priority
  //This outputs the DQ result, outcomes considered, and time taken
  Timer time;
  time.step("Start");
  List Result;
  Result=RankDQRankPriority( A,  N,  ENTRIES,  BREAK,  o1,  o2);
  time.step("End");
  if(!Result[0]){
    NumericVector TimeJumps(time);
    List R=List::create(false,Result[1],TimeJumps);
    return R;
  }
  NumericVector TimeJumps(time);
  List R=List::create(true,Result[1],TimeJumps);
  return R;
}
