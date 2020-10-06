library(Rcpp)
library(primes)
library(Rmpfr)

#CPTConvert takes a CPT(X) in R format and returns the corresponting entries vector (flattened, C++ version) as well as some other helpful vectors for utilising this
CPTConvert<-function(CPT){
  #Number of parents of X
  NPa=length(dim(CPT))-1
  if(NPa==0){
    #If X has 0 parents PA is the vector {0}
    PA=as.vector(c(0))
  }
  else{
    #If X has parents, PA is a vector giving the domain sizes of its parents
    PA=as.vector(dim(CPT)[1:NPa])
  }
  #dom - domain size of X
  dom=dim(CPT)[(NPa+1)]
  #initiate vector entries
  #this will be the CPT flattened as required
  ENTRIES=vector()
  if(NPa==0){
    #If X ha 0 parents, entries is just the single row of entries in CPT
    ENTRIES=as.vector(CPT)
  }
  else{
    #Asst is a matrix, each column gives a parental assignment to Pa(X)
    #These are given in lexicographic order
    Asst=matrix(ncol=NPa,nrow=prod(PA))
    for(i in 1:NPa){
      if(i==NPa){
        a=1
      }
      else{
        a=prod(PA[(i+1):NPa])
      }
      b=prod(PA[i:NPa])
      b=dim(Asst)[1]/b
      Asst[,i]=rep(1:PA[i],each=a,b)
    }
    #For each parental assignment
    for(j in 1:prod(PA)){
      #Asst[j,] - jth parental assignment lexicographically
      #cycle through the possible assignments to X to extract the row from CPT corresponding to this assignment
      #Append this row to entries
      for(k in 1:dom){
        ENTRIES=c(ENTRIES,CPT[matrix(c(Asst[j,],k),1)])
      }
    }
    ENTRIES=as.vector(ENTRIES)
    #Entries is now the rows of CPT end to end, listed in lexicographic order of the associated parental assignment
  }
  #return Entries (flattened CPT(X)), NPa (#parents of X), PA (vector of domain sizes of Pa(X)), and dom (domain size of X)
  return(list(ENTRIES,NPa,PA,dom))
}

#n.val takes a variable, n, and the CPTs of the associated CP-net and returns the domain size of n
n.val<-function(n,CPT){
  #if CPT(n) is one row (i.e. no parents), take the length of this row
  if(length(dim(CPT[n][[1]]))==0){
    length(CPT[n][[1]])
  }
  else{
    #if CPT(n) is a k+1 dimensional array (n has k parents), then by construction, the k+1'th dimension = domain size of n
    dim(CPT[n][[1]])[length(dim(CPT[n][[1]]))]
  }
}

#Takes a CP-net in R form and returns it in C++
CPNConvert<-function(CPN){
  #extract structure, A, and CPTs, CPT, from CPN
  A=CPN[[1]]
  CPT=CPN[[2]]
  #n - number of variables (#rows in A)
  n=dim(A)[1]
  #N- domain sizes of all variables
  N=c()
  for(i in 1:n){
    #append domain size of variable i:
    N=c(N,n.val(i,CPT))
  }
  #ENTRIES vector for full CP-net will be the flattened version of each CPT(X) (obtained via CPTConvert) in order
  ENTRIES=c()
  #Breaks vector starts with value 1
  BREAKS=c(1)
  for(j in 1:n){
    #For each variable, j:
    #ENT - flattened version of CPT(j)
    ENT=CPTConvert(CPT[[j]])[[1]]
    #Append ENT to entries
    ENTRIES=c(ENTRIES,ENT)
    #L=length(ENT)- size of the CPT just added to ENTRIES
    #next entry in breaks will be L greater than the last, indicating that the next break point will be L indices later
    l=length(BREAKS)
    BREAKS=c(BREAKS,(BREAKS[l]+length(ENT)))
  }
  #ENTRIES and BREAKS now complete
  #Return the CP-net in its new form
  A=as.matrix(A)
  N=as.vector(N)
  ENTRIES=as.vector(ENTRIES)
  BREAKS=as.vector(BREAKS)
  return(list(A,N,ENTRIES,BREAKS))
}


