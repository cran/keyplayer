# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

##################################################################################
### Node-level measure
### mreach.degree
#' Compute the M-reach Degree Centrality Score in a Netwrok
#'
#' \code{mreach.degree} computes the size of reachable set within M.
#' M-reach degree centrality generalizes the \code{\link[sna]{degree}} centrality
#' by delimiting specific neighborhoods.
#'
#' @param adj.matrix Matrix indicating the adjacency matrix of the network.
#'
#' @param node Integer indicating the column index of the chosen player
#' in the adjacenncy matrix. If not specified, scores for all nodes will be reported.
#'
#' @param M Number indicating the maximum geodistance between two nodes,
#' above which the two nodes are considered disconnected.
#' M hence defines the reachable set. The default is \code{Inf}.
#'
#' @param binary If \code{TRUE}, the adjacency matrix is binarized, and thus M
#' essentially means steps. If \code{FALSE}, the edge values are considered.
#'
#' @param cmode String indicating the type of centrality being evaluated.
#' \code{"outdegree"}, \code{"indegree"}, and \code{"total"} refer to
#' indegree, outdegree, and (total) degree respectively. \code{"all"} reports
#' all the above measures and is the default.
#'
#' @details The interprtation of the measure in binary and weighted adjacency matrix
#' are slightly different, while the function allows both. User can choose to binarize
#' the adjacency matrix, and in this case the reachable set is defined by nodes that
#' are reachable within M steps. Hence those more directly connected to the chosen set
#' form the reachable set. If weighted adjacency matrix is used, the reachable set
#' is defined by nodes that are reachable within geodistance M. Hence those with
#' shorter distance to the chosen set form the reachable set.
#'
#' It is true that the binarized setting is more general as only the network structure
#' matters there and the interpretation of weights is irrelevent.
#' Hence, by default \code{binary=TRUE}. See An and Liu (2015) for more details.
#'
#' @return A vector indicating the outdegree, indegree, or total-degree
#' mreach.degree score of the chosen node; or a data frame containing all
#' the above information.
#'
#'
#' @author Weihua An \email{weihuaan@@indiana.edu}; Yu-Hsin Liu \email{yuhsliu@@indiana.edu}
#'
#' @references
#' An, Weihua and Yu-Hsin Liu (2015). "keyplayer: An R Package for Locating Key Players in Social Networks."
#' \emph{Working Paper}, Indiana Univeristy.\cr
#'
#' @examples
#' # Create a 5x5 weighted and directed adjacency matrix,
#' # where edge values represent the strength of tie
#' W <- matrix(
#'   c(0,1,3,0,0,
#'     0,0,0,4,0,
#'     1,1,0,2,0,
#'     0,0,0,0,3,
#'     0,2,0,0,0),
#'     nrow=5, ncol=5, byrow = TRUE)
#'
#'
#' # List the 2-reach degree scores for every node where W is binarized
#' mreach.degree(W,M=2,cmode="all")
#'
#' @seealso
#' \code{\link[sna]{geodist}};
#' \code{\link{mreach.closeness}};
#' \code{\link{kpcent}};
#' \code{\link{kpset}}
#'
#' @export
#'
mreach.degree=function(adj.matrix, node, M = Inf, binary = TRUE, cmode = "all"){
  d <- data.matrix(adj.matrix, rownames.force = NA)
  rownames(d) <- c(colnames(d))

  # create geodistance matrix
  distances=sna::geodist(d,ignore.eval=binary)$gdist
  colnames(distances) <- c(colnames(d))
  rownames(distances) <- c(colnames(d))

  # set M constraint using distances matrix
  if(missing(M)){
    distances = distances
  }else{
    for(i in 1:ncol(distances)){
      for(j in 1:nrow(distances)){
        if(distances[i,j]>M) distances[i,j]=Inf
      }
    }
  }
  # calculate the score
  diag(distances)=Inf
  weights=1/distances
  weights[weights>0]=1
  outdegree=apply(weights,1,sum) # column vector, outdegree
  indegree=apply(weights,2,sum) # row vector, indegree
  total= outdegree + indegree

  if(missing(node)){
    kppos=cbind(outdegree, indegree, total)
  }else{
    kppos=c(outdegree[node], indegree[node], total[node])
    names(kppos) <- c("outdegree","indegree","total")
    kppos=as.data.frame(t(kppos))
    rownames(kppos) <- "score"
  }
  # reports
  if (cmode=="outdegree"){
    outdegree[node]
  }else if (cmode=="indegree"){
    indegree[node]
  }else if (cmode=="total"){
    total[node]
  }else{
    kppos
  }
}

###############################################################################
### mreach.closeness
#' Compute the M-reach Closeness Centrality Score in a Netwrok
#'
#' \code{mreach.closeness} refines the \code{\link{mreach.degree}} centrality by
#' using the (inverse) geodistance as weights.
#' The edge values should be properly interpreted as distances.
#'
#' @param adj.matrix Matrix indicating the adjacency matrix of the network.
#'
#' @param node Integer indicating the column index of the chosen player
#' in the adjacenncy matrix. If not specified, scores for all nodes will be reported.
#'
#' @param M Number indicating the maximum geodistance between two nodes,
#' above witch the two nodes are considered disconnected.
#' M hence defines the reachable set. The default is \code{Inf}.
#'
#' @param binary If \code{TRUE}, the adjacency matrix is binarized.
#' If \code{FALSE}, the edge values are considered.
#'
#' @param cmode String indicating the type of centrality being evaluated.
#' \code{"outdegree"}, \code{"indegree"}, and \code{"total"} refer to
#' indegree, outdegree, and (total) degree respectively. \code{"all"} reports
#' all the above measures and is the default.
#'
#' @details \code{mreach.closeness} refines the \code{\link{mreach.degree}} centrality
#' by using the (inverse) geodistance as weights, just as \code{\link[sna]{closeness}}
#' centrality refines \code{\link[sna]{degree}} centrality.
#' The delimiting definition of neighborhoods given by M induces several flexible
#' properties of the M-reach closeness centrality. It captures
#' the degree centrality when M is properly set (e.g. M=1 in a binarized network).
#' It captures the Gil-Schmidt power index (Gil and Schmidt, 1996)
#' and the cohesion centrality (Borgatti, 2006) when M is sufficiently large
#' (unconstrained). The normalization factor takes care of non-binary
#' edge values interpreted as distances. Also note that the geodistance matrix does
#' not necessarily to be symmetric, and thus the measure is directed.
#' see An and Liu (2015) for more details.
#'
#' @return A vector indicating the outdegree, indegree, or total-degree
#' cohesion score of the chosen player; or a data frame containing all
#' the above information. Note that we normalize the outdegree and indegree scores
#' to [0,1]. This means that the total-degree score is between [0,2].
#'
#' @author Weihua An \email{weihuaan@@indiana.edu}; Yu-Hsin Liu \email{yuhsliu@@indiana.edu}
#'
#' @references
#' An, Weihua and Yu-Hsin Liu (2015). "keyplayer: An R Package for Locating Key Players in Social Networks."
#' \emph{Working Paper}, Indiana Univeristy.\cr
#'
#' Borgatti, Stephen P. (2006). "Identifying Sets of Key Players in a Network."
#' \emph{Computational, Mathematical and Organizational Theory}, 12(1):21-34.\cr
#'
#' Gil, J and Schmidt, S (1996). "The Origin of the Mexican Network of Power."
#' Proceedings of the International Social Network Conference, Charleston, SC, 22-25.\cr
#'
#' @examples
#' # Create a 5x5 weighted and directed adjacency matrix, where edge values
#' # represent the strength of tie
#' W <- matrix(
#'   c(0,1,3,0,0,
#'     0,0,0,4,0,
#'     1,1,0,2,0,
#'     0,0,0,0,3,
#'     0,2,0,0,0),
#'     nrow=5, ncol=5, byrow = TRUE)
#'
#' # Transform the edge value to distance interpretaion
#' A <- W
#' A[W!=0] <- 1/W[W!=0]
#'
#' # List all types of 2-reach closeness scores for every node
#' mreach.closeness(A,M=2,cmode="all")
#'
#' @seealso
#' \code{\link[sna]{geodist}};
#' \code{\link{kpcent}};
#' \code{\link{kpset}}
#'
#'
#' @export

mreach.closeness=function(adj.matrix, node, M = Inf, binary = FALSE, cmode = "all"){
  d <- data.matrix(adj.matrix, rownames.force = NA)
  rownames(d) <- c(colnames(d))
  # create matrix of geodesic distances
  distances=sna::geodist(d,ignore.eval=binary)$gdist
  colnames(distances) <- c(colnames(d))
  rownames(distances) <- c(colnames(d))

  # set threshold constraint using distances matrix
  if(missing(M)){
    distances = distances
  }else{
    for(i in 1:ncol(distances)){
      for(j in 1:nrow(distances)){
        if(distances[i,j]>M) distances[i,j]=Inf
      }
    }
  }

  # calculate the score
  diag(distances)=Inf
  weights=1/distances
  sum.out=apply(weights,1,sum) # column vector, outdegree
  sum.in=apply(weights,2,sum) #  row vector, indegree
  outdegree = sum.out/((ncol(distances)-1)*max(weights)) # normalize to [0,1]
  indegree = sum.in/((ncol(distances)-1)*max(weights)) # normalize to [0,1]
  total=outdegree + indegree
  if(missing(node)){
    kppos=cbind(outdegree, indegree, total)
  }else{
    kppos=c(outdegree[node], indegree[node],total[node])
    names(kppos) <- c("outdegree","indegree","total")
    kppos=as.data.frame(t(kppos))
    rownames(kppos) <- "score"
  }
  # reports
  if (cmode=="outdegree"){
    outdegree[node]
  }else if (cmode=="indegree"){
    indegree[node]
  }else if (cmode=="total"){
    total[node]
  }else{
    kppos
  }
}

#####################################################################################
### fragment

#' Compute the Fragmentation Centrality Score in a Netwrok
#'
#' \code{fragment} measures the extent of fragmentation of a network after a
#' set of nodes is removed from the network. The more the fregmentation level of
#' the residual network is, the more central those nodes are.
#'
#' @param adj.matrix Matrix indicating the adjacency matrix of the network.
#'
#' @param nodes Integer indicating the column index of the chosen player
#' in the adjacenncy matrix. If there are multiple players,
#' use \code{c(index1,index2,...)}.
#' If not specified, scores for all nodes will be reported.
#'
#' @param M Number indicating the maximum geodistance between two nodes,
#' above witch the two nodes are considered disconnected.
#' M hence defines the reachable set. The default is \code{Inf}.
#'
#' @param binary If \code{TRUE}, the adjacency matrix is binarized.
#' If \code{FALSE}, the edge values are considered. The default is \code{FALSE}.
#'
#' @details A natural way to apply the fragmentation centrality is in the
#' context of counter-terrorism, as originally proposed in Borgatti (2006).
#' The measure uses geodistances to compute the fragmentation level of the
#' residual network, and thus edge values should be properly adjusted to
#' distance interpretation. The fragmentation centrality is not directional
#' as edge values are counted aggregately in a network level. \code{fragment}
#' keeps the flexible features of defining reachable set using M, discribed in
#' \code{\link{mreach.degree}} and \code{\link{mreach.closeness}}.
#
#'
#' @return Vector indicating fragment score(s) of the chosen player(s).
#' Score is normalized to [0,1].
#'
#' @author Weihua An \email{weihuaan@@indiana.edu}; Yu-Hsin Liu \email{yuhsliu@@indiana.edu}
#'
#' @references
#' Borgatti, Stephen P. 2006. "Identifying Sets of Key Players in a Network."
#' \emph{Computational, Mathematical and Organizational Theory}, 12(1):21-34.\cr
#'
#' @examples
#' # Create a 5x5 weighted and directed adjacency matrix, where edge values
#' # represent the strength of tie
#' W <- matrix(
#'   c(0,1,3,0,0,
#'     0,0,0,4,0,
#'     1,1,0,2,0,
#'     0,0,0,0,3,
#'     0,2,0,0,0),
#'     nrow=5, ncol=5, byrow = TRUE)
#'
#' # Transform the edge value to distance interpretaion
#' A <- W
#' A[W!=0] <- 1/W[W!=0]
#'
#' # List the fragmentation centrality scores for every node
#' fragment(A)
#'
#' @seealso
#' \code{\link[sna]{geodist}};
#' \code{\link{kpcent}};
#' \code{\link{kpset}}
#'
#'
#' @export

fragment=function(adj.matrix, nodes, M=Inf, binary=FALSE){
  d <- data.matrix(adj.matrix, rownames.force = NA)
  # geodesic distances
  distances=sna::geodist(d,ignore.eval = binary)$gdist
  diag(distances)=Inf # replace the diagonal with Infs
  weights=1/distances # take the reciprocal of distances
  m=max(weights)

  if(missing(nodes)){
    s <- c(1:ncol(adj.matrix))
    for(k in 1:ncol(adj.matrix)){
      d <- data.matrix(adj.matrix, rownames.force = NA)
      d <- d[-k,-k]
      # create matrix of geodesic distances
      distances=sna::geodist(d,ignore.eval = binary)$gdist
      for(i in 1:ncol(distances)){
        for(j in 1:nrow(distances)){
          if(distances[i,j]>M) distances[i,j]=Inf
        }
      }
      # cauculate the score
      diag(distances)=Inf
      weights=1/distances
      sum=sum(weights) # sum for the matrix elements, sum is a scalar
      s[k]=1-sum/(ncol(distances)*(ncol(distances)-1)*m)
    }
    s<-t(s)
    rownames(s)<-"fragment"
    t(s)
  }else{
    # remove objective nodes
    ColToDelete <- nodes
    d <- d[,-ColToDelete]
    RowToDelete <- t(nodes)
    d <- d[-RowToDelete,]

    # create matrix of geodesic distances
    distances=sna::geodist(d,ignore.eval = binary)$gdist

    # set threshold using distances matrix
    if(missing(M)){
      distances = distances
    }else{
      for(i in 1:ncol(distances)){
        for(j in 1:nrow(distances)){
          if(distances[i,j]>M) distances[i,j]=Inf
        }
      }
    }
    # cauculate the score
    diag(distances)=Inf
    weights=1/distances
    sum=sum(weights) # sum for the matrix elements, sum is a scalar
    1-sum/(ncol(distances)*(ncol(distances)-1)*m)
  }
}


################################################################################
### diffusion
#' Compute the Diffusion Centrality Score in a Network
#'
#' \code{diffusion} measures player's ability to disseminate information when
#' the process is modeled random. The adjacency matrix should carry this randomness,
#' where edge value p_{ij} represents the probability that information is passed from
#' i to j.
#'
#' @param adj.matrix Matrix indicating the adjacency matrix of the network.
#' It is suggested that the inputted adjacency matrix for the diffusion centrality
#' be properly transfomred to the probability interpretation.
#'
#' @param node Integer indicating the column index of the chosen player
#' in the adjacenncy matrix. If not specified, scores for all nodes will be reported.
#'
#' @param T Integer indicating the maximum number of iterations
#' of communication process. In the first iteration, the adjacency matrix
#' is as the input. In the nth iteration, the adjacency matrix becomes
#' the input adjacency matrix to the power of n. By default, T is the network size.
#'
#' @return A vector indicating the defusion centrality score(s) of
#' the chosen player(s).
#'
#' @details The diffusion centrality is developed by Banerjee et.al. (2013),
#' where it finds that this fairly general measure approximates the
#' empirical result of word-of-mouth diffusion well.
#' It is general because it can represents (or approximates) the degree centrality,
#' eigenvector centrality, and Katz-Bonacich centrality, in different settings.
#' See the Banerjee et.al. (2014) for proofs.
#'
#'
#' @author Weihua An \email{weihuaan@@indiana.edu}; Yu-Hsin Liu \email{yuhsliu@@indiana.edu}
#'
#' @references
#' Banerjee, A., A. Chandrasekhar, E. Duflo, and M. Jackson (2013):
#' "Diffusion of Microfinance," \emph{Science}, Vol. 341. p.363\cr
#'
#' Banerjee, A., A. Chandrasekhar, E. Duflo, and M. Jackson (2014):
#' "Gossip: Identifying Central Individuals in a Social Network,"
#' Working paper\cr
#'
#' @examples
#' # Create a 5x5 weighted and directed adjacency matrix, where edge values
#' # represent the strength of tie
#' W <- matrix(
#'   c(0,1,3,0,0,
#'     0,0,0,4,0,
#'     1,1,0,2,0,
#'     0,0,0,0,3,
#'     0,2,0,0,0),
#'     nrow=5, ncol=5, byrow = TRUE)
#'
#' # Transform the edge value to probability interpretaion
#' P <- W *0.2
#'
#' # List the diffusion centrality score for every node
#' diffusion(P, T = 2)
#'
#' @seealso
#' \code{\link{kpcent}};
#' \code{\link{kpset}}
#'
#'
#'
#' @export


diffusion=function(adj.matrix, node, T=ncol(adj.matrix)){
  d <- data.matrix(adj.matrix, rownames.force = NA)
  s <- 0
  for (t in 1:T){
    s=s+matpow::matpow(d,t)$prod1
  }
  # sum for each node (row), sum is a column vector, outdegree
  diffusion=apply(s,1,sum)
  if(missing(node)){
    score <- t(diffusion)
    rownames(score)<-"diffusion"
    t(score)
  }else{
    diffusion[node]
  }
}


##################################################################################
##################################################################################
### Group-level measure (also able to measure the node level)
### contract function
#' Group the Chosen Players in a Network
#'
#' \code{contract} is a fundmental tool for analyzing the group-level centrality.
#' It contracts the chosen nodes into one following a specified rule, and thus
#' the node-level method can be applied. There are four grouping criteria
#' in the \code{contract} function: \code{minimum}, \code{maximum}, \code{union},
#' and \code{addition},
#' each of them responses to a specific interpretation of the edge values
#' of the adjacency matrix. See details.
#'
#' @param adj.matrix Matrix indicating the adjacency matrix of the network.
#'
#' @param nodes Integer indicating the column index of the chosen player
#' in the adjacenncy matrix. Input \code{c(index1,index2,...)} for multiple players
#'
#' @param method Indication of which grouping criterion should be used.
#' \code{method="min"} indicates the "minimum" criterion (edge values as distances).
#' \code{method="max"} indicates the "maximum" criterion (edge values as non-cummulative strengths).
#' \code{method="add"} indicates the "addition" criterion (edge values as cummulative strengths).
#' \code{method="union"} indicates the "union" criterion (edge values as probability).
#' The default uses the "minimum" criterion as the distances interpretation is standard in \pkg{sna}.
#'
#' @details When the edge values are interpreted as distances, the contracting
#' follows a "minimum" rule that the distance (edge value) between an objective
#' and the contracting set is the minimun distance of each individual node in
#' the set. For example, suppose A to C has distance 2 and B to C has distance 1,
#' then the distance between C and the merged set AB is 1.
#'
#' When the edge values are interpreted as accumulative strength, the contracting follows
#' an "addition" rule that the strength (edge value) of tie between an objective
#' node and the contracting set is simply the summation of the edge value of
#' each node in the set. In the above example, when the strength interpretation
#' is used, the strangth of tie between C and the merged AB is 1+2=3. If the strength
#' is not accumulative, then the "maximum" rule should be applied.
#'
#' When the edge values are interpreted as probability, the contracting follows
#' an "union" rule that the link to the set can be formed when at least one of
#' the nodes in the set is connected. The rule assumes that the link to each
#' node in the set is independent to each other and thus computes the probability
#' using 1 minus the product of the probability of not being connected for each
#' node. For example, suppose A has probability 0.2 to reach C and B has probability
#' 0.5 to reach C, then C can be reached from merged AB with probability 1-(1-0.2)*(1-0.5)=0.6.
#'
#' @return A new adjacency matrix after contracting the chosen nodes (named
#' \code{set}).
#'
#' @author Weihua An \email{weihuaan@@indiana.edu}; Yu-Hsin Liu \email{yuhsliu@@indiana.edu}
#'
#' @references
#' Butts, Carter T. (2014). sna: Tools for Social Network Analysis. R package
#' version 2.3-2. \url{http://CRAN.R-project.org/package=sna}\cr
#'
#'
#' @keywords internal

contract=function(adj.matrix, nodes, method=c("min","max","union","add")){
  d <- data.matrix(adj.matrix, rownames.force = NA)
  if (missing(method)){
    # contract the objective nodes; row sum for outdegree; column sum for indegree
    if (length(nodes)==1){
      set <- d[,nodes]
      set.r <- d[nodes,]
    }else{
      max <- max(d)+1
      d <- replace(d,d==0,max)
      set <- d[,ncol(d)]
      for(i in 1:nrow(d)){
        set[i]=min(d[i,nodes])
      }
      set.r <- d[ncol(d),]
      for(j in 1:ncol(d)){
        set.r[j]=min(d[nodes,j])
      }
      d <- replace(d,d==max,0)
      set <- replace(set,set==max,0)
      set.r <- replace(set.r,set.r==max,0)
    }
  }else if (method=="add"){
    if (length(nodes)==1){
      set <- d[,nodes]
      set.r <- d[nodes,]
    }else{
      set <- rowSums(d[,nodes])
      set.r <- colSums(d[nodes,])
    }
  }else if (method=="union"){
    if (length(nodes)==1){
      set <- d[,nodes]
      set.r <- d[nodes,]
    }else{
      set <- rowSums(d[,nodes])
      for(i in 1:nrow(d)){
        set[i]=1-prod(1-d[i,nodes])
      }
      set.r <- colSums(d[nodes,])
      for(j in 1:ncol(d)){
        set.r[j]=1-prod(1-d[nodes,j])
      }
    }
  }else if(method=="max"){
    if (length(nodes)==1){
      set <- d[,nodes]
      set.r <- d[nodes,]
    }else{
      set <- d[,ncol(d)]
      for(i in 1:nrow(d)){
        set[i]=max(d[i,nodes])
      }

      set.r <- d[ncol(d),]
      for(j in 1:ncol(d)){
        set.r[j]=max(d[nodes,j])
      }
    }
  }else{
    if (length(nodes)==1){
      set <- d[,nodes]
      set.r <- d[nodes,]
    }else{
      max <- max(d)+1
      d <- replace(d,d==0,max)
      set <- d[,ncol(d)]
      for(i in 1:nrow(d)){
        set[i]=min(d[i,nodes])
      }
      set.r <- d[ncol(d),]
      for(j in 1:ncol(d)){
        set.r[j]=min(d[nodes,j])
      }
      d <- replace(d,d==max,0)
      set <- replace(set,set==max,0)
      set.r <- replace(set.r,set.r==max,0)
    }
  }

  # update the adj.matrix
  set.r <- t(set.r)
  set.r <- cbind(set.r,0)
  d <- cbind(d,set)
  d <- rbind(d,set.r)
  rownames(d)=c(colnames(d))
  ColToDelete <- nodes
  d <- d[,-ColToDelete] # if v is names, use d <- d[,!(names(x) %in% ColToDelete)]
  RowToDelete <- t(nodes)
  d <- d[-RowToDelete,] # if v is names, use d <- d[!(names(x) %in% RowToDelete),]
  d
}


##################################################################################
### group.mreach.degree
#' Compute the group-level mreach.degree Centrality Score in a Netwrok
#'
#' \code{group.mreach.degree} computes the size of the reachable set within M.
#'
#' @param adj.matrix Matrix indicating the adjacency matrix of the network.
#'
#' @param nodes Integer indicating the column index of the chosen player
#' in the adjacenncy matrix. If there are multiple players,
#' use \code{c(index1,index2,...)}
#'
#' @param M Number indicating the maximum distance between two nodes,
#' above which the two nodes are considered disconnected.
#' m hence defines a reachable set. The default is \code{Inf}.
#'
#' @param method Indication of which grouping criterion should be used.
#' \code{method="min"} indicates the "minimum" criterion (edge values as distances).
#' \code{method="max"} indicates the "maximum" criterion (edge values as non-cummulative strengths).
#' \code{method="add"} indicates the "addition" criterion (edge values as cummulative strengths).
#' \code{method="union"} indicates the "union" criterion (edge values as probability).
#' The default is the "minimum" criterion for mreach.degree centrality.
#'
#' @param binary If \code{TRUE}, the adjacency matrix is binarized.
#' If \code{FALSE}, the edge values are considered.
#'
#' @param cmode String indicating the type of centrality being evaluated.
#' The default is to report the total degree.
#' \code{"outdegree"} and \code{"indegree"} refer to indegree and outdegree
#' respectively. If \code{"all"}, all the three types are reported.
#'
#' @return A vector indicating the outdegree, indegree, or total-degree
#' mreach.degree score of the chosen player(s); or a data frame containing all
#' the above information.
#'
#' @author Weihua An \email{weihuaan@@indiana.edu}; Yu-Hsin Liu \email{yuhsliu@@indiana.edu}
#'
#' @keywords internal

group.mreach.degree=function(adj.matrix, nodes, M = Inf, method = "min", binary = TRUE, cmode = "total"){
  d <- data.matrix(adj.matrix, rownames.force = NA)
  rownames(d) <- c(colnames(d))

  if(NCOL(d[,nodes])==1){
    # create matrix of geodesic distances
    distances=sna::geodist(d,ignore.eval=binary)$gdist
    colnames(distances) <- c(colnames(d))
    rownames(distances) <- c(colnames(d))

    # set step constraint using distances matrix
    if(missing(M)){
      distances = distances
    }else{
      for(i in 1:ncol(distances)){
        for(j in 1:nrow(distances)){
          if(distances[i,j]>M) distances[i,j]=Inf
        }
      }
    }

    # calculate the score
    diag(distances)=Inf
    weights=1/distances
    weights[weights>0]=1
    kppos.out=apply(weights,1,sum) # column vector, outdegree
    kppos.in=apply(weights,2,sum) # row vector, indegree
    kppos.total=kppos.out + kppos.in
    kppos=c(kppos.out[nodes], kppos.in[nodes],kppos.total[nodes])
    names(kppos) <- c("outdegree","indegree","total")
    kppos=as.data.frame(t(kppos))
    rownames(kppos) <- "score"
    # reports
    if (missing(cmode)){
      kppos.total[nodes]
    }else if (cmode=="outdegree"){
      kppos.out[nodes]
    }else if (cmode=="indegree"){
      kppos.in[nodes]
    }else if (cmode=="all"){
      kppos
    }else{
      kppos.total[nodes]
    }
  }else{
    # geodesic distances
    distances=sna::geodist(d,ignore.eval=binary)$gdist
    diag(distances)=Inf
    weights=1/distances
    d=contract(d,nodes,method)

    # create matrix of geodesic distances
    distances=sna::geodist(d,ignore.eval=binary)$gdist
    colnames(distances) <- c(colnames(d))
    rownames(distances) <- c(colnames(d))

    # set step using distances matrix
    if(missing(M)){
      distances = distances
    }else{
      for(i in 1:ncol(distances)){
        for(j in 1:nrow(distances)){
          if(distances[i,j]>M) distances[i,j]=Inf
        }
      }
    }

    # calculate the score
    diag(distances)=Inf
    weights=1/distances
    weights[weights>0]=1
    kppos.out=apply(weights,1,sum) # column vector, outdegree
    kppos.in=apply(weights,2,sum) # row vector, indegree
    kppos.total=kppos.out + kppos.in
    kppos=c(kppos.out["set"], kppos.in["set"], kppos.total["set"])
    names(kppos) <- c("outdegree","indegree","total")
    kppos=as.data.frame(t(kppos))
    rownames(kppos) <- "score"
    if (missing(cmode)){
      kppos.total["set"]
    }else if (cmode=="outdegree"){
      kppos.out["set"]
    }else if (cmode=="indegree"){
      kppos.in["set"]
    }else if (cmode=="all"){
      kppos
    }else{
      kppos.total["set"]
    }
  }
}

##################################################################################
### group.mreach.closeness
#' Compute the Group-level mreach.closeness Centrality Score in a Netwrok
#'
#' \code{mreach.closeness} refines the \code{\link{mreach.degree}} centrality by
#' using the (inverse) geodistance as weights.
#' The edge values should be properly interpreted as distances.
#'
#' @param adj.matrix Matrix indicating the adjacency matrix of the network.
#'
#' @param nodes Integer indicating the column index of the chosen player
#' in the adjacenncy matrix. If there are multiple players,
#' use \code{c(index1,index2,...)}
#'
#' @param M Number indicating the maximum distance between two nodes,
#' above witch the two nodes are considered disconnected. The default is
#' \code{Inf}.
#'
#' @param method Indication of which grouping criterion should be used.
#' \code{method="min"} indicates the "minimum" criterion (edge values as distances).
#' \code{method="max"} indicates the "maximum" criterion (edge values as non-cummulative strengths).
#' \code{method="add"} indicates the "addition" criterion (edge values as cummulative strengths).
#' \code{method="union"} indicates the "union" criterion (edge values as probability).
#' The default is the "minimum" criterion for mrach.closeness centrality.
#'
#' @param binary If \code{TRUE}, the adjacency matrix is binarized.
#' If \code{FALSE}, the edge values are considered.
#'
#' @param cmode String indicating the type of centrality being evaluated.
#' The default is to report the total degree.
#' \code{"outdegree"} and \code{"indegree"} refer to indegree and outdegree
#' respectively. If \code{"all"}, all the three types are reported.
#'
#' @return A vector indicating the outdegree, indegree, or total-degree
#' cohesion score of the chosen players; or a data frame containing all
#' the above information. Note that we normalize the outdegree and indegree scores
#' to [0,1]. This means that the total-degree score is between [0,2].
#'
#' @author Weihua An \email{weihuaan@@indiana.edu}; Yu-Hsin Liu \email{yuhsliu@@indiana.edu}
#'
#' @keywords internal


group.mreach.closeness=function(adj.matrix, nodes, M = Inf, method = "min", binary = FALSE, cmode = "total"){
  d <- data.matrix(adj.matrix, rownames.force = NA)
  rownames(d) <- c(colnames(d))

  if(NCOL(d[,nodes])==1){
    # create matrix of geodesic distances
    distances=sna::geodist(d,ignore.eval=binary)$gdist
    colnames(distances) <- c(colnames(d))
    rownames(distances) <- c(colnames(d))

    # set threshold constraint using distances matrix
    if(missing(M)){
      distances = distances
    }else{
      for(i in 1:ncol(distances)){
        for(j in 1:nrow(distances)){
          if(distances[i,j]>M) distances[i,j]=Inf
        }
      }
    }

    # calculate the score
    diag(distances)=Inf
    weights=1/distances
    sum.out=apply(weights,1,sum) # column vector, outdegree
    sum.in=apply(weights,2,sum) # row vector, indegree
    kppos.out = sum.out/((ncol(distances)-1)*max(weights)) # normalize to [0,1]
    kppos.in = sum.in/((ncol(distances)-1)*max(weights)) # normalize to [0,1]
    kppos.total=kppos.out + kppos.in
    kppos=c(kppos.out[nodes], kppos.in[nodes],kppos.total[nodes])
    names(kppos) <- c("outdegree","indegree","total")
    kppos=as.data.frame(t(kppos))
    rownames(kppos) <- "score"
    # reports
    if (missing(cmode)){
      kppos.total[nodes]
    }else if (cmode=="outdegree"){
      kppos.out[nodes]
    }else if (cmode=="indegree"){
      kppos.in[nodes]
    }else if (cmode=="all"){
      kppos
    }else{
      kppos.total[nodes]
    }
  }else{
    # geodesic distances
    distances=sna::geodist(d,ignore.eval=binary)$gdist
    diag(distances)=Inf
    weights=1/distances
    m=max(weights)

    d=contract(d,nodes,method)

    # create matrix of geodesic distances
    distances=sna::geodist(d,ignore.eval=binary)$gdist
    colnames(distances) <- c(colnames(d))
    rownames(distances) <- c(colnames(d))

    # set threshold using distances matrix
    if(missing(M)){
      distances = distances
    }else{
      for(i in 1:ncol(distances)){
        for(j in 1:nrow(distances)){
          if(distances[i,j]>M) distances[i,j]=Inf
        }
      }
    }

    # calculate the score
    diag(distances)=Inf
    weights=1/distances
    sum.out=apply(weights,1,sum) # column vector, outdegree
    sum.in=apply(weights,2,sum) # row vector, indegree
    kppos.out = sum.out/((ncol(distances)-1)*m) # normalize to [0,1]
    kppos.in = sum.in/((ncol(distances)-1)*m) # normalize to [0,1]
    kppos.total=kppos.out + kppos.in
    kppos=c(kppos.out["set"], kppos.in["set"], kppos.total["set"])
    names(kppos) <- c("outdegree","indegree","total")
    kppos=as.data.frame(t(kppos))
    rownames(kppos) <- "score"
    if (missing(cmode)){
      kppos.total["set"]
    }else if (cmode=="outdegree"){
      kppos.out["set"]
    }else if (cmode=="indegree"){
      kppos.in["set"]
    }else if (cmode=="all"){
      kppos
    }else{
      kppos.total["set"]
    }
  }
}


###################################################################################
### kpcent

#' Compute Group Centraltiy in a Network
#'
#' \code{kpcent} reports the group-level centrality scores of the specified
#' type of centrality measure.
#'
#' @param adj.matrix Matrix indicating the adjacency matrix of the network.
#'
#' @param nodes Integer indicating the column index of the chosen player
#' in the adjacenncy matrix. If there are multiple players,
#' use \code{c(index1,index2,...)}
#'
#' @param type
#' \code{type="betweenness"} for \code{\link[sna]{betweenness}} centrality. \cr
#' \code{type="closeness"} for \code{\link[sna]{closeness}} centrality. \cr
#' \code{type="degree"} for \code{\link[sna]{degree}} centraslity. \cr
#' \code{type="diffusion"} for \code{\link{diffusion}} centrality. \cr
#' \code{type="evcent"} for \code{\link[sna]{evcent}} (eigenvector) centrality. \cr
#' \code{type="fragment"} for \code{\link{fragment}} centrality. \cr
#' \code{type="mreach.degree"} for \code{\link{mreach.degree}} centrality. \cr
#' \code{type="mreach.closeness"} for \code{\link{mreach.closeness}} centrality. \cr
#'
#'
#' @param method Indication of which grouping criterion should be used. \cr
#' \code{"min"} indicates the "minimum" criterion and is the default for
#' betweenness, closeness, fragmentation, and M-reach centralities. \cr
#' \code{"max"} indicates the "maximum" criterion and is the default for
#' degree and eigenvector centralities.\cr
#' \code{"add"} indicates the "addition" criterion.\cr
#' \code{"union"} indicates the "union" criterion and is the default for
#' diffusion centrality.\cr
#' See Details section for explanations on grouping method.
#'
#' @param M Positive number indicating the maximum geodistance between two nodes,
#' above witch the two nodes are considered disconnected. The default is
#' \code{Inf}. The option is applicable to mreach.degree, mreach.closeness,
#' and fragmentation centralities.
#'
#' @param T Integer indicating the maximum number of iterations
#' of communication process. For diffusion centrality only.
#' In the first iteration, the adjacency matrix
#' is as the input. In the nth iteration, the adjacency matrix becomes
#' the input adjacency matrix to the power of n. By default, T is the network size.
#'
#' @param binary If \code{TRUE}, the adjacency matrix is binarized.
#' If \code{FALSE}, the edge values are considered. By default, \code{binary=FALSE}
#'
#' @param cmode String indicating the type of centrality being evaluated.
#' The option is applicable to degree and M-reach centralities.
#' \code{"outdegree"}, \code{"indegree"}, and \code{"total"} refer to
#' indegree, outdegree, and (total) degree respectively. \code{"all"} reports
#' all the above measures. The default is to report the total degree.
#' The option also applies to closeness centrality, but with different options.
#' The default is to use the Gil-Schmidt power index as the closeness measure.
#' See \code{\link[sna]{closeness}} for complete options.
#'
#' @return A vector indicating the centrality score or a data frame containing
#' the summery information of the directed centrality scores.
#'
#' @details The basic idea of measuring the group-level centrality is to treat
#' a group of nodes as a large pseudo-node. We propose several methods to measure
#' the tie status between this pseudo node and other nodes, responding to several
#' common edge value interpretations (An and Liu, 2015).
#'
#' Minimum Criterion: the edge value between a group and an outside node
#' is measured as the minimal value among all the (nonzero) edge values between
#' any node in the group and the outside node. Suggested if edge values are
#' interpreted as distances.\cr
#' \emph{Example: suppose node A to C has distance 2 and B to C has distance 1,
#' then according to the minimum criterion, the distance between C and
#' the merged set AB is 1. Note that if B and C are not connected,
#' the algorithm takes the distance between A and C to describe
#' the distance between AB and C.}
#'
#' Maximun Criterion: the edge value between a group and an outside node
#' is measured as the maximal value among all the (nonzero) edge values between
#' any node in the group and the outside node. Suggested if edge values are
#' interpreted as non-cummulative strengths. \cr
#' \emph{Example: we keep using the above example, but the figure now indicates
#' the strength of tie. According to the maximum criterion, the strength of tie
#' between AB and C is 2.}
#'
#' Addition Criterion: the edge value between a group and an outside node
#' is measured as the sum of all the edge values between any node in the group
#' and the outside node. Suggested if edge values are as cummulative strengths. \cr
#' \emph{Example: according to the addition criterion, the strength of tie between
#' AB and C is 3}
#'
#' Union Criterion: the edge value between a group and an outside node is measured
#' as the probability that there is at least one path connecting the group with
#' the outside node. Suggested if edge values are as probability. \cr
#' \emph{Example: suppose A has probability 0.2 to reach C and B has probability
#' 0.5 to reach C, then C can be reached from merged AB with probability
#' 1-(1-0.2)*(1-0.5)=0.6 according to the union criterion.}
#'
#' @author Weihua An \email{weihuaan@@indiana.edu}; Yu-Hsin Liu \email{yuhsliu@@indiana.edu}
#'
#' @references
#' An, Weihua. (2015). "Multilevel Meta Network Analysis with Application to Studying Network Dynamics of Network Interventions." \emph{Social Networks} 43: 48-56.\cr
#'
#' An, Weihua and Yu-Hsin Liu (2015). "keyplayer: An R Package for Locating Key Players in Social Networks."
#' \emph{Working Paper}, Indiana Univeristy.\cr
#'
#' Banerjee, A., A. Chandrasekhar, E. Duflo, and M. Jackson (2013):
#' "Diffusion of Microfinance," \emph{Science}, Vol. 341. p.363\cr
#'
#' Banerjee, A., A. Chandrasekhar, E. Duflo, and M. Jackson (2014):
#' "Gossip: Identifying Central Individuals in a Social Network,"
#' Working paper\cr
#'
#' Borgatti, Stephen P. (2006). "Identifying Sets of Key Players in a Network."
#' \emph{Computational, Mathematical and Organizational Theory}, 12(1):21-34.\cr
#'
#' Butts, Carter T. (2014). sna: Tools for Social Network Analysis. R package
#' version 2.3-2. \url{http://CRAN.R-project.org/package=sna}\cr
#'
#' @seealso
#' \code{\link{kpset}}
#'
#' @examples
#' # Create a 5x5 weighted and directed adjacency matrix,
#' # where edge values represent the strength of tie
#' W <- matrix(
#'   c(0,1,3,0,0,
#'     0,0,0,4,0,
#'     1,1,0,2,0,
#'     0,0,0,0,3,
#'     0,2,0,0,0),
#'     nrow=5, ncol=5, byrow = TRUE)
#'
#' # List the degree centrality for group of node 2 and 3
#' kpcent(W,c(2,3),type="degree")
#'
#' # Transform the edge value to distance interpretaion
#' # Compute the fragmentation centrality for node 2
#' A <- W
#' A[W!=0] <- 1/W[W!=0]
#' kpcent(A,2,type="fragment")
#'
#' # Replicate the group-level degree centrality (normalized) when the weights
#' # are given by the inverse distances and report the outgoing score only
#' kpcent(A,c(2,3),type="mreach.closeness",binary=TRUE,M=1,cmode="outdegree")
#'
#' # Transform the edge value to probability interpretation
#' # Compute the diffusion centrality with number of iteration 20
#' P <- 0.1*W
#' kpcent(P,c(2,3),type="diffusion",T=20)
#'
#' @export


kpcent=function(adj.matrix, nodes, type, M=Inf, T=ncol(adj.matrix), method, binary=FALSE, cmode){
  d <- data.matrix(adj.matrix, rownames.force = NA)

  if(type == "betweenness"){
    d <- data.matrix(adj.matrix, rownames.force = NA)
    if (length(nodes)==1){
      betweenness = sna::betweenness(d, nodes = nodes, ignore.eval = binary, rescale = FALSE)
    }else{
      if(missing(method)){
        d <- contract(d, nodes, method="min")
      }else{
        d <- contract(d, nodes, method)
      }
      s = ncol(d)
      betweenness <- sna::betweenness(d, nodes = s, ignore.eval = binary, rescale = FALSE)
    }
    score <- betweenness

  }else if(type == "closeness"){
    if(missing(cmode)){
      if (length(nodes)==1){
        closeness = sna::closeness(d, nodes = nodes, cmode="suminvdir", ignore.eval = binary)
      }else{
        if(missing(method)){
          d <- contract(d, nodes, method="min")
        }else{
          d <- contract(d, nodes, method)
        }
        s = ncol(d)
        closeness <- sna::closeness(d, nodes = s, cmode="suminvdir", ignore.eval = binary)
      }
    }else{
      if (length(nodes)==1){
        closeness = sna::closeness(d, nodes = nodes, cmode, ignore.eval = binary)
      }else{
        if(missing(method)){
          d <- contract(d, nodes, method="min")
        }else{
          d <- contract(d, nodes, method)
        }
        s = ncol(d)
        closeness <- sna::closeness(d, nodes = s, cmode, ignore.eval = binary)
      }
    }
    score <- closeness

  }else if(type == "degree"){
    if(length(nodes)==1){
      indegree = sna::degree(d, nodes = nodes, cmode = "indegree", ignore.eval = binary)
      outdegree = sna::degree(d, nodes = nodes, cmode = "outdegree", ignore.eval = binary)
      total = sna::degree(d, nodes = nodes, cmode = "freeman", ignore.eval = binary)
    }else{
      if (missing(method)){
        d <- contract(d,nodes, method="max")
      }else{
        d <- contract(d, nodes, method)
      }
      s = ncol(d)
      indegree = sna::degree(d, nodes = s, cmode = "indegree", ignore.eval = binary)
      outdegree = sna::degree(d, nodes = s, cmode = "outdegree", ignore.eval = binary)
      total = sna::degree(d, nodes = s, cmode = "freeman", ignore.eval = binary)
    }

    set.degree=c(outdegree, indegree, total)
    names(set.degree) <- c("outdegree","indegree","total")
    set.degree=as.data.frame(t(set.degree))
    rownames(set.degree) <- "score"

    if (missing(cmode)){
      score <- total
    }else if (cmode=="outdegree"){
      score <- outdegree
    }else if (cmode=="indegree"){
      score <- indegree
    }else if (cmode=="all"){
      score <- set.degree
    }else{
      score <- total
    }
  }else if(type == "diffusion"){
    s <- 0
    if(NCOL(d[,nodes])==1){
      for (t in 1:T){
        s=s+matpow::matpow(d,t)$prod1
      }
      diffusion=apply(s,1,sum)
      score <- diffusion[nodes]
    }else{
      if(missing(method)){
        d <- contract(d, nodes, method="union")
      }else{
        d <- contract(d, nodes, method)
      }
      s=0
      for (t in 1:T){
        s=s+matpow::matpow(d,t)$prod1
      }
      diffusion=apply(s,1,sum)
      score <- diffusion["set"]
    }
  }else if(type == "evcent"){
    if(length(nodes)==1){
      evcent = sna::evcent(d, nodes = nodes, ignore.eval = binary)
    }else{
      if(missing(method)){
        d <- contract(d, nodes, method="max")
      }else{
        d <- contract(d, nodes, method)
      }
      s = ncol(d)
      evcent <- sna::evcent(d, nodes = s, ignore.eval = binary)
    }
    score <- evcent
  }else if(type == "fragment"){
    score <- fragment(adj.matrix, nodes, M, binary)
  }else if(type == "mreach.degree"){
    if(missing(method)){
      score <- group.mreach.degree(adj.matrix, nodes, M, method="min", binary, cmode)
    }else{
      score <- group.mreach.degree(adj.matrix, nodes, M, method, binary, cmode)
    }
  }else if(type == "mreach.closeness"){
    if(missing(method)){
      score <- group.mreach.closeness(adj.matrix, nodes, M, method="min", binary, cmode)
    }else{
      score <- group.mreach.closeness(adj.matrix, nodes, M, method, binary, cmode)
    }
  }
  score
}

####################################################################################
### kpset
### Step 1. create delta.score to compute the change from switching i to j

#' The change on objective function for greedy search implimentation
#'
#' \code{delta.score} calculates the base kp score using \code{kpcent} with
#' the specified chosen set of players called candidate. Then the function
#' replaces a player from another specified set of plyers called residual and
#' calculates the new kp score. The function finally reports the difference of
#' the two scores.
#' @param adj.matrix The adjacency matrix of the network. Input the matrix.
#' @param candidate A specified set of players which centrality is measured.
#' @param residual A specified set of players which the member can replace the
#' member in the candidate set.
#' @param i The specific member in candidate set to be replaced
#' @param j The specific member in residual set to replace \code{i} in candidate set
#' @param type Choose
#' \code{type="betweenness"} for betweenness centrality,
#' \code{type="closeness"} for closeness centrality,
#' \code{type="degree"} for degree centraslity,
#' \code{type="diffusion"} for diffusion centrality.
#' \code{type="evcent"} for eigenvector centrality,
#' \code{type="fragment"} for fragmentation centrality,
#' \code{type="mreach.degree"} for mreach.degree centrality, and
#' \code{type="mreach.closeness"} for mreach.closeness centrality.
#'
#' @param M Positive number indicating the maximum distance between two nodes,
#' above witch the two nodes are considered disconnected. The default is
#' \code{Inf}.
#' The option is applicable to mreach.degree, mreach.closeness, fragmentation,
#' and diffusion centralities.
#'
#' @param T Integer indicating the maximum number of iterations
#' of communication process. For diffusion centrality only.
#' In the first iteration, the adjacency matrix
#' is as the input. In the nth iteration, the adjacency matrix becomes
#' the input adjacency matrix to the power of n. By default, T is the network size.
#'
#' @param method Indication of which grouping criterion should be used.
#' \code{method="min"} indicates the "minimum" criterion (edge values as distances).
#' \code{method="max"} indicates the "maximum" criterion (edge values as non-cummulative strengths).
#' \code{method="add"} indicates the "addition" criterion (edge values as cummulative strengths).
#' \code{method="union"} indicates the "union" criterion (edge values as probability).
#' By default, the minimun criterion is used for betweenness, closeness, fragmentation,
#' mreach.degree, and mreach.closeness centralities.
#' The maximun criterion is used for degree and eigenvector centralities.
#' The union criterion is used for diffusion centrality.
#'
#' @param binary If \code{TRUE}, the adjacency matrix is binarized.
#' If \code{FALSE}, the edge values are considered. By default, \code{binary=FALSE}
#'
#' @param cmode String indicating the type of centrality being evaluated.
#' The option is applicable to degree, mreach.degree, and mreach.closeness centralities.
#' The default is to report the total degree.
#' \code{cmode="outdegree"} and \code{cmode="indegree"} refer to indegree and outdegree
#' respectively. If \code{cmode="all"}, all the three types are reported.
#' The option can also applicable to closeness centrality.
#' See \code{\link[sna]{closeness}} Details section.
#' The default is to use the Gil-Schmidt power index.
#'
#' @return The change on kp score
#' @keywords internal

delta.score=function(adj.matrix,candidate,residual,i,j,type,M=Inf,T=ncol(adj.matrix),method,binary=FALSE,cmode){

  s <- kpcent(adj.matrix,candidate,type,M,T,method,binary,cmode)
  candidate[i] <- residual[j]
  delta.s = kpcent(adj.matrix,candidate,type,M,T,method,binary,cmode)-s
  delta.s
}

################################################################################
### Step 2. keyplayers greedy search algorithm

#' Selecting the Most Central Group of Players in a Network
#'
#' \code{kpset} implements a greedy search algorithm to find the most
#' central players given the sepcified centraliy measure and the target group size.
#'
#' @param adj.matrix Matrix indicating the adjacency matrix of the network.
#'
#' @param size Integer indicating the target size of players.
#'
#' @param type
#' \code{type="betweenness"} for \code{\link[sna]{betweenness}} centrality. \cr
#' \code{type="closeness"} for \code{\link[sna]{closeness}} centrality. \cr
#' \code{type="degree"} for \code{\link[sna]{degree}} centraslity. \cr
#' \code{type="diffusion"} for \code{\link{diffusion}} centrality. \cr
#' \code{type="evcent"} for \code{\link[sna]{evcent}} (eigenvector) centrality. \cr
#' \code{type="fragment"} for \code{\link{fragment}} centrality. \cr
#' \code{type="mreach.degree"} for \code{\link{mreach.degree}} centrality. \cr
#' \code{type="mreach.closeness"} for \code{\link{mreach.closeness}} centrality. \cr
#'
#' @param method Indication of which grouping criterion should be used. \cr
#' \code{"min"} indicates the "minimum" criterion and is the default for
#' betweenness, closeness, fragmentation, and M-reach centralities. \cr
#' \code{"max"} indicates the "maximum" criterion and is the default for
#' degree and eigenvector centralities.\cr
#' \code{"add"} indicates the "addition" criterion.\cr
#' \code{"union"} indicates the "union" criterion and is the default for
#' diffusion centrality.\cr
#' See \code{\link{kpcent}} Details section for explanations on grouping method.
#'
#' @param M Positive number indicating the maximum geodistance between two nodes,
#' above witch the two nodes are considered disconnected. The default is
#' \code{Inf}. The option is applicable to mreach.degree, mreach.closeness,
#' and fragmentation centralities.
#'
#' @param T Integer indicating the maximum number of iterations
#' of communication process. For diffusion centrality only.
#' In the first iteration, the adjacency matrix
#' is as the input. In the nth iteration, the adjacency matrix becomes
#' the input adjacency matrix to the power of n. By default, T is the network size.
#'
#' @param binary If \code{TRUE}, the adjacency matrix is binarized.
#' If \code{FALSE}, the edge values are considered. By default, \code{binary=FALSE}
#'
#' @param cmode String indicating the type of centrality being evaluated.
#' The option is applicable to degree and M-reach centralities.
#' \code{"outdegree"}, \code{"indegree"}, and \code{"total"} refer to
#' indegree, outdegree, and (total) degree respectively. \code{"all"} reports
#' all the above measures. The default is to report the total degree.
#' The option also applies to closeness centrality, but with different options.
#' The default is to use the Gil-Schmidt power index as the closeness measure.
#' See \code{\link[sna]{closeness}} for complete options.
#'
#' @param iteration The upper boound of the number of search amount. See
#' Details section for the definition.
#'
#' @return \code{kpset} returns the column indices of the players who form
#' the most central set and its centrality score.
#'
#' @details
#' The most central individuals together do not necessarily form the most
#' central set in a netowrk as they may share overlapped connections. Hence,
#' a greedy search algorithm is implemented here to identify the key set of
#' players. The basic logic is as the following:
#' \enumerate{
#' \item Randomly select a set of nodes \emph{C}. The remaining set of nodes is denoted as \emph{R}.
#' \item Update the selected set \emph{C}.
#'   \enumerate{
#'   \item Replace node \emph{i} in \emph{C} with node \emph{j} in \emph{R}.
#'         Replace \emph{i} by \emph{j} if the \emph{j}th replacement improves
#'         the centrality of \emph{C} the most. (loop 1)
#'   \item Repeat loop 1 for each node in \emph{C} (loop 2)
#'   \item Stop if (a) the change in \emph{C}'s centrality score is smaller than
#'         a specified threshold or (b) the process reaches the specified number
#'         of iterations (i.e., the number of loop 2).
#'   }
#' \item Return the final candidate set and the centrality score.
#' }
#'
#' See empirical examples in An and Liu (2015).
#'
#' @author Weihua An \email{weihuaan@@indiana.edu}; Yu-Hsin Liu \email{yuhsliu@@indiana.edu}
#'
#' @references
#' An, Weihua. (2015). "Multilevel Meta Network Analysis with Application to Studying Network Dynamics of Network Interventions." \emph{Social Networks} 43: 48-56.\cr
#'
#' An, Weihua and Yu-Hsin Liu (2015). "keyplayer: An R Package for Locating Key Players in Social Networks."
#' \emph{Working Paper}, Indiana Univeristy.\cr
#'
#' Borgatti, Stephen P. (2006). "Identifying Sets of Key Players in a Network." \emph{Computational, Mathematical and Organizational Theory}, 12(1):21-34.\cr
#'
#' Butts, Carter T. (2014). sna: Tools for Social Network Analysis. R package
#' version 2.3-2. \url{http://CRAN.R-project.org/package=sna}\cr
#'
#' @examples
#' # Create a 5x5 weighted and directed adjacency matrix
#' W <- matrix(
#'   c(0,1,3,0,0,
#'     0,0,0,4,0,
#'     1,1,0,2,0,
#'     0,0,0,0,3,
#'     0,2,0,0,0),
#'     nrow=5, ncol=5, byrow = TRUE)
#'
#' # Find the most central player set sized 2 in terms of the degree centrality
#' kpset(W,size=2,type="degree")
#'
#' @export

kpset=function(adj.matrix, size, type, M=Inf, T=ncol(adj.matrix), method, binary=FALSE, cmode, iteration=1000){
  d <- data.matrix(adj.matrix, rownames.force = NA)

  N <- c(1:NCOL(d))
  C <- sample.int(NCOL(d),size)
  R <- N[!N %in% C]

  storage <- c(0:iteration)
  d1 = d
  d1[d1==0] = Inf
  m <- max(max(d),1/min(d1))
  for (t in 2:(iteration+1)){
    for (i in 1:NROW(C)){
      for (j in 1:NROW(R)){
        if(delta.score(adj.matrix,C,R,i,j,type,M,T,method,binary,cmode)>0) C[i]<-R[j]
      }
      R <- N[!N %in% C]
    }
    storage[t] <- kpcent(adj.matrix,C,type,M,T,method,binary,cmode)
    if (storage[t]-storage[t-1]<(1/(2*(NCOL(d)^2)*(m^3)))) break
  }

  call<-match.call()
  col1 <- C
  col1 <- sort(col1)
  col2 <- kpcent(adj.matrix,C,type,M,T,method,binary,cmode)
  names(col2) <- NULL
  result <- list(keyplayers=col1, centrality=col2)
  result
}
