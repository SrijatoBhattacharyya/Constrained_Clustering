
# rm(list=ls())
library(data.tree)
library(igraph)
# source("ComplexDomainFun.R")
source("ComplexDomain_utils.R")
library(vegan)
library(cccd)
library(salso)
# source("FEMFun.R")
library(igraph)





packages = c("deldir", "ggplot2", "igraph","cccd","fields",'plotly')

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = FALSE)
      library(x, character.only = TRUE)
    }else
      library(x, character.only = TRUE)
  }
)

### function for forest birth

forest_birth = function(graph0, subgraphs, trees, cluster){
  csize <- sapply(subgraphs, igraph::vcount)
  k <- length(csize)
  
  # Randomly select subgraph index for splitting wp proportional to cluster size
  
  clust_split = sample.int(k, 1, prob = csize - 1)
  submst = trees[[clust_split]]
  
  chosen_subgraph = subgraphs[[clust_split]]
  chosen_tree = submst
  
  edge_cutted = E(submst)[sample.int(csize[clust_split]-1, 1)]
  submst = delete_edges(submst, edge_cutted)

  connect_comp = components(submst)
  idx_new = (connect_comp$membership == 2)
  vid_new = V(submst)$vid[idx_new]
  vid_old = V(submst)$vid[!idx_new]
  
  trees[[clust_split]] <- igraph::subgraph(chosen_tree, which(!idx_new))
  trees[[k+1]]         <- igraph::subgraph(chosen_tree, which(idx_new))
  
  csize[clust_split] = length(vid_old)
  csize[k+1] = length(vid_new)
  
  cluster[vid_new] = k + 1
  
  edge_list <- ends(graph0, E(graph0), names = FALSE)
  
  # For each edge, check if one endpoint is in vid1 and the other in vid2.
  is_cross_edge <- (edge_list[, 1] %in% vid_old & edge_list[, 2] %in% vid_new) |
    (edge_list[, 1] %in% vid_new & edge_list[, 2] %in% vid_old)
  
  # Count the number of cross edges.
  num_cross_edges <- sum(is_cross_edge)
  
  birth_prior_ratio_term = log(num_cross_edges)
  subgraphs = create_subgraphs_from_clusters(graph0, cluster)
  
  
  return(list(subgraphs = subgraphs,
              cluster = cluster,
              csize = csize,
              trees = trees,
              chosen_subgraph = chosen_subgraph,
              birth_prior_ratio_term = birth_prior_ratio_term,
              vid1 = vid_old,
              vid2 = vid_new))
  
}

### function for forest death

forest_death <- function(graph0, subgraphs, trees, cluster) {
  # 0) ensure every edge has an eid
  if (is.null(E(graph0)$eid)) {
    E(graph0)$eid <- seq_len(ecount(graph0))
  }
  
  # 1) map internal index → your vid
  vid_map <- V(graph0)$vid
  
  # 2) get every graph0‐edge as (u,v) in internal indices
  el_idx = as_edgelist(graph0, names = FALSE)
  clu_u <- cluster[ el_idx[,1] ]
  clu_v <- cluster[ el_idx[,2] ]
  
  # 3) find edges that cross *any* two different clusters
  diff_edges <- which(clu_u != clu_v)
  if (length(diff_edges)==0) stop("No inter‑cluster edges found.")
  
  # 4) build the unique candidate cluster‐pairs and sample one
  cpairs <- t(apply(cbind(clu_u[diff_edges], clu_v[diff_edges]), 1, sort))
  uniq   <- unique(cpairs)
  if (is.null(nrow(uniq))) uniq <- matrix(uniq, nrow=1)
  chosen <- uniq[sample(nrow(uniq),1), ]
  cl_new <- chosen[1]  # will survive
  cl_old <- chosen[2]  # will be merged away
  
  # 5) grab vids in each cluster
  vids_new <- V(subgraphs[[cl_new]])$vid
  vids_old <- V(subgraphs[[cl_old]])$vid
  merged_vids <- c(vids_new, vids_old)
  
  # 6) full induced subgraph for MH bookkeeping
  merged_subgraph <- subgraph(graph0, which(vid_map %in% merged_vids))
  
  # 7) compute death‐prior ratio = –log(#cross‐edges between these two parts)
  idx_new <- which(vid_map %in% vids_new)
  idx_old <- which(vid_map %in% vids_old)
  is_cross_subset <- (el_idx[,1] %in% idx_old & el_idx[,2] %in% idx_new) |
    (el_idx[,1] %in% idx_new & el_idx[,2] %in% idx_old)
  cross_rows <- which(is_cross_subset)
  if (length(cross_rows)==0) stop("No crossing edges between these two clusters—shouldn’t happen!")
  death_prior_ratio_term <- -log(length(cross_rows))
  
  # 8) pick one of *those* specific crossing‐edges
  sel_global <- cross_rows[sample(length(cross_rows),1)]
  sel_eid    <- E(graph0)[sel_global]$eid
  
  # 9) collect the eids of the two old trees
  eids_old <- E(trees[[cl_old]])$eid
  eids_new <- E(trees[[cl_new]])$eid
  
  # 10) build the merged‐tree *only* from those three sets of edges
  merged_tree <- subgraph.edges(
    graph0,
    eids            = c(eids_old, eids_new, sel_eid),
    delete.vertices = TRUE
  )
  
  # 11) splice in the new tree, drop the old
  trees[[cl_new]] <- merged_tree
  trees[[cl_old]] <- NULL
  
  # 12) relabel your cluster‐vector
  cluster[cluster == cl_old] <- cl_new
  cluster[cluster >  cl_old] <- cluster[cluster > cl_old] - 1
  
  # 13) rebuild subgraphs & sizes
  subgraphs <- create_subgraphs_from_clusters(graph0, cluster)
  csize     <- sapply(subgraphs, vcount)
  
  # 14) return
  list(
    subgraphs              = subgraphs,
    cluster                = cluster,
    csize                  = csize,
    trees                  = trees,
    chosen_subgraph        = merged_subgraph,
    vid1                   = vids_old,
    vid2                   = vids_new,
    death_prior_ratio_term = death_prior_ratio_term
  )
}


library(igraph)

countCrossEdges <- function(graph0, vid1, vid2) {
  # Count edges in graph0 connecting a vertex in vid1 to a vertex in vid2.
  cross_edges <- E(graph0)[ (.from(vid1) & .to(vid2)) | (.from(vid2) & .to(vid1)) ]
  return(length(cross_edges))
}

library(igraph)

# Function to sample a uniform spanning tree (UST) from a given graph.
sampleUST <- function(g) {
  # Sample a spanning tree from g.
  tree_ust <- sample_spanning_tree(g)
  
  # If the result is an edge sequence (class "igraph.es"), convert to a graph.
  if (!is.igraph(tree_ust)) {
    tree_ust <- subgraph.edges(g, tree_ust, delete.vertices = FALSE)
  }
  return(tree_ust)
}

# Function to draw a UST from each subgraph in the subgraphs list.
drawUSTs <- function(subgraphs) {
  # Initialize an empty list to store the spanning trees.
  ust_list <- vector("list", length(subgraphs))
  
  # Loop over the list of subgraphs.
  for (i in seq_along(subgraphs)) {
    # Get the current subgraph.
    g <- subgraphs[[i]]
    
    # Sample a uniform spanning tree from the subgraph.
    ust_list[[i]] <- sampleUST(g)
  }
  
  return(ust_list)
}






#### function to compute log likelihood given clustering
cluster_log_density <- function(Y, sigmasq_mu, sigmasq_y, cluster) {
  # ensure same length
  if (length(Y) != length(cluster)) stop("Y and cluster must match.")
  
  total_ld <- 0
  # get unique cluster labels
  for (cl in unique(cluster)) {
    y_sub <- Y[cluster == cl]
    n_j   <- length(y_sub)
    
    # closed-form log det
    log_det <- (n_j - 1)*log(sigmasq_y) +
      log(sigmasq_y + n_j*sigmasq_mu)
    
    # quadratic form: y'y / sigmasq_y  - (sigmasq_mu/(sigmasq_y*(sigmasq_y + n_j*sigmasq_mu))) * (sum(y_sub))^2
    sum_y   <- sum(y_sub)
    sum_y2  <- sum(y_sub^2)
    quad    <- (sum_y2 / sigmasq_y) -
      (sigmasq_mu * sum_y^2) /
      (sigmasq_y * (sigmasq_y + n_j*sigmasq_mu))
    
    # accumulate log-density
    total_ld <- total_ld +
      ( -n_j/2 * log(2*pi)
        - 0.5 * log_det
        - 0.5 * quad )
  }
  
  total_ld
}


contract_to_quotient_graph <- function(graph, cluster) {
  # Load igraph (if not already loaded)
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("The 'igraph' package is required but not installed.")
  }
  
  # Check that the cluster vector has the same length as the number of vertices
  # if (length(cluster) != igraph::vcount(graph)) {
  #   stop("The length of the cluster vector must be equal to the number of vertices in the graph.")
  # }
  
  # Contract the graph using the cluster assignment vector.
  # This merges vertices based on cluster membership.
  graph_contracted <- igraph::contract(
    graph,
    mapping = cluster,
    vertex.attr.comb = list(name = function(x) paste(unique(x), collapse = ","))
  )
  
  # Simplify the contracted graph to remove self-loops and merge multiple edges.
  # For edge attributes such as 'weight', we can sum the weights.
  quotient_graph <- igraph::simplify(
    graph_contracted,
    remove.loops = TRUE,
    edge.attr.comb = list(weight = "sum", "ignore")
  )
  
  # Optionally, rename the vertices of the quotient graph using unique cluster labels.
  unique_clusters <- sort(unique(cluster))
  igraph::V(quotient_graph)$name <- as.character(unique_clusters)
  
  return(quotient_graph)
}



log_numSpanningTrees <- function(g) {
  
  # L = graph.laplacian(g, sparse = T)
  
  L <- laplacian_matrix(g, sparse = TRUE)
  
  L = as(L[-1,-1, drop = F], "symmetricMatrix") # make symmetric improves speed a lot
  
  out = as.numeric(Matrix::determinant(L, logarithm = T)$modulus) # Matrix::determinant  
  return(out)
}
# 


count_connected_pairs <- function(G, subgraphs) {
  # Create a membership vector for the vertices in G, 
  # indexed by the vertex IDs (vid) rather than names.
  membership <- rep(NA, vcount(G))
  # Assume each vertex in G has a "vid" attribute.
  names(membership) <- as.character(V(G)$vid)
  
  # Assign cluster labels to vertices based on their presence in the subgraphs.
  for (i in seq_along(subgraphs)) {
    # Extract the vertex ids ("vid") from this subgraph.
    vids_sub <- V(subgraphs[[i]])$vid  
    # Use these vids (converted to character for indexing) to set membership.
    membership[ as.character(vids_sub) ] <- i
  }
  
  # Initialize a vector to record unique inter-cluster pairs.
  unique_pairs <- c()
  
  # Define a small helper: Given an edge e in G, return the pair of vid values.
  get_edge_vids <- function(e) {
    # ends(G, e) returns a character vector of vertex indices if names are absent.
    # We use as.numeric() to convert these indices into numbers,
    # then use them to fetch the "vid" attribute from the graph.
    indices <- as.numeric(ends(G, e, names = FALSE))
    # Return the vid values.
    return(V(G)$vid[indices])
  }
  
  # Loop over all edges in G.
  for (e in E(G)) {
    vpair <- get_edge_vids(e)
    # Lookup cluster assignments using the "vid" as a key (converted to character).
    clust1 <- membership[ as.character(vpair[1]) ]
    clust2 <- membership[ as.character(vpair[2]) ]
    
    # If both endpoints have cluster assignments and they are in different clusters,
    # record this pair (using a sorted string "i_j" so that order doesn't matter).
    if (!is.na(clust1) && !is.na(clust2) && clust1 != clust2) {
      pair <- paste(sort(c(clust1, clust2)), collapse = "_")
      unique_pairs <- c(unique_pairs, pair)
    }
  }
  
  # Remove duplicate pairs.
  unique_pairs <- unique(unique_pairs)
  
  # Optionally, you could return both the count and the pairs.
  # For now, just return the count.
  return(length(unique_pairs))
}



library(igraph)

# Efficient edge‐classification for a single membership vector
get_cluster_edges <- function(graph, membership) {
  # Number of edges
  m <- ecount(graph)
  # “eids” is simply 1:m
  eids <- seq_len(m)
  # ends_mat[i,] = the two endpoints of edge i (as numeric vertex indices)
  ends_mat <- ends(graph, eids, names = FALSE)
  
  # A logical vector: TRUE if both ends of edge i are in the same cluster
  same_cluster <- membership[ends_mat[,1]] == membership[ends_mat[,2]]
  
  list(
    within  = eids[same_cluster],
    between = eids[!same_cluster]
  )
}

create_subgraphs_from_clusters <- function(graph0, cluster) {
  # 1. sanity‐check
  if (length(cluster) != igraph::vcount(graph0)) {
    stop("cluster vector length must equal number of vertices in graph0")
  }
  
  # 2. sorted labels
  clusters  <- sort(unique(cluster))
  
  # 3. preallocate named list
  subgraphs <- vector("list", length(clusters))
  names(subgraphs) <- as.character(clusters)
  
  # 4. fill in each slot by position
  for (cl in clusters) {
    vids <- which(cluster == cl)
    # <— use subgraph(), not induced_subgraph(), to dodge the impl arg
    subgraphs[[ as.character(cl) ]] <- igraph::subgraph(graph0, vids)
  }
  
  return(subgraphs)
}



sum_log_trees_subgraphs <- function(subgraphs_list) {
  # Initialize the sum of log number of spanning trees
  total_log_num_trees <- 0
  
  # Loop through each subgraph and sum their log number of spanning trees
  for (subgraph in subgraphs_list) {
    total_log_num_trees <- total_log_num_trees + log_numSpanningTrees(subgraph)
  }
  
  return(total_log_num_trees)
}





# function to get whether an edge is within a cluster or bewteen two clusters

getEdgeStatus <- function(membership, inc_mat) {
  membership_head = membership[inc_mat[, 1]]
  membership_tail = membership[inc_mat[, 2]]
  edge_status = rep('w', nrow(inc_mat))
  edge_status[membership_head != membership_tail] = 'b'
  return(edge_status)
}


# function to split an existing cluster given MST
# subgraphs: cluster_id -> subgraph
# cluster: vid -> cluster_id

splitCluster_common <- function(mstgraph, k, subgraphs, csize) { 
  clust_split = sample.int(k, 1, prob = csize - 1)
  mst_subgraph = subgraphs[[clust_split]]
  
  edge_cutted = E(mst_subgraph)[sample.int(csize[clust_split]-1, 1)]
  eid_cutted = edge_cutted$eid
  mst_subgraph = igraph::delete.edges(mst_subgraph, edge_cutted)
  connect_comp = igraph::components(mst_subgraph)
  idx_new = (connect_comp$membership == 2)
  vid_new = V(mst_subgraph)$vid[idx_new]
  vid_old = V(mst_subgraph)$vid[!idx_new]
  
  return(list(vid_old = vid_old, vid_new = vid_new, eid_cutted = eid_cutted,
              clust_old = clust_split, idx_new = idx_new))
}


# function to update if a split move is accepted
updateSplit_shannon <- function(split_res, subgraphs, k, csize, csize_z, eid_btw_mst, cluster, edge_status,
                                adj_list, adj_edge_list) {
  clust_split = split_res$clust_old
  vid_old = split_res$vid_old; vid_new = split_res$vid_new
  
  subgraph_split = subgraphs[[clust_split]]
  idx_new = split_res$idx_new
  subgraphs[[clust_split]] = igraph::induced_subgraph(subgraph_split, !idx_new)  # subgraph of old cluster
  subgraphs[[k+1]] = igraph::induced_subgraph(subgraph_split, idx_new) # subgraph of new cluster
  
  csize[clust_split] = length(vid_old)
  csize[k+1] = length(vid_new)
  
  cluster[vid_new] = k + 1
  eid_btw_mst = c(eid_btw_mst, split_res$eid_cutted)
  eid_btw_mst_new = split_res$eid_cutted
  # update edge status
  adj_vid_old = unlist(adj_list[vid_old])
  adj_eid_old = unlist(adj_edge_list[vid_old])
  clust_adj_old = cluster[adj_vid_old]
  idx_btw = which(clust_adj_old != clust_split)
  eid_btw = adj_eid_old[idx_btw]
  edge_status[eid_btw] = 'b'
  
  return(list(subgraphs = subgraphs, csize = csize, cluster = cluster,
              eid_btw_mst = eid_btw_mst, estatus = edge_status, eid_btw_mst_new = eid_btw_mst_new))
}



# function to merge two existing clusters
mergeCluster_shannon <- function(mstgraph, eid_btw_mst, subgraphs, csize, csize_z, cluster, edge_list, 
                                 change = F) {
  # edge for merging
  edge_merge = sample.int(length(eid_btw_mst), 1)
  # update cluster information
  # clusters of endpoints of edge_merge
  eid_merge = eid_btw_mst[edge_merge]
  # eid_merge = eid_btw_mst1
  
  clusters_merge = cluster[edge_list[eid_merge, ]]
  clusters_merge = sort(clusters_merge)
  c1 = clusters_merge[1]; c2 = clusters_merge[2] # note c1 < c2
  eid_btw_mst1 = eid_btw_mst[-edge_merge]
  # merge c2 to c1
  
  # vid of vertices in c2
  vid_old = igraph::V(subgraphs[[c2]])$vid
  # vid in merged cluster
  vid_new = c(igraph::V(subgraphs[[c1]])$vid, vid_old)
  
  csize_new = NULL; subgraphs_new = NULL
  if(change) {
    subgraphs_new = subgraphs
    subgraphs_new[[c1]] = igraph::induced_subgraph(mstgraph, vid_new)
    subgraphs_new[[c2]] = NULL
    
    csize_new = csize
    csize_new[c1] = length(vid_new)
    csize_new = csize_new[-c2]
  }
  
  # now drop c2
  return(list(vid_old = vid_old, vid_new = vid_new, clust_old = c2, clust_new = c1,
              edge_merge = edge_merge, subgraphs = subgraphs_new, csize = csize_new, 
              eid_btw_mst1 = eid_btw_mst1, eid_merge = eid_merge))
}


# function to update if a merge move is accepted
updateMerge_shannon <- function(res_merge, subgraphs, csize, csize_z, eid_btw_mst, cluster, edge_status,
                                adj_list, adj_edge_list, mstgraph) {
  clust_old = res_merge$clust_old; clust_new = res_merge$clust_new
  vid_old = igraph::V(subgraphs[[clust_old]])$vid
  vid_new = c(igraph::V(subgraphs[[clust_new]])$vid, vid_old)
  subgraphs[[clust_new]] = igraph::induced_subgraph(mstgraph, vid_new)
  subgraphs[[clust_old]] = NULL
  
  csize[clust_new] = length(vid_new)
  csize = csize[-clust_old]
  
  cluster[vid_old] = clust_new
  idx = which(cluster > clust_old)
  cluster[idx] = cluster[idx] - 1
  
  eid_btw_mst = eid_btw_mst[-res_merge$edge_merge]
  # update edge status
  adj_vid_old = unlist(adj_list[vid_old])
  adj_eid_old = unlist(adj_edge_list[vid_old])
  clust_adj_old = cluster[adj_vid_old]
  idx_within = which(clust_adj_old == clust_new)
  eid_within = adj_eid_old[idx_within]
  edge_status[eid_within] = 'w'
  
  return(list(subgraphs = subgraphs, csize = csize, cluster = cluster,
              eid_btw_mst = eid_btw_mst, estatus = edge_status))
}


#function for hyper-step

proposeMST <- function(graph0, edge_status, subgraphs) {
  nedge = length(edge_status)
  nb = sum(edge_status == 'b')
  nw = nedge - nb
  weight = numeric(nedge)
  weight[edge_status == 'w'] = runif(nw)
  weight[edge_status == 'b'] = runif(nb, 10, 20)
  E(graph0)$weight = weight
  mstgraph = mst(graph0, weights = E(graph0)$weight)
  
  # update subgraphs
  subgraphs_new = lapply(subgraphs, function(g, mstgraph) {induced_subgraph(mstgraph, V(g)$vid)},
                         mstgraph)
  # update eid_btw_mst
  eid_btw_mst = E(mstgraph)$eid[E(mstgraph)$weight >= 10]
  
  return(list(mstgraph = mstgraph, subgraphs = subgraphs_new, eid_btw_mst = eid_btw_mst))
}



# function to standardize Y
standardize <- function(x) {
  xmean = mean(x)
  x = x - xmean
  xscale = 2 * max(abs(x))
  x = x / xscale
  param = c('mean' = xmean, 'scale' = xscale)
  return(list(x = x, std_par = param))
}

# function to unstandardize Y
unstandardize <- function(x, std_par, nomean = F, s2 = F) {
  if(s2) {
    x = x * std_par['scale'] ^ 2
  } else {
    x = x * std_par['scale']
  }
  if(!nomean) x = x + std_par['mean']
  return(x)
}




### plot a clustered variable on a spatial graph 


plotSpatGraphClust=function(graph0,coords,membership,title="", index=NULL,plotly=F,remove.edge=T){
  ## plotly: whether to make it dynamic
  ## remove.edge: whether to remove the between-cluster edges
  
  inc_mat = igraph::get.edgelist(graph0)
  p=length(unique(membership))
  
  ##Delete the between-cluster edges
  if(remove.edge){
    edge_status=getEdgeStatus(membership, inc_mat)
    edge_cutted=which(edge_status=='b')
    graph0=igraph::delete.edges(graph0, edge_cutted)
  }
  
  edgelist<-igraph::get.edgelist(graph0) 
  edgedata <- data.frame(coords[edgelist[,1],],coords[edgelist[,2],])
  colnames(edgedata) <- c("X1","Y1","X2","Y2")
  
  
  if(missing(index)){
    g1=ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edgedata, size = 0.5, colour="grey") + 
      geom_point(data=data.frame(coords)%>%mutate(membership=membership), aes(lon, lat, colour=as.factor(membership)))+scale_colour_hue()+
      ggtitle(title)+
      theme(legend.title =element_blank(),plot.title = element_text(hjust = 0.5))+ labs(fill = "Cluster ID") + xlab("X") + ylab("Y")
  }else{
    mtsub=data.frame(coords[index,]);colnames(mtsub)=c('X1','Y1');
    g1=ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edgedata, size = 0.5, colour="grey") + 
      geom_point(data=data.frame(coords)%>%mutate(membership=membership), aes(lon, lat, colour = as.factor(membership)))+scale_colour_hue()+ggtitle(title)+
      theme(legend.title =element_blank(),plot.title = element_text(hjust = 0.5))+
      geom_text(data = mtsub, aes(x = X1,y=Y1, label =
                                    rownames(mtsub)), hjust = 0, size = 4) + labs(fill = "Cluster ID") + xlab("X") + ylab("Y")
    
  }
  
  if(plotly) {g1=ggplotly(g1)}
  g1
  
}



### function to compute rand index

rand_index_fast <- function(u, v) {
  # u, v : integer, factor or character vectors of cluster labels
  stopifnot(length(u) == length(v))
  
  # 1. contingency table
  tab <- table(u, v)
  
  # 2. sums of choose(count, 2)
  nij_sum <- sum(tab * (tab - 1) / 2)        # pairs in same cluster in both
  ai      <- rowSums(tab)
  bj      <- colSums(tab)
  ai_sum  <- sum(ai * (ai - 1) / 2)          # pairs in same cluster in u only
  bj_sum  <- sum(bj * (bj - 1) / 2)          # pairs in same cluster in v only
  
  # 3. total pairs
  n           <- length(u)
  total_pairs <- n * (n - 1) / 2
  
  # 4. true negatives: pairs in different clusters in both
  tn <- total_pairs - ai_sum - bj_sum + nij_sum
  
  # 5. Rand index
  (nij_sum + tn) / total_pairs
}







#' 
#' 
#' ### Main functions -----
#' 
#' #' MCMC for sampling from posterior with balancedness constraint
#'  (Remove likelihood ratio from acceptance ratio to sample from prior)
#' #'
#' #' @param Y Standardized vector of responses
#' #' @param z vector of balance variables 
#' #' @param graph0 Spatial graph. Should be an \code{igraph} object containing \code{length(Y)} vertices
#' #' @param init_val Named list of initial values. Should include the following items.
#' #'                 \code{'tree'}: A list of M (=1) initial spanning trees, where M is the number
#' #'                                of weak learners. Each should be an \code{igraph} object.
#' #'                 \code{'cluster'}: Initial cluster membership. Should be an integer matrix of size
#' #'                                   \code{length(Y) * M}. Each column should contain consecutive 
#' #'                                   integers from 1 to the number of clusters.
#' #'                 \code{'mu'}: A list of initial values of \eqn{\mu}. Each item should be a vector of
#' #'                              length equal to the number of clusters.
#' #'                 \code{'sigmasq_y'}: Initial value of noise variance.
#' #' @param hyperpar Named list of hyperparameters. Should include the following items.
#' #'                 \code{'M'}: Number of weak learners (in our case M = 1).
#' #'                 \code{'sigmasq_mu'}: Variance of the Gaussian prior for \eqn{\mu}.
#' #'                 \code{'lambda_s'}: Scale parameter of the prior for noise variance \eqn{\sigma^2_y}.
#' #'                 \code{'nu'}: Degree of freedom of the prior for noise variance \eqn{\sigma^2_y}.
#' #'                 \code{'q'}: parameter for controlling constraint attainment
#' #'                 \code{'c'}: hyperparameter of prior controlling number of clusters
#' #' @param MCMC Number of MCMC iterations
#' #' @param BURNIN Number of burnin iterations
#' #' @param THIN Length thinning intervals. Will retain samples every \code{THIN} iterations
#' #' @param seed Random seed
#' #'
#' #' @return A list of MCMC samples containing the following items. 
#' #'         \code{'cluster_out'}: A list of \code{length(Y) * M} cluster membership matrices.
#' #'         \code{'mu_out'}: A list of posterior sample of \eqn{\mu}. \code{mu_out[[i]][[j]][k]} 
#' #'                          is the constant value for the k-th cluster in the j-th weak learner
#' #'                          in the i-th posterior sample.
#' #'         \code{'sigmasq_y_out'}: A vector of posterior draws of noise variance.
#' #'         \code{'estimated_clus'}: A vector of estimated clustering.
cohesion_forest_balanced_BA <- function(Y, z, graph0, init_val, hyperpar, MCMC, BURNIN, THIN, seed = 1234) {
  # seed =202
  set.seed(seed)
  n = vcount(graph0)
  
  # hyper-parameter
  q = hyperpar['q']
  c =  hyperpar['c']
  M = hyperpar['M']
  sigmasq_mu = hyperpar['sigmasq_mu']
  lambda_s = hyperpar['lambda_s']
  nu = hyperpar['nu']
  hyper = c(sigmasq_mu, lambda_s, nu)
  
  cluster_size = c()
  
  if('name' %in% names(vertex_attr(graph0))) {
    graph0 = igraph::delete_vertex_attr(graph0, 'name')
  }
  
  # get edgelist matrix and adjacency matrix
  inc_mat = igraph::get.edgelist(graph0, names = F)
  adj_list = lapply(igraph::as_adj_list(graph0), FUN = function(x) {x$vid})
  adj_edge_list = lapply(igraph::as_adj_edge_list(graph0), FUN = function(x) {x$eid})
  
  
  # initial values
  mstgraph_lst = init_val[['trees']]
  mstgraph_lst_dummy = init_val[['trees']]
  mu = init_val[['mu']]   # vector of cluster-wise constant values
  cluster = init_val[['cluster']]  # n*1 matrix of 1s (intial cluster membership vector)
  cluster1 = init_val[['cluster']]
  sigmasq_y = init_val[['sigmasq_y']]
  k = as.numeric(apply(cluster, 2, max))  # initial number of clusters
  
  csize = list() # cluster size vector that would update in each iteration
  csize_z = list() # balance vector of z variables
  csize_z1 = list() ## dummy balance vector of z variables
  subgraphs = list()
  eid_btw_mst = list()  #vector of between cluster edge ids
  parent_vid = NULL
  
  subgraphs_dummy = list()
  csize_dummy = list()
  eid_btw_mst_dummy = list()
  cluster_dummy = init_val[['cluster']]
  csize_z_dummy = list()
  k_dummy = as.numeric(apply(cluster, 2, max))
  edge_status_dummy = as.matrix(rep("w", times = ecount(graph0)))
  
  subgraphs_dummy1 = list()
  csize_dummy1 = list()
  cluster_dummy1 = init_val[['cluster']]
  csize_z_dummy1 = list()
  k_dummy = as.numeric(apply(cluster, 2, max))
  edge_status_dummy = init_val[['edge_status']]
  g = matrix(0, nrow = n, ncol = 1)  # n*1 matrix of fitted mu's
  
  cluster_iter = cluster[,1]
  
  csize[[1]] = Rfast::Table(cluster_iter)
  csize_z[[1]] = as.vector(tapply(z, cluster[, 1], sum))
  mstgraph_iter = mstgraph_lst[[1]]
  clust_vid_iter = split(1:n, cluster_iter)
  subgraphs[[1]] = lapply(clust_vid_iter, function(vids, mstgraph) {
    igraph::induced_subgraph(mstgraph, vids)
  }, mstgraph_iter)
  inc_mat_mst = igraph::get.edgelist(mstgraph_iter, names = F)
  c1_m = cluster_iter[inc_mat_mst[, 1]]; c2_m = cluster_iter[inc_mat_mst[, 2]]
  idx_btw = which(c1_m != c2_m)
  idx_within = which(c1_m == c2_m)
  eid_btw_mst[[1]] = (igraph::E(mstgraph_iter)$eid)[idx_btw]
  
  
  # whether an edge in graph0 is within a cluster or bewteen two clusters
  # n*M matrix
  # edge_status = apply(cluster, 2, FUN = getEdgeStatus, inc_mat)
  # edge_status = as.matrix(rep("w", times = ecount(graph0)))
  edge_status[,1] = init_val[['edge_status']]
  
  ################# MCMC ####################
  M = 1
  ## MCMC results
  mu_out = list()
  sigmasq_y_out = numeric((MCMC-BURNIN)/THIN)
  cluster_out = array(0, dim = c((MCMC-BURNIN)/THIN, n, M))
  tree_out = list()
  log_post_out = numeric((MCMC-BURNIN)/THIN)
  shannon9 = c()
  shannon_index_vec = c()
  num = 1
  
  depth = 1
  acc = 0
  iter = 1
  m=1
  mat = matrix(0, nrow = n, ncol = (MCMC -BURNIN)/THIN)
  
  relative_entropy = c()
  # log_prior_prob_old_OCO = log(1 - psplit(1, a, b))
  shan_old = diversity(as.vector(csize_z[[1]]), index = "shannon") / log(length(as.vector(csize_z[[1]])))
  log_prior_prob_old_shan = (q * shan_old)
  acc_shan_vec = c()
  
  trees = list()
  trees[[1]] = subgraphs[[1]]
  iter = 1
  ## MCMC iteration
  subgraphs_iter = create_subgraphs_from_clusters(graph0, cluster_iter)
  subgraphs[[1]] = subgraphs_iter
  # print(length(cluster_iter))
  log_lik_old = cluster_log_density(Y, sigmasq_mu, sigmasq_y, cluster_iter)
  # print("done")
  
  split_accepted_count_full = 0
  merge_accepted_count_full = 0
  change_accepted_count_full = 0
  log_post_out_full = c()
  
  
  split_vec = rep(0, times = MCMC)
  merge_vec = rep(0, times = MCMC)
  change_vec = rep(0, times = MCMC)
  
  for(iter in 1:MCMC) {
    if(iter == BURNIN +1){
      split_accepted_count = 0
      merge_accepted_count = 0
      change_accepted_count = 0
    }
    # 
    k_iter = k[1]
    cluster_iter = cluster[, 1]
    csize_iter = csize[[1]]
    csize_z_iter = csize_z[[1]]
    subgraphs_iter = subgraphs[[1]]
    trees_iter = trees[[1]]
    
    e_iter = Y
    
    if(k_iter == 1) {rb = 1; rd = 0; rc = 0
    } else if(k_iter == n) {rb = 0; rd = 0.6; rc = 0.4
    } else {rb = 0.3; rd = 0.3; rc = 0.3; rh = 0.1}
    move = sample(4, 1, prob = c(rb, rd, rc, rh))
    # print(move)
    if(move == 1) { ## birth move
      ### split a cluster with by removing any edge with equal probability
      update_res2 =  forest_birth(graph0, subgraphs_iter, trees_iter, cluster_iter)
      cluster1 = as.vector(update_res2$cluster)
      
      csize_z1 = as.vector(tapply(z, cluster1, sum))
      new_shan = diversity(csize_z1, index = "shannon") / log(length(csize_z1))
      log_prior_prob_new_shan = (q * new_shan)
      birth_prior_ratio_term = update_res2$birth_prior_ratio_term
      
      # compute log-prior ratio
      log_A = (log_prior_prob_new_shan - log_prior_prob_old_shan) +log(1-c) 
      
      subgraphs1 = update_res2$subgraphs
      # chosen_subgraph = update_res2$chosen_subgraph
      
      ## log proposal ratio
      if(k_iter == n-1) {
        rd_new = 0.6
      } else {rd_new = 0.3}
      # log_P = log(rd_new) - log(rb) + log(chosen_cluster_size-1)+
      #   log_numSpanningTrees(chosen_subgraph)- numerator -
      #   log(count_connected_pairs(graph0, subgraphs1)) + log(nedge - 1) - log(vcount(chosen_subgraph)-1)
      # 
      # log_P = log(rd_new) - log(rb) +
      #   log_numSpanningTrees(chosen_subgraph)- numerator -
      #   log(count_connected_pairs(graph0, subgraphs1)) + log(vcount(graph0) - k_iter)
      # n_e = 
      
      log_P = log(rd_new) - log(rb) -
        log(ecount(contract_to_quotient_graph(graph0, cluster1))) - 
        log(countCrossEdges(graph0, update_res2$vid1, update_res2$vid2)) + 
        log(vcount(graph0) - k_iter) 
      
      log_lik_new = cluster_log_density(Y, sigmasq_mu, sigmasq_y, cluster1)
      log_L = log_lik_new - log_lik_old
      acc_prob = min(0, log_A + log_P + log_L)
      acc_prob = exp(acc_prob)
      if(runif(1) < acc_prob){
        # accept
        print("split is accepted")
        split_vec[iter] =1
        update_res = update_res2
        subgraphs[[1]] = update_res$subgraphs
        trees[[1]] = update_res$trees
        csize[[1]] = update_res$csize
        k[1] = k[1] + 1
        cluster[, 1] = update_res$cluster
        csize_z[[1]] = as.vector(tapply(z, cluster[, 1], sum))
        ## Compute and store the prior probability of the tree formed
        log_prior_prob_old_shan = log_prior_prob_new_shan
        old_shan = new_shan
        split_accepted_count_full = split_accepted_count_full + 1
        if(iter > BURNIN){
          split_accepted_count = split_accepted_count + 1
        }
        log_lik_old = log_lik_new
        # print(plotSpatGraphClust(graph0, coords, cluster[, 1]))
        
      }
    }
    if(move == 2) { ## death move
      update_res2 = forest_death(graph0, subgraphs_iter, trees_iter, cluster_iter)
      cluster1 = as.vector(update_res2$cluster)
      csize_z1 = as.vector(tapply(z, cluster1, sum))

      ### Log-prior ratio
      if(length(csize_z1) == 1){
        new_shan = 0
      }else{
        new_shan = diversity(csize_z1, index = "shannon") / log(length(csize_z1))
      }
      death_prior_ratio_term = update_res2$death_prior_ratio_term
      log_prior_prob_new_shan = (q * new_shan)
      log_A = (log_prior_prob_new_shan - log_prior_prob_old_shan) - log(1-c) 
      
      ## log proposal ratio
      if(k_iter == 2) {rb_new = 1
      }else {rb_new = 0.3}
      
      chosen_subgraph = update_res2$chosen_subgraph
      numerator = update_res2$count_inducing_spanning_tree
      
      log_P = log(rb_new) - log(rd) + 
        log(ecount(contract_to_quotient_graph(graph0, cluster_iter))) + 
        log(countCrossEdges(graph0, update_res2$vid1, update_res2$vid2)) - 
        log(vcount(graph0) - k_iter + 1) 
      
      
      # compute log-likelihood ratio
      log_lik_new = cluster_log_density(Y, sigmasq_mu, sigmasq_y, cluster1)
      log_L = log_lik_new - log_lik_old
      
      #   # acceptance probability
      acc_prob = min(0, log_A + log_P + log_L)
      acc_prob = exp(acc_prob)
      if(runif(1) < acc_prob){
        # accept
        print("merge is accepted")
        update_res = update_res2
        k[1] = k[1] - 1
        subgraphs[[1]] = update_res$subgraphs
        csize[[1]] = update_res$csize
        cluster[, 1] = update_res$cluster
        csize_z[[1]] = as.vector(tapply(z, cluster[, 1], sum))
        trees[[1]] = update_res$trees
        merge_accepted_count_full = merge_accepted_count_full + 1
        
        ## Compute and store the prior probability of the tree formed
        log_prior_prob_old_shan = log_prior_prob_new_shan
        old_shan = new_shan
        if(iter > BURNIN){
          merge_accepted_count = merge_accepted_count + 1
        }
        log_lik_old = log_lik_new
        merge_vec[iter] =1
        
        # print(plotSpatGraphClust(graph0, coords, cluster[, 1]))
        
      }
    }
    if(move == 3) { ## death move
      update_res1 = forest_death(graph0, subgraphs_iter, trees_iter, cluster_iter)
      cluster1 = as.vector(update_res1$cluster)
      csize_z1 = as.vector(tapply(z, cluster1, sum))
      # print(length(cluster1))
      
      ### Log-prior ratio
      if(length(csize_z1) == 1){
        new_shan = 0
      }else{
        new_shan = diversity(csize_z1, index = "shannon") / log(length(csize_z1))
      }
      death_prior_ratio_term = update_res1$death_prior_ratio_term
      log_prior_prob_new_shan = (q * new_shan)
      log_A1 = (log_prior_prob_new_shan - log_prior_prob_old_shan) - log(1-c) 
      
      ## log proposal ratio
      if(k_iter == 2) {rb_new = 1
      }else {rb_new = 0.3}
      
      
      log_P1 =log(rb_new) - log(rd) + 
        log(ecount(contract_to_quotient_graph(graph0, cluster_iter))) + 
        log(countCrossEdges(graph0, update_res1$vid1, update_res1$vid2)) - 
        log(vcount(graph0) - k_iter + 1) 
      
      
      
      subgraphs1 = update_res1$subgraphs
      trees1 = update_res1$trees
      update_res2 =  forest_birth(graph0, subgraphs1, trees1, cluster1)
      cluster1 = as.vector(update_res2$cluster)
      csize_z1 = as.vector(tapply(z, cluster1, sum))
      new_shan = diversity(csize_z1, index = "shannon") / log(length(csize_z1))
      log_prior_prob_new_shan = (q * new_shan)
      birth_prior_ratio_term = update_res2$birth_prior_ratio_term
      
      # compute log-prior ratio
      log_A = (log_prior_prob_new_shan - log_prior_prob_old_shan) 
      
      subgraphs1 = update_res2$subgraphs

      ## log proposal ratio
      if(k_iter == n-1) {
        rd_new = 0.6
      } else {rd_new = 0.3}
      
      
      log_P2 = log(rd_new) - log(rb) -
        log(ecount(contract_to_quotient_graph(graph0, cluster1))) - 
        log(countCrossEdges(graph0, update_res2$vid1, update_res2$vid2)) + 
        log(vcount(graph0) - k_iter +1) 
      
      log_lik_new = cluster_log_density(Y, sigmasq_mu, sigmasq_y, cluster1)
      
      # log_A = log_A1 + log_A2
      log_P = log_P1 + log_P2
      log_L = log_lik_new - log_lik_old
      
      # log_L = evalLogLikeRatio('merge', e_iter, update_res2$vid_old, update_res2$vid_new, sigmasq_y, sigmasq_mu)
      # print(log_L)
      #   # acceptance probability
      acc_prob = min(0, log_A + log_P + log_L)
      acc_prob = exp(acc_prob)
      if(runif(1) < acc_prob){
        # accept
        print("change is accepted")
        update_res = update_res2
        subgraphs[[1]] = update_res$subgraphs
        csize[[1]] = update_res$csize
        cluster[, 1] = update_res$cluster
        csize_z[[1]] = as.vector(tapply(z, cluster[, 1], sum))
        change_accepted_count_full = change_accepted_count_full + 1
        
        ## Compute and store the prior probability of the tree formed
        log_prior_prob_old_shan = log_prior_prob_new_shan
        old_shan = new_shan
        if(iter > BURNIN){
          change_accepted_count = change_accepted_count + 1
        }
        log_lik_old = log_lik_new
        trees[[1]] = update_res$trees
        change_vec[iter] =1
        
        
        # print(plotSpatGraphClust(graph0, coords, cluster[, 1]))
        
      }
    }
    if(move ==4){
      ## hyper
      print("hyper is accepted")
      trees[[1]] = drawUSTs(subgraphs_iter) 
    }
    
    # plotSpatGraphClust(mstgraph_m, coords, cluster[, m])
    
    # update mu_m
    k_iter = k[1]; cluster_iter = cluster[, 1]; csize_iter = csize[[1]]; 
    csize_z_iter = csize_z[[1]]
    Qinv_diag = 1 / (csize_iter/sigmasq_y + 1/sigmasq_mu)
    group_sum <- tapply(e_iter, cluster_iter, sum)
    b_t <- Qinv_diag * group_sum / sigmasq_y
    
    mu[[1]] = rnorm(k_iter, b_t, sqrt(Qinv_diag))
    g[, 1] = mu[[1]][cluster_iter]
    cluster_size[iter + 1] = k
    
    # update sigmasq_y
    Y_hat = g[, 1]
    rate = 0.5*(nu*lambda_s + sum((Y - Y_hat)^2))
    sigmasq_y = 1/rgamma(1, shape = (n+nu)/2, rate = rate)
    
    shannon_index_vec[num] = diversity(csize_z_iter, "shannon") / log(length(csize_z_iter))
    num = num + 1
    
    # if(iter %% THIN == 0){
    #   log_post_out_full[iter/THIN] = evalLogPost(graph0, mu, g, sigmasq_y, 
    #                                              k_iter, Y, hyper, q, l = shannon_index_vec[num -1])
    #   
    # }
    
    ## save result
    if(iter > BURNIN & (iter - BURNIN) %% THIN == 0) {
      mu_out[[(iter-BURNIN)/THIN]] = mu
      sigmasq_y_out[(iter-BURNIN)/THIN] = sigmasq_y
      cluster_out[(iter-BURNIN)/THIN, , ] = cluster
      acc_shan_vec[(iter-BURNIN)/THIN] = old_shan
      
    }
    
    if(iter %% 100 == 0){
      cat('Iteration', iter, 'done\n')
    }
    
  }
  
  mode(cluster_out) = 'integer'  # to save memory
  
  
  
  clus = as.matrix((cluster_out)[,,1])
  vec = as.vector(dlso(clus, loss=VI()))
  estimated_clus = vec
  acc_probability = acc /((MCMC - BURNIN)/THIN)
  acc_count = c(split_accepted_count, 
                merge_accepted_count, change_accepted_count)
  acc_vec = acc_count / (0.3 * (MCMC - BURNIN))
  
  return(list('mu_out' = mu_out,
              'sigmasq_y_out' = sigmasq_y_out,
              'cluster_out' = cluster_out, 
              'estimated_clus' = estimated_clus))
  
}




#' 
#' 
#' ### Main functions -----
#' 
#' #' MCMC for sampling from posterior with shrinkage constraint
#'  (Remove likelihood ratio from acceptance ratio to sample from prior)
#' #'
#' #' @param Y Standardized vector of responses
#' #' @param graph0 Spatial graph. Should be an \code{igraph} object containing \code{length(Y)} vertices
#' #' @param init_val Named list of initial values. Should include the following items.
#' #'                 \code{'tree'}: A list of M (=1) initial spanning trees, where M is the number
#' #'                                of weak learners. Each should be an \code{igraph} object.
#' #'                 \code{'cluster'}: Initial cluster membership. Should be an integer matrix of size
#' #'                                   \code{length(Y) * M}. Each column should contain consecutive 
#' #'                                   integers from 1 to the number of clusters.
#' #'                 \code{'mu'}: A list of initial values of \eqn{\mu}. Each item should be a vector of
#' #'                              length equal to the number of clusters.
#' #'                 \code{'sigmasq_y'}: Initial value of noise variance.
#' #' @param hyperpar Named list of hyperparameters. Should include the following items.
#' #'                 \code{'M'}: Number of weak learners (in our case M = 1).
#' #'                 \code{'sigmasq_mu'}: Variance of the Gaussian prior for \eqn{\mu}.
#' #'                 \code{'lambda_s'}: Scale parameter of the prior for noise variance \eqn{\sigma^2_y}.
#' #'                 \code{'nu'}: Degree of freedom of the prior for noise variance \eqn{\sigma^2_y}.
#' #'                 \code{'q'}: parameter for controlling constraint attainment
#' #'                 \code{'c'}: hyperparameter of prior controlling number of clusters
#' #' @param MCMC Number of MCMC iterations
#' #' @param BURNIN Number of burnin iterations
#' #' @param THIN Length thinning intervals. Will retain samples every \code{THIN} iterations
#' #' @param prior_cluster Pre-specified clustering
#' #' @param seed Random seed
#' #' @return A list of MCMC samples containing the following items. 
#' #'         \code{'cluster_out'}: A list of \code{length(Y) * M} cluster membership matrices.
#' #'         \code{'mu_out'}: A list of posterior sample of \eqn{\mu}. \code{mu_out[[i]][[j]][k]} 
#' #'                          is the constant value for the k-th cluster in the j-th weak learner
#' #'                          in the i-th posterior sample.
#' #'         \code{'sigmasq_y_out'}: A vector of posterior draws of noise variance.
#' #'         \code{'estimated_clus'}: A vector of estimated clustering.

cohesion_forest_shrinkage_BA <- function(Y, graph0, init_val, hyperpar, MCMC, BURNIN, THIN, prior_cluster, seed = 1234) {
  set.seed(seed)
  # set.seed(1234)
  n = vcount(graph0)
  # hyper-parameter
  q = hyperpar['q']
  c =  hyperpar['c']
  M = hyperpar['M']
  sigmasq_mu = hyperpar['sigmasq_mu']
  lambda_s = hyperpar['lambda_s']
  nu = hyperpar['nu']
  hyper = c(sigmasq_mu, lambda_s, nu)
  
  cluster_size = c()
  
  if('name' %in% names(vertex_attr(graph0))) {
    graph0 = igraph::delete_vertex_attr(graph0, 'name')
  }
  
  # get edgelist matrix and adjacency matrix
  inc_mat = igraph::get.edgelist(graph0, names = F)
  adj_list = lapply(igraph::as_adj_list(graph0), FUN = function(x) {x$vid})
  adj_edge_list = lapply(igraph::as_adj_edge_list(graph0), FUN = function(x) {x$eid})
  
  
  # initial values
  mstgraph_lst = init_val[['trees']]
  mstgraph_lst_dummy = init_val[['trees']]
  mu = init_val[['mu']]   # vector of cluster-wise constant values
  cluster = init_val[['cluster']]  # n*1 matrix of 1s (intial cluster membership vector)
  cluster1 = init_val[['cluster']]
  sigmasq_y = init_val[['sigmasq_y']]
  k = as.numeric(apply(cluster, 2, max))  # initial number of clusters
  
  csize = list() # cluster size vector that would update in each iteration
  csize_z = list() # balance vector of z variables
  csize_z1 = list() ## dummy balance vector of z variables
  subgraphs = list()
  eid_btw_mst = list()  #vector of between cluster edge ids
  parent_vid = NULL
  
  subgraphs_dummy = list()
  csize_dummy = list()
  eid_btw_mst_dummy = list()
  cluster_dummy = init_val[['cluster']]
  csize_z_dummy = list()
  k_dummy = as.numeric(apply(cluster, 2, max))
  edge_status_dummy = as.matrix(rep("w", times = ecount(graph0)))
  
  subgraphs_dummy1 = list()
  csize_dummy1 = list()
  cluster_dummy1 = init_val[['cluster']]
  csize_z_dummy1 = list()
  k_dummy = as.numeric(apply(cluster, 2, max))
  edge_status_dummy = init_val[['edge_status']]
  g = matrix(0, nrow = n, ncol = 1)  # n*1 matrix of fitted mu's
  
  cluster_iter = cluster[,1]
  # g[, 1] = mu[[1]][cluster_iter]
  
  
  csize[[1]] = Rfast::Table(cluster_iter)
  csize_z[[1]] = as.vector(tapply(z, cluster[, 1], sum))
  mstgraph_iter = mstgraph_lst[[1]]
  clust_vid_iter = split(1:n, cluster_iter)
  subgraphs[[1]] = lapply(clust_vid_iter, function(vids, mstgraph) {
    igraph::induced_subgraph(mstgraph, vids)
  }, mstgraph_iter)
  inc_mat_mst = igraph::get.edgelist(mstgraph_iter, names = F)
  c1_m = cluster_iter[inc_mat_mst[, 1]]; c2_m = cluster_iter[inc_mat_mst[, 2]]
  idx_btw = which(c1_m != c2_m)
  idx_within = which(c1_m == c2_m)
  eid_btw_mst[[1]] = (igraph::E(mstgraph_iter)$eid)[idx_btw]
  
  
  # whether an edge in graph0 is within a cluster or bewteen two clusters
  # n*M matrix
  # edge_status = apply(cluster, 2, FUN = getEdgeStatus, inc_mat)
  # edge_status = as.matrix(rep("w", times = ecount(graph0)))
  edge_status[,1] = init_val[['edge_status']]
  
  ################# MCMC ####################
  M = 1
  ## MCMC results
  mu_out = list()
  sigmasq_y_out = numeric((MCMC-BURNIN)/THIN)
  cluster_out = array(0, dim = c((MCMC-BURNIN)/THIN, n, M))
  tree_out = list()
  log_post_out = numeric((MCMC-BURNIN)/THIN)
  shannon9 = c()
  shannon_index_vec = c()
  num = 1
  
  depth = 1
  acc = 0
  iter = 1
  m=1
  mat = matrix(0, nrow = n, ncol = (MCMC -BURNIN)/THIN)
  
  relative_entropy = c()
  rand_old = rand_index_fast(cluster_iter, prior_cluster)
  log_prior_prob_old_rand = (q * rand_old)
  acc_shan_vec = c()
  
  trees = list()
  trees[[1]] = subgraphs[[1]]
  iter = 1
  ## MCMC iteration
  subgraphs_iter = create_subgraphs_from_clusters(graph0, cluster_iter)
  subgraphs[[1]] = subgraphs_iter
  log_lik_old = cluster_log_density(Y, sigmasq_mu, sigmasq_y, cluster_iter)

  split_accepted_count_full = 0
  merge_accepted_count_full = 0
  change_accepted_count_full = 0
  log_post_out_full = c()
  
  
  for(iter in 1:MCMC) {
    if(iter == BURNIN +1){
      
      split_accepted_count = 0
      merge_accepted_count = 0
      change_accepted_count = 0
    }
    # 
    k_iter = k[1]
    cluster_iter = cluster[, 1]
    csize_iter = csize[[1]]
    subgraphs_iter = subgraphs[[1]]
    trees_iter = trees[[1]]

    e_iter = Y
    
    if(k_iter == 1) {rb = 1; rd = 0; rc = 0
    } else if(k_iter == n) {rb = 0; rd = 0.6; rc = 0.4
    } else {rb = 0.3; rd = 0.3; rc = 0.3; rh = 0.1}
    move = sample(4, 1, prob = c(rb, rd, rc, rh))
    if(move == 1) { ## birth move
      ### split a cluster with by removing any edge with equal probability
      update_res2 =  forest_birth(graph0, subgraphs_iter, trees_iter, cluster_iter)
      cluster1 = as.vector(update_res2$cluster)
      
      rand_new = rand_index_fast(cluster1, prior_cluster)
      log_prior_prob_new_rand = (q * rand_new)
      # compute log-prior ratio
      log_A = (log_prior_prob_new_rand - log_prior_prob_old_rand) +log(1-c) 
      
      subgraphs1 = update_res2$subgraphs

      ## log proposal ratio
      if(k_iter == n-1) {
        rd_new = 0.6
      } else {rd_new = 0.3}
      
      log_P = log(rd_new) - log(rb) -
        log(ecount(contract_to_quotient_graph(graph0, cluster1))) - 
        log(countCrossEdges(graph0, update_res2$vid1, update_res2$vid2)) + 
        log(vcount(graph0) - k_iter) 
      
      # compute log-likelihood ratio
      log_lik_new = cluster_log_density(Y, sigmasq_mu, sigmasq_y, cluster1)
      log_L = log_lik_new - log_lik_old
      acc_prob = min(0, log_A + log_P + log_L)
      acc_prob = exp(acc_prob)
      if(runif(1) < acc_prob){
        # accept
        print("split is accepted")
        update_res = update_res2
        subgraphs[[1]] = update_res$subgraphs
        trees[[1]] = update_res$trees
        csize[[1]] = update_res$csize
        k[1] = k[1] + 1
        cluster[, 1] = update_res$cluster
        csize_z[[1]] = as.vector(tapply(z, cluster[, 1], sum))
        ## Compute and store the prior probability of the tree formed
        log_prior_prob_old_rand = log_prior_prob_new_rand
        rand_old = rand_new
        split_accepted_count_full = split_accepted_count_full + 1
        
        if(iter > BURNIN){
          split_accepted_count = split_accepted_count + 1
        }
        log_lik_old = log_lik_new
        # print(plotSpatGraphClust(graph0, coords, cluster[, 1]))
        
      }
    }
    if(move == 2) { ## death move
      update_res2 = forest_death(graph0, subgraphs_iter, trees_iter, cluster_iter)
      cluster1 = as.vector(update_res2$cluster)
      csize_z1 = as.vector(tapply(z, cluster1, sum))

      ### Log-prior ratio
      rand_new = rand_index_fast(cluster1, prior_cluster)
      log_prior_prob_new_rand = (q * rand_new)
      
      # compute log-prior ratio
      log_A = (log_prior_prob_new_rand - log_prior_prob_old_rand) -log(1-c) 
      
      ## log proposal ratio
      if(k_iter == 2) {rb_new = 1
      }else {rb_new = 0.3}
      
      
      log_P = log(rb_new) - log(rd) + 
        log(ecount(contract_to_quotient_graph(graph0, cluster_iter))) + 
        log(countCrossEdges(graph0, update_res2$vid1, update_res2$vid2)) - 
        log(vcount(graph0) - k_iter + 1) 
      
      # compute log-likelihood ratio
      log_lik_new = cluster_log_density(Y, sigmasq_mu, sigmasq_y, cluster1)
      log_L = log_lik_new - log_lik_old
      
      #   # acceptance probability
      acc_prob = min(0, log_A + log_P + log_L)
      acc_prob = exp(acc_prob)
      if(runif(1) < acc_prob){
        # accept
        print("merge is accepted")
        update_res = update_res2
        k[1] = k[1] - 1
        subgraphs[[1]] = update_res$subgraphs
        csize[[1]] = update_res$csize
        cluster[, 1] = update_res$cluster
        csize_z[[1]] = as.vector(tapply(z, cluster[, 1], sum))
        trees[[1]] = update_res$trees
        merge_accepted_count_full = merge_accepted_count_full + 1
        
        ## Compute and store the prior probability of the tree formed
        log_prior_prob_old_rand = log_prior_prob_new_rand
        rand_old = rand_new
        if(iter > BURNIN){
          merge_accepted_count = merge_accepted_count + 1
        }
        log_lik_old = log_lik_new
        
        # print(plotSpatGraphClust(graph0, coords, cluster[, 1]))
        
      }
    }
    if(move == 3) { ## death move
      update_res1 = forest_death(graph0, subgraphs_iter, trees_iter, cluster_iter)
      cluster1 = as.vector(update_res1$cluster)
      csize_z1 = as.vector(tapply(z, cluster1, sum))

      ### Log-prior ratio
      if(length(csize_z1) == 1){
        new_shan = 0
      }else{
        new_shan = diversity(csize_z1, index = "shannon") / log(length(csize_z1))
      }
      
      ## log proposal ratio
      if(k_iter == 2) {rb_new = 1
      }else {rb_new = 0.3}
      
      
      log_P1 =log(rb_new) - log(rd) + 
        log(ecount(contract_to_quotient_graph(graph0, cluster_iter))) + 
        log(countCrossEdges(graph0, update_res1$vid1, update_res1$vid2)) - 
        log(vcount(graph0) - k_iter + 1) 
      
      
      
      subgraphs1 = update_res1$subgraphs
      trees1 = update_res1$trees

      update_res2 =  forest_birth(graph0, subgraphs1, trees1, cluster1)
      cluster1 = as.vector(update_res2$cluster)
      rand_new = rand_index_fast(cluster1, prior_cluster)
      log_prior_prob_new_rand = (q * rand_new)
      
      # compute log-prior ratio
      log_A = (log_prior_prob_new_rand - log_prior_prob_old_rand) 
      
      subgraphs1 = update_res2$subgraphs
      # chosen_subgraph = update_res2$chosen_subgraph
      
      ## log proposal ratio
      if(k_iter == n-1) {
        rd_new = 0.6
      } else {rd_new = 0.3}
      
      
      log_P2 = log(rd_new) - log(rb) -
        log(ecount(contract_to_quotient_graph(graph0, cluster1))) - 
        log(countCrossEdges(graph0, update_res2$vid1, update_res2$vid2)) + 
        log(vcount(graph0) - k_iter +1) 
      
      log_lik_new = cluster_log_density(Y, sigmasq_mu, sigmasq_y, cluster1)
      
      log_P = log_P1 + log_P2
      log_L = log_lik_new - log_lik_old
      
      #   # acceptance probability
      acc_prob = min(0, log_A + log_P + log_L)
      acc_prob = exp(acc_prob)
      if(runif(1) < acc_prob){
        # accept
        print("change is accepted")
        update_res = update_res2
        subgraphs[[1]] = update_res$subgraphs
        csize[[1]] = update_res$csize
        cluster[, 1] = update_res$cluster
        csize_z[[1]] = as.vector(tapply(z, cluster[, 1], sum))
        change_accepted_count_full = change_accepted_count_full + 1
        
        ## Compute and store the prior probability of the tree formed
        log_prior_prob_old_rand = log_prior_prob_new_rand
        rand_old = rand_new
        if(iter > BURNIN){
          change_accepted_count = change_accepted_count + 1
        }
        trees[[1]] = update_res$trees
        log_lik_old = log_lik_new
        
        
        # print(plotSpatGraphClust(graph0, coords, cluster[, 1]))
        
      }
    }
    if(move ==4){
      ## hyper
      print("hyper is accepted")
      trees[[1]] = drawUSTs(subgraphs_iter) 
    }
    
    # plotSpatGraphClust(mstgraph_m, coords, cluster[, m])
    
    # update mu_m
    k_iter = k[1]; cluster_iter = cluster[, 1]; csize_iter = csize[[1]]; 
    Qinv_diag = 1 / (csize_iter/sigmasq_y + 1/sigmasq_mu)
    group_sum <- tapply(e_iter, cluster_iter, sum)
    b_t <- Qinv_diag * group_sum / sigmasq_y
    
    mu[[1]] = rnorm(k_iter, b_t, sqrt(Qinv_diag))
    g[, 1] = mu[[1]][cluster_iter]
    cluster_size[iter + 1] = k
  
    # update sigmasq_y
    Y_hat = g[, 1]
    rate = 0.5*(nu*lambda_s + sum((Y - Y_hat)^2))
    sigmasq_y = 1/rgamma(1, shape = (n+nu)/2, rate = rate)
    
     num = num + 1
    
    ## save result
    # if(iter %% THIN == 0){
    #   log_post_out_full[iter/THIN] = evalLogPost(graph0, mu, g, sigmasq_y, 
    #                                              k_iter, Y, hyper, q, l = rand_old)
    #   
    # }
    
    if(iter > BURNIN & (iter - BURNIN) %% THIN == 0) {
      mu_out[[(iter-BURNIN)/THIN]] = mu
      sigmasq_y_out[(iter-BURNIN)/THIN] = sigmasq_y
      cluster_out[(iter-BURNIN)/THIN, , ] = cluster
      acc_shan_vec[(iter-BURNIN)/THIN] = rand_old
      
    }
    
    if(iter %% 100 == 0){
      cat('Iteration', iter, 'done\n')
    }
    
  }
  
  mode(cluster_out) = 'integer'  # to save memory
  
  
  
  clus = as.matrix((cluster_out)[,,1])
  vec = as.vector(dlso(clus, loss=VI()))
  
  estimated_clus = vec
  acc_probability = acc /((MCMC - BURNIN)/THIN)
  acc_count = c(split_accepted_count, 
                merge_accepted_count, change_accepted_count)
  acc_vec = acc_count / (0.3 * (MCMC - BURNIN))
  return(list('mu_out' = mu_out,
              'sigmasq_y_out' = sigmasq_y_out,
              'cluster_out' = cluster_out, 
              'estimated_clus' = estimated_clus
  ))
  
}





