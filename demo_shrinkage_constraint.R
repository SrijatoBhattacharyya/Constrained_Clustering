rm(list=ls())
library(igraph)
library(fields)
library(ggplot2)
library(fossil)
library(caret)
library(MLmetrics)
library(mclust)
library(proxy)
library(dplyr)
# set working directory if necessary
# source('ComplexDomainFun.R')
# source("comp-dom-functions.R")
source("ComplexDomain_utils.R")
source('constraint_FUN.R')



## Ushaped
set.seed(200)
Ushape=gensfUbnd1(rate=1,rot_angle= 45,cluster.value=c(-1,0, 1))
sf_bnd=Ushape$bnd; sf_cluster=Ushape$cluster; ggplot(Ushape$cluster)+geom_sf(aes(fill=as.factor(beta)))
bnd = data.frame(sf::st_coordinates(sf_bnd)[, 1:2]);colnames(bnd)=c('x','y')
# plot(bnd)
#' #' #' ## Generate n uniform training locations
#' #' #'
n = 600
sf_coords=st_sample(sf_bnd,n);coords = st_coordinates(sf_coords);colnames(coords)=c('lon','lat')

# plot(coords) 


cluster_true=unlist(st_within(sf_coords,sf_cluster))

#' #' #' ## Generate piecewise functions
f_true=sf_cluster$beta[cluster_true]

clus_original = f_true + 2
# l1 = length(clus_original[clus_original == 2])
# l2 = length(clus_original[clus_original == 1])
# l3 = length(clus_original[clus_original == 3])
# diversity(c(l1, l2, l3))/log(3)
# c(l1, l2, l3)


Y = f_true + rnorm(n, 0, 0.2)

# rot_angle = 45
# 
# ggplot() +
#   # geom_sf(data=clusterD) +
#   geom_point(aes(x = lon, y = lat, col=Y), data = as.data.frame(coords)) +
#   scale_color_gradientn(colours = rainbow(8), name = 'Y')

# plot observed data
ggplot() + 
  geom_boundary(bnd) +
  geom_point(aes(x = lon, y = lat, col = Y), data = as.data.frame(coords)) +
  scale_color_gradientn(colours = rainbow(5), name = 'Y') +
  labs(x = 'Scaled Lon.', y = 'Scaled Lat.') + 
  ggtitle('Observed Data')

# print(plotSpatGraphClust(graph0, coords, clus_original))
plot(coords, col=clus_original, pch=16, main="true")


sf_bnd=Ushape$bnd1; sf_cluster=Ushape$cluster1; ggplot(Ushape$cluster1)+geom_sf(aes(fill=as.factor(beta)))
bnd = data.frame(sf::st_coordinates(sf_bnd)[, 1:2]);colnames(bnd)=c('x','y')
plot(bnd)
#' #' #' ## Generate n uniform training locations
#' #' #'
n = 600
# sf_coords=st_sample(sf_bnd,n);coords = st_coordinates(sf_coords);colnames(coords)=c('lon','lat')

# plot(coords) 
#' #' # generate n_ho uniform hold-out locations
# n_ho = 200
# sf_coords_ho=st_sample(sf_bnd,n_ho);coords_ho = st_coordinates(sf_coords_ho);colnames(coords_ho)=c('lon','lat')
# plot(coords_ho)

cluster_true=unlist(st_within(sf_coords,sf_cluster))
# cluster_ho_true=unlist(st_within(sf_coords_ho,sf_cluster))

#' #' #' ## Generate piecewise functions
f_true=sf_cluster$beta[cluster_true]
# f_ho_true=sf_cluster$beta[cluster_ho_true]
# Y_ho = f_ho_true + rnorm(n_ho, 0, 0.5)

prior_cluster = f_true + 2
# l1 = length(clus_original[clus_original == 2])
# l2 = length(clus_original[clus_original == 1])
# l3 = length(clus_original[clus_original == 3])
# diversity(c(l1, l2, l3))/log(3)
# c(l1, l2, l3)

Y2 = f_true + rnorm(n, 0, 0.2)

# rot_angle
# plot observed data
ggplot() + 
  geom_boundary(bnd) +
  geom_point(aes(x = lon, y = lat, col = Y2), data = as.data.frame(coords)) +
  scale_color_gradientn(colours = rainbow(5), name = 'Y2') +
  labs(x = 'Scaled Lon.', y = 'Scaled Lat.') + 
  ggtitle('Observed Data')


# print(plotSpatGraphClust(graph0, coords, prior_cluster))
plot(coords, col=prior_cluster, pch=16, main="prespecified cluster")
plot(coords, col=clus_original, pch=16, main="prespecified cluster")

rand.index(prior_cluster, clus_original)




# cluster_true - prior_cluster
### initializing parameters ----

# get mesh and triangulation
mesh = gen2dMesh(coords, bnd)
graph0 = constrainedDentri(n, mesh)
E(graph0)$eid = c(1:ecount(graph0))  # edge id
V(graph0)$vid = c(1:vcount(graph0))  # vertex id
mstgraph = mst(graph0)  # initial spanning tree
graph0 = delete_edge_attr(graph0, 'weight')
mstgraph0 = delete_edge_attr(mstgraph, 'weight')

# plot spatial graph
plotGraph(coords, graph0) + 
  geom_boundary(bnd) +
  labs(x = 'Scaled Lon.', y = 'Scaled Lat.') + 
  ggtitle('Spatial Graph')

# z = rep(1, times = vcount(graph0))

set.seed(200)
Ushape1=gensfUbnd(rate=1,rot_angle=45,cluster.value=c(100,100, 3))
# Ushape=gensfUbnd2(rate=1,rot_angle=45,cluster.value=c(-1,1, 0))

sf_bnd=Ushape1$bnd; sf_cluster1=Ushape1$cluster

f1_true=sf_cluster1$beta[cluster_true]
# f1_ho_true=sf_cluster1$beta[cluster_ho_true]
Y1 = f1_true + rnorm(n, 0, 0.001)
# Y_ho = f1_ho_true + rnorm(n_ho, 0, 0.001)

rot_angle = 45

# 
# ggplot() +
#   # geom_sf(data=clusterD) +
#   geom_point(aes(x = lon, y = lat, col=Y1), data = as.data.frame(coords)) +
#   scale_color_gradientn(colours = rainbow(8), name = 'z')

z = Y1


M = 1      # number of weak learners
k_max = 5   # maximum number of clusters per weak learner
mu = list() # initial values of mu
mstgraph_lst = list()  # initial spanning trees



cluster = matrix(1, nrow = n, ncol = M)  # initial cluster memberships
for(m in 1:M) {
  mu[[m]] = c(0)
  mstgraph_lst[[m]] = mstgraph0
}







clus_no = 3



# find lambda_s
nu = 3; q = 0.9
quant = qchisq(1-q, nu)
p = 5
q = 0
a= 1
b= 1
M = 1
hyperpar = c()
# hyperpar['sigmasq_mu'] = (0.5/(2*sqrt(M)))^2
hyperpar['sigmasq_mu'] = 0.5

# hyperpar['lambda_s'] = lambda_s
hyperpar['nu'] = nu
hyperpar['lambda_k'] = 4
hyperpar['M'] = M
# hyperpar['k_max'] = k_max
hyperpar['p'] = p
hyperpar['q'] = q
hyperpar['a'] = a
hyperpar['b'] = b




init_val = list()
init_val[['trees']] = mstgraph_lst
init_val[['mu']] = mu
init_val[['cluster']] = cluster
init_val[['sigmasq_y']] = 1



n = vcount(graph0)
# hyper-parameter
p = hyperpar['p']
q = hyperpar['q']
a = hyperpar['a']
b = hyperpar['b']
sigmasq_mu = hyperpar['sigmasq_mu']
lambda_s = hyperpar['lambda_s']
nu = hyperpar['nu']

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
eid_btw_mst_recent = list()
eid_cutted_lst = c()
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
# eid_btw_mst_dummy = list()
cluster_dummy1 = init_val[['cluster']]
csize_z_dummy1 = list()
k_dummy = as.numeric(apply(cluster, 2, max)) 
edge_status_dummy = as.matrix(rep("w", times = ecount(graph0)))

g = matrix(0, nrow = n, ncol = 1)  # n*1 matrix of fitted mu's

cluster_iter = cluster[,1]
g[, 1] = mu[[1]][cluster_iter]
csize[[1]] = Rfast::Table(cluster_iter)
# csize_z[[1]] = sum(z)
csize_z[[1]] = sum(z)
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
edge_status = as.matrix(rep("w", times = ecount(graph0)))

k_iter = k[1]
mstgraph_iter = mstgraph_lst[[1]]
edge_status_iter = edge_status[, 1]
cluster_iter = cluster[, 1]
csize_iter = csize[[1]]
csize_z_iter = csize_z[[1]]
subgraphs_iter = subgraphs[[1]]
eid_btw_mst_iter = eid_btw_mst[[1]]
e_iter = Y 
for (i in 1: (clus_no - 1)) {
  k_iter = k[1]
  mstgraph_iter = mstgraph_lst[[1]]
  edge_status_iter = edge_status[, 1]
  cluster_iter = cluster[, 1]
  csize_iter = csize[[1]]
  csize_z_iter = csize_z[[1]]
  subgraphs_iter = subgraphs[[1]]
  eid_btw_mst_iter = eid_btw_mst[[1]]
  e_iter = Y 
  
  split_res = splitCluster_common( mstgraph_iter, k_iter, subgraphs_iter, csize_iter)
  
  vid_new = split_res$vid_new; vid_old = split_res$vid_old
  parent_vid = split_res$parent_vid
  
  # compute log-prior ratio
  rd_new = 0.25
  
  ## output following the split
  update_res = updateSplit_shannon(split_res, subgraphs_iter, k_iter, csize_iter, csize_z_iter, 
                                   eid_btw_mst_iter, 
                                   cluster_iter, edge_status_iter, adj_list, adj_edge_list)
  
  print("split is accepted")
  # accept
  
  subgraphs[[1]] = update_res$subgraphs
  csize[[1]] = update_res$csize
  eid_btw_mst[[1]] = update_res$eid_btw_mst
  cluster[, 1] = update_res$cluster
  csize_z[[1]] = as.vector(tapply(z, cluster[, 1], sum))
  k[1] = k[1] + 1
  edge_status[, 1] = update_res$estatus
  # depth = depth + 1
  
  k_iter = k[1]
  mstgraph_iter = mstgraph_lst[[1]]
  edge_status_iter = edge_status[, 1]
  cluster_iter = cluster[, 1]
  csize_iter = csize[[1]]
  csize_z_iter = csize_z[[1]]
  subgraphs_iter = subgraphs[[1]]
  eid_btw_mst_iter = eid_btw_mst[[1]]
  e_iter = Y
  
}



M = 1    # number of weak learners
k_max = 5   # maximum number of clusters per weak learner
mu = list() # initial values of mu
for(m in 1:M) {
  mu[[m]] = c(0)
  mstgraph_lst[[m]] = mstgraph0
}




init_val = list()
init_val[['trees']] = mstgraph_lst
init_val[['mu']] = mu
init_val[['cluster']] = cluster
init_val[['sigmasq_y']] = 1
init_val[['edge_status']] = edge_status[, 1]
plotSpatGraphClust(mstgraph_lst[[1]], coords, cluster[, 1])


# standardize Y
std_res = standardize(Y)
Y_std = std_res$x; std_par = std_res$std_par

# find lambda_s
nu = 3; q = 0.9
quant = qchisq(1-q, nu)
lambda_s = quant * var(Y_std) / nu

hyperpar = c()

MCMC = 100000# MCMC iterations
BURNIN = 50000# burnin period length
THIN = 5   # thinning intervals

c = 0.85
q = 300
hyperpar = c()
# hyperpar['sigmasq_mu'] = (0.5/(2*sqrt(M)))^2
hyperpar['sigmasq_mu'] = var(Y_std)
hyperpar['lambda_s'] = lambda_s
hyperpar['nu'] = nu
hyperpar['M'] = M
hyperpar['q'] = q
hyperpar['c'] = c

#
MCMC_res = cohesion_forest_shrinkage_BA(Y_std, graph0, init_val, hyperpar, MCMC, BURNIN, THIN, prior_cluster, seed = 2024)
vec = MCMC_res$estimated_clus
plotSpatGraphClust(graph0, coords, vec)

