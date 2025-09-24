### Functions related to complex domains ###

library(mgcv)
library(fdaPDE)
library(igraph)
library(sf)
library(ggplot2)
library(GGally)
# library(ggnet)
library(network)
library(fields)

# library(rgl)
### Functions for 2-d constrained domains -----



# function to generate equally spaced points along a given line
refineLine <- function(start, end, n) {
  grids_x = seq(start[1], end[1], length.out = n)
  grids_y = seq(start[2], end[2], length.out = n)
  grids = cbind(grids_x, grids_y)
  return(grids)
}

### Constrained KNN -----

# function to get a KNN graph given a distance matrix
KNNGraph <- function(dist_mat, k_nn = 5, cross_dist = F, return_graph = T) {
  
  n1 = nrow(dist_mat); n2 = ncol(dist_mat)
  if(cross_dist) {
    # dist_mat is cross distance matrix
    adj_list = apply(dist_mat, 1, function(x) order(x)[1:k_nn])
  } else {
    adj_list = apply(dist_mat, 1, function(x) order(x)[2:(k_nn+1)])
  }
  adj_list = t(adj_list)
  
  if(return_graph & n1 == n2) {
    require(igraph)
    require(Matrix)
    
    i_all = c(); j_all = c(); x_all = c()
    for(i in 1:n1) {
      for(cidx in 1:k_nn) {
        i_all = c(i_all, i)
        j_all = c(j_all, adj_list[i, cidx])
        x_all = c(x_all, dist_mat[i, adj_list[i, cidx] ])
      }
    }
    adj_mat = sparseMatrix(i = i_all, j = j_all, x = x_all, dims = c(n1, n2))
    # get knn graph
    knngraph = graph_from_adjacency_matrix(adj_mat, mode = 'directed', weighted = T)
    knngraph = as.undirected(knngraph, mode = 'collapse', edge.attr.comb = 'first')
    return(knngraph)
  } else {
    return(lapply(1:n1, function(i) adj_list[i, ]))
  }
}





###############################







### -- ARTIFICIAL: DATA GENERATION --

### Functions related to complex domains ###
rm(list=ls())
library(igraph)
library(fields)
library(ggplot2)

library(fdaPDE)
library(igraph)
library(sf)
library(ggplot2)
library(GGally)
# library(ggnet)
library(network)



### Functions for 2-d constrained domain -----

#' function to generate sf U shape
#' genMesh(sfc_bnd, coords = NULL, n_ref = 100, n_bnd = NULL, ...)


#' function to generate bnd (sfc) and cluster (sfc format)
#' @param  rate  expansion rate
#' @param  rot_angle rotation angle
#' @param  cluster whether to generate partitioned Ubnd
#' @examples
#'
#' ##### Generate U-shape
#' # rate is the expansion rate
#' Ushape=gensfUbnd(rate=1,rot_angle=45,cluster.value=c(-1,0.5,1))
#' sf_bnd=Ushape$bnd; sf_cluster=Ushape$cluster; ggplot(Ushape$cluster)+geom_sf(aes(fill=as.factor(beta)))
#' bnd = data.frame(sf::st_coordinates(sf_bnd)[, 1:2]);colnames(bnd)=c('x','y')
#' plot(bnd)
#' #' #' #' ## Generate n uniform training locations
#' #' #' #'
#' n = 500
#' sf_coords=st_sample(sf_bnd,n);coords = st_coordinates(sf_coords);colnames(coords)=c('lon','lat')
#' 
#' plot(coords)
#' #' #' # generate n_ho uniform hold-out locations
#' n_ho = 200
#' sf_coords_ho=st_sample(sf_bnd,n_ho);coords_ho = st_coordinates(sf_coords_ho);colnames(coords_ho)=c('lon','lat')
#' plot(coords_ho)
#' 
#' cluster_true=unlist(st_within(sf_coords,sf_cluster))
#' cluster_ho_true=unlist(st_within(sf_coords_ho,sf_cluster))
#' 
#' #' #' #' ## Generate piecewise functions
#' f_true=sf_cluster$beta[cluster_true]
#' f_ho_true=sf_cluster$beta[cluster_ho_true]
#' Y = f_true + rnorm(n, 0, 0.1)
#' Y_ho = f_ho_true + rnorm(n_ho, 0, 0.1)
#' 
#' rot_angle = 45
#' #' #' # Generate piecewise smooth true functions
#' coords2 = rotate(coords[, 1], coords[, 2], angle = -rot_angle);coords2=data.frame(lon=coords2$x,lat=coords2$y)
#' coords2_ho = rotate(coords_ho[, 1], coords_ho[, 2], angle = -rot_angle);coords2_ho=data.frame(lon=coords2_ho$x,lat=coords2_ho$y)
#' 
#' f_true = mgcv::fs.test(coords2[,1], coords2[,2], b = 1.5, exclude = F)
#' f_ho_true = mgcv::fs.test(coords2_ho[,1], coords2_ho[,2], b = 1.5, exclude = F)
#' 
#' attributes(f_true) = NULL
#' attributes(f_ho_true) = NULL
#' 
#' f_true[cluster_true == 1] = -f_true[cluster_true == 1]
#' f_true[cluster_true == 2] = 2 * f_true[cluster_true == 2]
#' 
#' f_ho_true[cluster_ho_true == 1] = -f_ho_true[cluster_ho_true == 1]
#' f_ho_true[cluster_ho_true == 2] = 2 * f_ho_true[cluster_ho_true == 2]
#' ############
#' Y = f_true + rnorm(n, 0, 0.1)
#' Y_ho = f_ho_true + rnorm(n_ho, 0, 0.1)
#' 
#' ggplot() +
#'   # geom_sf(data=clusterD) +
#'   geom_point(aes(x = lon, y = lat, col=Y), data = as.data.frame(coords)) +
#'   scale_color_gradientn(colours = rainbow(5), name = 'Y')
#' 
#' #' #' #' ## Generate graphs on meshes
#' #' #' #'
#' # mesh=genMesh(sf_bnd, n_ref = 100,graph=TRUE,type='hexagonal') # type=hexagonal or regular
#' # sf_mesh=mesh$mesh;
#' # g0=mesh$g
#' # plotGraph(coords,g0)
#' 
#' #' #' ## Generate graphs on points
#' mesh = gen2dMesh(coords, bnd)
#' g0 = constrainedDentri(n, mesh, gaurantee_connected = T)
#' E(g0)$eid = as.integer(1:ecount(g0))  # edge id
#' V(g0)$vid = as.integer(1:vcount(g0))  # vertex id
#' E(g0)$weight=runif(ecount(g0))
#' plotGraph(coords, g0)
#' mstgraph = mst(g0)
#' plotGraph(coords, mstgraph)
#' vcount(mstgraph)

gensfUbnd=function(rate=1,rot_angle=45,cluster.value=NULL){
  require(sf)
  ubnd = genURegion(angle = rot_angle)
  outer = rate*as.matrix(cbind(ubnd$x, ubnd$y))
  sfc_ubnd = st_polygon(list(as.matrix(rbind(outer, outer[1,])))) %>% st_geometry()
  if(!is.null(cluster.value)){
    centroid = st_point(c(0, 0))
    sfc_circle = centroid %>% st_buffer(rate* 0.9) %>% st_geometry()
    sfc_clust12 = st_cast(st_difference(sfc_ubnd, sfc_circle), 'POLYGON')
    sfc_clust3 = st_intersection(sfc_ubnd, sfc_circle)
    beta_true_uniq=cluster.value;
    sf_uclust = st_sf(beta = beta_true_uniq,geometry = c(sfc_clust12, sfc_clust3))
    return(list(bnd=sfc_ubnd, cluster=sf_uclust))
  }else{
    return(list(bnd=sfc_ubnd))
  }
}

gensfUbnd1=function(rate=1,rot_angle=45,cluster.value=NULL){
  require(sf)
  ubnd = genURegion(angle = rot_angle)
  outer = rate*as.matrix(cbind(ubnd$x, ubnd$y))
  sfc_ubnd = st_polygon(list(as.matrix(rbind(outer, outer[1,])))) %>% st_geometry()
  if(!is.null(cluster.value)){
    centroid = st_point(c(0, 0))
    sfc_circle = centroid %>% st_buffer(rate* 0.9) %>% st_geometry()
    sfc_clust12 = st_cast(st_difference(sfc_ubnd, sfc_circle), 'POLYGON')
    sfc_clust3 = st_intersection(sfc_ubnd, sfc_circle)
    beta_true_uniq=cluster.value;
    sf_uclust = st_sf(beta = beta_true_uniq,geometry = c(sfc_clust12, sfc_clust3))
    
    centroid1 = st_point(c(0, 0))
    sfc_circle1 = centroid1 %>% st_buffer(rate* 1.3) %>% st_geometry()
    sfc_clust121 = st_cast(st_difference(sfc_ubnd, sfc_circle1), 'POLYGON')
    sfc_clust31 = st_intersection(sfc_ubnd, sfc_circle1)
    beta_true_uniq1=cluster.value;
    sf_uclust1 = st_sf(beta = beta_true_uniq1, geometry = c(sfc_clust121, sfc_clust31))
    
    return(list(bnd=sfc_ubnd, cluster=sf_uclust, bnd1 = sfc_ubnd, cluster1 = sf_uclust1))
  }else{
    return(list(bnd=sfc_ubnd))
  }
}


gensfUbnd2=function(rate=1,rot_angle=45,cluster.value=NULL){
  require(sf)
  ubnd = genURegion(angle = rot_angle)
  outer = rate*as.matrix(cbind(ubnd$x, ubnd$y))
  sfc_ubnd = st_polygon(list(as.matrix(rbind(outer, outer[1,])))) %>% st_geometry()
  if(!is.null(cluster.value)){
    centroid = st_point(c(0, 0))
    sfc_circle = centroid %>% st_buffer(rate* 1.3) %>% st_geometry()
    sfc_clust12 = st_cast(st_difference(sfc_ubnd, sfc_circle), 'POLYGON')
    sfc_clust3 = st_intersection(sfc_ubnd, sfc_circle)
    beta_true_uniq=cluster.value;
    sf_uclust = st_sf(beta = beta_true_uniq,geometry = c(sfc_clust12, sfc_clust3))
    
    centroid1 = st_point(c(0, 0))
    sfc_circle1 = centroid1 %>% st_buffer(rate* 1.3) %>% st_geometry()
    sfc_clust121 = st_cast(st_difference(sfc_ubnd, sfc_circle1), 'POLYGON')
    sfc_clust31 = st_intersection(sfc_ubnd, sfc_circle1)
    beta_true_uniq1=cluster.value;
    sf_uclust1 = st_sf(beta = beta_true_uniq1, geometry = c(sfc_clust121, sfc_clust31))
    
    return(list(bnd=sfc_ubnd, cluster=sf_uclust, bnd1 = sfc_ubnd, cluster1 = sf_uclust1))
  }else{
    return(list(bnd=sfc_ubnd))
  }
}


#' function to rotate locations counterclockwise (in degrees)
#' @param x input x.
#' @param y input y
#' @returns roate
#' @examples
#' rotate(x,y, 1)
#' rotate(x,y,45)
rotate <- function(x, y, angle = 0) {
  angle = angle / 180 * pi
  x_rot = x * cos(angle) - y * sin(angle)
  y_rot = x * sin(angle) + y * cos(angle)
  return(list(x = x_rot, y = y_rot))
}

#' function to rotate locations counterclockwise (in degrees)
#' @param coords input matrix or data.frame.
#' @returns roate
#' @examples
#' rotate2(coords, 1)
#' rotate2(coodds,45)
rotate2 <- function(coords, angle = 0) {
  res = rotate(coords[, 1], coords[, 2], angle)
  coords[, 1] = res$x; coords[, 2] = res$y
  return(coords)
}

#' function to check if points are inside a polygon
insidePolygon <- function(bnd, x, y) {
  bnd$x = round(bnd$x, 9)
  bnd$y = round(bnd$y, 9)
  pg = st_polygon(list( cbind(bnd$x, bnd$y) ))
  
  coords = as.data.frame(cbind(x, y))
  pt = st_as_sf(coords, coords = c(1, 2))
  return( st_intersects(pt, pg, sparse = F)[, 1] )
}


genURegion <- function (r0 = 0.1, r = 0.3, l = 3, n_theta = 20, n_line = 30, angle = 0) {
  
  rr = r + (r - r0)
  theta <- seq(pi, pi/2, length.out = n_theta)
  x = rr * cos(theta)
  y = rr * sin(theta)
  
  x = c(x, seq(l/n_line, l - l/n_line, length.out = n_line))
  y = c(y, rep(rr, n_line))
  
  theta = seq(pi/2, -pi/2, length = n_theta)
  x = c(x, (r - r0) * cos(theta) + l)
  y = c(y, (r - r0) * sin(theta) + r)
  
  x = c(x, seq(l - l/n_line, l/n_line, length.out = n_line))
  y = c(y, rep(r0, n_line))
  
  theta = seq(pi/2, pi, length.out = round(n_theta * r0 / (2*r - r0)))
  x = c(x, r0 * cos(theta))
  y = c(y, r0 * sin(theta))
  n = length(x)
  x = c(x, x[n:1])
  y = c(y, -y[n:1])
  
  # rotate
  rot = rotate(x, y, angle)
  return(rot)
}

#' function to generate uniform locations in a region
#' @returns  (regular df) random coordinates in bnd
#' @examples
#' coords = genLocations(n, ubnd)
genLocations <- function(n, bnd = NULL) {
  
  if(is.null(bnd)) {
    # generate locations in a unit square
    coords = cbind(runif(n), runif(n))
  } else {
    x_min = min(bnd$x); x_max = max(bnd$x)
    y_min = min(bnd$y); y_max = max(bnd$y)
    
    x = runif(3*n, x_min, x_max); y = runif(3*n, y_min, y_max)
    idx_in = insidePolygon(bnd, x, y)
    coords = cbind(x[idx_in], y[idx_in])
    
    if(nrow(coords) >= n) {
      coords = coords[1:n, ]
    } else {
      for(i in 1:(n - nrow(coords))) {
        x = runif(1, x_min, x_max); y = runif(1, y_min, y_max)
        while(!insidePolygon(bnd, x, y))
          x = runif(1, x_min, x_max); y = runif(1, y_min, y_max)
          coords = rbind(coords, c(x, y))
      }
    }
  }
  colnames(coords) = c('lon', 'lat')
  return(coords)
}

#' function to generate grids in a region
#' Generate regular grid in a region
#' @examples genGrids(5,5,bnd)
genGrids <- function(nx, ny, bnd = NULL) {
  
  if(is.null(bnd)) {
    # generate locations in a unit square
    coords_grid = expand.grid(seq(0, 1, length.out = nx), seq(0, 1, length.out = ny))
    coords_grid = as.matrix(coords_grid)
  } else {
    x_min = min(bnd$x); x_max = max(bnd$x)
    y_min = min(bnd$y); y_max = max(bnd$y)
    
    # to avoid points on boundaries
    x_buffer = (x_max - x_min) * 1e-5
    y_buffer = (y_max - y_min) * 1e-5
    
    grid_x = seq(x_min + x_buffer, x_max - x_buffer, length.out = nx)
    grid_y = seq(y_min + y_buffer, y_max - y_buffer, length.out = ny)
    coords_grid = as.matrix(expand.grid(grid_x, grid_y))
    
    x = coords_grid[, 1]; y = coords_grid[, 2]
    idx_in = insidePolygon(bnd, x, y)
    coords_grid = coords_grid[idx_in, ]
  }
  colnames(coords_grid) = c('lon', 'lat')
  return(coords_grid)
}

#' function to create a 2d mesh
#' note the first and last boundary nodes are the same
#' @examples
#' mesh = gen2dMesh(coords, ubnd)
#' g0 = constrainedDentri(n, mesh, gaurantee_connected = T)
#' E(g0)$eid = as.integer(1:ecount(g0))  # edge id
#' V(g0)$vid = as.integer(1:vcount(g0))  # vertex id
#' E(g0)$weight=runif(ecount(g0))

gen2dMesh <- function(coords, bnd, ...) {
  
  n = nrow(coords); n_bnd = length(bnd$x) - 1
  coords_all = rbind(coords, cbind(bnd$x, bnd$y)[1:n_bnd, ])
  
  # get boundary segments
  segments = cbind( (n+1):(n+n_bnd), c((n+2):(n+n_bnd), n+1) )
  
  mesh = create.mesh.2D(coords_all, segments = segments, ...)
  mesh$n_int = n  # number of interior nodes
  # edges that connect boundary nodes
  mesh$bnd_edges = apply(mesh$edges, 1, FUN = function(x) any(x > n))
  return(mesh)
}



#' function to connect each component to its nearest neighbor such that we can obtain a connected graph from a disconnected one
connectGraph <- function(graph, dist_mat) {
  n = vcount(graph)
  connected_comp = components(graph)
  memberships = connected_comp$membership
  n_comp = connected_comp$no
  if (n_comp == 1)
    return(graph)
  
  edge_list = as_edgelist(graph)
  for (i in 1:(n_comp - 1)) {
    # merge one component to its nearest component
    idx_1 = memberships == 1
    idx_2 = memberships != 1
    cdist_mat = dist_mat[idx_1, idx_2, drop = F]
    idx_min = which(cdist_mat == min(cdist_mat), arr.ind = TRUE)[1, ]
    
    # find the pair of vertices that has minimum distance
    vid_1 = c(1:n)[idx_1][ idx_min[1] ]
    vid_2 = c(1:n)[idx_2][ idx_min[2] ]
    
    # connect vid_1 with vid_2
    edge_new = sort(c(vid_1, vid_2))
    edge_list = rbind(edge_list, edge_new)
    
    # merge two components where vid_1 and vid_2 lie in
    comp_id_merged = memberships[vid_2]
    memberships[memberships == comp_id_merged] = 1
  }
  
  return(graph_from_edgelist(edge_list, directed = F))
}

# function to generate equally spaced points along a given line
refineLine <- function(start, end, n) {
  grids_x = seq(start[1], end[1], length.out = n)
  grids_y = seq(start[2], end[2], length.out = n)
  grids = cbind(grids_x, grids_y)
  return(grids)
}


#' function to estimate geodesic distance on a 2-d constrained domain
#' by using shortest path length between two vertices on a dense graph
#' Calculate geodesic distance
#' @examples coords_all = rbind(coords, coords_ho)
#' dist_res = gdist(coords_all, ubnd)
#' dist_res = gdist(coords=coords,coords2=coords_ho,bnd=ubnd)
gdist <- function(coords, bnd, nx = 50, ny = 50, k_nn = 8,
                  coords2 = NULL, crs = NULL, return_both = F) {
  
  require(igraph)
  require(fields)
  require(FNN)
  
  # generate grids
  coords_grids = genGrids(nx, ny, bnd)
  n_grids = nrow(coords_grids)
  
  # get nearest neighbor graph
  coords_all = rbind(coords_grids, coords)
  if(!is.null(coords2))
    coords_all = rbind(coords_all, coords2)
  if (is.null(crs)) {
    nn_res_grid = get.knn(coords_grids, k = k_nn)
    nn_res_nongrid = get.knnx(coords_grids, coords_all[-(1:n_grids), ], k = k_nn)
  } else {
    coords_all_sf = st_as_sf(as.data.frame(coords_all), coords = c('lon', 'lat'), crs = crs)
    coords_all_sf = st_geometry(coords_all_sf)
    nn_res_grid = st_knn(coords_all_sf[1:n_grids], k = k_nn)
    nn_res_nongrid = st_knn(coords_all_sf[1:n_grids], coords_all_sf[-(1:n_grids)], k = k_nn)
  }
  knngraph = constrainedKNN(nn_res_grid, nn_res_nongrid, coords_all, bnd)
  
  # check if the graph is connected
  if(igraph::components(knngraph)$no != 1)
    stop("Disconnected kNN graph; try larger 'k'")
  
  # compute geodesic distance matrix
  from = (n_grids + 1):(n_grids + nrow(coords))
  n1 = ifelse(is.null(coords2), 0, nrow(coords))
  to = (n_grids + n1 + 1):nrow(coords_all)
  cdist_mat = igraph::distances(knngraph, v = from, to = to)
  
  if(return_both & !is.null(coords2)) {
    # return distance matrix of coords as well
    from = (n_grids + 1):(n_grids + nrow(coords))
    dist_mat = igraph::distances(knngraph, v = from, to = from)
    return(list('dist' = dist_mat, 'cdist' = cdist_mat))
  } else {
    # return cross distance matrix only
    return(cdist_mat)
  }
}

#' function to obtain intrinsic coordinates on a rotated U-shape domain
#' adapted from R package "mgcv"
intrinsicCoordsU <- function(coords, r0 = 0.1, r = 0.5, l = 3, angle = 0) {
  # rotate to horizontal
  coords = rotate2(coords, -angle)
  
  n = nrow(coords)
  q = pi * r/2 # 1/2 length of semi-circle part of centre curve
  a = rep(0, n); d = rep(0, n) # along and distance to arrays
  x = coords[, 1]; y = coords[, 2]
  
  ## convert x,y to along curve and distance to curve (a,d)
  ## co-ordinates. 0 distance along is at (x = -r, y = 0)
  
  ind = (x >= 0) & (y > 0)
  a[ind] = q + x[ind]
  d[ind] = y[ind]-r
  
  ind = (x >= 0) & (y <= 0)
  a[ind] = -q - x[ind]
  d[ind] = -r - y[ind]
  
  ind = (x < 0)
  a[ind] = -atan(y[ind] / x[ind]) * r
  d[ind] = sqrt(x[ind]^2 + y[ind]^2) - r
  
  return(cbind(a, d))
}

#' function to simulate a GP on a U-shape domain using intrinsic coordinates
#' get intrinsic coordinate (a, d)
#' Generate GP process on U-shape domain
#' @examples X_all = matrix(0, nrow = n + n_ho, ncol = p)
#' X_all = simGpU(p, coords_all, angle = rot_angle)
simGpU <- function(n, coords, r0 = 0.1, r = 0.5, l = 3, angle = 0) {
  
  coords_in = intrinsicCoordsU(coords, r0, r, l, angle)
  
  # length scale
  phi_a = 1; phi_d = 1
  coords_in[, 1] = coords_in[, 1] / phi_a
  coords_in[, 2] = coords_in[, 2] / phi_d
  # get covariance matrix based on (a, d)
  dist_mat = as.matrix(dist(coords_in))
  corr = exp(-dist_mat)
  varcov = corr
  chol_varcov = t(chol(varcov))
  
  # simulate GP
  X = matrix(0, nrow = nrow(coords_in), ncol = n)
  for (j in 1:n)
    X[, j] = chol_varcov %*% rnorm(nrow(coords_in))
  return(X)
}

#' evaluate test function on a horizontal U-shape domain
#'
evalFunU <- function(coords, X, r0 = 0.1, r = 0.5, l = 3, angle = 0) {
  coords_in = intrinsicCoordsU(coords, r0, r, l, angle)
  a = coords_in[, 1]; d = coords_in[, 2]
  f = a * X[, 1] + d ^ 2
  return(f)
}

### Constrained KNN/RNN on 2-d constrained domain -----

# helper function to check if 3 points are in counterclockwise order
ccw <- function(A, B, C) {
  return( (C[2] - A[2]) * (B[1] - A[1]) > (B[2] - A[2]) * (C[1] - A[1]) )
}

# function to ckeck if two segments AB and CD intersect
segmentIntersect <- function(A, B, C, D) {
  return( ccw(A, C, D) != ccw(B, C, D) & ccw(A, B, C) != ccw(A, B, D) )
}

# helper function to get RNN list for one location
# return idx in d that is less than r
# if no idx is found, return which.min(d)
getRNNList <- function(d, r) {
  idx = which(d <= r)
  if(length(idx) == 0)
    idx = which.min(d)
  return(idx)
}

#' function to obtain constrained RNN given geodesic distances
#' Create igraph RNN
#' @examples
#' coords_all = rbind(coords, coords_ho)
#' dist_res = gdist(coords_all, ubnd)
#' dist_res = gdist(coords=coords,coords2=coords_ho,bnd=ubnd)
#' constrainedRNN(dist_res,bnd=ubnd)
constrainedRNN <- function(dist_mat, r = 1, bnd = NULL, return_graph = F) {
  
  n1 = nrow(dist_mat); n2 = ncol(dist_mat)
  nn_list = apply(dist_mat, 1, getRNNList, r = r)
  
  if (return_graph & n1 == n2) {
    require(igraph)
    require(Matrix)
    i_all = c(); j_all = c(); x_all = c()
    
    for(i in 1:n1) {
      for(idx in 1:length(nn_list[[i]])) {
        j = nn_list[[i]][idx]
        if (i == j) next
        i_all = c(i_all, i)
        j_all = c(j_all, j)
        x_all = c(x_all, dist_mat[i, j])
      }
    }
    adj_mat = sparseMatrix(i = i_all, j = j_all, x = x_all, dims = c(n1, n2))
    # get rnn graph
    rnngraph = graph_from_adjacency_matrix(adj_mat, mode = 'directed', weighted = T)
    rnngraph = as.undirected(rnngraph, mode = 'collapse', edge.attr.comb = 'first')
    
    # check if the graph is connected
    if(components(rnngraph)$no > 1)
      warning("Disconnected RNN graph; try larger 'r'")
    
    return(rnngraph)
  } else {
    return(nn_list)
  }
}

#' function to get a KNN graph given a distance matrix
#' Create igraph KNN
#' @examples
#' coords_all = rbind(coords, coords_ho)
#' dist_res = gdist(coords_all, ubnd)
#' dist_res = gdist(coords=coords,coords2=coords_ho,bnd=ubnd)
#' KNNGraph(dist_res,k_nn=5)
#' or use cccd function nng()
KNNGraph <- function(dist_mat, k_nn = 5, cross_dist = F, return_graph = T) {
  
  n1 = nrow(dist_mat); n2 = ncol(dist_mat)
  if(cross_dist) {
    # dist_mat is cross distance matrix
    adj_list = apply(dist_mat, 1, function(x) order(x)[1:k_nn])
  } else {
    adj_list = apply(dist_mat, 1, function(x) order(x)[2:(k_nn+1)])
  }
  adj_list = t(adj_list)
  
  if(return_graph & n1 == n2) {
    require(igraph)
    require(Matrix)
    
    i_all = c(); j_all = c(); x_all = c()
    for(i in 1:n1) {
      for(cidx in 1:k_nn) {
        i_all = c(i_all, i)
        j_all = c(j_all, adj_list[i, cidx])
        x_all = c(x_all, dist_mat[i, adj_list[i, cidx] ])
      }
    }
    adj_mat = sparseMatrix(i = i_all, j = j_all, x = x_all, dims = c(n1, n2))
    # get knn graph
    knngraph = graph_from_adjacency_matrix(adj_mat, mode = 'directed', weighted = T)
    knngraph = as.undirected(knngraph, mode = 'collapse', edge.attr.comb = 'first')
    return(knngraph)
  } else {
    return(lapply(1:n1, function(i) adj_list[i, ]))
  }
}

#' function to build constrained kNN graph on a 2-d constrained domain given NNs
#' remove NN pairs if the corresponding edge across the boundary
#' nn_res: list returned by FNN::get.knn()
constrainedKNN <- function(nn_res_grid, nn_res_nongrid, coords, bnd) {
  
  ##
  ##
  require(igraph)
  require(Matrix)
  
  # get boundary segments
  n_bnd = length(bnd$x)
  bnd_segments = cbind(bnd$x[-n_bnd], bnd$y[-n_bnd], bnd$x[-1], bnd$y[-1])
  
  # combine nn_res
  nn_res = list(nn.index = rbind(nn_res_grid$nn.index, nn_res_nongrid$nn.index),
                nn.dist = rbind(nn_res_grid$nn.dist, nn_res_nongrid$nn.dist))
  
  # get sparse adjacency matrix
  # and remove edges across the boundary
  i_all = c(); j_all = c(); x_all = c()
  n = nrow(nn_res$nn.index)
  k_nn = ncol(nn_res$nn.index)  # number of nn
  for(i in 1:n) {
    coords_i = coords[i, ]
    cnt_nn = 0  # count how many valid nn (i.e., not crossing boundary)
    for(cidx in 1:k_nn) {
      j = nn_res$nn.index[i, cidx]
      coords_j = coords[j, ]
      
      # check if segment (i, j) crosses domain boundary
      intersect_bnd = apply(bnd_segments, 1, function(seg) segmentIntersect(
        coords_i, coords_j, seg[1:2], seg[3:4]
      ))
      
      if(!any(intersect_bnd)) {
        # add to adj matrix
        i_all = c(i_all, i); j_all = c(j_all, j)
        x_all = c(x_all, nn_res$nn.dist[i, cidx])
        cnt_nn = cnt_nn + 1
      }
    }
    
    # check if the current node has any neighbor
    if(cnt_nn == 0)
      stop("Disconnected kNN graph; try larger 'k'")
  }
  
  adj_mat = sparseMatrix(i = i_all, j = j_all, x = x_all, dims = c(n, n))
  
  # get knn graph
  knngraph = graph_from_adjacency_matrix(adj_mat, mode = 'directed', weighted = T)
  knngraph = as.undirected(knngraph, mode = 'collapse', edge.attr.comb = 'first')
  return(knngraph)
}


### Plotting functions -----

#' function to plot complex domain
#' @examples
#' geom_boundary(ubnd)
geom_boundary <- function(bnd, ...) {
  
  n = length(bnd$x)
  segments = cbind(bnd$x[-n], bnd$y[-n], bnd$x[-1], bnd$y[-1])
  segments = data.frame(segments)
  names(segments) = c('x1', 'y1', 'x2', 'y2')
  return(geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = segments, ...))
}

#' function to plot ellipse
geom_ellipse <- function(rx, ry, xc, yc, color = "red", size = 0.3, ...) {
  x = xc + rx * cos(seq(0, pi, length.out=200))
  ymax = yc + ry * sin(seq(0, pi, length.out=200))
  ymin = yc + ry * sin(seq(0, -pi, length.out=200))
  annotate("ribbon", x = x, ymin = ymin, ymax = ymax, color = color, size = size, fill = NA, ...)
}



### Util functions -----

# function to get knn of sfc objects
st_knn <- function(data, query = NULL, k = 1) {
  require(sf)
  if (is.null(query))
    dist_mat = st_distance(data)
  else
    dist_mat = st_distance(query, data)
  dist_mat = as.matrix(dist_mat)
  nn_index = t(apply(dist_mat, 1, WhichNMin, n = k))
  nn_dist = t(sapply(1:nrow(nn_index), function(i) dist_mat[i, nn_index[i, ]]))
  return(list('nn.index' = nn_index, 'nn.dist' = nn_dist))
}


# function to obtain indices of n smallest values in a vector
WhichNMin <- function(x, n = 1) {
  n = min(n, length(x))
  threshold = sort(x, partial = n)[n]
  which_n_min = which(x <= threshold)
  # deal with ties
  if(length(which_n_min) > n) {
    x = x[which_n_min]
    which_n_min = which_n_min[order(x)[1:n]]
  }
  return(which_n_min)
}




### Plotting functions -----

# function to plot complex domain


# plotting functions from SCC
plotGraph_old <- function(coords, graph, title = NULL){
  require(ggplot2)
  edgelist = get.edgelist(graph) 
  edgedata = data.frame(coords[edgelist[,1 ], ], coords[edgelist[, 2], ])
  colnames(edgedata) = c("X1", "Y1", "X2", "Y2")
  
  ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data = edgedata, size = 0.5, colour = "red") +
    labs(title = title, x = "lon", y = "lat")+
    theme(plot.title = element_text(hjust = 0.5))
}

# function to plot spatial field
plotField <- function(coords, Data, col_lim = NULL, legend_name = NULL, title = NULL, colors = rainbow(10)){
  if(missing(col_lim)) {col_lim = range(Data)}
  ggplot() + 
    geom_tile(data = data.frame(coords), aes(lon, lat, fill = Data)) +
    scale_fill_gradientn(colours = colors, limits = col_lim, name = legend_name, na.value = 'white') +
    ggtitle(title) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          # legend.title=element_blank(), 
          plot.title = element_text(hjust = 0.5, size = 10),
          legend.text = element_text(size = 9, hjust = 0, margin = margin(l = 3)))
}








#' plot Graph function
#' convert igraph to network
#' @examples
#' plotGraph(coords,g0,label=1:60,label.size=3)
#' plotGraph(coords,g0,directed=TRUE, label=1:60, label.size=2)
#' help(ggnet2) for more details

plotGraph=function(coords,g0,directed=FALSE,...){
  library(GGally)
  
  net = network(get.edgelist(g0), directed = directed,matrix.type='edgelist')
  ## assign coordinates as node attributes
  
  net%v%'x'=coords[,1];net%v%'y'=coords[,2];
  if(!directed){
    ggnet2(net, mode = c("x", "y"),color='green',size=2, ...)
  }else{
    ggnet2(net, mode = c("x", "y"),color='green',size=2,arrow.size=2,arrow.gap=0.025,...)
  }
}

#' plot Attributed Graph function
#' convert igraph to network
#' @examples
#' plotGraphData(coords,g0,Data)

plotGraphData=function(coords, graph, Data, title = NULL, col_lim = NULL, index=NULL){
  edgelist<-get.edgelist(graph)
  edgedata <- data.frame(coords[edgelist[,1],], coords[edgelist[,2],])
  colnames(edgedata) <- c("X1","Y1","X2","Y2")
  if(missing(col_lim)) {col_lim = NA}
  if(missing(index)){
    ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edgedata, size = 0.5, colour="grey") +
      geom_point(data=data.frame(coords), aes(lon, lat, color = Data))+
      scale_color_gradientn(colours = rainbow(10), limits = col_lim) +
      ggtitle(title) +
      theme(legend.title =element_blank(),plot.title = element_text(hjust = 0.5))
  }
  else{
    mtsub=data.frame(coords[index,]);colnames(mtsub)=c('X1','Y1');
    ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edgedata, size = 0.5, colour="grey") +
      geom_point(data=data.frame(coords), aes(lon, lat, color = Data))+scale_color_gradientn(colours
                                                                                             = rainbow(5)) +
      ggtitle(title)+
      theme(legend.title =element_blank(),plot.title = element_text(hjust = 0.5))+
      geom_text(data = mtsub, aes(x = X1,y=Y1, label =
                                    rownames(mtsub)), hjust = 0, size = 4)
    
  }
}




### Redefine functions related to FEM ###
### code adapted from fdaPDE ###

# function to returns the index of the triangle containing the point
# returns multiple triangles for points on edges if all_tri == TRUE
# code adapted from fdaPDE
R_insideIndex2 <- function (mesh, location, all_tri = FALSE)
{
  #  insideIndex returns the index of the triangle containing the point
  # (X,Y) if such a triangle exists, and NaN otherwise.
  #  TRICOEF may have already been calculated for efficiency,
  #  but if the function is called with four arguments, it is calculated.
  
  
  eps=2.2204e-016
  small = 10000*eps
  
  nodes = mesh$nodes
  triangles = mesh$triangles
  X = location[1]
  Y = location[2]
  
  ntri   = dim(triangles)[[1]]
  indtri   = matrix(1:ntri,ncol=1)
  
  #  compute coefficients for computing barycentric coordinates if needed
  
  tricoef = fdaPDE:::R_tricoefCal(mesh)
  
  #  compute barycentric coordinates
  r3 = X - nodes[triangles[,3],1]
  s3 = Y - nodes[triangles[,3],2]
  lam1 = ( tricoef[,4]*r3 - tricoef[,2]*s3)
  lam2 = (-tricoef[,3]*r3 + tricoef[,1]*s3)
  lam3 = 1 - lam1 - lam2
  
  #  test these coordinates for a triple that are all between 0 and 1
  int  = (-small <= lam1 & lam1 <= 1+small) & 
    (-small <= lam2 & lam2 <= 1+small) & 
    (-small <= lam3 & lam3 <= 1+small)
  
  #  return the index of this triple, or NaN if it doesn't exist
  indi = indtri[int]
  if (length(indi)<1)
  {
    ind = NA
  }else if(all_tri) {
    ind = indi
  } else{
    ind = min(indi)
  }
  
  ind
}



# function to obtain a constrained Delaunay triangulation graph from a mesh
constrainedDentri <- function(n, mesh, threshold = 1000) {
  coords = mesh$nodes[1:n, ]
  
  # drop edges that connect boundary nodes
  rid_drop = mesh$bnd_edges
  edge_list = mesh$edges[!rid_drop, ]
  
  # compute edge length
  distance = sqrt( rowSums((coords[edge_list[, 1], ] - coords[edge_list[, 2], ]) ^ 2) )
  
  rid_drop = distance > threshold
  edge_list = edge_list[!rid_drop, ]
  distance = distance[!rid_drop]
  
  graph0 = graph_from_edgelist(edge_list, directed = F)
  E(graph0)$weight = distance
  return(graph0)
}


