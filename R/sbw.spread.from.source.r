########################################################################################################################
## Criteria that influence the SBW spreading to the neighbors:
## A. The species of the source cell
## B. The distance between source and target cells
## C. The position of the target cell with respect with the main wind direction
##########################################################################################

#' @example:: load("data/landscape.rda"); load("data/mask.rda"); source("R/default.params.r"); params = default.params();
#' land=landscape;nc=ncol(mask); wind_dir=params$wind_dir; radius=5; cell.size=2000

sbw.spread.from.source = function(land, nc, params, wind_dir, radius, cell.size, side = 1){  ##From Est (works for different angle between 0 and 360

  ## Select all the cells that are currently defoliated  
  source.cells = land[land$ny.def>0,]
  
  ## Find the neighbors of all source cells in a radius
  # 'radius' parameter in nn2 should be in meters, as the coordinates
  # the maximum number of neighbors to search for (k) is number of cells of a square of size 2* radius,
  # so this maximum number of cells is (2*num_of_cells_in_the_radius)^2
  k=((radius*1000/cell.size)*2)^2
  neighs_r = nn2(land[,c("x", "y")], land[land$cell.id %in% source.cells$cell.id, c("x", "y")],
                 k=k, radius=radius*1000, searchtype='radius')  
  nn.indx = neighs_r[[1]] # position of cells within the dataframe 'land'
  nn.dist = neighs_r[[2]]
  nn.dist[nn.indx==0] = NA
  nn.indx[nn.indx==0] = NA
  
  ## We need to find out 
  ## 1. the species of the source cell associated to each neighbor (spps)
  ## 2. if the neighbors are defoliated or not (defol)
  ## 3. the distance between the target and the source cells (dist)
  ## As the nn.indx may contain NAs we should do
    # # option A)
    # ids_vector = NULL
    # for(i in 1:nrow(nn.indx)){
    #   if(i==1){
    #     ids_vector = land$cell.id[nn.indx[i,]]
    #   } else{
    #     ids_vector = rbind(ids_vector, land$cell.id[nn.indx[i,]])
    #   }
    # }
    # # option B)
    # ids_values = land$cell.id[nn.indx]
    # ids_matrix = nn.indx
    # ids_matrix[!is.na(ids_matrix)] = ids_values[!is.na(ids_values)]

  # Choose option B for all the data frames
  ids.values = land$cell.id[nn.indx] # the real id of each cell (it is different from its position within the 'land' dataframe)
  spps.values = land$spp[nn.indx]
  defol.values = land$ny.def[nn.indx]
  ids = spps = defol = nn.indx  
  ids[!is.na(nn.indx)] = ids.values[!is.na(ids.values)]
  spps[!is.na(nn.indx)] = spps.values[!is.na(spps.values)]
  defol[!is.na(nn.indx)] = defol.values[!is.na(defol.values)]
  dist = matrix(nn.dist, nrow=nrow(source.cells))
  
  ## Maximum radius for each target cell and the corresponding weighted and linear distances
  ## We use the weight_dist  # EF: NOT USED. REMOVE?
    # r = dist[,radius]
    # weight_dist = round(exp(-dist^2/(r/2)^2),3)
    # linear_dist = round((r-dist)/r,3)
  
  ## ---> Criteria A
  ## Matrix with the cell.id of the neighbor cells that are not defoliated yet
  neigh_nodefol = ids * (defol==0) 
  
  ## Copy in this data frame the cell.id of the source cells (1st column)
  neigh_nodefol[,1] = ids[,1]
  neigh_nodefol = data.frame(neigh_nodefol)
  names(neigh_nodefol)[1] = "source"  ## EF: AND THE OTHER ARE THE TARGET CELLS?
  
  ## Relation between the source cell, the target cell and the species of the source cell
  source_wspp = pivot_longer(neigh_nodefol, cols=2:k, values_to="target") %>% dplyr::select(-name) %>% 
    filter(target!=0) %>% left_join(dplyr::select(land, cell.id, spp), by=c("source" = "cell.id")) %>% 
    rename(source_spp = spp) 
  
  ## Apply the weight corresponding to these species
  source_wspp$w_spp = ifelse(source_wspp$source_spp %in% c("EPN", "SAB"), 0.7, 0.3)
  
  ## ---> Criteria B
  ## Matrix with the distance of the neighbor cells that are not defoliated yet
  neigh_dist = dist * (defol==0) 
  
  ## Copy in this data frame the cell.id of the source cells (1st column)
  neigh_dist[,1] =  ids[,1]
  neigh_dist = data.frame(neigh_dist)
  names(neigh_dist)[1] = "source" 
  
  ## Relation between the source cell and the distance with the target cells (not explicitly indicated in this data frame,
  ## but we keep the order of them)
  source_dist =  pivot_longer(neigh_dist, cols=2:k, values_to="dist") %>% dplyr::select(-name) %>% filter(dist!=0)
  r = max(source_dist$dist)-2000
  
  ## Apply the weight corresponding to the distance
  source_dist$w_dist = exp(-(source_dist$dist-2000)^2/(r/2)^2)
    
  ## ---> Criteria C  
  x <- c(k:-k)
  G <- expand.grid(x,x)
  names(G) <- c("x","y")
  # Cell id
  G$z <- nc*G$y + G$x #c(1:length(x)^2)
  # Distance from focal cell
  G$w <- sqrt((side*G$x)^2 + (side*G$y)^2)
  # Each cell of the neighborhood is given an angle according to the convention: N=0, E=90, S=180, W=270
  G$theta <- atan2(G$x, G$y) + 2*pi*abs(sign(G$x))*(1/2)*(1-sign(G$x)) #in rads
  G$angle <- (180/pi)*G$theta #in degrees

  # Weight wind: 
  # cells in the direction of the prevailing wind (wind_dir) are given higher weights
  # Presenting wind_dir is a user-define number identical for all cells.
  # However, it could be a raster with the prevailing wind in each cell (verification that the code still works would be needed)
  # Weights decrease on cells away from the preferred wind direction
  G$wind_left <- (G$angle - wind_dir) %% 360
  G$wind_right <- (wind_dir - G$angle) %% 360
  G$wind <- G$wind_left*(G$wind_left<=180) + G$wind_right*(G$wind_right<180)
  G$w_wind <- (1/180)*G$wind  # test: matrix(G$w_wind, nrow = 5, byrow = TRUE)
  # Remove focal cell as w_wind is not correct at the center
  G <- filter(G, z!=0)
  
  
  ## The position of each cell are in: G$z. And G$z is like position = ids[i,] - ids[i,1] (for the source cell 'i')
  source_wind = source_wspp %>% dplyr::select(source, target) %>% mutate(position=target-source) %>% 
    left_join(dplyr::select(G, z, w_wind), by=c("position"="z"))
    
  ## Merge the three criteria and compute the final weight 
  res = cbind(source_wspp, source_dist[,-1], source_wind[,"w_wind"]) %>% 
    mutate(w_add = params$val_w_spp*w_spp + params$val_w_dist*w_dist + params$val_w_wind*w_wind, 
           w_multi=params$val_w_spp*w_spp * params$val_w_dist*w_dist * params$val_w_wind*w_wind) %>%  
    group_by(target) %>% summarise(final_w_add=sum(w_add), final_w_multi=sum(w_multi)) 
  
  ## Rescaling the final weight to [0,1]
  ## Mathieu 31/05/2024 > Do not rescaling, this is an artifact. Keep the "ecological" probability
  ## and fix a upper threshold above everything is defoliated
  res$spread_potential_add = (res$final_w_add-min(res$final_w_add))/(max(res$final_w_add)-min(res$final_w_add)) 
  res$spread_potential_multi = (res$final_w_multi-min(res$final_w_multi))/(max(res$final_w_multi)-min(res$final_w_multi)) 
  
  ## Questions:  
  ## Rescale the variables to the range [0,1] before applying any weight ??
  ## Do we want to apply a weight to each of these factors?
  ## Do we use a multiplicative or an additive formula?
  ## weight_factor_spp*w_spp + weight_factor_distance*w_dist + weight_factor_angel*w_wind
  
  return(res)  
}


