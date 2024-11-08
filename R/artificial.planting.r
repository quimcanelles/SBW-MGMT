#' Artificial planting
#' 
#' Based on function "Forest transitions"
#'
#' Determines the planted species of a cell after a clear cutting / SBW mortality.
#' 
#' @param land A \code{landscape} data frame with forest stand records in rows
#' @param target.cell.ids A vector of \code{cell.id} codes for those cells that may change stand composition
#' @param suitab A data frame with climatic suitability and soil suitability for potential colonizing species 
#' for each forest stand. It is the result of the function \code{suitability()}
#' @param scn A character indicating the type of transition according artificial planting scenario: \code{A} for "scn1", etc.
#' 
#' @return A vector with the name of the new species / group of species
#' 
#' @export
#' 
#' @examples
#' source("R/default.params.r")
#' source("R/suitability.r")
#' 
#' load("data/landscape.rda")
#' load("data/default.tables.rda")
#' params = default.params()
#' suitab = suitability(landscape, params, default.tables) 
#' 
#' land=landscape; target.cells=landscape$cell.id[runif(20,1,nrow(landscape))]
#' params$ap.rate=0.5; params$epn.rege.rate=0.8; params$pet.rege.rate=0.2
#' artificial.planting(landscape, landscape$cell.id[runif(20,1,nrow(landscape))], suitab, params)
#' 



artificial.planting = function(land, target.cells, suitab, params){   
  
  # 1. If target data.frame is empty
  if(length(target.cells)==0)
    return(numeric())
  
  
  # 2. Artificial planting according to SCN
  ap.target.cells = sample(target.cells, round(length(target.cells)*params$ap.rate))
  no.ap.cells <- target.cells[!(target.cells %in% ap.target.cells)]
  
  ## If target data.frame is empty
  if(length(ap.target.cells)==0)
    return(numeric())
  
  
  # 3. Probability of regeneration according to SCN
  prob.reg <- data.frame(spp=NA, potential.spp=NA, ptrans=NA)
  for(ispp in levels(land$spp)){
    aux <- data.frame("spp"=c(ispp,ispp),
                         "potential.spp"=c("EPN","PET"), 
                         "ptrans"= c(params$epn.rege.rate, params$pet.rege.rate))
                       
    prob.reg <- rbind(prob.reg, aux)
  
  }
  prob.reg <- prob.reg[-1,]
  
  
  # 4. Merge probability and suitability
  ## Base data frame
  subland = filter(land, cell.id %in% ap.target.cells)
  subland$spp = as.character(subland$spp)
  current.spp = subland$spp
  
  ## Join to the subland data frame the probability of transition to potential.spp (according to initial spp)
  ## Finally join the climatic-soil suitability index (i.e. modifier) of the potential spp
  subland1 = left_join(subland, prob.reg, by="spp") %>%
             left_join(suitab, by=c("cell.id", "potential.spp"))
  
  ## Stability criteria: if the species is present in the target location, then
  ## soil conditions are assumed to be optimal (not limiting)
  subland1$suit.soil[subland1$spp == subland1$potential.spp] = 1

  
  # 5. Final probability
  subland1$p = subland1$ptrans * pmin(subland1$suit.clim, subland1$suit.soil)
  
  ## Reshape the data frame, so we have a column for each potential species with 
  ## the corresponding transition probability (one row per target cell)
  ## Substitute dcast by "gather" or "spread" from tidyverse
  aux = reshape2::dcast(subland1, formula = cell.id ~ potential.spp, value.var = "p")
  
  ## Function that returns spp id according to probability x
  select.spp = function(x){
    if(sum(x)==0)
      return(0)
    id.spp = sample(1:length(x), 1, replace=FALSE, prob=x)
    return(id.spp)
  }
  
  
  # 6. Select a new spp according to the final probabilities and assign the corresponding species name
  ## If after all filters, p for all potential.spp is 0, the current species remains
  spp.names = names(aux)[2:ncol(aux)]
  id.spp = apply(aux[,2:ncol(aux)], 1, select.spp)
  new.spp = numeric(length=length(id.spp))
  new.spp[id.spp!=0] = spp.names[id.spp[id.spp!=0]]
  new.spp[id.spp==0] = as.character(current.spp[id.spp==0])
  
  
  subland$new.spp = new.spp

  
  ## Return the vector with the name of the new spp
  return(subland)
  
}

