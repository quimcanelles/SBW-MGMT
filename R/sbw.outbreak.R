
#@ exemple:: load("data/landscape.rda"); land=landscape; source("R/default.params.r"); preoutbreak=1; outbreak=1; calm=1; collapse=1; 
#duration.last.outbreak=1; params = default.params(); load("data/default.tables.rda"); tbls=default.tables;load("data/mask.rda")

sbw.outbreak = function(land, params, tbls, preoutbreak=1, outbreak=1, calm=1, collapse=1, 
                        duration.last.outbreak=1, mask=NA){

  # 0.  Fix function  
  `%notin%` = Negate(`%in%`)

  ## Determine cell resolution in meters
  cell.size = (mask@extent[2] - mask@extent[1])/ncol(mask) 
  
  ## 1. Defliation in the different phases of the SBW outbreak
  ## 08/04/24: In the pre-epidemic phase, the epicenters should be only SAB or EPN. 
  ## However, cells in the bubble around an epicenter could be any deciduous. 
  ## The intensity should be, of course, less in these non-host cells.
  cat("   Defoliation in the ")
  if(preoutbreak>0){
    cat("pre-epidemic phase ", "\n")

    ## "potential_source" are the cells where sbw can start.
    ## All the function applies to "potential_source" cells instead of "land" cells.
    ## That means, that sbw outbreak only can start in SAB and EPN location with optimal climatic conditions for sbw
    ## (0.5<temp<2.8), with time_since_last_outbread >= 30, and in cells that are not currently defoliated
    potential_source = filter(land, ny.def0>=5, tssbw>=30, spp %in% c("SAB", "EPN"), temp>0.5, temp<2.8)

    ## First epicenters
    epicenter = sample(potential_source$cell.id, size=rdunif(1,3,6-preoutbreak), replace=F, #QC: size = the later preoutbreak the more new epicenters. Callibrated with historical data.
                          prob=potential_source$ny.def0*(1200-potential_source$elev)/100)  # elevation threshold from Bouchard et al. 20xx
    
    ## ID (or index or indicator ;-)) of the focus cells of the epicenters
    sbw.new.sprd = epicenter
    
    ## Find between 20 to 40 SAB or EPN neighs for each epicenter and add to sbw.new.sprd.
    ## Each epicenter core has a different size indicated by the kmin and kmax parameters!
    for(i in 1:length(epicenter)){
      neighs = nn2(land[,c("x", "y")], land[land$cell.id==epicenter[i], c("x", "y")],
                   k=rdunif(1, params$kmin.bubble, params$kmax.bubble) , searchtype='priority')  
      nn.indx = neighs[[1]]
      
      cells_filtered <- land$cell.id[nn.indx][land$spp[nn.indx] %in% c("SAB", "EPN")] #QC: epicenters only in SAB and EPN
      num_cells <- rdunif(1, 20, 40)                                                  #QC: only between 20 and 40. Accroding HistData should be between 8 and 20.
      nn.cellid <- if (length(cells_filtered) < num_cells) {cells_filtered} else {sample(cells_filtered, num_cells)} #QC: ensuring to have enough cells
      
      sbw.new.sprd = unique(c(sbw.new.sprd, nn.cellid))
    }

 
    ## Ask Mathieu:
    sbw.bubble.sprd.host = land$cell.id[land$cell.id %in% sbw.new.sprd & land$ny.def0>=5 & land$tssbw>=30 &
                                    land$spp %in% c("SAB", "EPN") & land$temp>0.5 & land$temp<2.8]
    sbw.bubble.sprd.nonhost = land$cell.id[land$cell.id %in% sbw.new.sprd & land$ny.def0>=5 & land$tssbw>=30 &
                                     !(land$spp %in% c("SAB", "EPN", "NonFor")) & land$temp>0.5 & land$temp<2.8]
    
    # and finally assign intensity to all of them (rewrite intensity just assigned to epicenter cores)
    # before it was prob=c(0.2,0.4,0.3,0.1) for intensities from 0 to 3, where 0 means non actual defoliation
    # and this prob applied only to host cells
    land$curr.intens.def[land$cell.id %in% sbw.bubble.sprd.host] = 
      sample(0:3, size=length(sbw.bubble.sprd.host), replace=T, prob=c(0.2,0.4,0.3,0.1))
    land$curr.intens.def[land$cell.id %in% sbw.bubble.sprd.nonhost] = 
      sample(0:1, size=length(sbw.bubble.sprd.nonhost), replace=T, prob=c(0.2,0.8)) #prob=c(1,0))
    # Number of cells by intensity per tree species (to be checked a bit) 08/04/2024
    # table(land$curr.intens.def[land$cell.id%in%sbw.new.sprd], land$spp[land$cell.id%in%sbw.new.sprd])
  }

  if(outbreak>0){

    cat("epidemic phase ", "\n")
    ## Spatial spreading of the current outbreak to cells not yet defoliated, that is, cells with ny.def0>=5 & tssbw>=30
    ## The function 'spread.tonew' returns cell.ids
    ## The radius has to be variable to allow spreading further away or limit outbreak
    radius = rdunif(1, params$radius.outbreak.mid-params$radius.outbreak.range, params$radius.outbreak.mid+params$radius.outbreak.range) 
    radius = pmax(radius, 1)
    sbw.new.sprd.df = sbw.spread.from.source(land, nc=ncol(mask), params, wind_dir=params$wind_dir, radius=radius, cell.size=cell.size)
    sbw.new.sprd.df$effective_spread = runif(nrow(sbw.new.sprd.df), 0 , 1) <= sbw.new.sprd.df$spread_potential_add
    #get sbw.new.sprd. QC: differentiate sbw.new.sprd from  sbw.new.sprd.df, otherwise step 2, 3 and 5 will not work correctly.
    sbw.new.sprd=sbw.new.sprd.df$target[sbw.new.sprd.df$effective_spread==TRUE]
    ## Add spread potential (final_w_multi) to 'land' df
    land$spread.weight.multi[land$cell.id %in% sbw.new.sprd.df$target] = sbw.new.sprd.df$final_w_multi
    
    ## Only if some new cells are defoliated, assign level of defoliation
    if(length(sbw.new.sprd)>0){
      ## Level of defoliation of the cells recently integrated in the outbreak (the sbw.new.spread cells)
      ## It can be 0 (no-defoliation), 1, 2 or 3!
      ## Probabilities are different according to the spp
      ## Non forested cells are excluded
      sbw.new.sprd.host = land$cell.id[land$cell.id %in% sbw.new.sprd & land$spp %in% c("SAB", "EPN")]
      sbw.new.sprd.nonhost = land$cell.id[land$cell.id %in% sbw.new.sprd & !(land$spp %in% c("SAB", "EPN", "NonFor"))]
      
      land$curr.intens.def[land$cell.id %in% sbw.new.sprd.host] = 
        sample(1:3, size=length(sbw.new.sprd.host), replace=T, prob=c(0.5,0.3,0.2))
      land$curr.intens.def[land$cell.id %in% sbw.new.sprd.nonhost] = 
        sample(0:2, size=length(sbw.new.sprd.nonhost), replace=T, prob=c(0.2,0.7,0.1))
    }
  }

  if(collapse>0){ 
    cat("collapse phase ", "\n")
    ## if collapse, reduce number of new cells, only spontaneously
    potential = filter(land, ny.def0>=5, tssbw>=30, spp %in% c("SAB", "EPN"),temp>0.5, temp<2.8)
    sbw.new.sprd = sample(potential$cell.id, size=max(2,round(rlnorm(1, log(8), 0.5))), replace=F,  #QC:between 0:15 cells and some are 50cells. Calibrated with HistData
                          prob=potential$ny.def0*(1200-potential$elev)/100)
    
    #Add intensity
    land$curr.intens.def[land$cell.id %in% sbw.new.sprd] = 
      sample(0:2, size=length(sbw.new.sprd), replace=T, prob=c(0.5,0.3,0.2))
  }

  if(calm>0){    #calm>0 | collapse>0
    cat("calm phase ", "\n")

    ## not add new locations to the current outbreak if it is collapsing
    ## when calm, few spontaneously cells will have to be here and there defoliated
    potential = filter(land, ny.def0>=5, tssbw>=30, spp %in% c("SAB", "EPN"),temp>0.5, temp<2.8)
    sbw.new.sprd = sample(potential$cell.id, size=max(2,round(rlnorm(1, log(8), 0.5))), replace=F,  #QC:between 0:15 cells and some are 50cells. Calibrated with HistData
                          prob=potential$ny.def0*(1200-potential$elev)/100)
    #Add intensity
    land$curr.intens.def[land$cell.id %in% sbw.new.sprd] = 
      sample(0:2, size=length(sbw.new.sprd), replace=T, prob=c(0.5,0.3,0.2))
  }
  
  
  ## 2. Update SBW tracking variables for newly defoliated cells
  land$ny.def[land$cell.id %in% sbw.new.sprd] = 
    ifelse(land$curr.intens.def[land$cell.id %in% sbw.new.sprd]==0, 0, 
           land$ny.def[land$cell.id %in% sbw.new.sprd]+1)
  land$ny.def0[land$cell.id %in% sbw.new.sprd] = 
    ifelse(land$curr.intens.def[land$cell.id %in% sbw.new.sprd]>0, 0, 
           land$ny.def0[land$cell.id %in% sbw.new.sprd]+1)
  land$cum.intens.def[land$cell.id %in% sbw.new.sprd] = 
    land$curr.intens.def[land$cell.id %in% sbw.new.sprd]
  
  
  ## 3. Update level of defoliation intensity of already defoliated cells, that is, 
  ## cells with ny.def0<5 & tssbw>=30, and ny.def<=18 --> otherwise, the defoliation is forced to stop.

  sbw.curr.sprd = land$cell.id[land$ny.def0<5 & land$tssbw>=30 & land$ny.def<=18 &
                                 land$cell.id %notin% sbw.new.sprd]
  # sbw.curr.sprd = land$cell.id[land$cum.intens.def>0 & land$cell.id %notin% sbw.new.sprd] #QC test
  curr.outbreak = filter(land, cell.id %in% sbw.curr.sprd)
  ## Level of defoliation of these cells. It can be 0 (no-defol), 1, 2 or 3!
  land$curr.intens.def[land$cell.id %in% sbw.curr.sprd] = 
    intens.def.curr(filter(land, cell.id %in% sbw.curr.sprd), params, tbls, preoutbreak, outbreak, collapse, calm)
  ## Update SBW tracking variables
  land$ny.def[land$cell.id %in% sbw.curr.sprd] = land$ny.def[land$cell.id %in% sbw.curr.sprd]+
    ifelse(land$curr.intens.def[land$cell.id %in% sbw.curr.sprd]==0, 0, 1)
  land$ny.def0[land$cell.id %in% sbw.curr.sprd] = 
    ifelse(land$curr.intens.def[land$cell.id %in% sbw.curr.sprd]>0, 0, 
           land$ny.def0[land$cell.id %in% sbw.curr.sprd]+1)
  land$cum.intens.def[land$cell.id %in% sbw.curr.sprd] = 
    land$cum.intens.def[land$cell.id %in% sbw.curr.sprd]+land$curr.intens.def[land$cell.id %in% sbw.curr.sprd]  
  
  
  ## 4. Stop defoilating cells that have been defoliated for more than 18 years
  ## Warning, not consider cells that have just been defoliated
  ## And avoid recurrence of defoliation in this oubreak by making tssbw==0
  sbw.stop = land$cell.id[land$ny.def0==0 & land$ny.def>18 & land$cell.id %notin% sbw.curr.sprd]
  land$curr.intens.def[land$cell.id %in% sbw.stop] = 0
  land$ny.def0[land$cell.id %in% sbw.stop] = 1
  land$tssbw[land$cell.id %in% sbw.stop] = 0
  if(calm==8){land$curr.intens.def[land$cell.id %notin% sbw.new.sprd] = 0} #QC test to solve problems
  
  
  ## 5. Increase number of year of non-defoliation for the remaing cells
  land$ny.def0[land$cell.id %notin% c(sbw.new.sprd, sbw.curr.sprd, sbw.stop)] = 
    land$ny.def0[land$cell.id %notin% c(sbw.new.sprd, sbw.curr.sprd, sbw.stop)]+1
  
  
  ## 6. For those cells that have just accumulated 5 years of non-defoliation, 
  ## reset the counter of number of years of defoliation and 
  ## cumulative intensity of defoliaion
  land$ny.def[land$ny.def0==5] = 0
  land$cum.intens.def[land$ny.def0==5] = 0
  
  
  ## 7. Kill cells according to number of years of defoliation and species composition
  ## with probability as cumulative*current intensity of defoliation
  if(outbreak>0 | collapse>0){
    kill.cells = forest.mortality(land)
  }
  else{
    kill.cells = integer()
  }
  
  
  ## 8. Set outbreak parameters
  if(preoutbreak>0){ 
    preoutbreak = preoutbreak-1
    phase = "preoutbreak"
    if(preoutbreak==0)
      duration.last.outbreak = outbreak = rdunif(1,0,2)+params$duration.last.outbreak  # Gray2008 
  }
  else if(outbreak>0){ 
    outbreak = outbreak-1
    phase = "outbreak"
    if(outbreak==0)
      collapse = rdunif(1,4,6)
  }
  else if(collapse>0){  
    phase = "collapse"
    collapse = collapse-1
    if(collapse==0) 
      calm = round(rnorm(1, 15, 2))
  }
  else if(calm>0){
    phase = "calm"
    calm = calm-1    
    if(calm==0)
      preoutbreak = rdunif(1,3,4)
  }
  
  res = list(kill.cells=kill.cells, land.sbw=land, preoutbreak=preoutbreak, outbreak=outbreak, collapse=collapse, calm=calm)
  return(res) 
} 

