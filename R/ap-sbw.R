#' Quebec Landscape Dynamic Model: Artificial planting and SBW
#'
#' Run the landscape dynamic model QLDM that includes the processes of spruce budworm outbreaks, clear-cuts,
#' artificial planting, post-disturbance regeneration, and natural forest succession.
#'
#' @param scn Indicates the scn (setted in default.tables) that defines harvesting and artificial planting rates. 'scn0':'scn8'
#' @param is.sbw A flag to indicate that spruce budworm outbreaks are a simulated process of the model
#' @param is.harvesting A flag to indicate that harvesting by clear and partial cuts is a simulated process of the model
#' @param custom.params List with the model paramaters and default and/or user-defined values 
#' @param rcp Climate projection, either \code{NA} (default), 'rcp45' or 'rcp85' 
#' @param nrun Number of replicates to run the model
#' @param out.path String with the directory path to save the \code{landscape} data frame at each time step indicated in \code{out.seq}
#'  
#' @example 
#' params = default.params()
#' res = ap.sbw(scn="test", is.sbw=T, is.harvesting=T, custom.params=params, 
#' rcp='rcp45', nrun=3, out.path=NA)

#' @case.example
#' scn="scn0"; is.sbw=T; is.harvesting=T; custom.params=default.params(); rcp='rcp45'; nrun=3; out.path="outputs/test"; 
#' custom.params$preoutbreak = 2; custom.params$harv.rate=0.0025; custom.params$ap.rate=0.5; custom.params$epn.rege.rate=0.8; custom.params$pet.rege.rate=0.2

ap.sbw = function(scn, is.sbw = FALSE, is.harvesting = FALSE, custom.params = NA, 
                  rcp = NA, nrun = 1, out.path = NA, landscape = NULL){
  
  
  ###########################################################################-
  #################        0. LOAD INITIAL ARGUMENTS        #################
  
  #Library
  library(stringr)
  library(RANN)
  library(sp)
  library(raster)
  library(scales)
  library(tidyverse) 
  options(dplyr.summarise.inform=F)
  
  #Functions
  source("R/artificial.planting.r")
  source("R/buffer.mig.r")
  source("R/forest.transition.r")
  source("R/data.r")
  source("R/default.params.r")
  source("R/forest.mortality.r")
  source("R/harvest.area.r")
  source("R/intens.def.curr.r")
  source("R/intensity.defoliation.r")
  source("R/sbw.outbreak.r")
  source("R/sbw.spread.from.source.r")
  source("R/suitability.r")
  ## Function to select items not in a vector
  `%notin%` = Negate(`%in%`)
  
  #Data
  load("data/mask.rda")
  load("data/prec_rcp45_ouranos.rda")
  load("data/temp_rcp45_ouranos.rda")
  load("data/default.tables.rda")  # as v2 but without scn.df table
  if(is.null(landscape)){
    load("data/landscape.rda")
  }
  
  #####################################################################-
  #################        1. DATA PREPARATION        #################
  
  cat("A. Data preparation ...\n")
  
  ### 1.1. PARAMETERS ###
  ## Get the list of default parameters and update user-initialized parameters
  params = default.params()
  if(!any(is.na(custom.params))){
    # Check class of custom.params
    if((!inherits(custom.params, "list"))) {
      stop("'custom.params' must be a named list")
    }
    ## Check that the names of the customized parameters are correct
    if(!all(names(custom.params) %in% names(params)))
      stop("Wrong custom parameters names")
    params = custom.params
  }
  
  ### 1.2. TIME STEPS ###
  ## Build the discrete time sequence according to time.step
  time.seq = seq(1, params$time.horizon, params$time.step) 
  
  ## Set out.seq to save function outputs
  if(params$save.land.df){
      out.seq = seq(0, params$time.horizon, params$freq.save)
      out.seq = out.seq[out.seq!=0]
      if(!all(out.seq %in% time.seq)){warning("Not all time steps in the output sequence provided are simulated.", call.=F)}
      if(is.na(out.path)) stop("Directory path to save outputs not provided") } 
  if(!(params$save.land.df)){
    out.seq=0
  }

  ### 1.3. CELL RESOLUTION ###
  cell.res = res(mask)[1] ## cell resolution in meters
  km2.pixel = raster::res(mask)[1] * raster::res(mask)[2] / 10^6
  
  ### 1.4. SET CLIMATIC PROJECTIONS ###
  ## Load precipitation and temperature projections provided with the package according to the climatic scenario.
  prec.chg = F; temp.chg = F
  if(!is.na(rcp)){
    prec.proj = get(paste0("prec_", rcp)); prec.chg = T
    temp.proj = get(paste0("temp_", rcp)); temp.chg = T
  }
  
  ### 1.5. INITIALIZE TRACKING DATA FRAMES ###
  breaks = c(0,20,40,60,80,100,999)
  tags = c("C10","C30", "C50", "C70", "C90", "OLD")
  track.spp.age.class = data.frame(run=NA, year=NA, mgmt.unit=NA, spp=NA, age.class=NA, area=NA)
  track.suit.class = data.frame(run=NA, year=NA, bioclim.domain=NA, potential.spp=NA, poor=NA, med=NA, good=NA)
  track.sbw.defol.intens = data.frame(run=NA, year=NA, phase=NA, cell.id=NA, spp=NA, curr.intens.def=NA, bioclim.domain=NA)
  track.sbw.defol.intens.sm = data.frame(run=NA, year=NA, phase=NA, curr.intens.def=NA, ncell=NA, pct=NA)
  track.sbw.kill = data.frame(run=NA, year=NA, phase=NA, cell.id=NA, spp=NA, ny.def=NA, curr.intens.def=NA, bioclim.domain=NA)
  track.sbw.kill.sm = data.frame(run=NA, year=NA, spp=NA, ny.def=NA, curr.intens.def=NA, area=NA)
  track.sbw.kill.sm2 = data.frame(run=NA, year=NA, area=NA)
  track.cut = data.frame(run=NA, year=NA, mgmt.unit=NA, spp=NA, age=NA)
  track.cut.sm = data.frame(run=NA, year=NA, a.age=NA, ncells=NA)
  track.ap = data.frame(run=NA, year=NA, initial.spp=NA, ap.spp=NA, mgmt.unit=NA)
  track.ap.sm = data.frame(run=NA, year=NA, ncells=NA)
  track.change = data.frame(cell.id=NA, run=NA, year=NA, spp=NA, trans=NA, bioclim.domain=NA)
  track.change.sm = data.frame(run=NA, year=NA, trans=NA, n_obs=NA)
  track.spp = data.frame(run=NA, year=NA, BOJ=NA, EPN=NA, ERS=NA, NonFor=NA, OTH.FEU.N=NA, OTH.FEU.S=NA, OTH.RES.N=NA, OTH.RES.S=NA, PET=NA, SAB=NA)


  ######################################################################-
  #################        2. START SIMULATIONS        #################
  
  cat("\n") 
  cat(paste("B. Simulations ...\n"))
  
  # irun=1
  for(irun in 1:nrun){

    ## Main landscape data frame 
    tbls = default.tables
    land = landscape
    land$year.change=NA 
    
    ## Mark the initial sbw phase
    duration.last.outbreak = params$outbreak + params$current.duration
    preoutbreak = params$preoutbreak 
    outbreak = params$outbreak    
    collapse = params$collapse    
    calm = params$calm
    current.phase = ifelse(params$preoutbreak>0, "preoutbreak",
                   ifelse(params$outbreak>0, "outbreak", 
                          ifelse(params$calm>0, "calm",
                                 ifelse(params$collapse>0, "collapse", 
                                        stop("At least one phase has to have an initial duration")))))

    ## Record initial distributions:
    aux =  group_by(land, curr.intens.def) %>% summarize(ncell=length(cell.id)) %>% mutate(pct=NA)
    aux$pct[aux$curr.intens.def==0] = NA
    aux$pct[aux$curr.intens.def!=0] = aux$ncell[aux$curr.intens.def!=0] / sum(na.omit(land$curr.intens.def>0))
    track.sbw.defol.intens.sm = rbind(track.sbw.defol.intens.sm, data.frame(run=irun, year=params$year.ini, phase=current.phase, aux))
    if(any(land$cum.intens.def>0)){
      track.sbw.defol.intens = rbind(track.sbw.defol.intens, data.frame(run=irun, year=params$year.ini, phase=current.phase,
                                                                        filter(land, curr.intens.def>0) %>% dplyr::select(cell.id, spp, curr.intens.def, bioclim.domain)))      
    }
    track.spp = rbind(track.spp, cbind(data.frame(run=irun, year=params$year.ini),
                                       pivot_wider(as.data.frame.table(table(land$spp)), 
                                                   names_from = Var1, values_from = Freq)))

    ## Simulation per time step ##
    # t=1
    for(t in time.seq){
      
      ## Print replicate and time step
      cat("\n") 
      cat(paste0("Replicate ", irun, "/", nrun,". Time step ", t, "/", params$time.horizon, ": ", 
                 t+params$year.ini, ". Scn: ",scn,"\n"))
      
      ### 2.1. SET LAND TRANSITION AND LAND NEWSPP 
      land$transition=NA
      land$new.spp = NA
      land$spread.weight.multi = NA
      
      ### 2.2. UPDATE CLIMATIC VARIABLES ###
      ## Update climatic variables at each time step if climate change is activated
      ## Column 1 is cell.index, the following columns are temp (precip) in 2020-2025, 2025-2030, 2030-2035, ... etc.
      ## The last column (temp95) then corresponds to the period 2095-2100
      ## The first time step (t=5) we start with climate 2020-2025
      if(temp.chg & t < params$time.horizon){
        cat("  Update temperature projections\n")
        aux = temp.proj[,c(1,t+74)] 
        names(aux) = c("cell.id", "temp")
        aux[is.na(aux)] = mean(na.omit(aux[,c(2)])) #*To fix the 2 rows with NA causing troubles.
        land = dplyr::select(land, -temp) %>% left_join(aux, by="cell.id")
      }
      if(prec.chg & t < params$time.horizon){
        cat("  Update precipitation projections\n")
        aux = prec.proj[,c(1,t+74)] 
        names(aux) = c("cell.id", "prec")
        aux[is.na(aux)] = mean(na.omit(aux[,c(2)])) #*To fix the 2 rows with NA cousing troubles.
        land = dplyr::select(land, -prec) %>% left_join(aux, by="cell.id")
      }
      
      
      ### 2.3. UPDATE SUITABILITY
      suitab = suitability(land, params, tbls)
      
      
      ##################################### PROCESSES OF CHANGE #####################################
      #### 1. SBW module ####
      if(is.sbw){
        cat ("  B.1. SBW \n")
        sbw.out = integer()
        sbw.out = sbw.outbreak(land, params, tbls, preoutbreak, outbreak, calm, collapse, duration.last.outbreak, mask)
        kill.cells = sbw.out$kill.cells
        land.sbw = sbw.out$land.sbw    
        preoutbreak = sbw.out$preoutbreak 
        outbreak = sbw.out$outbreak    
        collapse = sbw.out$collapse    
        calm = sbw.out$calm
        phase = ifelse(preoutbreak>0, "preoutbreak",  ifelse(outbreak>0, "outbreak",
                       ifelse(calm>0, "calm","collapse")))
        
        ## Tracking sbw.defol.intens
        if(sum(na.omit(land.sbw$curr.intens.def))>0){
          track.sbw.defol.intens = rbind(track.sbw.defol.intens, data.frame(run=irun, year=t+params$year.ini, phase=phase,
                                                                            filter(land.sbw, curr.intens.def>0) %>% dplyr::select(cell.id, spp, curr.intens.def,bioclim.domain)))
          ncell.def = c(NA, rep(sum(na.omit(land.sbw$curr.intens.def>0)),length(unique(land.sbw$curr.intens.def))-1))
          track.sbw.defol.intens.sm = rbind(track.sbw.defol.intens.sm,
                     data.frame(run=irun, year=t+params$year.ini, phase=phase, group_by(land.sbw, curr.intens.def) %>%
                     summarize(ncell=length(cell.id)) %>% mutate(pct=ncell/ncell.def)))
        }
        else{
          track.sbw.defol.intens = track.sbw.defol.intens
          track.sbw.defol.intens.sm = rbind(track.sbw.defol.intens.sm,
                     data.frame(run=irun, year=t+params$year.ini, phase=phase, curr.intens.def=0, ncell=nrow(land.sbw), pct=NA))
        }
        
        ## Tracking sbw.defol.kill
        if(length(kill.cells)>0){
          track.sbw.kill = rbind(track.sbw.kill, data.frame(run=irun, year=t+params$year.ini, phase=phase, 
                                                            filter(land, cell.id%in%kill.cells) %>% dplyr::select(cell.id, spp, ny.def, curr.intens.def, bioclim.domain)))
          aux = filter(land, cell.id %in% kill.cells) %>% group_by(spp, ny.def, curr.intens.def) %>% 
            summarise(area=length(spp)*km2.pixel) 
          names(aux) = c("spp", "ny.def", "curr.intens.def", "area")
          track.sbw.kill.sm = rbind(track.sbw.kill.sm, data.frame(run=irun, year=t+params$year.ini, aux))  
          area.sbw = sum(aux$area)
          track.sbw.kill.sm2 = rbind(track.sbw.kill.sm2, data.frame(run=irun, year=t+params$year.ini, area=area.sbw))
        }


        # Modify land data frame for killed cells
        land$age[land$cell.id %in% kill.cells] = 0
        land$tssbw = land$tssbw + params$time.step   
        land$tssbw[land$cell.id %in% kill.cells] = 0

        # Copy the sbw variables to the main data frame
        land$ny.def[land$cell.id %in% land.sbw$cell.id] = land.sbw$ny.def
        land$ny.def[land$cell.id %in% kill.cells] = 0
        land$ny.def0[land$cell.id %in% land.sbw$cell.id] = land.sbw$ny.def0
        land$ny.def0[land$cell.id %in% kill.cells] = 0
        land$cum.intens.def[land$cell.id %in% land.sbw$cell.id] = land.sbw$cum.intens.def
        land$cum.intens.def[land$cell.id %in% kill.cells] = 0
        land$curr.intens.def[land$cell.id %in% land.sbw$cell.id] = land.sbw$curr.intens.def
        land$curr.intens.def[land$cell.id %in% kill.cells] = 0
        ## Add spread potential of the epidimic phase
        land$spread.weight.multi[land$cell.id %in% land.sbw$cell.id] = land.sbw$spread.weight.multi
        
      }#end is.sbw
    
       
      #### 2. Harvest area ####
      if(is.harvesting){
        cat ("  B.2. Harvest area \n")
        harv.out = integer()
        harv.out = harvest.area(land, params, params$is.harvprem) 
      
        ## Tracking
        if(nrow(harv.out)>0){
          track.cut = rbind(track.cut, data.frame(run=irun, year=t+params$year.ini, mgmt.unit=harv.out$mgmt.unit, spp=harv.out$spp, age=harv.out$age))
          track.cut.sm = rbind(track.cut.sm, data.frame(run=irun, year=t+params$year.ini, a.age=trunc(mean(harv.out$age)), ncells=nrow(harv.out)))
        }
        if(nrow(harv.out)==0){
          track.cut = rbind(track.cut, data.frame(run=irun, year=t+params$year.ini, mgmt.unit=NA, spp=NA, age=NA))
          track.cut.sm = rbind(track.cut.sm, data.frame(run=irun, year=t+params$year.ini, a.age=NA, ncells=0))
          track.ap = rbind(track.ap, data.frame(run=irun, year=t+params$year.ini, initial.spp=NA, ap.spp=NA, mgmt.unit=NA))
          track.ap.sm = rbind(track.ap.sm, data.frame(run=irun, year=t+params$year.ini, ncells=NA))
        }
      
        ## Modifying land data frame
        land$age[land$cell.id %in% harv.out$cell.id] = 0
        land$age.class[land$cell.id %in% harv.out$cell.id] = "C10" ##**** Class C10??
        land$new.spp[land$cell.id %in% harv.out$cell.id] = NA
        land$transition[land$cell.id %in% harv.out$cell.id] = "Harv"
        land$year.change[land$cell.id %in% harv.out$cell.id] = params$year.ini+t
        no.ap.cells = filter(land, land$cell.id %in% harv.out$cell.id)
      
        ##### 3. Artificial planting ####
        ap.out = integer()
        if(nrow(harv.out)>0){
          cat ("  B.3. Artificial planting \n")
          ap.out = artificial.planting(land, harv.out$cell.id, suitab, params)
      
          ## Tracking
          if(length(ap.out)==0){
            track.ap = rbind(track.ap, data.frame(run=irun, year=t+params$year.ini, initial.spp=NA, ap.spp=NA, mgmt.unit=NA))
            track.ap.sm = rbind(track.ap.sm, data.frame(run=irun, year=t+params$year.ini, ncells=NA))
          }
          if(length(ap.out)>0){
            ## Tracking
            track.ap = rbind(track.ap, data.frame(run=irun, year=t+params$year.ini, initial.spp=ap.out$spp, ap.spp=ap.out$new.spp, mgmt.unit=ap.out$mgmt.unit))
            track.ap.sm = rbind(track.ap.sm, data.frame(run=irun, year=t+params$year.ini, ncells=nrow(ap.out)))
         
            ## Modifying land data frame
            land$new.spp[land$cell.id %in% ap.out$cell.id] = ap.out$new.spp
            land$transition[land$cell.id %in% ap.out$cell.id] = "Harv_AP"
            land$year.change[land$cell.id %in% ap.out$cell.id] = params$year.ini+t
            ## Clear cutted cells without artificial planting
            no.ap.cells = filter(land, land$cell.id %in% harv.out$cell.id & !land$cell.id %in% ap.out$cell.id)
          }
        } #end of AP: nrow(harv.out)>0
      } #end is.harvesting
      
      
      #### 4. Forest transition ####
      cat ("  B.4. Forest transition \n")
      
      ## After harvesting
      if(exists("no.ap.cells")){ #!!! or nrow(no.ap.cells)>0??
        land$new.spp[land$cell.id %in% no.ap.cells$cell.id] = forest.transition(land, no.ap.cells$cell.id, suitab, params, default.tables, type.trans="C")
        land$transition[land$cell.id %in% no.ap.cells$cell.id] = "Harv_ForTrans"
        land$year.change[land$cell.id %in% no.ap.cells$cell.id] = params$year.ini+t
        land$tscomp[land$cell.id %in% no.ap.cells$cell.id] = 0
      }
        
      ## After SBW outbreak
      if(is.sbw){
        if(exists("kill.cells")){
        land$new.spp[land$cell.id %in% kill.cells] = forest.transition(land, kill.cells, suitab, params, tbls, type.trans="O")
        land$transition[land$cell.id %in% kill.cells] = "SBW_ForTrans"
        land$year.change[land$cell.id %in% kill.cells] = params$year.ini+t
        land$tscomp[land$cell.id %in% kill.cells] = 0}
      }
        
        
      ## Natural succession of tree spp at every 40 years starting at Time Since Chancge of Species Composition = 70
      chg.comp.cells = filter(land, (age-age.matu) %in% seq(40,400,40) & tscomp>=70) %>% dplyr::select(cell.id)
      if(length(unlist(chg.comp.cells))>0){
        land$new.spp[land$cell.id %in% unlist(chg.comp.cells)] = forest.transition(land, unlist(chg.comp.cells), suitab, params,default.tables, type.trans="S")  
        land$transition[land$cell.id %in% unlist(chg.comp.cells)] = "NatTrans"
        land$year.change[land$cell.id %in% unlist(chg.comp.cells)] = params$year.ini+t
      }
      
        
      #### 5. Update land.df ####
      land$spp = as.factor(ifelse(is.na(land$new.spp),as.character(land$spp), land$new.spp))
      land$age = land$age + params$time.step
      land$tscomp = land$tscomp + params$time.step
      
        
      ##################################### TRACKING AND SPATIAL OUTS #####################################
      cat ("  B.5. Traking spatial outputs \n")
      ## Age classes distribution per species and management unit **different than age.class???      
      land$age.class = cut(land$age, breaks=breaks, include.lowest=TRUE, right=TRUE, labels=tags)
      track.spp.age.class = rbind(track.spp.age.class, data.frame(run=irun, year=t+params$year.ini, 
                                                                   group_by(land, mgmt.unit, spp) %>% count(age.class) %>%  
                                                                     mutate(area=n*km2.pixel)) %>% dplyr::select(-n))
      ## Suitability classes distribution per bioclim.domain   
      suitab = suitability(land, params, tbls) 
      aux = left_join(suitab, dplyr::select(land, cell.id, bioclim.domain), by="cell.id") %>%
        group_by(bioclim.domain, potential.spp) %>% summarise(poor=sum(suit.clim==0)*km2.pixel, 
                                                              med=sum(suit.clim==0.5)*km2.pixel, good=sum(suit.clim==1)*km2.pixel) 
      track.suit.class = rbind(track.suit.class, data.frame(run=irun, year=t+params$year.ini, aux))
      rm(suitab)
      
      ## Track change
      track.change = rbind(track.change, na.omit(data.frame(cell.id=land$cell.id, run=irun, year=t+params$year.ini, spp=land$spp, trans=land$transition, bioclim.domain=land$bioclim.domain)))
      aux = na.omit(track.change %>% group_by(run, year, trans) %>% summarise(n_obs = n()))
      track.change.sm = rbind(track.change.sm,aux)
      
      
      ## Track change of species
      track.spp = rbind(track.spp, cbind(data.frame(run=irun, year=t+params$year.ini),
                                         pivot_wider(as.data.frame.table(table(land$spp)), names_from = Var1, values_from = Freq)))
      
      ## If required, save landscape data frame at the time steps required
      if(params$save.land.df & t %in% out.seq){
        if(!file.exists(out.path))
          dir.create(file.path(out.path), showWarnings = T, recursive = T) 
        saveRDS(land, file=paste0(out.path, "/landscape_run", irun, "_step", t, ".rds"))
      }
      
      ### Cleaning global environment 
      rm(sbw.out, land.sbw, kill.cells, chg.comp.cells)
      if(is.harvesting) 
        rm(ap.out, harv.out, no.ap.cells)
      
      ## Stop the current run if phase changes
      if(params$stop.end.phase & current.phase!=phase){
       break()
      }
      current.phase = phase
  
    }#end t
  }#end irun
  
  
  ######################################################################-
  #################        3. FINALIZATION        #################
  cat("\n", "C. Build outputs...\n")
  
  res = list(SppByAgeClass = track.spp.age.class[-1,],
              SuitabilityClass = track.suit.class[-1,],
              LandChange = track.change[-1,],
              LandChange.sm = track.change.sm[-1,],
              ArtificialPlanting = track.ap[-1,],
              ArtificialPlanting.sm = track.ap.sm[-1,],
              Cuts = track.cut[-1,],
              Cuts.sm = track.cut.sm[-1,],
              SBWDefol = track.sbw.defol.intens[-1,],
              SBWDefol.sm = track.sbw.defol.intens.sm[-1,],
              SBWKill = track.sbw.kill[-1,],
              SBWKill.sm = track.sbw.kill.sm[-1,],
              SBWKill.sm2=track.sbw.kill.sm2[-1,],
              SppChange=track.spp[-1,])
  
  
  ## If required, save res
  if(params$save.land.df){
    if(!file.exists(out.path))
      dir.create(file.path(out.path), showWarnings = T) 
    saveRDS(res, file=paste0(out.path, "/ap.sbw_results.rds"))
    fileConn=file(paste0(out.path, "/conditions.txt"))
    writeLines(c(paste0("scn: ",scn),paste0("is.sbw: ",is.sbw), paste0("is.harvesting: ",is.harvesting), 
                 paste0("custom.params: ",custom.params), paste0("rcp: ",rcp),paste0("nrun: ",nrun),paste0("time.horizon: ",params$time.horizon),paste0("freq.save: ",params$freq.save)), fileConn)
    close(fileConn)
  }
  
  cat("\n", "END OF SIMULATION", "\n")
  return(res)
}#end ap.sbw function

