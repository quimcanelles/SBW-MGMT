#############################################################################################################
#### Script for runninng ap-sbw model #######################################################################
#### 

rm(list=ls())
source("R/ap-sbw.r")

scn="scn0"; #scn="scn2", scn="scn3", scn="scn4", scn="scn5", scn="scn6", scn="scn7", scn="scn8"
nrun=2; 
time.step=1; #only 5 
time.horizon=80; #0:80 & n*5
custom.params=NA; 
rcp='rcp45'#'rcp85',NA
is.sbw=T; 
out.path="outputs/test1"


############################################ RUN ap-sbw() ##################################################

## 1 basic run
r = ap.sbw(scn="scn0", is.sbw=T, is.harvesting=T, custom.params=NA, rcp='rcp45', nrun=1, out.path=NA )


## 2 Save outputs
source("R/default.params.r")  
custom.params = default.params()
custom.params$save.land.df=T; custom.params$outbreak=12
r = ap.sbw(scn="scnDefault", is.sbw=T, is.harvesting=T, custom.params=custom.params, rcp='rcp45', nrun=1, out.path="outputs/test_scnDefault" )


## 3 Change default parameters according to scn
source("R/default.params.r")  
custom.params = default.params()
custom.params$save.land.df=T; custom.params$outbreak=12
custom.params$harv.rate=0.0025; custom.params$ap.rate=0.5; custom.params$epn.rege.rate=0.2: custom.params$pet.rege.rate=0.8
r = ap.sbw(scn="scn2", is.sbw=F, is.harvesting=T, custom.params=custom.params, rcp='rcp45', nrun=3, out.path="outputs/test_scn2" )
  
  
## 4 Change values of an input table, eg. soil.suitability of SAB
#*** incloure custom tables a la funció principal
data(default.tables)
soil.suitability = default.tables[["soil.suitability"]]
soil.suitability[soil.suitability$spp=="SAB",2:6] = 0
default.tables[["soil.suitability"]] = soil.suitability
r = ap.sbw(scn="scn1", is.sbw=F, is.harvesting=T, custom.params=custom.params, rcp='rcp85', nrun=3, out.path=NA )

sbw.outbreak(custom.params = NULL, custom.tables = tbl, rcp = NA, prec.proj = NA, temp.proj = NA,  
                time.horizon = 10, nrun = 1, time.save=5, out.path = NA) 




############################################ RUN ap-sbw() per scn ##################################################
## Run the AP-SBW model for a set of scenarios --------------------------------------------------------------------------
## Read an excel file with a list of parameters values for each testing scenario,
## then, pass these values to the corresponding elements in the "params" list.
rm(list=ls())
library(readxl)
source("R/default.params.r"); source("R/ap-sbw.r")
scenario.params = read_excel("data/params_scenarios.xlsx", sheet = "scenario_params_paper")

scenarios = scenario.params$scenario[2:nrow(scenario.params)]
scn = "scn0" #scn0, scn1, scn2, ...,  scn12
for (scn in (scenarios)){
  custom.params = default.params()  # to have a named list of the parameters 
  for(i in 1:length(custom.params)){
    val <- scenario.params[scenario.params$scenario == scn, names(custom.params)[i]]
    custom.params[[i]] = val[[1]]
  }
  # transform class of the parameters
  custom.params$save.land.df = T
  custom.params$stop.end.phase = ifelse(custom.params$stop.end.phase %in% c("FALSE", "F"), F, T)
  custom.params$enable.succ = ifelse(custom.params$enable.succ %in% c("FALSE", "F"), F, T)
  custom.params$is.harvprem = ifelse(custom.params$is.harvprem %in% c("FALSE", "F"), F, T)
  custom.params$freq.save = 5
  
  # run the model
  out.path=paste0("outputs/test3/", scn)
  res = ap.sbw(scn=scn, is.sbw=T, is.harvesting=T, custom.params=custom.params, rcp='rcp45', nrun=10, out.path)
  
}










############################################ RUN ap-sbw() as Nu ##################################################  
## Run the AP-SBW model for a set of scenarios --------------------------------------------------------------------------
## Read an excel file with a list of parameters values for each testing scenario,
## then, pass these values to the corresponding elements in the "params" list.

source("R/default.params.r"); source("R/ap-sbw.r")

scenario.params = read_excel("data/params_scenarios.xlsx", sheet = "scenario_params")
scenarios = c("scn01", "scn02"); scn = "scn01"
for(scn in scenarios){
  custom.params = default.params()  # to have a named list of the parameters 
  for(i in 1:length(custom.params))
    custom.params[[i]] = scenario.params[scenario.params$scenario==scn, names(custom.params)[i]]
  # transform class of the parameters
  custom.params$save.land.df = ifelse(custom.params$save.land.df %in% c("FALSE", "F"), F, T)
  custom.params$stop.end.phase = ifelse(custom.params$stop.end.phase %in% c("FALSE", "F"), F, T)
  custom.params$enable.succ = ifelse(custom.params$enable.succ %in% c("FALSE", "F"), F, T)
  custom.params$is.harvprem = ifelse(custom.params$is.harvprem %in% c("FALSE", "F"), F, T)
  # run the model
  out.path=paste0("outputs/", scn)
  res = ap.sbw(scn=scn, is.sbw=T, is.harvesting=T, custom.params=custom.params, rcp='rcp45', nrun=1, out.path)
  if(!file.exists(out.path))
    dir.create(file.path(out.path), showWarnings = T) 
  saveRDS(res, paste0(out.path, "/", "ap_sbw.rds"))
}


## Run the AP-SBW model from scratch --------------------------------------------------------------------------
## Load initial conditions and remove the current outbreak from the landscape to start from scratch
load("C:/Users/nuria.aquilue/OneDrive - ctfc.cat/QBCMOD/AP-SBW/AP-SBW/data/landscape.rda")
summary(landscape)
landscape$ny.def = 0
landscape$ny.def0 = 30
landscape$cum.intens.def = 0
landscape$curr.intens.def  = 0
## Run 
r = ap.sbw(scn="scn1", is.sbw=T, is.harvesting=T, custom.params=params, rcp='rcp45', #is.harvprem=F,
           nrun=2, out.path="outputs/test", landscape = landscape)
