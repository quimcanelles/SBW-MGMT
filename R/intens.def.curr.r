
intens.def.curr = function(curr.outbreak, params, tbls, preoutbreak, outbreak, collapse, calm){
  
  `%notin%` = Negate(`%in%`)
  
  ## Intensity of defoliation for those cells that were not defoliated in the previous year
  no.defol = filter(curr.outbreak, curr.intens.def==0)
  if(nrow(no.defol)>0){
    if(collapse>0 | calm>0)
      {no.defol$new.intens = 0}
    else
      {no.defol.host = no.defol[no.defol$spp%in%c("SAB", "EPN"),]
      no.defol.host$new.intens = ifelse(abs(rnorm(nrow(no.defol.host), 0, 2))<=no.defol.host$cum.intens.def*no.defol.host$ny.def0, 0, 
                                 sample(1:3, prob=c(1,0.25,0.05), replace=T))
    
      no.defol.nonhost = no.defol[no.defol$spp%notin%c("SAB", "EPN"),]
      no.defol.nonhost$new.intens = ifelse(abs(rnorm(nrow(no.defol.nonhost), 0, 2))<=no.defol.nonhost$cum.intens.def*no.defol.nonhost$ny.def0, 0, 
                                    sample(0:2, prob=c(0.5,0.4,0.1), replace=T))
      
      no.defol=rbind(no.defol.host,no.defol.nonhost)}
    # table(no.defol$new.intens)
  }
  
  ## Intensity of defoliation for those cells that were defoliated in the previous year
  active.defol = filter(curr.outbreak, curr.intens.def>0)
  if(collapse>=4){ 
    to.collapse = rdunif(1, 0.05*nrow(active.defol), 0.15*nrow(active.defol))
    to.collapse = sample(active.defol$cell.id, size=to.collapse, replace=F, prob=(1/active.defol$ny.def)^2)
    active.defol$new.intens[active.defol$cell.id %in% to.collapse] =  0
    active.defol$new.intens[active.defol$cell.id %notin% to.collapse] = 
      pmax(active.defol$curr.intens.def[active.defol$cell.id %notin% to.collapse] -
             rdunif(nrow(active.defol)-length(to.collapse), 0, 1) , 0)
  }
  if(collapse>1 & collapse<4){
    to.collapse = rdunif(1, 0.2*nrow(active.defol), 0.6*nrow(active.defol))
    to.collapse = sample(active.defol$cell.id, size=to.collapse, replace=F, prob=(1/active.defol$ny.def)^2)
    active.defol$new.intens[active.defol$cell.id %in% to.collapse] =  0
    active.defol$new.intens[active.defol$cell.id %notin% to.collapse] = 
      pmax(active.defol$curr.intens.def[active.defol$cell.id %notin% to.collapse] -
             rdunif(nrow(active.defol)-length(to.collapse), 0, 1) , 0)
  }
  else if(collapse==1){
    active.defol$new.intens = 0
  }
  else if(outbreak>0 | preoutbreak>0){  
    active.defol$new.intens = ifelse(runif(nrow(active.defol),0,1)<0.75, active.defol$curr.intens.def,
                                      intensity.defoliation(active.defol, params, tbls))
    active.defol$new.intens[active.defol$new.intens==1 & active.defol$ny.def>1 & active.defol$spp%in%c("SAB","EPN")] =
      sample(1:3, size=length(active.defol$new.intens[active.defol$new.intens==1 & active.defol$ny.def>1] & active.defol$spp%in%c("SAB","EPN")),
             replace=T, prob=c(0.7,0.2,0.1))
  }
  else if(calm>1){
    to.collapse = rdunif(1, 0.2*nrow(active.defol), 0.6*nrow(active.defol))
    to.collapse = sample(active.defol$cell.id, size=to.collapse, replace=F)
    active.defol$new.intens[active.defol$cell.id %in% to.collapse] =  0
    active.defol$new.intens[active.defol$cell.id %notin% to.collapse] = 
      rdunif(nrow(active.defol)-length(to.collapse), 1, 3) 
  } 
  else if(calm==1){
    active.defol$new.intens = 0
  }
  # table(active.defol$ny.def, active.defol$new.intens )
  # table(active.defol$curr.intens, active.defol$new.intens)
  if(nrow(no.defol)==0)
    kk = active.defol
  else if(nrow(active.defol)==0)
    kk = no.defol
  else
    kk = rbind(no.defol, active.defol)
  kk = left_join(dplyr::select(curr.outbreak, cell.id), dplyr::select(kk, cell.id, new.intens), by="cell.id")
  # sum(kk$new.intens>0)/length(kk$new.intens)
  # round(table(kk$new.intens)[-1]/sum(kk$new.intens>0),2)
  
  return(kk$new.intens)
  
}

