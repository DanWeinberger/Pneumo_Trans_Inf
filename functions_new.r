
## logLL for carriage
calc_Multinomial_log_LL = function(observed_count_of_carriage, modeled_prevalence_of_carriers){
  LL = sum(observed_count_of_carriage * log(modeled_prevalence_of_carriers)) 
  if(any(modeled_prevalence_of_carriers<=0)) LL=-1000001
  return(LL)
}

### model

SIS_ode <- function(times, state, parms) {
  # parameters
  betaVT <- parms$betaVT
  betaNVT <- parms$betaNVT
  clearVT = parms$clearVT
  clearNVT = parms$clearNVT 
  comp =parms$comp
  no.agegps=parms$no.agegps
  age.out=parms$ageout
  age.in= parms$agein
  pop=df$Population
  
  # states
  S <- state[1:no.agegps]
  VT <- state[(1:no.agegps)+no.agegps]
  NVT <- state[(1:no.agegps)+2*no.agegps]
  B <- state[(1:no.agegps)+3*no.agegps]
  N <- S + VT + NVT + B
  
  FOI.VT = (betaVT %*% ((VT+B)/N*pop)) %>% as.vector()
  FOI.NVT = (betaNVT %*% ((NVT+B)/N*pop)) %>% as.vector() 
 
  dS <- -FOI.VT*S - FOI.NVT*S +clearVT*VT +clearNVT*NVT -age.out*S + age.in*c(sum(N),S[-no.agegps])
  dVT <- FOI.VT*S -comp*FOI.NVT*VT  -clearVT*VT +clearNVT*B -age.out*VT + age.in*c(0,VT[-no.agegps])
  dNVT <- FOI.NVT*S -comp*FOI.VT*NVT -clearNVT*NVT +clearVT*B - age.out*NVT + age.in*c(0,NVT[-no.agegps])    
  dB <- comp*FOI.VT*NVT +comp*FOI.NVT*VT -clearVT*B -clearNVT*B - age.out*B + age.in*c(0,B[-no.agegps])
  
  return(list(cbind(dS, dVT, dNVT, dB)))
}

sample_mixing <- function(no.samples,df_mixing_by_id_curr,df_mixing_by_id){
  L=dim(df_mixing_by_id_curr)[1]
  df_mixing_by_id_prop=subset(df_mixing_by_id_curr, p_agegp=="agegp1") #intitialise for while loop
  
  # get proposaed new set of participants through resampling
  if (no.samples>0){
    while(length(unique(df_mixing_by_id_prop$p_agegp)) != length(unique(df_mixing_by_id_curr$p_agegp))){ #prevent from havin no participants of one age group
      smpl=sample(1:L,no.samples,replace=F)
      df_mixing_by_id_prop=df_mixing_by_id_curr[-smpl,]
      df_mixing_by_id_prop = rbind (df_mixing_by_id[sample(1:L,no.samples,replace=T),] , df_mixing_by_id_prop)
    }
  }else df_mixing_by_id_prop = df_mixing_by_id_curr
  return(df_mixing_by_id_prop)
}
 
get_mixing_mat <- function(df_mixing_by_id_prop,df){
  # crerate proposed contact matrix
  pop=df_mixing_by_id_prop %>% group_by(p_agegp) %>% summarise(count=length(local_id))
  Contact_matrix_raw = df_mixing_by_id_prop %>% group_by(p_agegp) %>% 
        summarise_at(paste0("agegp",1:length(unique(df_mixing_by_id_prop$p_agegp))), funs(sum(.,na.rm=T)))
  Contact_matrix_raw = Contact_matrix_raw[,-1] %>% as.matrix() %>% t()
  M = t(t(Contact_matrix_raw)/pop$count) # col=participants; row=contacts
  names(dimnames(M))=list("contact","participant")
  
  M=t(M)*df$Population; M=(M+t(M))/2; M=M/df$Population; Mixing_mat=t(M)/df$Population
  rownames(Mixing_mat)=df$Age_groups
  colnames(Mixing_mat)=rownames(Mixing_mat)
  return(Mixing_mat) 
}

