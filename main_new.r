
################################
# code for Estimating the contribution of different population age strata 
# to vaccine serotype pneumococcal transmission in the pre vaccine era, a modelling study
# written by Stefan Flasche 
################################

#### setup
require("MASS")
require("deSolve")
require("binom")
require("coda")
require("ggplot2")
require("dplyr")
require("tidyr")
require("rootSolve")
require("Rcpp")
require("RcppArmadillo")

setwd("C:/Users/Stefan/Documents/github/Pneumo_Trans_Inf_private")
source("functions_new.r")
sourceCpp("functions_fast.cpp")

# user input
setting="EW"
MCMC.length=1000000
Thin=MCMC.length/10000 #only save every "thin" posterior
Save_out=T

# sensitivity analyses
Fit_type="Susceptibility" #"Susceptibility" ,"Transmissibility"
Dual_carr="VT" #VT , NVT
NTwaningSens = F #alternative waning assumption for Nha trang?
lowComp = F #alternative assumption on competition parameter. 
nonPhysicalContacts = F #also include nonphysical contacts (no such recoreded in Kilifi)

# adjust accordingly
if(Fit_type=="Susceptibility") {Fit_type_num = 1}else if(Fit_type=="Transmissibility"){Fit_type_num = 2}
sim.load=paste0("sim_",setting,"_",substr(Fit_type,1,3))
sim=paste0("Sim_",setting,"_",substr(Fit_type,1,3))
if(Dual_carr=="VT"){
  dc = 1
}else{
  dc=0
  sim=paste0(sim,"_NVT")}
if(lowComp) sim=paste0(sim,"_lowComp") 
if(nonPhysicalContacts) sim=paste0(sim,"_allContacts") 

# set up data
df_Kilifi=data.frame(Setting="Kilifi",
              Age_groups=c("<1y","1-5y","6-14y","15-20y","21-49y","50+"),
              Age_group_upper = c(1,6,15,21,50,62),
              Population=c(9617,45170,68547,33289,72143,24214),
              VT.prev=c(0.41,0.338,.146,.142,.072,.038),
              NVT.prev=c(.46,.44,.39,.25,.22,.20),
              N.carr=c(39,127,82,56,97,104),
              VT.clear=c(.062,.143,.343,.343,.343,.343)/7, 
              NVT.clear=c(.086,.159,.349,.349,.349,.349)/7)
df_EW=data.frame(Setting="EW",
               Age_groups=c("<1y","1-5y","6-14y","15-20y","21-49y","50+"),
               Age_group_upper = c(1,6,15,21,50,82),
               Population=c(629200,3009200,5999000,4132900,21368700,17905900), 
               VT.prev=c(0.431,0.254,.077,.0,.0039,.0), 
               NVT.prev=c(.118,.145,.096,.125,.022,.0),
               N.carr=c(51,138,52,8,232,2),
               VT.clear=1/c(72,34.8,18,17.5,17,17), 
               NVT.clear=1/c(72,34.8,18,17.5,17,17))
df_NT=data.frame(Setting="NT",
                 Age_groups=c("<1y","1-5y","6-17y","18-49y","50+"),
                 Age_group_upper = c(1,6,18,50,76),
                 Population=c(2094,15239,40324,115538,37544),
                 VT.prev=c(0.266,.245,.018,.007,.008), 
                 NVT.prev=c(.074,.202,.092,.012,.013),
                 N.carr=c(41,63,55,262,98),
                 VT.clear=1/c(72,34.8,18,17,17), #from E&W - sensitivity analyses: use Kilifi
                 NVT.clear=1/c(72,34.8,18,17,17))
df_Fin=data.frame(Setting="Fin",
                  Age_groups=c("<1y","1-5y","6-17y","18+"),
                  Age_group_upper = c(1,6,18,81),
                  Population=c(56748,284828,764192,4130843),
                  VT.prev=c(203/1207,422/1235,30/258,25/789)*(1-.319), 
                  NVT.prev=c(203/1207,422/1235,30/258,25/789)*(.319), 
                  N.carr=(c(1200,1235,258,789)/3) %>% round(0), 
                  VT.clear=c(.67,.67,1,1)/30,
                  NVT.clear=c(.67,.67,1,1)/30)

if(NTwaningSens){
  sim = "Sim_NT_Sus_SensWaning"
  df_NT$VT.clear = c(.062,.143,.343,.343,.343)/7
  df_NT$NVT.clear= c(.086,.159,.349,.349,.349)/7
}

#### set parameters
if(setting=="Kilifi"){df=df_Kilifi
}else if(setting=="EW"){df=df_EW;
}else if(setting=="NT"){df=df_NT
}else if(setting=="Fin"){df=df_Fin}
competition=0.1 
if(lowComp) competition=0.3
no.agegps=dim(df)[1]
agegp.l=df$Age_group_upper - c(0,df$Age_group_upper[-no.agegps])
state_prop=c(agegp.l*100-2,rep(1,no.agegps),rep(1,no.agegps),rep(0,no.agegps))

#### Mixing
df_mixing_by_id = read.csv(paste0("Data/ContactsbyID_",df$Setting[1],".csv"),header=T)
if (nonPhysicalContacts) df_mixing_by_id = read.csv(paste0("Data/ContactsbyID_",df$Setting[1],"_all.csv"),header=T)
no.samples=no.samples = 5 #1% of samples resampled every step

df_mixing_by_id_curr = df_mixing_by_id
Mix_mat = get_mixing_mat(df_mixing_by_id_curr,df)

#number of infant contacts
df_mixing_by_id_curr %>% replace(is.na(.), 0) %>%
  filter(p_agegp=="agegp1") %>% 
  dplyr::select(starts_with("age")) %>%
  summarise_all(mean, na.rm = T) %>%
  sum

#### MCMC
acceptance=0
LL_curr=-10000000
p_curr=rep(0.001,no.agegps*2); names(p_curr)=paste(rep(c("VT","NVT"),each=no.agegps),rep(df$Age_groups,2)); d=length(p_curr)
covariance_matrix0=diag(p_curr)*0.01; covariance_matrix=covariance_matrix0
scaling_max=10; scaling=1

#overwrite parameter initialisation and re-create most recent mixing matrix
VT.prev.model.curr=NA; NVT.prev.model.curr=NA
U1_FOI_prop_curr=NA; all_FOI_prop_curr=NA
load(file=paste0("Sims\\",sim.load,"\\R_out.Rdata")); covariance_matrix0=covariance_matrix
if(NTwaningSens) df = df_NT
if(!nonPhysicalContacts){
  mix_id_counts = mixingIDsav[rowSums(mixingIDsav)>0,df_mixing_by_id$local_id %>% as.character()] %>% tail(1) %>% c()
  df_mixing_by_id_curr = df_mixing_by_id[rep(1:dim(df_mixing_by_id)[1], mix_id_counts),]
  Mix_mat = get_mixing_mat(df_mixing_by_id_curr,df)
}

#create variables to store results
psav=matrix(NA,nrow=MCMC.length/Thin,ncol=2*no.agegps); rownames(psav)=paste("MCMC",1:(MCMC.length/Thin)); colnames(psav)=paste(rep(c("VT","NVT"),each=no.agegps),df$Age_groups)
prevsav=array(NA,dim=c(MCMC.length/Thin,no.agegps,2),dimnames=list(MCMC.sample=1:(MCMC.length/Thin), Age_groups=df$Age_groups, type=c("VT.prev","NVT.prev")))
FOIcontributionsav=array(NA, dim=c(MCMC.length/Thin,no.agegps,2),dimnames=list(MCMC.sample=1:(MCMC.length/Thin), Age_groups=df$Age_groups,FOI=c("U1","All")))
mixingIDsav=matrix(0,nrow=MCMC.length/Thin,ncol=dim(df_mixing_by_id)[1]); colnames(mixingIDsav) = df_mixing_by_id$local_id
maxeigvsav=rep(NA,MCMC.length/Thin)
InfectionMatrix = array(NA, dim=c(MCMC.length/Thin,no.agegps,no.agegps),dimnames=list(MCMC.sample=1:(MCMC.length/Thin), infectee=df$Age_groups,Infector=df$Age_groups))

LL_curr=-10000000

for (i in 1:MCMC.length){
  p_prop = 0.05 * (mvrnorm(1,(p_curr),scaling*(2.38^2)/d*covariance_matrix0)) + (1-0.05) * (mvrnorm(1,(p_curr),scaling*(2.38^2)/d*covariance_matrix))
  if(setting %in% c("Kilifi","EW")) {p_prop[c("VT 15-20y","VT 50+", "NVT 15-20y","NVT 50+")] = p_prop[c("VT 6-14y","VT 21-49y","NVT 6-14y","NVT 21-49y")]
  }else if(setting == "NT") {p_prop[c("VT 18-49y", "NVT 18-49y")] = p_prop[c("VT 50+","NVT 50+")]
  }else if(setting == "Fin") {p_prop[] = p_prop[]}
  df_mixing_by_id_prop = sample_mixing(no.samples,df_mixing_by_id_curr,df_mixing_by_id)
  Mix_mat = get_mixing_mat(df_mixing_by_id_prop,df)
  
  if(any(p_prop<0) | any(p_prop>1)){
    LL_prop=-100000
  } else{
  
    ageing=1/(agegp.l*365)
    susVT=p_prop[1:no.agegps]
    susNVT=p_prop[(1:no.agegps)+no.agegps]
    param.list=list(betaVT=sweep(Mix_mat,Fit_type_num,susVT,"*"),
                    betaNVT=sweep(Mix_mat,Fit_type_num,susNVT,"*") ,
                    clearVT=df$VT.clear ,
                    clearNVT=df$NVT.clear ,
                    comp=competition ,
                    no.agegps=no.agegps, 
                    ageout=ageing,
                    agein=c(1/sum(1/ageing) , ageing[-no.agegps]),
                    Population=df$Population)
    
    # solve model steady state
    res = runsteady(y=state_prop,
                    fun=SIS_ode_cpp,
                    parms=param.list,
                    times=c(0,1e5))
    
    steady.state=res$y %>% matrix(nrow=no.agegps); rownames(steady.state)=df$Age_groups; colnames(steady.state)=c("S","VT","NVT","B")
    VT.prev.model = (steady.state[,"VT"] + .5*steady.state[,"B"] )/rowSums(steady.state) # assume 50% of dual carriage is VT as observed in Malawi
    NVT.prev.model= (steady.state[,"NVT"] + .5*steady.state[,"B"])/rowSums(steady.state)
    
    #calc LL
    LL_prop=0
    for (age in 1:no.agegps)
      LL_prop = LL_prop + calc_Multinomial_log_LL(c(df$VT.prev[age], df$NVT.prev[age], 1-df$NVT.prev[age]-df$VT.prev[age])*df$N.carr[age],
                                                  c(VT.prev.model[age],NVT.prev.model[age],1-VT.prev.model[age]-NVT.prev.model[age]) )
      
  }
  rand_num=runif(1)
  log_LL_ratio= LL_prop - LL_curr
  
  if(log(rand_num)<(log_LL_ratio)){
    
    # calculate where the U1s get their VT infection from at steady state
    VT.infected = VT.prev.model*df$Population
    U1_FOI_prop=as.matrix(param.list$betaVT)[1,]*VT.infected # the contributions by age to the U1 row of the FOI 
    all_FOI_prop = colSums(as.matrix(param.list$betaVT)*VT.infected) # the contributions by age to the overall FOI 
    # calculate more generally infection flow matric
    FOI_mat = sweep(as.matrix(param.list$betaVT),2,VT.infected,"*")
    Susecptible4VTinf = (steady.state[,"S"]/rowSums(steady.state) + steady.state[,"NVT"]/rowSums(steady.state)*competition)*df$Population
    Infection_mat = sweep(FOI_mat,1,Susecptible4VTinf,"*")
    
    # overwrite stuff
    p_curr=p_prop
    LL_curr=LL_prop
    acceptance=acceptance+1
    state_prop=res$y
    VT.prev.model.curr=   VT.prev.model
    NVT.prev.model.curr=  NVT.prev.model
    U1_FOI_prop_curr= U1_FOI_prop
    all_FOI_prop_curr= all_FOI_prop
    Infection_mat_curr = Infection_mat
    df_mixing_by_id_curr = df_mixing_by_id_prop
  }

  # thininng
  if((i%%Thin)==0){
    it=i/Thin
    
    # save stuff
    psav[it,]=p_curr
    prevsav[it,,"VT.prev"]=   VT.prev.model.curr
    prevsav[it,,"NVT.prev"]=  NVT.prev.model.curr
    FOIcontributionsav[it,,"U1"]=U1_FOI_prop_curr
    FOIcontributionsav[it,,"All"]=all_FOI_prop_curr
    InfectionMatrix[it,,]=Infection_mat_curr
    df_mixing_by_id_curr_idtab = df_mixing_by_id_curr$local_id %>% table
    mixingIDsav[it,names(df_mixing_by_id_curr_idtab)]=df_mixing_by_id_curr_idtab
    maxeigvsav[it]=(Mix_mat %>% eigen)$value %>% abs() %>% max()
  }
  
  #adaptive MCMC
  if(i>Thin*2){
    covariance_matrix=cov(psav[1:it,])
    scaling = exp(log(scaling) + (0.9999)^(i)*(acceptance/(i)-0.234))
    scaling = min(scaling, scaling_max)
  }
  
  if((i %% Thin) == 0) print(round(c(i,acceptance/i, LL_curr, LL_prop, scaling),2))
 
}



if(Save_out) {
  save(maxeigvsav,p_curr,LL_curr,covariance_matrix, mixingIDsav, InfectionMatrix, psav, prevsav, df, FOIcontributionsav, sim, file=paste0("Sims\\",sim,"\\R_out.Rdata"))
  source("plots.r")
  setting %>% print()
  }
