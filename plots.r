
#### setup
#setwd("D:\\Github\\Pneumo_Trans_Inf")

require("ggplot2")
require("coda")
require("dplyr")
require("tidyr")
require("MultinomialCI")
require("binom")
require("GGally")

#setting="NT"
#Fit_type="Susceptibility" #"Susceptibility" ,"Transmissibility"
sim=paste0("sim_",setting,"_",substr(Fit_type,1,3))

load(file=paste0("Sims\\",sim,"\\R_out.Rdata"))
no.agegps=dim(df)[1]

infs_VT_ann = (df$Population* df$VT.prev * df$VT.clear * 365) %>% sum() #annual number of new VT infections

# plot p chain convergence
df_psav = psav %>% cbind(maxeigvsav) %>% as.data.frame() 
df_psav = df_psav[!(df_psav %>% rowSums() %>% is.na() %>% as.vector()),]
df_psav$MCMC = rownames(df_psav)
df_psav = df_psav %>% gather(p,value,-MCMC)
df_psav$MCMC = gsub("MCMC ","",df_psav$MCMC) %>% as.numeric
p1 = ggplot(df_psav, aes(x=MCMC, y=value, group=p)) +
  geom_line() +
  facet_wrap(~p, scales="free_y")
ggsave(paste0("Sims\\",sim,"\\convergence.tiff"), p1, compression="lzw", unit="cm", width=23,height=15)

#plot correlations
df_psav %>% spread(p,value) %>% ggpairs %>%
  ggsave(paste0("Sims\\",sim,"\\correlation.tiff"), . , compression="lzw", unit="cm", width=35,height=25)

#gweke convergence assessment, any Z>2 means indication of lack of converagence
any(abs(((df_psav %>% spread(p,value))[,-1] %>% as.mcmc() %>% geweke.diag)$z)>2)
 
#median parameter estimates       
df_psav %>% group_by(p) %>% dplyr::summarise(median_ps=median(value))

#plot Mixing matrix convergence
df_mix = mixingIDsav[rowSums(mixingIDsav)>0,] %>% apply(2,quantile, probs=c(0.025,.5,.975)) %>% t() %>% as_tibble()
df_mix$N = 1:dim(df_mix)[1]
names(df_mix) = c("low","mid","high","N")
p_mix = ggplot(df_mix, aes(x=N, y=mid, ymin=low, ymax=high)) +
  geom_point() +
  geom_linerange() +
  ylab("number of inclusions of the participant") + xlab("contact survey participant number")
ggsave(paste0("Sims\\",sim,"\\convergence_mixing.tiff"), p_mix, compression="lzw", unit="cm", width=50,height=8)


#plot fit to data
require(plyr)
df_prevsav = prevsav %>% adply(1:3) 
names(df_prevsav)[4]="prevalence"
detach("package:plyr", unload=TRUE)
df_prevsav = df_prevsav %>% group_by(Age_groups, type) %>% dplyr::summarise(prev_med = median(prevalence, na.rm=T),
                                                   prev_lo = quantile(prevalence, probs=0.025, na.rm=T),
                                                   prev_hi = quantile(prevalence, probs=0.975, na.rm=T))
df_prevsav$data="model"
df_dat = df[,c("Age_groups","VT.prev","NVT.prev")] %>% gather(type,prev_med,-Age_groups); df_dat$prev_lo=NA; df_dat$prev_hi=NA; df_dat$data="data"
df_dat[c(1:no.agegps),c("prev_lo","prev_hi")] =  binom.confint(df$VT.prev*df$N.carr,df$N.carr, method="exact")[,c("lower","upper")]
df_dat[(1:no.agegps)+no.agegps,c("prev_lo","prev_hi")] =  binom.confint(df$NVT.prev*df$N.carr,df$N.carr, method="exact")[,c("lower","upper")]
df_fit = rbind(df_prevsav %>% as.data.frame(), df_dat)

p2 = ggplot(df_fit, aes(x=data, y=prev_med, ymin=prev_lo, ymax=prev_hi)) +
  geom_errorbar(width=0) +
  geom_point() +
  facet_grid(type~Age_groups) +
  xlab("") +
  theme_bw() + scale_y_continuous("carriage prevalence") 
ggsave(paste0("Sims\\",sim,"\\fit.tiff"), p2, compression="lzw", unit="cm", width=23,height=15)
