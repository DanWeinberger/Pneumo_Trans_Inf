  
  #### setup

  Fit_type="Susceptibility" #"Susceptibility" ,"Transmissibility"=="Transmissibility"){Fit_type_num = 2}
  resolution=1000
  
  require("ggplot2")
  require("dplyr")
  require("tidyr")
  require("MultinomialCI")
  require("binom")
  require("gridExtra")
  require("cowplot")
  require("car")
  require("networkD3")
  require("scales")
  
  # set up data
  df_Kilifi=data.frame(Setting="Kilifi",
                       Age_groups=c("<1y","1-5y","6-14y","15-20y","21-49y","50+"),
                       Age_group_upper = c(1,6,15,21,50,62),
                       Population=c(9617,45170,68547,33289,72143,24214),
                       VT.prev=c(0.41,0.338,.146,.142,.072,.038),#PCV10
                       NVT.prev=c(.46,.44,.39,.25,.22,.20),
                       N.carr=c(39,127,82,56,97,104),
                       VT.clear=c(.062,.143,.343,.343,.343,.343)/7, #weekly conversion to daily. 1-5y = (1*0-2y +1.5*2-3.5y +1.5*35-5y)/4
                       NVT.clear=c(.086,.159,.349,.349,.349,.349)/7)
  df_EW=data.frame(Setting="EW",
                   Age_groups=c("<1y","1-5y","6-14y","15-20y","21-49y","50+"),
                   Age_group_upper = c(1,6,15,21,50,82),
                   Population=c(629200,3009200,5999000,4132900,21368700,17905900), #2004
                   VT.prev=c(0.431,0.254,.077,.0,.0039,.0), # PCV13 #first sample of each person in the longitudinal study 2001/2
                   NVT.prev=c(.118,.145,.096,.125,.022,.0),
                   N.carr=c(51,138,52,8,232,2),
                   VT.clear=1/c(72,34.8,18,17.5,17,17), #conversion from duration in days to clearance rate from melegaro et al
                   NVT.clear=1/c(72,34.8,18,17.5,17,17))
  df_NT=data.frame(Setting="NT",
                   Age_groups=c("<1y","1-5y","6-17y","18-49y","50+"),
                   Age_group_upper = c(1,6,18,50,76),
                   Population=c(2094,15239,40324,115538,37544),
                   VT.prev=c(0.266,.245,.018,.007,.008), #PCV13 from PhD thesis from Carla rather than Lay's team!
                   NVT.prev=c(.074,.202,.092,.012,.013),
                   N.carr=c(41,63,55,262,98),
                   VT.clear=1/c(72,34.8,18,17,17), #from E&W
                   NVT.clear=1/c(72,34.8,18,17,17))
  df_Fin=data.frame(Setting="Fin",
                    Age_groups=c("<1y","1-5y","6-17y","18+"),
                    Age_group_upper = c(1,6,18,81),
                    Population=c(56748,284828,764192,4130843),
                    VT.prev=c(203/1207,422/1235,30/258,25/789)*(1-.319), #PCV10 + 6A from Nurhonen model 
                    NVT.prev=c(203/1207,422/1235,30/258,25/789)*(.319), 
                    N.carr=(c(1200,1235,258,789)/3) %>% round(0), #conservatively accounting for random effects / the resuts are pooled across 3 surveys
                    VT.clear=c(.67,.67,1,1)/30, #convert to daily from monthly 
                    NVT.clear=c(.67,.67,1,1)/30)
  
  df_all_set = list(Kilifi= df_Kilifi, EW=df_EW, NT=df_NT, Fin=df_Fin)
  
  ###################################
  # compare demogrphics -  distributions
  ###################################  
  df_all=NULL
  for(i in c("Kilifi","EW","NT","Fin")){
    
    sim=paste0("sim_",i,"_",substr(Fit_type,1,3))
    load(file=paste0("Sims\\",sim,"\\R_out.Rdata"))
    df=df_all_set[[i]]
    df$Pop_prop=cumsum(df$Population/sum(df$Population))
    df$VT.carriers=df$VT.prev *df$Population
    df$prop.VT.carriers=cumsum(df$VT.carriers/sum(df$VT.carriers))
    df_all = df_all %>% rbind(df)
  }
  df_all$ab_hi=gsub("<1y","1",df_all$Age_groups)
  df_all$ab_hi=gsub("1-5y","5",df_all$ab_hi)
  df_all$ab_hi=gsub("6-17y","17",df_all$ab_hi)
  df_all$ab_hi=gsub("18-49y","49",df_all$ab_hi)
  df_all$ab_hi=gsub("50\\+","80",df_all$ab_hi)
  df_all$ab_hi=gsub("6-14y","14",df_all$ab_hi)
  df_all$ab_hi=gsub("15-20y","20",df_all$ab_hi)
  df_all$ab_hi=gsub("21-49y","49",df_all$ab_hi)
  df_all$ab_hi=gsub("18\\+","80",df_all$ab_hi)
  df_all$ab_hi = df_all$ab_hi %>% as.numeric()
  
  
  df_all = df_all[,c("ab_hi","Pop_prop","VT.prev","NVT.prev","Setting","prop.VT.carriers")] %>% rbind(data.frame(ab_hi=0, Pop_prop=0, VT.prev=c(.41,.266,.431,.121), NVT.prev=c(0.46,.074,.118,.048), Setting=c("Kilifi","NT","EW","Fin"), prop.VT.carriers=0))
  
  p_demographics <- ggplot(df_all, aes(x=ab_hi, y=Pop_prop, group=Setting, col=Setting)) +
    geom_step(direction="vh")+
    theme_bw() +
    xlab("Age (in years)") + ylab("Proportion of the\npopulation")
  ggsave(paste("Sims/Combined/",Fit_type,"/demographics.tiff",sep=""),p_demographics, compression="lzw", unit="cm", width=12,height=6, dpi=resolution)

  p_carriers <- ggplot(df_all, aes(x=ab_hi, y=prop.VT.carriers, group=Setting, col=Setting)) +
    geom_step(direction="vh") +
    theme_bw() +
    xlab("Age (in years)") + ylab("Proportion of the\ncarrying population")
  ggsave(paste("Sims/Combined/",Fit_type,"/Carriers.tiff",sep=""),p_carriers, compression="lzw", unit="cm", width=12,height=6, dpi=resolution)
  
  df_all_car = df_all %>% gather(Type, prevalence, -Setting, -Pop_prop, -ab_hi, -prop.VT.carriers)
  p_carr <- ggplot(df_all_car, aes(x=ab_hi, y=prevalence, group=Setting, col=Setting)) +
    geom_step(direction="vh", alpha=0.6) +
    facet_grid(.~Type) +
    theme_bw() +
    xlab("Age (in years)") + ylab("prevalence")
  ggsave(paste("Sims/Combined/",Fit_type,"/Carriage.tiff",sep=""),p_carr, compression="lzw", unit="cm", width=18,height=6, dpi=resolution)
  
  ###################################
  # show fit
  ###################################  
  df_all=NULL
  for(i in c("Kilifi","EW","NT","Fin")){
    sim=paste0("sim_",i,"_",substr(Fit_type,1,3))
    load(file=paste0("Sims\\",sim,"\\R_out.Rdata"))
    no.agegps=dim(df)[1]
    
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
    df_fit = rbind(df_prevsav %>% data.frame(), df_dat)
    df_fit$setting=i
    df_all = df_all %>% rbind(df_fit)
  }
  
  df_all$Age_groups = factor(df_all$Age_groups, levels = c("<1y", "1-5y", "6-14y", "6-17y", "15-20y", "18-49y","18+", "21-49y", "50+"))
  
 
  p_fit = ggplot(data=subset(df_all, data=="data"), aes(x=Age_groups, y=prev_med, ymax=prev_hi,ymin=prev_lo, group=setting)) +
    geom_ribbon(data=subset(df_all, data=="model"), alpha=0.2) +
    geom_line(data=subset(df_all, data=="model"), alpha=0.8) +
    geom_pointrange()+
    facet_grid(type~setting, scales="free_x") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  ggsave(paste("Sims/Combined/",Fit_type,"/Fit.tiff",sep=""),p_fit, compression="lzw", unit="cm", width=30,height=14, dpi=resolution)
  
  
  ###################################
  # contact matrices
  ###################################  
  df_all=NULL
  for (i in c("EW","Kilifi","NT","Fin")){
    df=read.csv(file=paste0("Data/AvNumberOfContacts_",i,"_new.csv"))
    names(df)[1]="contact"
    # get observed age stratified probability of two specific persons contacting each other
    df[,-1] = df[,-1] / df_all_set[[i]][,"Population"]
    df_s = df %>% gather(source,value,-contact)
    df_s$setting = i
    df_s$value = df_s$value / mean(df_s$value) # standardise to make surveys compareable
    df_all = df_all %>% rbind(df_s)
  }
  df_all$source =  gsub("X","",df_all$source)
  df_all$source =  gsub("\\.","-",df_all$source)
  df_all$source =  gsub("50-","50+",df_all$source)
  df_all$source =  gsub("18-","18+",df_all$source)
  df_all$source =  gsub("18\\+49y","18-49y",df_all$source)
  df_all$source =  gsub("-1y","<1y",df_all$source)
  
  df_all$contact = factor(df_all$contact, levels = c("<1y", "1-5y", "6-14y", "6-17y", "15-20y", "18-49y","18+", "21-49y", "50+"))
  df_all$source  = factor(df_all$source, levels = c("<1y", "1-5y", "6-14y", "6-17y", "15-20y", "18-49y", "18+", "21-49y", "50+"))
  
  p_contacts = ggplot(df_all, aes(x=source, y=contact, fill=value)) +
    geom_tile(color="white") +
    facet_wrap(~setting, scales="free") +
    scale_fill_gradient2(trans = "log", breaks=c(.25,.5,1,2,4), low="darkred", high="darkblue",
                        name = "Probability of contact for\ntwo specific individuals\nstandardised by average") +
    theme(panel.background = element_rect(colour = "white", fill="white"), 
          axis.line = element_line(color="white"),
          axis.ticks = element_line(color="white"),
          text = element_text(size=8),
          axis.text.x = element_text(angle =90, hjust = .5))+
    xlab("age of participant") + ylab("age of contact")
    
  ggsave(paste("Sims/Combined/",Fit_type,"/Contacts.tiff",sep=""),p_contacts, compression="lzw", unit="cm", width=19,height=17, dpi=resolution)
  
  
  ###################################
  # bar-chart who infects whom 
  ################################### 
  df_all=NULL
  for(i in c("Kilifi","EW","NT","Fin")){
    sim=paste0("sim_",i,"_",substr(Fit_type,1,3))
    load(file=paste0("Sims\\",sim,"\\R_out.Rdata"))
    InfectionMatrix_summary = InfectionMatrix %>% apply(2:3,quantile, probs=c(.025,.5,.975))
    InfectionMatrix_median = InfectionMatrix_summary["50%",,]
    dd=dim(InfectionMatrix_median)[1]
    if(i=="NT") InfectionMatrix_median[1,1]=0.01 # get rid of 0
    
    InfectionMatrix_median_std = InfectionMatrix_median %>% apply(1,function(x) x/sum(x)) %>% t()
    df = InfectionMatrix_median_std %>% as.data.frame() %>% gather(key="Age_infector", value="value")
    df$Age_infectee = rep(colnames(InfectionMatrix_median_std),dd)
    df$Setting = i
    
    df$Age_infectee = factor(df$Age_infectee, levels = unique(df$Age_infectee))
    df$Age_infector = factor(df$Age_infector, levels = unique(df$Age_infector))
    
    df_all = rbind(df_all,df)
    
  } 
  df_all$Age_infectee = factor(  df_all$Age_infectee, levels = c("<1y", "1-5y", "6-14y", "6-17y", "15-20y", "18-49y","18+", "21-49y", "50+"))
  df_all$Age_infector = df_all$Age_infector %>% recode('"<1y"="Infant"; "1-5y"="Pre-school";"6-14y"="School"; "6-17y"="School"; "15-20y"="School"; "18-49y"="Adult";"18+"="Adult"; "21-49y"="Adult"; "50+"="Adult"')
  df_all$Age_infector = factor(df_all$Age_infector, levels = c("Infant","Pre-school","School","Adult"))
  
  p_bar = ggplot(df_all,aes(x=Age_infectee, y=value, fill=Age_infector)) +
    geom_bar(stat="identity", color="white", alpha=0.7) +
    facet_grid(Setting~., scale="free") +
    scale_y_continuous(labels=percent) +
    coord_flip() +
    labs(y = "proportion of infections", x="Age of infected person", fill="Age of\ninfecting person") 
  ggsave(paste("Sims/Combined/",Fit_type,"/infections_bar.tiff",sep=""),p_bar, compression="lzw", unit="cm", width=15,height=15, dpi=resolution)
  

  ###################################
  # infections in U1 
  ###################################  
  df.all.annual.infections = NULL
  for(i in c("Kilifi","EW","NT","Fin")){
    sim=paste0("sim_",i,"_",substr(Fit_type,1,3))
    load(file=paste0("Sims\\",sim,"\\R_out.Rdata"))
    df.annual.infections <- (InfectionMatrix*365) %>% as.data.frame.table(responseName = "infections") 
    df.annual.infections$setting=i
    df.annual.infections <- merge(df.annual.infections, df[,c("Age_groups","Population")], by.x="Infector", by.y="Age_groups") 
    df.all.annual.infections = rbind(df.all.annual.infections,df.annual.infections)
  }
  
  #calc # of annual infections caused per average person
  df.all.annual.infections %>% subset(infectee == "<1y") %>%
    group_by(setting,MCMC.sample) %>% summarise(cum.infections = sum(infections)/sum(Population)) %>%
    group_by(setting) %>% summarise(mid = median(cum.infections),
                                    low = quantile(cum.infections, probs=.025),
                                    hi = quantile(cum.infections, probs=.975))
  
  # how many children infect infants (infants = infectee) per year?
  df_infs <- df.all.annual.infections %>% subset(infectee == "<1y") %>% 
    mutate(cum.inf.pp=infections/Population) %>% 
    group_by(setting, Infector) %>% summarise(mid = median(cum.inf.pp),
                                              low = quantile(cum.inf.pp, probs=.025),
                                              hi = quantile(cum.inf.pp, probs=.975)) 
  df_infs %>% subset(Infector=="1-5y")
  
  # what do pre & school children contribute to transmission to infants?
  df.all.annual.infections.cumInfectees <- df.all.annual.infections %>% 
    subset(infectee == "<1y")
  df.all.annual.infections.cumInfectees$Agegp <- df.all.annual.infections.cumInfectees$Infector %>% 
    recode('"<1y"="Infant"; "1-5y"="Pre-school";"6-14y"="School"; "6-17y"="School"; "15-20y"="School"; "18-49y"="Adult";"18+"="Adult"; "21-49y"="Adult"; "50+"="Adult"')
  dfInfsProp <- df.all.annual.infections.cumInfectees %>% 
    group_by(setting, MCMC.sample) %>% mutate(cumInf = sum(infections)) %>%
    mutate(infProp=infections/cumInf) 
  dfInfsProp %>% subset(Agegp %in% c("Pre-school","School")) %>% 
    group_by(setting, MCMC.sample) %>% summarise(infPropSum = sum(infProp)) %>% 
    group_by(setting) %>% summarise(mid = median(infPropSum),
                                    low = quantile(infPropSum, probs=.025),
                                    hi = quantile(infPropSum, probs=.975))
  dfInfsProp %>% subset(Agegp %in% c("Pre-school")) %>% 
    group_by(setting, MCMC.sample) %>% summarise(infPropSum = sum(infProp)) %>% 
    group_by(setting) %>% summarise(mid = median(infPropSum),
                                    low = quantile(infPropSum, probs=.025),
                                    hi = quantile(infPropSum, probs=.975))
  dfInfsProp %>% subset(Agegp %in% c("School")) %>% 
    group_by(setting, MCMC.sample) %>% summarise(infPropSum = sum(infProp)) %>% 
    group_by(setting) %>% summarise(mid = median(infPropSum),
                                    low = quantile(infPropSum, probs=.025),
                                    hi = quantile(infPropSum, probs=.975))
  dfInfsProp %>% subset(Agegp %in% c("Infant")) %>% 
    group_by(setting, MCMC.sample) %>% summarise(infPropSum = sum(infProp)) %>% 
    group_by(setting) %>% summarise(mid = median(infPropSum),
                                    low = quantile(infPropSum, probs=.025),
                                    hi = quantile(infPropSum, probs=.975))
  
  # plot
  dfInfsProp %<>% group_by(setting, Infector) %>% summarise(mid = median(infProp),
                                                            low = quantile(infProp, probs=.025),
                                                            hi = quantile(infProp, probs=.975)) %>%
    mutate(outcome="proportion of all infections")
  df_infs %<>% mutate(outcome="annual infections per person") %>% rbind(dfInfsProp)
  df_infs$Agegp <- df_infs$Infector %>% recode('"<1y"="Infant"; "1-5y"="Pre-school";"6-14y"="School"
                                               ; "6-17y"="School"; "15-20y"="School"; "18-49y"="Adult"
                                               ;"18+"="Adult"; "21-49y"="Adult"; "50+"="Adult"')
  df_infs$Infector = factor(df_infs$Infector, levels = c("<1y", "1-5y", "6-14y", "6-17y", "15-20y", "18-49y","18+", "21-49y", "50+"))
  df_infs$Agegp = factor(df_infs$Agegp, levels = c("Infant","Pre-school","School","Adult"))
  p.u1 <- df_infs %>% ggplot(aes(x=Infector, y=mid, ymax=hi, ymin=low, color=Agegp)) +
    geom_linerange() +
    geom_point(size=1.5) +
    facet_grid(outcome~setting, scales="free") + theme_bw() +
    ylab("") + xlab("") +
    scale_y_log10(breaks=c(0.1,.5,1,5,10,50,100), minor_breaks=c(seq(.1,1,by=0.1),seq(2,10,by=1),seq(20,100,by=10))) +
    theme(text = element_text(size=8), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
    scale_color_discrete(name="Age group")
  ggsave(paste("Sims/Combined/",Fit_type,"/infections_U1.tiff",sep=""),p.u1, compression="lzw", unit="cm", width=19,height=10, dpi=resolution)
  
  
  
  ###################################
  # infections in all ages 
  ###################################  
  df.all.annual.infections = NULL
  for(i in c("Kilifi","EW","NT","Fin")){
    sim=paste0("sim_",i,"_",substr(Fit_type,1,3))
    load(file=paste0("Sims\\",sim,"\\R_out.Rdata"))
    df.annual.infections <- (InfectionMatrix*365) %>% as.data.frame.table(responseName = "infections") 
    df.annual.infections$setting=i
    df.annual.infections <- merge(df.annual.infections, df[,c("Age_groups","Population")], by.x="Infector", by.y="Age_groups") 
    df.all.annual.infections = rbind(df.all.annual.infections,df.annual.infections)
  }
  
  #calc # of annual infections caused per average person
  df.all.annual.infections %>% group_by(setting,MCMC.sample, Infector, Population) %>% summarise(cum.infections = sum(infections)) %>%
    group_by(setting,MCMC.sample) %>% summarise(cum.infections = sum(cum.infections)/sum(Population)) %>%
    group_by(setting) %>% summarise(mid = median(cum.infections),
                                    low = quantile(cum.infections, probs=.025),
                                    hi = quantile(cum.infections, probs=.975))
  
  # how many people infect infants (infants = infector) per year?
  df_infs <- df.all.annual.infections %>% 
    group_by(setting,MCMC.sample, Infector, Population) %>% summarise(cum.infections = sum(infections)) %>% 
    mutate(cum.inf.pp=cum.infections/Population) %>% 
    group_by(setting, Infector) %>% summarise(mid = median(cum.inf.pp),
                                    low = quantile(cum.inf.pp, probs=.025),
                                    hi = quantile(cum.inf.pp, probs=.975)) 
  df_infs %>% subset(Infector=="<1y")
    
  # what do pre & school children contribute to transmission?
  df.all.annual.infections.cumInfectees <- df.all.annual.infections %>% group_by(setting,MCMC.sample, Infector, Population) %>% summarise(cum.infections = sum(infections))
  df.all.annual.infections.cumInfectees$Agegp <- df.all.annual.infections.cumInfectees$Infector %>% recode('"<1y"="Infant"; "1-5y"="Pre-school";"6-14y"="School"; "6-17y"="School"; "15-20y"="School"; "18-49y"="Adult";"18+"="Adult"; "21-49y"="Adult"; "50+"="Adult"')
  dfInfsProp <- df.all.annual.infections.cumInfectees %>% 
    group_by(setting, MCMC.sample) %>% mutate(cumInf = sum(cum.infections)) %>%
    mutate(infProp=cum.infections/cumInf) 
  dfInfsProp %>% subset(Agegp %in% c("Pre-school","School")) %>% 
    group_by(setting, MCMC.sample) %>% summarise(infPropSum = sum(infProp)) %>% 
    group_by(setting) %>% summarise(mid = median(infPropSum),
                                    low = quantile(infPropSum, probs=.025),
                                    hi = quantile(infPropSum, probs=.975))
  # what do pre school children contribute to transmission?
  dfInfsProp %>% subset(Agegp %in% c("Pre-school")) %>% 
    group_by(setting, MCMC.sample) %>% summarise(infPropSum = sum(infProp)) %>% 
    group_by(setting) %>% summarise(mid = median(infPropSum),
                                    low = quantile(infPropSum, probs=.025),
                                    hi = quantile(infPropSum, probs=.975))
  # what do school children contribute to transmission?
  dfInfsProp %>% subset(Agegp %in% c("School")) %>% 
    group_by(setting, MCMC.sample) %>% summarise(infPropSum = sum(infProp)) %>% 
    group_by(setting) %>% summarise(mid = median(infPropSum),
                                    low = quantile(infPropSum, probs=.025),
                                    hi = quantile(infPropSum, probs=.975))
  # what do infants contribute to transmission?
  dfInfsProp %>% subset(Agegp %in% c("Infant")) %>% 
    group_by(setting, MCMC.sample) %>% summarise(infPropSum = sum(infProp)) %>% 
    group_by(setting) %>% summarise(mid = median(infPropSum),
                                    low = quantile(infPropSum, probs=.025),
                                    hi = quantile(infPropSum, probs=.975))

  # plot
  dfInfsProp %<>% group_by(setting, Infector) %>% summarise(mid = median(infProp),
                                    low = quantile(infProp, probs=.025),
                                    hi = quantile(infProp, probs=.975)) %>%
    mutate(outcome="proportion of all infections")
  df_infs %<>% mutate(outcome="annual infections per person") %>% rbind(dfInfsProp)
  df_infs$Agegp <- df_infs$Infector %>% recode('"<1y"="Infant"; "1-5y"="Pre-school";"6-14y"="School"
                                               ; "6-17y"="School"; "15-20y"="School"; "18-49y"="Adult"
                                               ;"18+"="Adult"; "21-49y"="Adult"; "50+"="Adult"')
  df_infs$Infector = factor(df_infs$Infector, levels = c("<1y", "1-5y", "6-14y", "6-17y", "15-20y", "18-49y","18+", "21-49y", "50+"))
  df_infs$Agegp = factor(df_infs$Agegp, levels = c("Infant","Pre-school","School","Adult"))
  p.all <- df_infs %>% ggplot(aes(x=Infector, y=mid, ymax=hi, ymin=low, color=Agegp)) +
    geom_linerange() +
    geom_point(size=1.5) +
    facet_grid(outcome~setting, scales="free") + theme_bw() +
    ylab("") + xlab("") +
    scale_y_log10(breaks=c(0.1,.5,1,5,10,50,100), minor_breaks=c(seq(.1,1,by=0.1),seq(2,10,by=1),seq(20,100,by=10))) +
    theme(text = element_text(size=8), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_color_discrete(name="Age group")
  ggsave(paste("Sims/Combined/",Fit_type,"/infections_all.tiff",sep=""),p.all, compression="lzw", unit="cm", width=19,height=10, dpi=resolution)
    
 
