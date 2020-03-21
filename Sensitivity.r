
##########
# additional sensitivity analyses
##########

require("tidyverse")
require("cowplot")
require("scales")

# pick scenario
setting="NT"
Fit_type="Susceptibility" #"Susceptibility" ,"Transmissibility"
NTwaningSens = T #alternative waning assumption for Nha trang?
lowComp = F #alternative assumption on competition parameter. 
nonPhysicalContacts = F #also include nonphysical contacts (no such recoreded in Kilifi)
resolution=300 # for plots


sim=paste0("Sim_",setting,"_",substr(Fit_type,1,3))
if(lowComp) sim=paste0(sim,"_lowComp") 
if(nonPhysicalContacts) sim=paste0(sim,"_allContacts")
if(NTwaningSens) sim=paste0(sim,"_SensWaning")

load(file=paste0("Sims\\",sim,"\\R_out.Rdata"))

InfectionMatrix_summary = InfectionMatrix %>% apply(2:3,quantile, probs=c(.025,.5,.975))
InfectionMatrix_median = InfectionMatrix_summary["50%",,]
dd=dim(InfectionMatrix_median)[1]
if(setting=="NT") InfectionMatrix_median[1,1]=0.01 # get rid of 0

InfectionMatrix_median_std = InfectionMatrix_median %>% apply(1,function(x) x/sum(x)) %>% t()
df = InfectionMatrix_median_std %>% as.data.frame() %>% gather(key="Age_infector", value="value")
df$Age_infectee = rep(colnames(InfectionMatrix_median_std),dd)
df$Setting = setting

df$Age_infectee = factor(df$Age_infectee, levels = unique(df$Age_infectee))
df$Age_infector = factor(df$Age_infector, levels = unique(df$Age_infector))

df$Age_infectee = factor(  df$Age_infectee, levels = c("<1y", "1-5y", "6-14y", "6-17y", "15-20y", "18-49y","18+", "21-49y", "50+"))
df$Age_infector = df$Age_infector %>% car::recode('"<1y"="Infant"; "1-5y"="Pre-school";"6-14y"="School"; "6-17y"="School"; "15-20y"="School"; "18-49y"="Adult";"18+"="Adult"; "21-49y"="Adult"; "50+"="Adult"')
df$Age_infector = factor(df$Age_infector, levels = c("Infant","Pre-school","School","Adult"))

p_bar = ggplot(df,aes(x=Age_infectee, y=value, fill=Age_infector)) +
  geom_bar(stat="identity", color="white", alpha=0.7) +
  facet_grid(Setting~., scale="free") +
  scale_y_continuous(labels=percent) +
  coord_flip() +
  labs(y = "proportion of infections", x="Age of infected person", fill="Age of\ninfecting person") 
ggsave(paste("Sims/",sim,"/infections_barplot.tiff",sep=""),p_bar, compression="lzw", unit="cm", width=15,height=7, dpi=resolution)

