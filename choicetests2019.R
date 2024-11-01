library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
setwd("~/MyDocs/MEGA/UCI/Schiedea/Analysis/scent/input/moth_choices")
moths <- read_csv("Oahu scent choices 2019 - Scent additions (3).csv") %>% mutate_if(is.character, factor)
moths$treat <- factor(ifelse(moths$trt %in% c("KAAL2","KAAL3"), "KAAL", as.character(moths$trt)))
moths$trial <- factor(with(moths, paste(date, observer,stick)))
moths$total.trt <- moths$visit.trt + moths$approach.trt 
moths$total.ctrl <- moths$visit.ctrl + moths$approach.ctrl
moths$prop.visit <- with(moths, visit.trt / (visit.trt + visit.ctrl))
moths$prop.approach <- with(moths, approach.trt / (approach.trt + approach.ctrl))
moths$prop.total <- with(moths, total.trt / (total.trt + total.ctrl))
moths$duration <- with(moths, as.numeric(stop - start)/3600)
moths$rate.total <- with(moths, ((total.ctrl + total.trt) / duration)/3600)

(obs.hrs <- moths %>% group_by(treat, base) %>% 
    summarise_at(vars(duration), funs(sum)) %>%
    unite(base.treat, c(base,treat), remove=F))

#RATE: average approaches or visits to either treatments or control
rates <- moths %>% group_by(base, treat) %>% 
  summarise_at(vars(total.trt, total.ctrl, rate.total), funs(mean))
#COUNTS
counts <- moths %>% group_by(base, treat) %>% 
  summarise_at(vars(approach.trt, approach.ctrl, visit.trt, visit.ctrl, total.trt, total.ctrl), funs(sum)) %>%
  mutate(approach=approach.trt+approach.ctrl, visit=visit.trt+visit.ctrl, total=total.trt+total.ctrl)
counts.samp <- full_join(counts, moths %>% group_by(base, treat) %>% summarize(n=n()))
#PROPORTIONS
props <- moths %>% group_by(base, treat) %>% 
  summarise_at(vars(prop.approach, prop.visit, prop.total), funs(mean(., na.rm = TRUE)))
#N TESTS
sampsize <- moths %>% group_by(base, treat) %>% 
  summarize_at(vars(total.trt, total.ctrl, duration), funs(sum, length)) %>%
  unite(base.treat, c(base,treat), remove=F)

write.csv(cbind(counts, props, sampsize), "choices2019.csv")

mothsK <- moths %>% filter(treat =="KAAL")
mothsE <- moths %>% filter(treat !="KAAL")

###Raw proportion plot###
library(tidyr)
mothsK.l <- mothsK %>% gather("interaction", "prop", prop.visit:prop.total)  %>% filter(interaction != "prop.approach") %>% mutate(n=ifelse(interaction=="prop.total",total.trt+total.ctrl,visit.trt+visit.ctrl))
ggplot(mothsK.l, aes(x=base, y=prop, fill=interaction)) + 
  geom_violin() +
  geom_point(aes(size=n), alpha=0.8, position=position_jitterdodge(dodge.width=0.9, jitter.width=0.2)) + 
  geom_hline(yintercept=0.5) + 
  scale_x_discrete("Base", labels=rev(c("S. kealiae", "S. globosa"))) +
  scale_y_continuous("Percent of interactions with S. kaalae side", limits=c(0,1), labels = scales::percent) +
  scale_fill_grey(end=0.6, labels=c("Approach","Visit")) +
  scale_size_area()+
  theme_minimal()  + theme(legend.position = "bottom")

###visitation plots###
##all these visitation rate plots are misleading - should be done with LSM percentages to each side!

#wide to long
moths.l <- moths %>% 
  dplyr::select(-starts_with("prop")) %>% dplyr::select(-starts_with("total")) %>%
  gather("int.side", "count", approach.trt:visit.ctrl, factor_key=T) %>% 
  separate(int.side, c("interaction","side"), remove=F) %>%
  unite(base.treat, c(base,treat), remove=F) %>%
  mutate(rate=count/obs.hrs$duration[match(base.treat, obs.hrs$base.treat)]) %>% #this is the counts divided by the total hours for that base and treatment, so when added up in ggplot it is the average visitation rate - hacky!
  mutate_if(is.character, as.factor)

library(ggplot2)
moths.l.plot <- filter(moths.l, treat %in% c("KAAL", "KAMIX2", "HOMIX2") & base.treat != "GLOB_KAMIX2")
moths.sampsize.plot <- sampsize %>% 
  filter(treat %in% c("KAAL", "KAMIX2", "HOMIX2") & base.treat != "GLOB_KAMIX2") %>%
  mutate(label = paste0("N=",total.trt_length, "  I=", total.trt_sum+total.ctrl_sum))

moths.l.plot$int.side <- factor(moths.l.plot$int.side, levels=c("approach.ctrl","visit.ctrl","visit.trt","approach.trt"))

ggplot(moths.l.plot, aes(x=base.treat, y=rate*(2*as.integer(as.factor(side))-3), fill=int.side, group=interaction)) + 
  geom_col() + coord_flip() + geom_hline(yintercept=0) + 
  #geom_text(mapping=aes(x=base.treat, label=label), y=-4.5, data=moths.sampsize.plot, inherit.aes=F) +
  scale_x_discrete("Treatment", labels=rev(c("stick + bagged mix K","stick + bagged mix H", "S. kealiae + bagged S. kaalae", "S. globosa + bagged S. kaalae"))) +
  scale_y_continuous("Interactions per hour", labels=abs, breaks=-5:5, limits=c(-5,5)) +
  scale_fill_brewer("", type="div", palette="PuOr", direction = -1, labels=c("Approach: control", "Visit: control", "Visit: treatment",  "Approach: treatment")) +
  #ggtitle("Moth choice tests 2019", subtitle="P. brevipalpis interactions with either bags or inflorescences") + 
  theme_minimal()  + theme(legend.position = "top")
ggsave("choices2019.png", height=3, width=7, dpi=300)

#plot with se
moths.tot <- moths %>% 
  dplyr::select(-starts_with("prop")) %>%
  gather(int.side, count, c(approach.trt:visit.ctrl, total.trt:total.ctrl), factor_key=T) %>% 
  separate(int.side, c("interaction","side"), remove=F) %>%
  unite(base.treat, c(base,treat), remove=F) %>%
  mutate(vrate=count/duration) %>% #this is the visitation rate/hr for that choice test
  mutate_if(is.character, as.factor) %>%
  filter(treat %in% c("KAAL", "KAMIX2", "HOMIX2") & base.treat != "GLOB_KAMIX2")

moths.tot.agg <- moths.tot %>% 
  group_by(base.treat, int.side, interaction, side) %>% 
  summarize(vratem=mean(vrate), se=sd(vrate)/sqrt(n()), n=n()) %>%
  filter(interaction !="approach") %>% 
  arrange(-row_number()) #flip for plotting bars on top of each other

moths.tot.agg$int.side <- factor(moths.tot.agg$int.side, levels=c("total.ctrl","visit.ctrl","visit.trt","total.trt"))
#moths.tot.agg$interaction <- factor(moths.tot.agg$interaction, levels=c("total","visit"))

ggplot(moths.tot.agg, aes(x=base.treat, y=vratem*(2*as.integer(as.factor(side))-3), fill=int.side, group=interaction)) + 
  geom_col(position="identity") + coord_flip() + geom_hline(yintercept=0) + 
  geom_errorbar(aes(ymin=(vratem-se)*(2*as.integer(as.factor(side))-3), ymax=(vratem+se)*(2*as.integer(as.factor(side))-3)), width=0.2, position=position_dodge(width=0.1)) +
  scale_x_discrete("Treatment", labels=rev(c("S. kealiae + bagged S. kaalae", "S. globosa + bagged S. kaalae"))) +
  scale_y_continuous("Interactions per hour", labels=abs, breaks=-9:9, limits=c(-9,9)) +
  scale_fill_brewer("", type="div", palette="PuOr", direction = -1, labels=c("Approach: control", "Visit: control", "Visit: treatment",  "Approach: treatment")) +
  theme_minimal()  + theme(legend.position = "top")

moths.tot.agg$int.side <- factor(moths.tot.agg$int.side, levels=c("total.ctrl","total.trt","visit.ctrl","visit.trt"))
ggplot(moths.tot.agg, aes(x=base.treat, y=vratem, fill=int.side)) + 
  geom_col(position="dodge") + geom_hline(yintercept=0) + 
  geom_errorbar(aes(ymin=(vratem-se), ymax=(vratem+se)), width=0.2, position=position_dodge(width=0.9)) +
  scale_x_discrete("Treatment", labels=rev(c("S. kealiae", "S. globosa"))) +
  scale_y_continuous("Interactions per hour", labels=abs, breaks=1:9) +
  scale_fill_brewer("", type="div", palette="PuOr", direction = -1, labels=c("Approach: control", "Approach: S. kaalae scent", "Visit: control", "Visit: S. kaalae scent")) +
  theme_minimal()  + theme(legend.position = "bottom")
    
#####t-tests####
t.test(asin(sqrt(mothsK$prop.total))/pi, mu = 0.5)

t.test(asin(sqrt(mothsK$prop.approach[mothsK$base=="GLOB"]))/pi, mu = 0.5)
t.test(asin(sqrt(mothsK$prop.approach[mothsK$base=="KEAL"]))/pi, mu = 0.5)

t.test(asin(sqrt(mothsK$prop.visit[mothsK$base=="GLOB"]))/pi, mu = 0.5)
t.test(asin(sqrt(mothsK$prop.visit[mothsK$base=="KEAL"]))/pi, mu = 0.5)

t.test(asin(sqrt(mothsE$prop.approach[mothsE$base=="stick" & mothsE$treat=="HOMIX2"]))/pi, mu = 0.5)
t.test(asin(sqrt(mothsE$prop.approach[mothsE$base=="stick" & mothsE$treat=="KAMIX2"]))/pi, mu = 0.5)


######GLMs#####
library(glmmTMB)

visit.betabin  <- glmmTMB(cbind(visit.trt, visit.ctrl)~base, data=mothsK, family="betabinomial") 
app.betabin    <- glmmTMB(cbind(approach.trt, approach.ctrl)~base, data=mothsK, family="betabinomial") 
total.betabin  <- glmmTMB(cbind(total.trt, total.ctrl)~base, data=mothsK, family="betabinomial") 

summary(total.betabin)
summary(visit.betabin)

#combine totals and visits - actually, doesn't make sense to model them together since dispersions are different --- unless dispersion allowed to vary! with dispformula
mothsK.test <- 
  bind_rows(list(total = mothsK %>% dplyr::select(trt, base, total.trt, total.ctrl) %>% rename(treat=total.trt, ctrl=total.ctrl), visit = mothsK %>%  dplyr::select(trt, base, visit.trt, visit.ctrl) %>% rename(treat=visit.trt, ctrl=visit.ctrl)), .id="interaction") %>% unite("base.int", base, interaction, remove=F)
  
(all.betabin <- glmmTMB(cbind(treat,ctrl)~base.int, data=mothsK.test, family="betabinomial", dispformula = ~interaction))
summary(all.betabin)

#test if means = 50% (I think)
library(emmeans)
total.betabin %>% ref_grid() %>% emmeans(~ base) %>% test()
visit.betabin %>% ref_grid() %>% emmeans(~ base) %>% test()
all.betabin %>% ref_grid() %>% emmeans(~ base.int) %>% test()

#separate tests for each combination
for(i in unique(mothsK.test$base.int)) {
  print(i)
  mod <- glmmTMB(cbind(treat, ctrl)~1, data=mothsK.test %>% filter(base.int==i), family="betabinomial") 
  #print(mod)
print(summary(mod))
}

mothsK.test %>% group_by(base.int) %>% summarize_if(is.numeric, sum)

#####EMMEANS########
get.emmeans <- . %>% ref_grid() %>% emmeans(~ base) %>% summary %>%
  mutate(response=plogis(emmean), uSE = plogis(emmean+SE), lSE = plogis(emmean-SE))

emm.both <- data.frame(bind_rows(list(total=get.emmeans(total.betabin), visit=get.emmeans(visit.betabin)), .id="interaction"))
emm.all <- all.betabin %>% ref_grid() %>% emmeans(~ base.int) %>% test %>%
  mutate(response=plogis(emmean), uSE = plogis(emmean+SE), lSE = plogis(emmean-SE)) %>%
  separate(base.int, c("base", "interaction"))

####GLM proportion plot#####
countsK <- counts.samp %>% filter(treat=="KAAL") %>% 
  gather("interaction", "count", approach:total) %>% filter(interaction !="approach") %>% 
  mutate(fraction = paste(ifelse(interaction=="visit", visit.trt, total.trt), "/", count)) %>%
  left_join(obs.hrs)

(glm.plot <- ggplot(emm.all, aes(x=base, y=response, fill=interaction)) + 
  geom_col(position="dodge") + 
  geom_errorbar(aes(ymin=lSE, ymax=uSE), width=0.2, position=position_dodge(width=0.9)) +
  geom_hline(yintercept=0.5, linetype=2) + 
  geom_text(aes(label=paste("P =",round(p.value,3),ifelse(p.value < 0.05,"*","")), y=0.95), position=position_dodge(width=0.9))+
  scale_x_discrete("Inflorescence presented", labels=rev(c("S. kealiae", "S. globosa"))) +
  scale_y_continuous(expression("Interactions on side with bagged"~italic("S. kaalae")), limits=c(0,1), labels = scales::percent, expand=c(0,0)) +
  scale_fill_grey("Interaction",start=0.7,end=0.4, labels=c("Approach","Visit")) + 
  theme_classic()  + theme(legend.position = "top", axis.text=element_text(color="black", size=11), 
                           axis.text.x=element_text(face="italic"), axis.ticks.x = element_blank()))

#for slides
glm.plot +  guides(fill=F) +
  geom_text(data=countsK, aes(label=paste0(fraction,"\n",gsub("total","approache",interaction),"s"), y=0.22),  position=position_dodge(width=0.9))+
  geom_text(data=countsK%>%filter(interaction=="total"), aes(label=paste(n,"choice tests","\n",round(duration,0),"observer-h"), y=0.08))


#for manuscript
glm.plot  +
  geom_text(data=countsK, aes(label=fraction, y=0.05),  position=position_dodge(width=0.9))

ggsave("choices2019_prop_manu.pdf",height = 5, width=5)
ggsave("choices2019_prop_manu.png",height = 5, width=5, dpi=300) 

