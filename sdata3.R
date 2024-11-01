#Read in inventory and GCMS data, process and reshape
#John Powers
#Jan 15 2019

#####Packages#####

library(knitr)
library(shiny)
library(rgl)
library(scales)
library(XML)
library(png)

library(googlesheets)
library(readr)
library(plyr); library(dplyr)
library(reshape2)
library(stringr)
library(maptools)

#library(ade4) #also contains CCA, load before vegan!
#library(candisc) #also contains scores, load before vegan!
library(vegan3d)
library(indicspecies)
library(phytools)
library(geomorph) #kmult
library(nlme)
library(MASS) #lda
library(corrplot)
library(cluster)
library(varSelRF)
library(beanplot)

#####Summarize#####
sampl <-gs_title("Schiedea Volatiles Sampling")
s <- gs_read(sampl, col_types = cols(.default = col_character(),
                                     Flow = col_double(),
                                     Mass = col_double(),
                                     RunDate = col_date(format = ""),
                                     SampleDate = col_date(format = ""),
                                     Start = col_time(format = ""),
                                     Stop = col_time(format = ""),
                                     Equi = col_double(),
                                     Duration = col_double(),
                                     Total = col_double()
))
#save(sampl, s, file="sch_vol_samp.Rdata")
load("sch_vol_samp.Rdata")
s$SampleDateStart <- paste(s$SampleDate, s$Start)
is.na(s$SampleDateStart) <- (is.na(s$SampleDate) + is.na(s$Start))>0
s$SampleDateStart <- as.POSIXct(s$SampleDateStart)
s$SampleDateStop <- paste(s$SampleDate, s$Stop)
is.na(s$SampleDateStop) <- (is.na(s$SampleDate) + is.na(s$Stop))>0
s$SampleDateStop <- as.POSIXct(s$SampleDateStop)

greenhouse <- matrix(c(-117.8475746, 33.6472914), nrow=1)
sunsets <- do.call("c", lapply(lapply(s$SampleDateStart, sunriset, crds=greenhouse, direction="sunset", POSIXct.out=T), function(x) x$time))
s$StartSunset <- difftime(s$SampleDateStart, sunsets,units="hours")
s$StopSunset <- difftime(s$SampleDateStop,sunsets,units="hours")
solarnoons <- do.call("c", lapply(lapply(s$SampleDateStart, solarnoon, crds=greenhouse, POSIXct.out=T), function(x) x$time))
s$StartNoon <- difftime(s$SampleDateStart,solarnoons,units="hours")

#load("inventory.Rdata")
s$Speciespop <- paste(s$Species,s$Population)
sf <- s[s$GC %in% c("G"),]
st <- count(sf, Speciespop,Sex)
stt <- data.frame(do.call('rbind', strsplit(as.character(st$Speciespop),' ',fixed=TRUE)), st$Sex, st$n)
names(stt) <- c("Species","Population","Sex","Samples")
sttw <- dcast(stt, Species+Population~Sex, sum, value.var="Samples")
nambi <- sttw[sttw$Species=="ambient", "-"]
sttw <- sttw[sttw$Species!="ambient",-which(names(sttw) %in% "-")]
# put it back
#gs <- sampl %>% gs_ws_new(ws_title = "Summary", input = sttw, trim = TRUE, verbose = FALSE)

datafile <- "svhyb_ambient_blanks.csv"
data <- read.table(datafile, header=TRUE, sep="\t", as.is=T)
data <- data[!duplicated(data),]#generated some reports twice
data <- data[,c(1:14,25,27)]#data[,c(1,4,25,2,3,10)]
colnames(data)[c(1,10,15)] <- c("Path","Fraction","Match")
data$Path <- substr(data$Path,52,nchar(data$Path)-4)
data$SampleDate <- gsub(".*(\\\\|/)(.*?)", "\\2", data$Path)
pad_func <- function(num) { formatC(as.integer(num), width=3, format="d", flag="0") }
data$Sample <- factor(ifelse(substr(data$SampleDate,1,5)=="BLANK", data$SampleDate, 
                      gsub("^KAHO", "SKAHO", str_replace(gsub("(.*?)_.*", "\\1", data$SampleDate), "(\\d+)", pad_func))))
data$SampleDate <- factor(data$SampleDate)
data$Path <- factor(data$Path)
data$Model <- factor(data$Model)
data$Fraction <- as.numeric(sub("%","",data$Fraction))/100
data$Purity <- as.numeric(sub("%","",data$Purity))/100
data$Min..Abund. <- as.numeric(sub("%","",data$Min..Abund.))/100
data$Name <- factor(sub(">","",data$Name))
data$Width <- as.numeric(sub(">","",data$Width))#important?
data$CAS <- factor(gsub("(\\d{2,7})(\\d{2})(\\d{1})$","\\1-\\2-\\3",data$CAS))#add dashes

#Multiply Fractions by total peak areas from sum_chromatograms.R
tics <- read.table("svhyb_TICsums.csv", header=TRUE, sep="\t", as.is=T)
tics$SampleDate <- gsub("(.*?).CDF", "\\1", tics$file)
tics$Sample <- factor(ifelse(substr(tics$SampleDate,1,5)=="BLANK", tics$SampleDate, gsub("^KAHO", "SKAHO", str_replace(gsub("(.*?)_.*", "\\1", tics$SampleDate), "(\\d+)", pad_func))))
data$Area <- data$Fraction * tics$sumTIC[match(data$Sample, tics$Sample)]
str(data)

#####Filtering######

  CASsubs <- read.delim("CASsubs3.csv", header=T, stringsAsFactors = F)
  data$CAS <- as.character(data$CAS)
   for(i in 1:nrow(CASsubs)){
     data$CAS[data$CAS==CASsubs[i,"CAS"]] <- CASsubs[i,"newCAS"]
   }
  data$CAS <- as.factor(data$CAS)

load("chemspc3.Rdata")
data$Name <- as.character(data$Name)
for(i in 1:nrow(chems.pc)){
  data$Name[data$CAS==as.character(chems.pc[i,"CAS"])] <- gsub("~\\{(.*?)\\}", "\\1", chems.pc[i,"IUPAC.Name"])
}

subs <- read.delim("subs3.csv", header=T, stringsAsFactors = F)
for(i in 1:nrow(subs)){
  data$Name[data$CAS==subs[i,"CAS"]] <- subs[i,"newName"]
}
data$Name <- as.factor(data$Name)

vol.all <- dcast(data, Sample~CAS, sum, value.var="Area")
rownames(vol.all) <- vol.all[,1]
vol.all[,1] <- NULL

#mds.all <- metaMDS(vol.all)
#ordiplot(mds.all, type="n")
#type.all <- as.factor(substr(as.character(levels(data$Sample)),1,4))
#text(mds.all, what="sites", col=as.integer(type.all), cex=0.5)

vol.df <- data.frame(Sample=rownames(vol.all))
vol.df$Samp <- substr(vol.df$Sample,0,1)=="S" & substr(vol.df$Sample,0,2)!="SI"
vol.df$Ambi <- substr(vol.df$Sample,0,1)=="A"
vol.df$Blank <-substr(vol.df$Sample,0,1)=="B"
vol.df$Index <- as.integer(substr(vol.df$Sample,nchar(as.character(vol.df$Sample))-2,nchar(as.character(vol.df$Sample))))
is.na(vol.df$Index) <- vol.df$Blank
vol.df$Samp <- vol.df$Samp & vol.df$Index %in% as.integer(s$Sample)[s$GC %in% c("G")]
specs <- c("KAAL", "KAHO", "HOOK")
vol.df$Focal <- substr(vol.df$Sample, 2,5) %in% specs
vol.df$Total <- s$Total[match(vol.df$Index,s$Sample)]
vol.df$Type <- factor(substr(vol.df$Sample,0,1))

chems <- data.frame(CAS=factor(colnames(vol.all)), Name=data$Name[match(colnames(vol.all),data$CAS)])
chems$RT <- sapply(chems$CAS, function(x) {median(data$RT[data$CAS==x])})
chems$Match <- sapply(chems$CAS, function(x) {median(data$Match[data$CAS==x])})
chems$RTV <- sapply(chems$CAS, function(x) {var(data$RT[data$CAS==x])})
chems$MatchV <- sapply(chems$CAS, function(x) {var(data$Match[data$CAS==x])})
chems$Siloxane <- grepl("sil", chems$Name, ignore.case=TRUE)
cont <- c("Benzene, 1-chloro-2-(trifluoromethyl)-", "Benzene, 1-chloro-3-(trifluoromethyl)-")
          #,"(Z)-hex-3-en-1-ol","oct-1-en-3-ol") #Chlorofluorcarbons!!
#c("Benzene","Caprolactam","2-Pentanone, 4-hydroxy-4-methyl-","2-Bromo dodecane","Ethanol, 2-butoxy-","Styrene","Ethanol, 2-(2-ethoxyethoxy)-","3-Penten-2-one, 4-methyl-","Nonanal","Decanal", "5-Hepten-2-one, 6-methyl-", "Homosalate", "Oxybenzone", "Toluene", "Disulfide, dimethyl","2-(oxiran-2-yl)oxirane","(E)-hex-3-en-1-ol", "(Z)-hex-3-en-1-ol")#some of these don't work with IUPAC names! also consider wound volatiles  "oct-1-en-3-ol"
chems$Cont <- chems$Name %in% cont
chems$SampC <- sapply(vol.all[vol.df$Samp,], function(x) {sum(x>0)})
chems$SampFC <- sapply(vol.all[vol.df$Focal,], function(x) {sum(x>0)})
chems$AmbiC <- sapply(vol.all[vol.df$Ambi,], function(x) {sum(x>0)})
chems$BlankC <- sapply(vol.all[vol.df$Blank,], function(x) {sum(x>0)})
chems$SampMax <- sapply(vol.all[vol.df$Samp,], max)
chems$SampFMax <- sapply(vol.all[vol.df$Focal,], max)
chems$SampMean <- sapply(vol.all[vol.df$Samp,], mean)
chems$SampFMean <- sapply(vol.all[vol.df$Focal,], mean)
chems$SampMedianPos <- sapply(vol.all[vol.df$Samp,], function(x) {median(x[x>0])})
chems$SampFMedianPos <- sapply(vol.all[vol.df$Focal,], function(x) {median(x[x>0])})
chems$SampP <- sapply(vol.all[vol.df$Samp,], function(x) {sum(x>0)/length(x)})
chems$SampFP <- sapply(vol.all[vol.df$Focal,], function(x) {sum(x>0)/length(x)})
chems$AmbiP <- sapply(vol.all[vol.df$Ambi,], function(x) {sum(x>0)/length(x)})
chems$BlankP <- sapply(vol.all[vol.df$Blank,], function(x) {sum(x>0)/length(x)})

chems$AmbiT <- sapply(chems$CAS, function(x) { t.test(vol.all[vol.df$Samp,x], vol.all[vol.df$Ambi,x], alternative="greater")$p.value })
chems$BlankT <- sapply(chems$CAS, function(x) { t.test(vol.all[vol.df$Samp,x], vol.all[vol.df$Blank,x], alternative="greater")$p.value })
chems$AmbiRatio <- sapply(chems$CAS, function(x) { 
  SampMean <- mean(vol.all[vol.df$Samp,x]);
  AmbiMean <- mean(vol.all[vol.df$Ambi,x]);
  return(ifelse(SampMean==0, 0, ifelse(AmbiMean==0, Inf, SampMean/AmbiMean)))  })
#chems$Ambi90Ratio <- sapply(chems$CAS, function(x) { 
#  SampMean <- quantile(vol.all[vol.df$Samp,x],probs=0.9);
#  AmbiMean <- quantile(vol.all[vol.df$Ambi,x],probs=0.9);
#  return(ifelse(SampMean==0, 0, ifelse(AmbiMean==0, Inf, SampMean/AmbiMean)))  })

chems$AmbiFT <- sapply(chems$CAS, function(x) { t.test(vol.all[vol.df$Focal,x], vol.all[vol.df$Ambi,x], alternative="greater")$p.value })
chems$BlankFT <- sapply(chems$CAS, function(x) { t.test(vol.all[vol.df$Focal,x], vol.all[vol.df$Blank,x], alternative="greater")$p.value })
chems$AmbiFRatio <- sapply(chems$CAS, function(x) { 
  SampFMean <- mean(vol.all[vol.df$Focal,x]);
  AmbiMean <- mean(vol.all[vol.df$Ambi,x]);
  return(ifelse(SampFMean==0, 0, ifelse(AmbiMean==0, Inf, SampFMean/AmbiMean)))  })

#criteria
#is.trueorna <- function(x) { is.na(x) | x }
RTMin <- 3
RTMax <- 18
MatchMin <- 75 #was 70
SampMaxMin <- 120000 #was 20000
SampCMin <- 1#(2 or above)
AmbiRatioMin <- 4 #was 3
Alpha <- 0.05
SampCMinNoAmbi <- 4
BlankMaxNoAmbi <- 4

chems$Real <- chems$Candidate <- with(chems,  RT > RTMin & RT < RTMax & Match > MatchMin & SampMax > SampMaxMin &  !Siloxane & !Cont & SampC > SampCMin)

#Filter the data set to examine compounds found at levels at least 3x the mean in FOCAL samples than ambients. Then apply false discovery rate to t test results based on that set of compounds.
chems$AmbiTadj <- 1
chems[chems$Real & chems$AmbiRatio > AmbiRatioMin,]$AmbiTadj <- with(chems, p.adjust(chems[Real & AmbiRatio > AmbiRatioMin,"AmbiT"], method="fdr"))
chems$BlankTadj <- 1
chems[chems$Real & chems$AmbiRatio > AmbiRatioMin,]$BlankTadj <- with(chems, p.adjust(chems[Real & AmbiRatio > AmbiRatioMin,"BlankT"], method="fdr"))
chems$Real <- with(chems, Real & ((AmbiRatio > AmbiRatioMin & AmbiTadj < Alpha & BlankTadj < Alpha)   #was AmbiTadj, BlankTadj !! 
                                  | (AmbiC == 0 & BlankC < BlankMaxNoAmbi & SampC > SampCMinNoAmbi)) )
chems$Real[chems$CAS%in%c("100-52-7","96-04-8")] <- TRUE #benzaldehyde Ambi90Ratio > 2.6 #heptan-2,3-dione - real in other species

chems$AmbiFTadj <- 1
chems[chems$Real & chems$AmbiFRatio > AmbiRatioMin,]$AmbiFTadj <- with(chems, p.adjust(chems[Real & AmbiFRatio > AmbiRatioMin,"AmbiFT"], method="fdr"))
chems$BlankFTadj <- 1
chems[chems$Real & chems$AmbiFRatio > AmbiRatioMin,]$BlankFTadj <- with(chems, p.adjust(chems[Real & AmbiFRatio > AmbiRatioMin,"BlankFT"], method="fdr"))
chems$FReal <- with(chems, Real & ((AmbiFRatio > AmbiRatioMin & AmbiFTadj < Alpha & BlankFTadj < Alpha) 
                                  | (AmbiC == 0 & SampFC > SampCMinNoAmbi)) )

######Flower volatiles table######
scales <- setNames(c("vol","fvol","favol","svol"), c("Emission/hr","Emission/hr/flower","Emission/hr/flower/area","Relative"))
vol <- vol.all[vol.df$Samp, chems$Real]#get RealLit from sfiltering.R
colnames(vol) <- chems$Name[chems$Real]
spec <- factor(substr(row.names(vol), 2,5))
species1 <-read.table("species.csv",header=TRUE,sep="\t",row.names=2) #load species data
species1 <- species1[levels(spec),]
species1$SpName <- factor(species1$SpName)
species <- species1[as.character(spec),]
sv <- as.data.frame(s) #load sample data
rownames(sv) <- sv$Filename
sv <- sv[rownames(vol),]
sv$Pop <- sv$Population
sv$Population <- mapvalues(sv$Population, from = c("794","866","891","899","879","881","892","904","3587"), to = c("WK","WK","WK","WK","879WKG","881KHV","892WKG","904WPG","3587WP"))
sv$Leaves <- as.integer(sv$Leaves)
sv$Inflor <- as.integer(sv$Inflor)
sv$Mphase <- as.integer(sv$Mphase)
sv$Fphase <- as.integer(sv$Fphase)
sv$Closed <- as.integer(sv$Closed)
sv$Buds <- as.integer(sv$Buds)
sv$Flrs <- as.integer(sv$Flrs)
#test effects of pumping/equil time
#summary(lme(sqrt(Abundance)~Total*as.integer(Flrs), random=~ 1|spec, data=sv))
#caprovol <- vol / vol$Caprolactam # / as.integer(sv$Flrs)[row(vol)]
vol <- vol / sv$Duration[row(vol)] #normalize by just pumping time !!!
svol <- as.data.frame(prop.table(as.matrix(vol),1)) #make samples sum to 1
ssvol <- as.data.frame(scale(svol,center=FALSE))#AND THEN scale each compound 0-1 so compounds are equally weighted
fvol <- vol / as.integer(sv$Flrs)[row(vol)] #normalize by number of flowers
sv$Abundance <- rowSums(fvol)
favol <- fvol / species$Sepal[row(fvol)]^2
pavol <- vol #convert to presence/absence
pavol[pavol>0] <-1

#use the the Schiedea Willyard Fig2 tree
schtree <- read.nexus("fig2.nex")
schtree$edge.length=schtree$edge.length/1000
sch<-multi2di(schtree, random=FALSE)
sch$edge.length<-sch$edge.length+0.00001
sch$tip.label <- toupper(paste("",substr(sch$tip.label,1,4),sep=""))
sch$tip.label[33] <- "SILEN"
#exclude hybrids to avoid tree conflicts
hybrids <- c("SALM","KAHO")
volnh <- vol[!spec %in% hybrids,]
fvolnh <- fvol[!spec %in% hybrids,]
pavolnh <- pavol[!spec %in% hybrids,]
speciesnh <- species[!spec %in% hybrids,]
species1nh <- species1[!rownames(species1) %in% hybrids,]
specnh <- factor(spec[!spec %in% hybrids])
svnh <- sv[!spec %in% hybrids,]
sch.p<-drop.tip(sch,sch$tip.label[-match(levels(specnh), sch$tip.label)])

#species averages - no hybrids
plantsi <- factor(paste(svnh$Population, svnh$Plant, svnh$Cutting, sep="-"))
volplants <- aggregate(x = volnh, by=list(plantsi),FUN=mean)
rownames(volplants) <-  volplants[,1]
volplants[,1] <- NULL
specplants <- aggregate(x = specnh, by=list(plantsi),FUN=unique)$x
vol1 <-aggregate(x = volplants, by=list(specplants),FUN=mean)
rownames(vol1) <-  vol1[,1]
vol1[,1] <- NULL
spec1 <- factor(rownames(vol1))
vol1 <- vol1[order(match(spec1,sch.p$tip.label)),]
fvol1 <-aggregate(x = fvolnh, by=list(specnh),FUN=mean)
rownames(fvol1) <-  fvol1[,1]
fvol1[,1] <- NULL
fvol1 <- fvol1[order(match(spec1,sch.p$tip.label)),]
svol1 <- as.data.frame(prop.table(as.matrix(vol1),1)) #make samples sum to 1
ssvol1 <- as.data.frame(scale(svol1, center=FALSE))
species1t <- species1nh[order(match(spec1,sch.p$tip.label)),] 
favol1 <- fvol1 / species1t$Sepal[row(fvol1)]^2
#convert to presence/absence
pavol1 <- vol1
pavol1[pavol1>0] <-1

schu<-multi2di(schtree, random=FALSE)
schu$tip.label <- toupper(paste("",substr(schu$tip.label,1,4),sep=""))
schu.p<-drop.tip(schu,schu$tip.label[-match(levels(specnh), schu$tip.label)])

set.seed(1)
scol <- sample(rainbow(nlevels(spec)))
set.seed(1)
ccol <- sample(rainbow(length(vol)))
#ccol[sapply(vol, max, na.rm=TRUE)<8*params$cutoff] <- "grey"
pcol <- setNames(c("green","blue","red"), c("Moth","Selfing","Wind"))

vnames <- paste(spec, sv$Population, sv$Sex, sep=" ")
save.image("sdata3_hexoct.Rdata")
