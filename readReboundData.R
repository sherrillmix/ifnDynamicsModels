acuteCut<-30
chronicCut<-300

dat<-read.csv('out/allLongitudinal.csv',stringsAsFactors=FALSE)
acuteChronicMM<-dat[dat$qvoa|(dat$time<acuteCut|dat$time>=chronicCut),]
acuteChronicMM$study<-'MM'
acuteChronicMM$virus<-acuteChronicMM$id
acuteChronicMM$ic50_IFNa2<-acuteChronicMM$ic50
acuteChronicMM$ic50_IFNb<-acuteChronicMM$beta
acuteChronicMM$class<-ifelse(acuteChronicMM$qvoa,'Outgrowth',ifelse(acuteChronicMM$time<=acuteCut,'Acute','Chronic'))
#acuteChronicMM$displayClass<-ifelse(acuteChronicMM$qvoa,'Outgrowth',ifelse(acuteChronicMM$isNadir,'Nadir',ifelse(acuteChronicMM$isLast,'Last',ifelse(acuteChronicMM$time<=acuteCut,'Acute',NA))))
acuteChronicMM$displayClass<-ifelse(acuteChronicMM$qvoa,'Outgrowth',ifelse(acuteChronicMM$time<=acuteCut,'Acute','Chronic'))
acuteChronicMM$virus<-acuteChronicMM$id
acuteChronicMM$repCap<-acuteChronicMM$replication
acuteChronicMM$infectivity<-acuteChronicMM$p24Release<-NA
#distinct "patient" for each time point (doesn't work because acute)
#acuteChronicMM$pat<-acuteChronicMM$sample

#voaMM<-dat[dat$qvoa,]
#voaMM$class<-'QVOA'
#voa$study<-'MM'
#voa$virus<-voa$id
#rebound2<-read.csv('data/Table S4.05.08.2020.csv',stringsAsFactors=FALSE,skip=2)
rebound2<-read.csv('data/TableS4_071420.csv',stringsAsFactors=FALSE,skip=2)
rebound2<-rebound2[!is.na(rebound2$Isolate.ID2)&rebound2$Isolate.ID2!='',]
rebound2<-rebound2[,!apply(is.na(rebound2),2,all)]
studyLookup<-c('MNU-0628'='RESERVOIR','NCT00051818'='INTERRUPT','NCT02227277'='IFNa2b treatment','NCT02463227'='VRC01','NCT02588586'='3BNC117','NCT02825797'='3BNC117/10-1074')
rebound2$study<-studyLookup[dnar::fillDown(rebound2$Study.number)]
rebound2$pat<-sub('\\..*','',rebound2$Isolate.ID2)
imc<-rebound2[rebound2$Type1=='IMC',]
rebound2<-rebound2[rebound2$Type1!='IMC',]
rebound2$displayClass<-rebound2$class<-ifelse(rebound2$Type1=='Rebound','Rebound',ifelse(grepl('Post',rebound2$Type1),'Post-ATI',ifelse(grepl('Pre|week',rebound2$Type1),'Pre-ATI','Outgrowth')))
rebound2$virus<-rebound2$Isolate.ID2
rebound2$ic50_IFNa2<-rebound2$IFNa2.IC50..pg.ml.5
rebound2$ic50_IFNb<-rebound2$IFNb.IC50..pg.ml.5
rebound2$repCap<-rebound2$Replicative.capacity.........ng.p24.ml.3
rebound2$infectivity<-NA #rebound2$Infectivity..IU.pg.RT.8
rebound2$p24Release<-NA #as.numeric(sub('%$','',ifelse(rebound2$p24.Particle.release....9=='',NA,rebound2$p24.Particle.release....9)))/100


rebound<-read.csv('data/Data Master 2020_ReboundandQVOA_20200504.csv',stringsAsFactors=FALSE)
rebound<-rebound[!grepl('_BE$',rebound$ID),]
table(rebound[!is.na(rebound$ic50_IFNa2),'Study'])
rebound$study<-sub(' / ','/',rebound$Study)
rebound$displayClass<-rebound$class<-ifelse(rebound$Type=='Rebound','Rebound',ifelse(grepl('Post',rebound$Type),'Post-ATI',ifelse(grepl('Pre',rebound$Type),'Pre-ATI','Outgrowth')))
rebound$pat<-sub(' \\(R-?[0-9]+\\)','',rebound$Patient)
rebound<-rebound[!is.na(rebound$Type),]
rebound[rebound$study=='OUTGROWTH','study']<-'MM'
rebound$virus<-rebound$ID
#table(comboA$study)

pair<-read.csv('rebound/donorRecipient.csv',stringsAsFactors=FALSE)
pair$class<-ifelse(pair$donor,'Donor','Recipient')
pair$type<-'CHAVI cohort'
#pair$label<-sprintf('CHAVI %s',pair$class)
pair$displayClass<-sprintf('%s %ss',ifelse(pair$class=='Donor','Chronic','Acute'),pair$class)
pair$ic50_IFNa2<-pair$IFNa2.Pooled.Donor.cells.IC50..pg..ml
pair$ic50_IFNb<-pair$IFNbeta.Pooled.Donor.cells.IC50..pg.ml
#pair[,'class']<-c('Recipient'='Acute Recipient','Donor'='Chronic Donor')[pair$class]
pair$pat<-pair$sample
pair$repCap<-pair$Replicative.capacity.Pooled.Donor.cells.p24.d7
pair[,colnames(dat)[!colnames(dat) %in% colnames(pair)]]<-NA
pair$virus<-pair$Renamed
pair$source<-'shilpa'
pair$study<-'Transmission'
pair$repCap<-pair$infectivity<-pair$p24Release<-NA

#pat,simpleClass,study,ic50,class
targetCols<-c('pat','class','study','ic50_IFNa2','ic50_IFNb','displayClass','virus','repCap','infectivity','p24Release')
combined<-rbind(acuteChronicMM[,targetCols],rebound2[,targetCols],pair[,targetCols])
combined$simpleClass<-ifelse(combined$class %in% c('Chronic Donor','Donor','1 Year','Nadir','Last'),'Chronic',combined$class)
combined$simpleClass[combined$simpleClass=='Recipient']<-'Acute'
combined$simpleClass[combined$simpleClass=='Pre-ATI']<-'Outgrowth' #may want own class
#combined$simpleClass[combined$simpleClass=='Pre-ATI']<-'QVOA' #may want own class
combined$simplePat<-ifelse(combined$simpleClass %in% c('Acute','Chronic','Donor','Recipient'),'',combined$pat)
combined$label<-ifelse(combined$simplePat=='',combined$displayClass,sprintf('%s%s%s',combined$displayClass,ifelse(combined$simplePat!='',' ',''),combined$simplePat))
combined$voaVsRebound<-ifelse(combined$class %in% c('Post-ATI','Pre-ATI','Outgrowth'),'Outgrowth',ifelse(combined$class %in% c('Rebound'),'Rebound',NA))
combined<-combined[!is.na(combined$ic50_IFNb)|!is.na(combined$ic50_IFNa2),]

fastRegex<-'MM15|WEAU'
slowRegex<-'MM55|MM62'
standardRegex<-'MM14|MM23|MM33|MM34|MM39|MM40'
combined$speed<-ifelse(grepl(fastRegex,combined$virus),'Fast',ifelse(grepl(slowRegex,combined$virus),'Slow',ifelse(grepl(standardRegex,combined$virus),'Standard','Other')))


if(FALSE){
zz<-table(comboA$pat)
zz2<-table(combined$pat[!is.na(combined$ic50_IFNa2)])
zz<-table(comboB$pat)
zz2<-table(combined$pat[!is.na(combined$ic50_IFNb)])
all(names(zz2) == names(zz))
cbind(zz,zz2)
all(zz2>=zz)
}
