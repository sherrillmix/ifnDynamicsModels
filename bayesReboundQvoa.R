library('rstan')
library(dnar)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source('functions.R')
source('rebound/functions.R')

stanCode4_withMixture<-'
  data {
    int<lower=0> nVirus;
    int<lower=0> nPatient;
    int<lower=0> nState;
    int<lower=0> nStudy;
    real ic50[nVirus];
    int<lower=1,upper=nPatient> patients[nVirus];
    int<lower=1,upper=nState+1> states[nVirus];
    int<lower=1,upper=nStudy> studies[nPatient];
    int<lower=0,upper=1> isTreat[nPatient];
  }
  parameters {
    matrix[nPatient,nState] baseIc50Raw;
    vector[nState] stateMeans;
    vector[nStudy-1] studyMeansRaw;
    vector<lower=0>[nState] stateSds;
    vector<lower=0>[nState] stateIsoSds;
    real<lower=0,upper=1> postProp;
    real<lower=0,upper=1> preProp;
    real treatMean;
  }
  transformed parameters{
    matrix[nPatient,nState] expectedIC50;
    vector[nStudy] studyMeans;
    studyMeans[1]=0;
    studyMeans[2:nStudy]=studyMeansRaw;
    for(ii in 1:nPatient){
      for(jj in 1:(nState)){
        expectedIC50[ii,jj]=studyMeans[studies[ii]]+stateMeans[1]+baseIc50Raw[ii,1]*stateSds[1];
        if(jj>1)expectedIC50[ii,jj]=expectedIC50[ii,jj]+stateMeans[jj]+baseIc50Raw[ii,jj]*stateSds[jj];
        if(jj==2 && isTreat[ii])expectedIC50[ii,jj]=expectedIC50[ii,jj]+treatMean;
      }
    }
  }
  model {
    for(ii in 1:nVirus){
      if(states[ii]<nState)ic50[ii]~normal(expectedIC50[patients[ii],states[ii]],stateIsoSds[states[ii]]);
      if(states[ii]==nState+1){
        target += log_sum_exp(
          log(postProp)+normal_lpdf(ic50[ii]|expectedIC50[patients[ii],2],stateIsoSds[2]),
          log(1-postProp)+normal_lpdf(ic50[ii]|expectedIC50[patients[ii],nState],stateIsoSds[nState])
        );
      }
      if(states[ii]==nState){
        target += log_sum_exp(
          log(preProp)+normal_lpdf(ic50[ii]|expectedIC50[patients[ii],2],stateIsoSds[2]),
          log(1-preProp)+normal_lpdf(ic50[ii]|expectedIC50[patients[ii],nState],stateIsoSds[nState])
        );
      }
    }
    postProp~beta(1,1);
    preProp~beta(1,1);
    stateIsoSds~gamma(1,1);
    stateSds~gamma(1,1);
    for(ii in 1:nState){
      baseIc50Raw[,ii]~normal(0,1);
    }
    treatMean~normal(0,10);
    studyMeansRaw~normal(0,10);
  }
'
mod4 <- stan_model(model_code = stanCode4_withMixture)


source('readReboundData.R')
if(!exists('fitA_withMix')){
fitBayes2<-function(model,patient,states,studies,treats,ic50,baseState='Acute',mixState='Post-ATI',chains=50,baseStudy='Other',logFunc=log,stateId=NULL,iter=6000,...){
  patientId<-structure(1:length(unique(patient)),.Names=sort(unique(patient)))
  studyId<-structure(1:length(unique(studies)),.Names=unique(studies[order(studies!=baseStudy)]))
  if(is.null(stateId))stateId=structure(1:length(unique(states)),.Names=unique(states[order(states!=baseState,states==mixState)]))
  if(any(!states %in% names(stateId)))stop('Missing state IDs')
  patientStudy<-tapply(studies,patient,unique)[names(patientId)]
  isTreat<-tapply(treats,patient,unique)[names(patientId)]
  if(any(sapply(patientStudy,length)!=1))stop('Multiple studies for a single patient')
  dat=list(
    nVirus=length(ic50),
    nPatient=max(patientId),
    nState=max(stateId)-1,
    nStudy=max(studyId),
    ic50=logFunc(ic50),
    patients=patientId[patient],
    states=stateId[states],
    studies=studyId[patientStudy[names(patientId)]],
    isTreat=isTreat
  )
  fit <- sampling(model, data = dat, iter=iter, chains=chains,thin=2,control=list(adapt_delta=.99,max_treedepth=15),...)
  return(list('fit'=fit,pats=patientId,states=stateId,studies=studyId,patientStudy=patientStudy,dat=dat))
}
fitA_withMix<-withAs(combined=combined[!is.na(combined$ic50_IFNa2),],fitBayes2(mod4,combined$pat,combined$simpleClass,ifelse(combined$study %in% c('Transmission'),combined$study,'Other'),combined$study %in% c('BEAT','IFNa2b treatment'),combined$ic50_IFNa2,chains=50,stateId=structure(1:5,.Names=c('Acute','Rebound','Chronic','Outgrowth','Post-ATI')),iter=5000))
fitB_withMix<-withAs(combined=combined[!is.na(combined$ic50_IFNb),],fitBayes2(mod4,combined$pat,combined$simpleClass,ifelse(combined$study %in% c('Transmission'),combined$study,'Other'),combined$study %in% c('BEAT','IFNa2b treatment'),combined$ic50_IFNb,chains=50,stateId=structure(1:5,.Names=c('Acute','Rebound','Chronic','Outgrowth','Post-ATI')),iter=5000))
alphaAdjust<-exp(mean(as.matrix(fitA_withMix$fit)[,sprintf('studyMeans[%d]',fitA_withMix$studies['Transmission'])]))
betaAdjust<-exp(mean(as.matrix(fitB_withMix$fit)[,sprintf('studyMeans[%d]',fitB_withMix$studies['Transmission'])]))
#save(fitA_withMix,fitB_withMix,file='work/mixFits_2020-05-14.Rdat')
#load('work/mixFits_2020-05-14.Rdat')
#save(fitA_withMix,fitB_withMix,file='work/mixFitsUpdate_2020-07-20.Rdat')

}

#pdf('test.pdf');plotSummary(fitA_withMix,addAcute=FALSE);plotSummary(fitB_withMix,addAcute=FALSE);dev.off()


zz<-as.matrix(fitB_withMix$fit)
meanCrI<-function(xx)c(quantile(xx,c(.025,.975)),mean(xx))[c(1,3,2)]
exp(-apply(zz[,grepl('stateMeans\\[4\\]',colnames(zz)),drop=FALSE],2,meanCrI))
exp(-apply(zz[,'postProp',drop=FALSE]*zz[,'stateMeans[2]']+(1-zz[,'postProp',drop=FALSE])*zz[,'stateMeans[4]'],2,meanCrI))
exp(-apply(zz[,'preProp',drop=FALSE]*zz[,'stateMeans[2]']+(1-zz[,'preProp',drop=FALSE])*zz[,'stateMeans[4]'],2,meanCrI))
mean(zz[,'preProp',drop=FALSE]*zz[,'stateMeans[2]']+(1-zz[,'preProp',drop=FALSE])*zz[,'stateMeans[4]']>zz[,'postProp',drop=FALSE]*zz[,'stateMeans[2]']+(1-zz[,'postProp',drop=FALSE])*zz[,'stateMeans[4]'])


#ordering<-c('Pre-ATI A06','Post-ATI A06','Pre-ATI A08','Rebound A08','Post-ATI A08','Pre-ATI A09','Rebound A09','Post-ATI A09','Pre-ATI 9201','Rebound 9201','Pre-ATI 9202','Rebound 9202','Pre-ATI 9203','Rebound 9203','Pre-ATI 9207','Rebound 9207','Rebound A08','Rebound S-22','Rebound S-23','Rebound S-30','Rebound 601','Rebound BEAT-004','Rebound BEAT-030','Rebound BEAT-044','Outgrowth B106','Outgrowth B199','Outgrowth MM14','Outgrowth MM15','Outgrowth MM23','Outgrowth MM34','Outgrowth MM40','Outgrowth MM34','Acute','1 Year','Nadir','Last','Acute Recipients','Chronic Donors')
ordering<-c('Pre-ATI A06','Post-ATI A06','Pre-ATI A08','Rebound A08','Post-ATI A08','Pre-ATI A09','Rebound A09','Post-ATI A09','Pre-ATI 9241','Rebound 9241','Pre-ATI 9242','Rebound 9242','Pre-ATI 9243','Rebound 9243','Pre-ATI 9244','Rebound 9244','Rebound A08','Rebound S22','Rebound S23','Rebound S30','Rebound 601','Rebound 004','Rebound 030','Rebound 044','Outgrowth B106','Outgrowth B199','Outgrowth MM14','Outgrowth MM15','Outgrowth MM23','Outgrowth MM34','Outgrowth MM40','Outgrowth MM34','Acute','Chronic','Acute Recipients','Chronic Donors')
pos<-structure(1:length(unique(combined$label[!is.na(combined$label)])),.Names=unique(combined$label[!is.na(combined$label)][orderIn(combined$label[!is.na(combined$label)],ordering)]))
posStudy<-sapply(names(pos),function(xx)combined[combined$label==xx&!is.na(combined$label),'study'][1])
posPat<-sapply(names(pos),function(xx)combined[combined$label==xx&!is.na(combined$label),'pat'][1])
posPat['Chronic']<-'XXX'
studySpace<-.5
patSpace<--.25
pos<-pos+cumsum(c(0,posStudy[-length(posStudy)]!=posStudy[-1]))*studySpace
pos<-pos+cumsum(c(0,posPat[-length(posPat)]==posPat[-1]))*patSpace
plotSummary<-function(fit,ylab='IFNa2 IC50 (pg/ml)',mar=c(6.9,4,.1,.1),xWidth=.4,addAcute=TRUE,xaxis=TRUE,logYAxis=TRUE,reps=2,cols=structure(rep('#00000033',nStates),.Names=stateNames),combine24=FALSE,cols2=cols,quantRange=c(.025,.975),subtract12=FALSE){
  if(combine24){
    extraState<-tail(names(fit$states),1)
    fit$states<-fit$states[-length(fit$states)]
  }
  mat<-as.matrix(fit$fit)
  states<-mat[,sprintf('stateMeans[%d]',fit$states)]
  stateSds<-mat[,sprintf('stateSds[%d]',fit$states)]
  stateIsoSds<-mat[,sprintf('stateIsoSds[%d]',fit$states)]
  nStates<-ncol(states)
  if(addAcute){
    states[,2:nStates]<-states[,2:nStates]+states[,1]
    stateSds[,2:nStates]<-sqrt(stateSds[,2:nStates]^2+stateSds[,1]^2)
    stateIsoSds[,2:nStates]<-sqrt(stateSds[,2:nStates]^2+stateSds[,1]^2)
    ylim<-range(exp(fit$dat$ic50))
    stateNames<-names(fit$states)
  }else if(subtract12){
    stateNames<-names(fit$states)[-2]
  }else{
    states<-states[,-1,drop=FALSE]
    stateSds<-stateSds[,-1,drop=FALSE]
    stateIsoSds<-stateIsoSds[,-1,drop=FALSE]
    nStates<-nStates-1
    stateNames<-names(fit$states[-1,drop=FALSE])
  }
  stateNames[stateNames=='QVOA']<-'Outgrowth'
  statesPat<-do.call(cbind,lapply(1:nStates,function(ii)rnorm(nrow(states)*reps,states[,ii],stateSds[,ii])))
  statesIso<-do.call(cbind,lapply(1:nStates,function(ii)rnorm(nrow(states)*reps,states[,ii],sqrt(stateSds[,ii]^2+stateIsoSds[,ii]^2))))
  if(combine24){
    postProp<-mat[,'postProp']
    preProp<-mat[,'preProp']
    propRand<-runif(nrow(states)*reps,0,1)
    propRand2<-runif(nrow(states)*reps,0,1)
    states4<-states[,3]
    statesPat4<-statesPat[,3]
    statesIso4<-statesIso[,3]
    states[,3]<-ifelse(runif(nrow(states),0,1)<preProp,states[,1],states4)
    statesPat[,3]<-ifelse(propRand<preProp,statesPat[,1],statesPat4)
    statesIso[,3]<-ifelse(propRand<preProp,statesIso[,1],statesIso4)
    states<-cbind(states,states[,3])
    statesPat<-cbind(statesPat,statesPat[,3])
    statesIso<-cbind(statesPat,statesIso[,3])
    states[,4]<-ifelse(runif(nrow(states),0,1)<postProp,states[,1],states4)
    statesPat[,4]<-ifelse(propRand<postProp,statesPat[,1],statesPat4)
    statesIso[,4]<-ifelse(propRand<postProp,statesIso[,1],statesIso4)
    nStates<-4
    stateNames<-c(stateNames,extraState)
  }
  if(subtract12){
    states[,1]<-states[,2]-states[,1]
    statesPat[,1]<-statesPat[,2]-statesPat[,1]
    statesIso[,1]<-statesIso[,2]-statesIso[,1]
    colName<-colnames(states)[2]
    states<-states[,-2,drop=FALSE]
    statesPat<-statesPat[,-2,drop=FALSE]
    statesIso<-statesIso[,-2,drop=FALSE]
    nStates<-nStates-1
  }
  if(!addAcute)ylim<-exp(range(c(-1,1,apply(statesIso,2,quantile,c(quantRange[1],quantRange[2])))))
  makeDense<-function(xx)density(xx,from=quantile(xx,quantRange[1]),to=quantile(xx,quantRange[2]))
  denses<-apply(states,2,makeDense)
  denses2<-apply(statesPat,2,makeDense)
  denses3<-apply(statesIso,2,makeDense)
  #denses2<-lapply(1:ncol(expects),function(ii)makeDense(rnorm(n*10,expects[,ii],sds[,fit$dat$states[ii]])))
  par(mar=mar)
  plot(1,1,log='y',yaxt='n',ylab=ylab,xlab='',xaxt='n',type='n',cex.lab=1.2,ylim=ylim,xlim=c(.5,nStates+.5),las=1,mgp=c(2.45,1,0))
  if(logYAxis)logAxis(las=1,mgp=c(3,.6,0))
  else axis(2,c(.2,.5,1,2,5),c('0.2','0.5','1','2','5'),las=1)
  for(ii in 1:nStates){
    polygon(ii+c(denses[[ii]]$y,-rev(denses[[ii]]$y))/max(denses[[ii]]$y)*xWidth,exp(c(denses[[ii]]$x,rev(denses[[ii]]$x))),col=cols[stateNames[ii]])
    #polygon(ii+c(denses2[[ii]]$y,-rev(denses2[[ii]]$y))/max(denses2[[ii]]$y)*xWidth,exp(c(denses2[[ii]]$x,rev(denses2[[ii]]$x))),col=cols[stateNames[ii]])
    polygon(ii+c(denses3[[ii]]$y,-rev(denses3[[ii]]$y))/max(denses3[[ii]]$y)*xWidth,exp(c(denses3[[ii]]$x,rev(denses3[[ii]]$x))),col=cols2[stateNames[ii]])
  }
  if(xaxis)for(ii in 1:nStates)axis(1,ii,stateNames[ii])
  if(!addAcute)abline(h=1,lty=2)
  return(stateNames)
}
pdf('out/voaRebound_bayesSummary.pdf',height=6,width=7)
  layout(rbind(1:2,0,3:4),width=c(8,2.2),height=c(1,.01,1))
  #classCols<-structure(c(rep("#9EC0E1E6",3),rep("#77BCA9B3",2), rep("#9FB755B3",2), "#84C47DB3", "#B99A4BB3", "#C77C62B3", "#E581A0E6"), .Names = c("Outgrowth","Pre-ATI","Post-ATI","Acute Recipients","Acute", "Chronic Donors", "Chronic","1 Year", "Nadir", "Last", "Rebound")) 
  classCols<-structure(c(rep(classCol['qvoa'],2),classCol['postQvoa'],rep("#aaaaaa",2), rep("#DDDDDD",2), "#aaaaaa", "#aaaaaa", "#aaaaaa", classCol['rebound'],patCols), .Names = c("Outgrowth","Pre-ATI","Post-ATI","Acute Recipients","Acute", "Chronic Donors", "Chronic","1 Year", "Nadir", "Last", "Rebound",names(patCols))) 
  stCols<-sprintf('%s99',substring(classCols,1,7))
  stCols2<-sprintf('%s33',substring(classCols,1,7))
  names(stCols)<-names(stCols2)<-names(classCols)
  names(stCols2)[names(stCols2)=='Donor']<-names(stCols)[names(stCols)=='Donor']<-'Chronic'
  plotFunc<-function(combined,ic50Col,fit,ylab='IFNa2 IC50 (pg/ml)',letters=LETTERS[1:2]){
    cex.axis<-1
    out<-withAs(combined=combined[!is.na(combined[,ic50Col])&!is.na(combined$label),],
      plotQvoa2(combined[,ic50Col],combined$label,pos,ifelse(combined$displayClass %in% c(),combined$pat,combined$displayClass),combined$study,combined$speed,ylab=ylab,mar=c(5.5,3.7,.1,.1),cex.axis=cex.axis,startDown=TRUE,pats=ifelse(combined$study %in% c('Transmission','MM'),NA,combined$pat),classCols=classCols,labelXAxis=FALSE)
    )
    patPos<-tapply(out$pos,sub('BEAT-','',sub('.* ','',names(out$pos))),mean)
    patPos<-patPos[!grepl('^Acute|^Month|^Recipient|^Donor|^Nadir|^Year|^Last|^Chronic',names(patPos))]
    slantAxis(1,patPos,names(patPos),cex=cex.axis,location=.8)
    text(grconvertX(par('fig')[1]+diff(par('fig')[1:2])*.0025,from='ndc'),grconvertY(par('fig')[4]-diff(par('fig')[3:4])*.005,from='ndc'),letters[1],xpd=NA,adj=c(0,1),cex=2)
    #plotSummary(fitA)
    states<-plotSummary(fit,addAcute=FALSE,ylab='Fold change from acute',mar=c(4.5,4,.1,1.35),cols=stCols,cols2=stCols2,xaxis=FALSE,combine24=TRUE)
    slantAxis(1,1:length(states),sub('Outgrowth','Pre-ATI',states),textOffsets=c(-.2,-.2,-.2,-.2),location=.7,axisArgs=list(tcl=-.4),srt=-45,cex=cex.axis)
    text(grconvertX(par('fig')[1]+diff(par('fig')[1:2])*.03,from='ndc'),grconvertY(par('fig')[4]-diff(par('fig')[3:4])*.005,from='ndc'),letters[2],xpd=NA,adj=c(0,1),cex=2)
  }
  tmp<-combined
  tmp$alphaAdjust<-tmp$ic50_IFNa2/ifelse(tmp$study=='Transmission',alphaAdjust,1)
  tmp$betaAdjust<-tmp$ic50_IFNb/ifelse(tmp$study=='Transmission',betaAdjust,1)
  plotFunc(tmp,'alphaAdjust',fitA_withMix,ylab=expression('IFN'*alpha*'2 IC'[50]*' (pg/ml)'))
  plotFunc(tmp,'betaAdjust',fitB_withMix,ylab=expression('IFN'*beta*' IC'[50]*' (pg/ml)'),letters=LETTERS[3:4])
  #remove adjustment based on single IC50 and use bayesian estimated
  #combined$ic50_IFNb[combo$study=='Transmission']<-combo$beta[combo$study=='Transmission']/6386*2230
  #plotSummary(fitB,ylab='IFNb IC50 (pg/ml)')
dev.off()
file.copy('out/voaRebound_bayesSummary.pdf','out/Fig._4.pdf',overwrite=TRUE)
#system('pdfjam out/voaRebound_bayesSummary.pdf --nup 1x2 --outfile tmp.pdf;pdfcrop tmp.pdf out/Fig._4.pdf')



mat<-as.matrix(fitA_withMix$fit)
matB<-as.matrix(fitB_withMix$fit)
message('Chronic')
print(meanCrI(exp(-mat[,sprintf('stateMeans[%d]',fitA_withMix$states['Chronic'])])))
print(meanCrI(exp(-matB[,sprintf('stateMeans[%d]',fitB_withMix$states['Chronic'])])))
message('Rebound')
print(meanCrI(exp(mat[,sprintf('stateMeans[%d]',fitA_withMix$states['Rebound'])])))
print(meanCrI(exp(matB[,sprintf('stateMeans[%d]',fitB_withMix$states['Rebound'])])))
message('Outgrowth (without mixture)')
print(meanCrI(exp(mat[,sprintf('stateMeans[%d]',fitA_withMix$states['Outgrowth'])])))
print(meanCrI(exp(matB[,sprintf('stateMeans[%d]',fitB_withMix$states['Outgrowth'])])))
message('Postprop (with mixture)')
print(meanCrI(exp(-(mat[,'postProp']*mat[,sprintf('stateMeans[%d]',fitA_withMix$states['Rebound'])]+(1-mat[,'postProp'])*mat[,sprintf('stateMeans[%d]',fitA_withMix$states['Outgrowth'])]))))
print(meanCrI(exp(-(matB[,'postProp']*matB[,sprintf('stateMeans[%d]',fitB_withMix$states['Rebound'])]+(1-matB[,'postProp'])*matB[,sprintf('stateMeans[%d]',fitB_withMix$states['Outgrowth'])]))))
message('Preprop (with mixture)')
print(meanCrI(exp(-(mat[,'preProp']*mat[,sprintf('stateMeans[%d]',fitA_withMix$states['Rebound'])]+(1-mat[,'preProp'])*mat[,sprintf('stateMeans[%d]',fitA_withMix$states['Outgrowth'])]))))
print(meanCrI(exp(-(matB[,'preProp']*matB[,sprintf('stateMeans[%d]',fitB_withMix$states['Rebound'])]+(1-matB[,'preProp'])*matB[,sprintf('stateMeans[%d]',fitB_withMix$states['Outgrowth'])]))))
message('Treatment effect')
print(meanCrI(exp(mat[,'treatMean'])))
print(meanCrI(exp(matB[,'treatMean'])))
message('Post prop')
print(meanCrI(mat[,'postProp']))
print(meanCrI(matB[,'postProp']))
message('Pre prop')
print(meanCrI(mat[,'preProp']))
print(meanCrI(matB[,'preProp']))
message('Pre vs post prop')
print(meanCrI(mat[,'preProp']>mat[,'postProp']))
print(meanCrI(matB[,'preProp']>matB[,'postProp']))

mean(mat[,'preProp',drop=FALSE]*mat[,'stateMeans[2]']+(1-mat[,'preProp',drop=FALSE])*mat[,'stateMeans[4]']>mat[,'postProp',drop=FALSE]*mat[,'stateMeans[2]']+(1-mat[,'postProp',drop=FALSE])*mat[,'stateMeans[4]'])
mean(matB[,'preProp',drop=FALSE]*matB[,'stateMeans[2]']+(1-matB[,'preProp',drop=FALSE])*matB[,'stateMeans[4]']>matB[,'postProp',drop=FALSE]*matB[,'stateMeans[2]']+(1-matB[,'postProp',drop=FALSE])*matB[,'stateMeans[4]'])
