library('rstan')
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

dat<-read.csv('out/allLongitudinal.csv',stringsAsFactors=FALSE)
meta<-read.csv('out/allLongitudinalMeta.csv')
maxDates<-tapply(meta$time,meta$mm,max)
artStarts<-tapply(meta$daysBeforeArt+meta$time,meta$mm,unique)
founders<-read.csv('out/founders.csv',stringsAsFactors=FALSE,row.names=1)
superTimes<-structure(founders$superTime,.Names=rownames(founders))

if(!exists('fit')){
ic50_bp_mar<-'
  data {
    int<lower=0> nVirus;
    int<lower=0> nPatient;
    int<lower=0> nArt;
    real cd4[nVirus];
    real vl[nVirus];
    real ic50[nVirus];
    int<lower=0,upper=nPatient> patients[nVirus];
    vector<lower=0,upper=1>[nPatient] isFast;
    vector<lower=0,upper=1>[nPatient] isNon;
    int<lower=0> artIds[nVirus];
    int<lower=0> days[nVirus]; real weeks[nVirus];
    //real daysBeforeArt[nVirus];
    real artStart[nArt];
    real<lower=0> hasArt[nVirus];
    int<lower=0> nSample;
    int<lower=0,upper=nSample> sample[nVirus];
    int<lower=0> tMax;
    matrix[nPatient,tMax] vlMat;
    matrix[nPatient,tMax] cd4Mat;
  }
  parameters {
    vector[nPatient] acuteRaw;
    real acuteMean;
    real<lower=0> acuteSD;
    //vector[nPatient] nadirTimeRaw;
    real<lower=0,upper=tMax> nadirTimeMean;
    real nadirTimeFast;
    real nadirTimeNon;
    real<lower=0> nadirTimePhi;
    vector[nPatient] nadirChangeRaw;
    real nadirChangeMean;
    real<lower=0> nadirChangeSD;
    real<lower=0> sigma;
    vector[nPatient] cd4BetaRaw;
    real cd4BetaMean[3];
    real<lower=0> cd4BetaSD;
    real fastAcute;
    real fastChangeMean;
    real nonAcute;
    real nonChangeMean;
    //vector[nPatient] vlBetaRaw;
    //real vlBetaMean[2];
    //real<lower=0> vlBetaSD;
  }
  transformed parameters{
    vector[nPatient] acute;
    //real expectedIC50[nVirus];
    //vector<lower=0,upper=tMax>[nPatient] nadirTime;
    vector[nPatient] nadirChange;
    vector[nPatient] cd4Beta;
    //vector[nPatient] vlBeta;
    matrix[nPatient,tMax] lp;
    //acute=acuteMean+acuteRaw*acuteSD;
    acute=acuteMean+acuteRaw*acuteSD+isNon*nonAcute+isFast*fastAcute;
    nadirChange=nadirChangeMean+nadirChangeRaw*nadirChangeSD+fastChangeMean*isFast+nonChangeMean*isNon;
    //cd4Beta=(cd4BetaMean[1]*(1-isFast)+cd4BetaMean[2]*isFast)+cd4BetaRaw*cd4BetaSD;
    //cd4Beta=(cd4BetaMean[1]+cd4BetaMean[2]*isFast+cd4BetaMean[3]*isNon)+cd4BetaRaw*cd4BetaSD;
    cd4Beta=(cd4BetaMean[1]*(1-isFast).*(1-isNon)+cd4BetaMean[2]*isFast+cd4BetaMean[3]*isNon)+cd4BetaRaw*cd4BetaSD;
    for (s in 1:tMax){
      //lp[,s]=rep_vector(neg_binomial_2_lpmf(s|nadirTimeMean,nadirTimePhi),nPatient);
      for(jj in 1:nPatient)lp[jj,s]=neg_binomial_2_lpmf(s|nadirTimeMean*exp(isFast[jj]*nadirTimeFast)*exp(isNon[jj]*nadirTimeNon),nadirTimePhi);
      //lp[,s]=rep_vector(log(normal_cdf(s+.5,nadirTimeMean,nadirTimeSD)-normal_cdf(s-.5,nadirTimeMean,nadirTimeSD)),nPatient);
      for (ii in 1:nVirus){
        lp[patients[ii],s] = lp[patients[ii],s] + normal_lpdf(ic50[ii] | weeks[ii] < s ? 
          acute[patients[ii]]+nadirChange[patients[ii]]*weeks[ii]/s  : 
          acute[patients[ii]]+nadirChange[patients[ii]] +(cd4Beta[patients[ii]])*(cd4[ii]-cd4Mat[patients[ii],s])
          ,sigma);
      }
    }
  }
  model {
    for(ii in 1:nPatient)target += log_sum_exp(lp[ii,]);
    //ic50~normal(expectedIC50,sigma);
    sigma~gamma(1,.1);
    nadirTimePhi~cauchy(0,10);//gamma(1,.01);
    //nadirTimeRaw~normal(0,1);
    nadirChangeSD~gamma(1,.1);
    nadirChangeRaw~normal(0,1);
    acuteSD~gamma(1,.1);
    acuteRaw~normal(0,1);
    nadirChangeMean~normal(0,10);
    fastChangeMean~normal(0,10);
    fastAcute~normal(0,10);
    cd4BetaMean~normal(0,10);
    nonChangeMean~normal(0,10);
    nonAcute~normal(0,10);
    cd4BetaSD~gamma(1,.1);
    cd4BetaRaw~normal(0,1);
    nadirTimeFast~normal(0,10);
    nadirTimeNon~normal(0,10);
    //vlBetaMean~normal(0,10);
    //vlBetaSD~gamma(1,.1);
    //vlBetaRaw~normal(0,1);
  }\n
'
ic50Mod <- stan_model(model_code = ic50_bp_mar)

bayesIC50_3<-function(mod,ic50,time,timePreArt,patient,chains=50,nIter=30000,fastProgressors=c(),cd4=c(),vl=c(),baseVl=c(),meta=NULL,superTimes=NULL,nonProgressors=FALSE,...){
  #treating all postART   as one group (currently only WEAU)
  patientId<-structure(1:length(unique(patient)),.Names=sort(unique(patient)))
  artId<-structure(1:length(unique(patient[!is.na(timePreArt)])),.Names=sort(unique(patient[!is.na(timePreArt)])))
  artStart<-tapply(time+timePreArt,patient,unique)
  sample<-paste(patient,time,sep='_')
  sampleId<-structure(1:length(unique(sample)),.Names=unique(sample[order(patient,time)]))
  ####patient-week matrix of inferred CD4 and VL
  weekMax<-tapply(ceiling(time/7),patient,max)+1
  newDat<-do.call(rbind,mapply(function(xx,yy){data.frame('pat'=yy,'week'=min(round(time/7)):xx,stringsAsFactors=FALSE)},weekMax,names(weekMax),SIMPLIFY=FALSE))
  newDat$day<-newDat$week*7-3.5
  newDat$vl<-newDat$cd4<-NA
  for(ii in unique(newDat$pat)){
    #thisDat<-meta[meta$mm==ii&meta$time<=weekMax[ii]*7,]
    #newDat[newDat$pat==ii,'vl']<-approx(thisDat$time,thisDat$vl,newDat[newDat$pat==ii,'day'],rule=2)$y
    #newDat[newDat$pat==ii,'cd4']<-approx(thisDat$time,thisDat$cd4,newDat[newDat$pat==ii,'day'],rule=2)$y
    thisDat<-unique(dat[dat$pat==ii&dat$time<=weekMax[ii]*7,c('time','CD4','vl')])
    newDat[newDat$pat==ii,'vl']<-approx(thisDat$time,thisDat$vl,newDat[newDat$pat==ii,'day'],rule=2)$y
    newDat[newDat$pat==ii,'cd4']<-approx(thisDat$time,thisDat$CD4,newDat[newDat$pat==ii,'day'],rule=2)$y
  }
  cd4Mat<-tapply(newDat$cd4,list(newDat$pat,newDat$week),c)
  vlMat<-tapply(newDat$vl,list(newDat$pat,newDat$week),c)
  # no data beyond last point so presumably shouldn't affect (careful if CD4 data ends far before IC50)
  cd4Mat<-t(apply(cd4Mat,1,dnar::fillDown))
  vlMat<-t(apply(vlMat,1,dnar::fillDown))
  weekOffset<-1-min(round(time/7))
  tMax<-150 #MM15 ends at 103
  superIds<-structure(1:(sum(!is.na(superTimes))+1),.Names=c('NA',names(superTimes)[!is.na(superTimes)]))
  superId<-ifelse(time<superTimes[patient],NA,superIds[patient])
  superId[is.na(superId)]<-1
  newDat$isSuper<-newDat$day>=ifelse(is.na(superTimes),Inf,superTimes)[newDat$pat]
  dat=list(
    nVirus=length(ic50),
    nArt=max(artId),
    nPatient=max(patientId),
    ic50=log(ic50),
    patients=patientId[patient],
    artIds=ifelse(is.na(artId[patient]),9999,artId[patient]),
    days=ifelse(time<0,0,time),
    artStart=artStart[names(artId)],
    nSample=max(sampleId),
    sample=sampleId[sample],
    daysBeforeArt=ifelse(is.na(timePreArt),9999,timePreArt),
    hasArt=!is.na(timePreArt),
    isFast=names(patientId) %in% fastProgressors,
    isNon=names(patientId) %in% nonProgressors,
    cd4=as.numeric(cd4),
    vl=as.numeric(vl),
    baseVl=baseVl[names(patientId)],
    weeks=round(time/7),
    vlMat=log(vlMat[names(patientId),weekOffset+1:tMax]),
    cd4Mat=cd4Mat[names(patientId),weekOffset+1:tMax]/100,
    tMax=tMax,
    superId=superId,
    isSuper=superId>1,
    nSuper=sum(!is.na(superTimes))
  )
  fit <- sampling(mod, data = dat, iter=nIter, chains=chains,thin=2,control=list(adapt_delta=.98,max_treedepth=12),...)
  return(list('fit'=fit,pats=patientId,arts=artId,artStart=artStart,sample=sampleId,dat=dat,superId=superIds,simDat=newDat))
}
baseVl<-dnar::withAs(xx=unique(meta[meta$time<365*2&meta$time>180,c('time','vl','mm')]),tapply(xx$vl,xx$mm,function(xx)mean(log(xx),na.rm=TRUE)))
dat$week<-as.integer(round(dat$time/7))

fit<-dnar::withAs(xx=dat[!is.na(dat$ic50)&!dat$qvoa,],bayesIC50_3(ic50Mod,xx$ic50,xx$time,xx$timeBeforeArt,xx$pat,fastProgressors=c('MM15','WEAU'),nonProgressors=c('MM55','MM62'),cd4=(xx$fillCD4)/100,vl=log(xx$fillVl),baseVl=baseVl,nIter=5000,chains=50,meta=meta,superTimes=superTimes<50))
fitB<-dnar::withAs(xx=dat[!is.na(dat$beta)&!dat$qvoa,],bayesIC50_3(ic50Mod,xx$beta,xx$time,xx$timeBeforeArt,xx$pat,fastProgressors=c('MM15','WEAU'),nonProgressors=c('MM55','MM62'),cd4=(xx$fillCD4)/100,vl=log(xx$fillVl),baseVl=baseVl,nIter=5000,chains=50,meta=meta,superTimes=superTimes))
}
#fitR<-dnar::withAs(xx=dat[!is.na(dat$replication)&!dat$qvoa,],bayesIC50_3(ic50Mod,xx$replication,xx$time,xx$timeBeforeArt,xx$pat,fastProgressors=c('MM15','WEAU'),nonProgressors=c('MM55','MM62'),cd4=(xx$fillCD4)/100,vl=log(xx$fillVl),baseVl=baseVl,nIter=5000,chains=50,meta=meta,superTimes=superTimes))
#save(fit,file='out/bayesFit_20200602.Rdat') #CD4 marginalized
#save(fit,file='out/bayesFit_20200602.Rdat') #CD4 marginalized
#save(fit,file='out/bayesFit_20200602c.Rdat')# CD4+VL marginalized 2000 rep, 125 weeks
#save(fit,file='out/bayesFit_20200604.Rdat')# CD4+VL marginalized 2000 rep, 100 weeks, setpoint
#save(fit,file='out/bayesFit_20200605.Rdat')# CD4+VL (split out fast), marginalized 1000 rep, super, 100 weeks, setpoint
#save(fit,fitB,fitR,file='out/bayesFit_20200607.Rdat')# non acute, CD4 (split out fast), marginalized 1000 rep, super, 100 weeks, setpoint
#save(fit,fitB,fitR,file='out/bayesFit_20200609.Rdat')# non acute, CD4 (split out fast), marginalized 5000 rep, super, 100 weeks, setpoint, separate lp, heirarchical nadir
#save(fit,fitB,fitR,file='out/bayesFit_20200612.Rdat')# remove non-progressor acute
#save(fit,fitB,fitR,file='out/bayesFit_20200615.Rdat')#  non-progressor acute and nadirChange
#save(fit,fitB,file='out/bayesFit_20200623.Rdat')#  neg binom time back in 
#save(fit,fitB,file='out/bayesFit_20200711.Rdat')# seperate cd4 beta for typ and non
#save(fit,fitB,file='out/bayesFit_20200714.Rdat')# seperate nadir times
#load(file='out/bayesFit_20200623.Rdat')
#
print(fit$fit,pars=c('nadirTimeRaw','nadirChangeRaw','expectedIC50','acuteRaw','lp','vlBetaRaw','cd4BetaRaw','superBetaRaw'),include=FALSE)
print(fitB$fit,pars=c('nadirTimeRaw','nadirChangeRaw','expectedIC50','acuteRaw','lp','vlBetaRaw','cd4BetaRaw','superBetaRaw'),include=FALSE)
#print(fitR$fit,pars=c('nadirTimeRaw','nadirChangeRaw','expectedIC50','acuteRaw','lp','vlBetaRaw','cd4BetaRaw','superBetaRaw'),include=FALSE)



softMax<-function(xx)exp(xx)/sum(exp(xx))
logsumexp <- function (x) {
    y = max(x)
  y + log(sum(exp(x - y)))
}
softmax <- function (x) {
    exp(x - logsumexp(x))
}

meanCrI<-function(xx)c(mean(xx,na.rm=TRUE),quantile(xx,c(.025,.975),na.rm=TRUE))
calcPreds<-function(mat,dat,newDat){
  newDat$scaleCd4<-newDat$cd4/100
  newDat$scaleVl<-log(newDat$vl)
  weekMax<-tapply(round(dat$days/7),names(dat$patients),max)
  vlCd4<-unique(data.frame('pat'=names(dat$patients),'day'=dat$days,'vl'=dat$vl,'cd4'=dat$cd4))
  newDat$patId<-sapply(newDat$pat,function(xx)dat$patients[names(dat$patients)==xx][1])
  weeks<-outer(newDat$pat,1:dat$tMax,function(xx,yy)yy)
  isNadir<-newDat$week<weeks
  propNadir<-newDat$week/weeks
  vlMat<-dat$vlMat[newDat$pat,]
  cd4Mat<-dat$cd4Mat[newDat$pat,]
  preds<-do.call(cbind,parallel::mclapply(split(mat,sort(rep(1:30,length.out=nrow(mat)))),function(xx,...){
    thisMat<-matrix(xx,ncol=ncol(mat))
    colnames(thisMat)<-colnames(mat)
    apply(thisMat,1,function(xx){
      pS<-do.call(rbind,lapply(1:dat$nPatient,function(pat)softmax(xx[grep(sprintf('lp\\[%d,',pat),names(xx))])))
      #cd4Bs<-xx[sprintf('cd4Beta[%d]',dat$patients)]
      #vlBs<-xx[sprintf('vlBeta[%d]',dat$patients)]
      #acutes<-xx[sprintf('acute[%d]',dat$patients)]
      #nadirs<-xx[sprintf('nadirChange[%d]',dat$patients)]
      cd4Bs<-xx[sprintf('cd4Beta[%d]',newDat$patId)]
      acutes<-xx[sprintf('acute[%d]',newDat$patId)]
      nadirs<-xx[sprintf('nadirChange[%d]',newDat$patId)]
      pred<-acutes+ifelse(isNadir,propNadir*nadirs,nadirs+(newDat$scaleCd4-cd4Mat)*cd4Bs)
      out<-apply(pred*pS[newDat$patId,],1,sum)
      out
    })
  },mc.cores=40))
  sigmas<-mat[,'sigma']
  isos<-apply(rbind(mat[,'sigma'],preds),2,function(xx)rnorm(length(xx[-1]),xx[-1],xx[1]))
  crI<-cbind(t(apply(preds,1,meanCrI)),t(apply(isos,1,meanCrI))[,-1])
  colnames(crI)<-c('mean','lower','upper','predL','predU')
  list('preds'=preds,'isos'=isos,'simData'=newDat,'crI'=crI)
}
predIc50<-calcPreds(as.matrix(fit$fit),fit$dat,fit$simDat)
predIc50B<-calcPreds(as.matrix(fitB$fit),fitB$dat,fitB$simDat)
#predIc50R<-calcPreds(as.matrix(fitR$fit),fitR$dat,fitR$simDat)


source('functions.R')
plotPointsLine<-function(dat,ic50,ii,ylab,addTitle=TRUE,sims=NULL,addFit=TRUE,filterAfter=TRUE,superTimes=NULL,artStart=NULL,ylim=range(ic50,na.rm=TRUE)){
  plot(dat$time/7,ic50,yaxt='n',log='y',bg=patCols[dat$pat],pch=21,type='n',xlab='',ylab=ylab,xaxt='n',cex=1.4,ylim=ylim)
  if(addTitle)title(ii,line=-1)
  thisDat<-dat[dat$pat==ii,]
  thisIc50<-ic50[dat$pat==ii]
  if(sum(!is.na(thisIc50))==0)return()
  if(addFit){
    if(is.null(sims)){
      thisFit<-lm(I(log(thisIc50))~time+time2,dat=thisDat)
      fakeDays<-(min(thisDat$time)):(max(thisDat$time)+50)
      fakeDf<-data.frame('time'=fakeDays,'time2'=fakeDays^2,'logTime'=log(fakeDays),'logTime2'=log(fakeDays)^2,'logTime3'=log(fakeDays)^3,'logTime4'=log(fakeDays)^4)
      predIc50<-predict(thisFit,fakeDf,interval='confidence')
      lines(fakeDays/7,exp(predIc50[,'fit']),col=patCols[ii])
      polygon(c(fakeDays/7,rev(fakeDays)/7),c(exp(predIc50[,'lwr']),rev(exp(predIc50[,'upr']))),col=patCols2[ii],border=NA)
      predIc50<-predict(thisFit,fakeDf,interval='prediction')
      polygon(c(fakeDays/7,rev(fakeDays)/7),c(exp(predIc50[,'lwr']),rev(exp(predIc50[,'upr']))),col=patCols3[ii],border=NA)
    }else{
      if(is.list(sims)&&all(names(sims)==c('preds','isos','simData','crI'))){
        selector2<-sims$simData$pat==ii&sims$simData$day>0
        lines(sims$simData[selector2,'day']/7,exp(sims$crI[selector2,'mean']),col=patCols[ii])
        polygon(c(sims$simData[selector2,'day'],rev(sims$simData[selector2,'day']))/7,exp(c(sims$crI[selector2,'lower'],rev(sims$crI[selector2,'upper']))),,col=patCols2[ii],border=NA)
        polygon(c(sims$simData[selector2,'day'],rev(sims$simData[selector2,'day']))/7,exp(c(sims$crI[selector2,'predL'],rev(sims$crI[selector2,'predU']))),,col=patCols3[ii],border=NA)
        #if(!is.null(artStarts)&&!is.na(artStarts[ii]))abline(v=artStarts[ii]/7,lty=ifelse(ii=='WEAU',2,1))
        if(!is.null(artStarts)&&!is.na(artStarts[ii])){
          if(ii=='WEAU'){
            nRects<-31
            vertBreaks<-seq(par('usr')[3],par('usr')[4],length.out=nRects+1)
            rect(artStarts[ii]/7,10^vertBreaks[seq(1,nRects,2)],maxDates[ii]/7,10^vertBreaks[seq(2,nRects+1,2)],col='#00000008',border=NA)
            rect(artStarts[ii]/7,10^vertBreaks[1],maxDates[ii]/7,10^vertBreaks[nRects+1],col='#00000004',border=NA)
            #abline(v=artStarts[ii]/7,lty=2)
            #points(rep(artStarts[1]/7,30),10^seq(par('usr')[3],par('usr')[4],length.out=30),pch=21,col=NA,bg='#00000055',cex=.5) 
          }else{
            rect(artStarts[ii]/7,10^par('usr')[3],maxDates[ii]/7,10^par('usr')[4],col='#00000011',border=NA)
          }
        }
      }else{
        sim<-sims[[ii]][,sims[[ii]]['time',]<ifelse(is.na(artStart[ii])||!filterAfter,Inf,artStart[ii])]
        polygon(c(sim['time',],rev(sim['time',]))/7,exp(c(sim['lowCI',],rev(sim['highCI',]))),border=NA,col=patCols2[ii])
        polygon(c(sim['time',],rev(sim['time',]))/7,exp(c(sim['lowPred',],rev(sim['highPred',]))),border=NA,col=patCols3[ii])
        lines(sim['time',]/7,exp(sim['mean',]),col=patCols[ii])
        if(ii!='WEAU'||!filterAfter)abline(v=artStart[ii]/7,lty=2) #MAGIC NUMBER. SUPPRESSING WEAU VERTICAL LINE SINCE DIFFERS FROM OTHERS
      }
    }
  }
  if(!is.null(superTimes) && !is.na(superTimes[ii])){
    baseY<-10^(par('usr')[3]+diff(par('usr')[3:4])*.02)
    topY<-10^(par('usr')[3]+diff(par('usr')[3:4])*.09)
    segments(rep(superTimes[ii]/7,3)+c(-7,0,7),c(topY,baseY,topY),rep(superTimes[ii]/7,3)+c(0,7,-7),c(baseY,topY,topY),lwd=1.2)
  }
  #points(thisDat$time/7,thisIc50,pch=21+thisDat$bulk,bg=patCols[ii])
  points(thisDat$time/7,thisIc50,pch=21+thisDat$bulk,bg=c(NA,ifelse(any(!thisDat$bulk&!thisDat$qvoa),'#ffff3399',NA),sprintf('%scc',classCol['qvoa']))[1+thisDat$bulk+2*thisDat$qvoa],cex=1+.2*thisDat$qvoa)
  normCol<-patCols[ii]
  for(isBulk in sort(unique(thisDat$bulk)))dnar::withAs(xx=cbind(thisDat,thisIc50)[thisDat$bulk==isBulk,],points(xx$time/7,xx$thisIc50,pch=21+xx$bulk,bg=c(normCol,ifelse(any(!thisDat$bulk&!thisDat$qvoa),'#ffff33',normCol),sprintf('%s',classCol['qvoa']))[1+xx$bulk+2*xx$qvoa],cex=1+.1*xx$qvoa))
}
plotCondenseIfn<-function(dat,ic50,ylab,showLegend=TRUE,sims=NULL,addFit=TRUE,filterAfter=TRUE,subplotLetters=LETTERS[1:3],superTimes=superTimes,artStarts=NULL,ylimExpand=c(1,1)){
  par(mar=c(0,0,0,0))
  layout(lay2,width=c(.42,rep(1,2),.01),height=c(.16,c(1,1,1,.2,1,.2,1),ifelse(showLegend,1.3,.32)))
  counter<-1
  for(ii in patOrder){
    plotPointsLine(dat,ic50,ii,ylab,sims=sims,addFit=addFit,filterAfter=filterAfter,superTimes=superTimes,ylim=range(ic50,na.rm=TRUE)*ylimExpand)
    if(counter>4)axis(1,seq(0,6,2)*100,cex.axis=1.2,mgp=c(2.75,.4,0),tcl=-.3)
    if(counter>4)axis(1,seq(1:3)*100,rep('',3),cex.axis=1.2,mgp=c(2.75,.7,0),tcl=-.3)
    if(counter%%2==1)logAxis(2,las=1,cex.axis=1.1,mgp=c(3,.7,0),axisMax=max(ic50,na.rm=TRUE)*1.3)
    labCex<-1.7
    if(counter==5)text(par('usr')[1]-.35*diff(par('usr')[1:2]),10^mean(par('usr')[3:4]),ylab,srt=90,xpd=NA,cex=labCex)
    if(counter==9)text(max(par('usr')[1:2]),10^(par('usr')[3]-.25*diff(par('usr')[3:4])),'Weeks from onset of symptoms',xpd=NA,cex=labCex)
    if(counter==9&showLegend)legend(par('usr')[1]-diff(par('usr')[1:2])*.3,10^(par('usr')[3]-diff(par('usr')[3:4])*.45),c(ifelse(is.null(sims),'Quadratic regression','Bayesian model'),ifelse(is.null(sims),'95% confidence interval','95% credible interval'),'95% prediction interval','Limiting dilution isolate','Bulk isolate'),col=c(patCols[1],NA,NA,'black','black'),pt.bg=c(NA,patCols2[1],patCols3[1],patCols[1],patCols[1]),lty=c(1,NA,NA,NA,NA),pch=c(NA,22,22,21,22),border=NA,pt.cex=c(3.2,3.2,3.2,1.4,1.4),cex=1.1,xjust=0,yjust=1,xpd=NA)
    if(counter==1)text(grconvertX(-.29,from='npc'),grconvertY(1.11,from='npc'),subplotLetters[1],xpd=NA,adj=c(0,1),cex=2.5) 
    if(counter==7)text(grconvertX(-.29,from='npc'),grconvertY(1.13,from='npc'),subplotLetters[2],xpd=NA,adj=c(0,1),cex=2.5)
    if(counter==9)text(grconvertX(-.29,from='npc'),grconvertY(1.13,from='npc'),subplotLetters[3],xpd=NA,adj=c(0,1),cex=2.5)
    counter<-counter+1
  }
}
for(showLegend in c(TRUE,FALSE)){
pdf(sprintf('out/bayesFit2%s.pdf',ifelse(showLegend,'_legend','')),width=3,height=7)
noQvoa<-dat#dat[!dat$qvoa,]
plotCondenseIfn(noQvoa,noQvoa$ic50,ylab=expression('IFN'*alpha*'2 IC'[50]*' (pg/ml)'),sims=predIc50,filterAfter=TRUE,showLegend=showLegend,superTimes=superTimes)
  #text(grconvertX(.01,from='ndc'),grconvertY(.99,from='ndc'),'A',xpd=NA,adj=c(0,1),cex=2.5)
plotCondenseIfn(noQvoa,noQvoa$beta,ylab=expression('IFN'*beta*' IC'[50]*' (pg/ml)'),sims=predIc50B,filterAfter=TRUE,subplotLetters=LETTERS[4:6],showLegend=showLegend,superTimes=superTimes,ylimExpand=c(1,2))
  #text(grconvertX(.01,from='ndc'),grconvertY(.99,from='ndc'),'B',xpd=NA,adj=c(0,1),cex=2.5)
#plotCondenseIfn(noQvoa,noQvoa$replication,ylab='Replicative capacity (day 7 p24 ng/ml)',sims=simsR)
#plotCondenseIfn(noQvoa,noQvoa$replication,ylab='Replicative capacity (day 7 p24 ng/ml)',sims=predIc50R,filterAfter=TRUE,subplotLetters=LETTERS[1:3],showLegend=showLegend,superTimes=superTimes)
dev.off()
}
system('pdftk out/bayesFit2.pdf cat 1 2 output tmp.pdf')
system('pdfcrop --margin "7 0 0 0" tmp.pdf tmp2.pdf')
system('pdfjam tmp2.pdf --nup 2x1 --outfile tmp3.pdf')
system('pdfcrop tmp3.pdf out/Fig._2.pdf')

ps<-dnar::withAs(mat=as.matrix(fit$fit),do.call(rbind,parallel::mclapply(1:fit$dat$nPatient,function(pat)apply(apply(mat[,grep(sprintf('lp\\[%d,',pat),colnames(mat))],1,softmax),1,mean),mc.cores=10)))
psB<-dnar::withAs(mat=as.matrix(fitB$fit),do.call(rbind,parallel::mclapply(1:fitB$dat$nPatient,function(pat)apply(apply(mat[,grep(sprintf('lp\\[%d,',pat),colnames(mat))],1,softmax),1,mean),mc.cores=10)))
pdf('out/nadirTimes.pdf',width=14,height=4)
par(mfrow=c(2,5))
for(ii in 1:nrow(ps)){
  plot(1:150,ps[ii,],main=names(fit$pats)[ii],type='l',xlab='Week',ylab='Probability of nadir')
}
for(ii in 1:nrow(ps)){
  plot(1:150,psB[ii,],main=names(fit$pats)[ii],type='l',xlab='Week',ylab='Probability of nadir')
}
dev.off()

maxNadirTimes<-apply(ps,1,function(xx)sum(1:150*xx))
maxNadirTimesB<-apply(psB,1,function(xx)sum(1:150*xx))
names(maxNadirTimes)<-names(fit$pats)
names(maxNadirTimesB)<-names(fitB$pats)
pdf('out/cd4_vs_ic50_test.pdf',width=8,height=4)
for(jj in c('ic50','beta')){
if(jj=='ic50')maxTimes<-maxNadirTimes
else maxTimes<-maxNadirTimesB
par(mfrow=c(1,2),mar=c(3,3,1,3))
for(ii in sort(unique(dat$pat))){
  thisDat<-dat[dat$pat==ii&!dat$qvoa,]
  thisTime<-ifelse(jj=='beta',maxNadirTimesB[ii],maxTimes[ii])*7
  plot(1,1,type='n',xlab='CD4 count (cells/ul)',ylab=ifelse(jj=='beta','IFNb IC50','IFNa2 IC50'),log='y',yaxt='n',xlim=range(thisDat$fillCD4),ylim=range(dat[,jj],na.rm=TRUE),main=ii)
  xx<-data.frame('ic50'=tapply(thisDat[,jj],thisDat$time,function(xx)exp(mean(log(xx),na.rm=TRUE))),
    'cd4'=tapply(thisDat$fillCD4,thisDat$time,unique))
  xx$time<-as.numeric(rownames(xx))
  points(xx$cd4,xx$ic50,pch=21,bg=ifelse(xx$time<=thisTime,'#FF000033','blue'),col=ifelse(xx$time<=thisTime,'#00000033','#000000'))
  segments(xx$cd4[-nrow(xx)],xx$ic50[-nrow(xx)],xx$cd4[-1],xx$ic50[-1])
  dnar::logAxis(las=1)
  #
  plot(1,1,type='n',xlab='Time',ylab=ifelse(jj=='beta','IFNb IC50','IFNa2 IC50'),log='y',yaxt='n',xlim=range(dat$time),ylim=range(dat[,jj],na.rm=TRUE),main=ii,mgp=c(2.5,1,0))
  xx$time<-as.numeric(rownames(xx))
  points(xx$time,xx$ic50,pch=21,bg=ifelse(xx$time<=thisTime,'#FF000033','blue'),col=ifelse(xx$time<=thisTime,'#00000033','#000000'))
  segments(xx$time[-nrow(xx)],xx$ic50[-nrow(xx)],xx$time[-1],xx$ic50[-1])
  dnar::logAxis(las=1)
  cd4Range<-c(0,1500)
  xx$convertCd4<-10^(par('usr')[3]+xx$cd4/1500*diff(par('usr')[3:4]))
  lines(xx$time,xx$convertCd4)
  abline(v=thisTime,lty=2)
}
par(mfrow=c(3,2),mar=c(0,0,0,0))
patExclude<-c("WEAU",'MM55','MM62','MM15')
for(ii in sort(unique(dat$pat[!dat$pat %in% patExclude]))){
  plot(1,1,type='n',xlab='CD4 count (cells/ul)',ylab=ifelse(jj=='beta','IFNb IC50','IFNa2 IC50'),log='y',yaxt='n',xlim=c(-550,350),ylim=10^c(range(log10(dat[dat$time>300&!dat$pat %in% patExclude,jj]),na.rm=TRUE)))
  thisDat<-dat[dat$pat==ii&!dat$qvoa,]
  thisDat$cd4<-thisDat$fillCD4-fit$dat$cd4Mat[ii,maxTimes[ii]]*100
  xx<-data.frame('ic50'=tapply(thisDat[,jj],thisDat$time,function(xx)exp(mean(log(xx),na.rm=TRUE))),
    'cd4'=tapply(thisDat$cd4,thisDat$time,unique))
  xx$time<-as.numeric(rownames(xx))
  #xx$cd4<-xx$cd4-fit$dat$cd4Mat[ii,maxTimes[ii]]*100
  points(thisDat$cd4,thisDat[,jj],pch=21,bg=ifelse(thisDat$time<=maxTimes[ii]*7,patCols2[ii],patCols[ii]),col=ifelse(thisDat$time<=maxTimes[ii]*7,'#00000033','#000000'))
  segments(xx$cd4[-nrow(xx)],xx$ic50[-nrow(xx)],xx$cd4[-1],xx$ic50[-1],col=ifelse(xx$time[-nrow(xx)]<maxTimes[ii]*7,'#00000022','#000000'))
  dnar::logAxis(las=1)
  #print(dnar::withAs(xx=xx[xx$time>maxTimes[ii],],cor.test(log(thisDat$ic50),thisDat$cd4)))
}
}
dev.off()

pdf('out/cd4_vs_ic50.pdf',width=4,height=8)
for(jj in c('ic50','beta')){
par(mfrow=c(5,2),mar=c(3,3.5,1,4.5))
for(ii in patOrder){
  thisDat<-dat[dat$pat==ii&!dat$qvoa,]
  thisTime<-ifelse(jj=='beta',maxNadirTimesB[ii],maxNadirTimes[ii])
  plot(1,1,type='n',xlab='',ylab=ifelse(jj=='beta','IFNb IC50','IFNa2 IC50'),log='y',yaxt='n',xlim=range(thisDat$time)/7,ylim=range(dat[,jj],na.rm=TRUE),main=ii,mgp=c(2.6,1,0))
  xx<-data.frame('ic50'=tapply(thisDat[,jj],thisDat$time,function(xx)exp(mean(log(xx),na.rm=TRUE))),
    'cd4'=tapply(thisDat$fillCD4,thisDat$time,unique))
  xx$time<-as.numeric(rownames(xx))/7
  #thisDat$cd4<-thisDat$CD4-fit$dat$cd4Mat[ii,thisTime]*100
  #xx$cd4<-xx$cd4-fit$dat$cd4Mat[ii,thisTime]*100
  thisDat$cd4<-thisDat$CD4
  #points(xx$time,xx$ic50,pch=21,bg=ifelse(xx$time<=thisTime,patCols3[ii],patCols[ii]),col=ifelse(xx$time<=thisTime,'#00000033','#000000'))
  points(thisDat$time/7,thisDat[,jj],pch=21,bg=ifelse(thisDat$time/7<=thisTime,patCols3[ii],patCols[ii]),col=ifelse(thisDat$time/7<=thisTime,'#00000033','#000000'))
  segments(xx$time[-nrow(xx)],xx$ic50[-nrow(xx)],xx$time[-1],xx$ic50[-1])
  dnar::logAxis(las=1)
  #cd4Range<-c(0,max(xx$cd4)*1.1)
  cd4Range<-c(0,1200)
  convertY<-function(xx,cd4,usr)10^(usr[1]+(xx-cd4[1])/diff(cd4)*diff(usr))
  xx$convertCd4<-convertY(xx$cd4,cd4Range,par('usr')[3:4])
  axis(4,convertY(pretty(cd4Range),cd4Range,par('usr')[3:4]),pretty(cd4Range),las=1,col='blue',col.axis='blue')
  lines(xx$time,xx$convertCd4/1.01,col='white')
  lines(xx$time,xx$convertCd4,col='blue')
  text(par('usr')[2]+.52*diff(par('usr')[1:2]),10^mean(par('usr')[3:4]),'CD4 count',srt=-90,xpd=NA,col='blue',cex=1.2)
  box()
  abline(v=thisTime,h=convertY(0,cd4Range,par('usr')[3:4]),lty=2)
}
}
dev.off()


exampleCurve<-function(fit,type='typical',times=1:600,confInt=.05){
  mat<-as.matrix(fit)
  getParams<-function(mat,cols,sdCol){
    out<-data.frame('mean'=apply(mat[,cols,drop=FALSE],1,sum))
    out$pat<-rnorm(nrow(out),out$mean,mat[,sdCol])
    return(out)
  }
  if(type=='fast'){
    acute<-getParams(mat,c('acuteMean','fastAcute'),'acuteSD')
    change<-getParams(mat,c('nadirChangeMean','fastChangeMean'),'nadirChangeSD')
  }else if (type=='non'){
    acute<-getParams(mat,c('acuteMean','nonAcute'),'acuteSD')
    change<-getParams(mat,c('nadirChangeMean','nonChangeMean'),'nadirChangeSD')
  }else{
    acute<-getParams(mat,'acuteMean','acuteSD')
    change<-getParams(mat,c('nadirChangeMean'),'nadirChangeSD')
  }
  if(type=='fast')  nadirTime<-rnbinom(nrow(mat),mu=mat[,'nadirTimeMean']*exp(mat[,'nadirTimeFast']),size=mat[,'nadirTimePhi'])
  else if(type=='non') nadirTime<-rnbinom(nrow(mat),mu=mat[,'nadirTimeMean']*exp(mat[,'nadirTimeNon']),size=mat[,'nadirTimePhi'])
  else nadirTime<-rnbinom(nrow(mat),mu=mat[,'nadirTimeMean'],size=mat[,'nadirTimePhi'])
  sigma<-mat[,'sigma']
  timeMat<-outer(nadirTime^-1,times/7)
  timeMat[timeMat>1]<-1
  out<-list(mean=acute$mean+timeMat*change$mean,pat=acute$pat+timeMat*change$pat)
  out$iso<-out$pat
  out$iso[,]<-rnorm(prod(dim(out$pat)),unlist(out$pat),sigma)
  quants<-lapply(out,function(xx){out<-t(apply(xx,2,meanCrI));colnames(out)<-c('mean','lower','upper');out})
  return(quants)
}
times<-1:700
example<-exampleCurve(fit$fit,times=times)
exampleFast<-exampleCurve(fit$fit,type='fast',times=times)
exampleSlow<-exampleCurve(fit$fit,type='non',times=times)
exampleB<-exampleCurve(fitB$fit,times=times)
exampleFastB<-exampleCurve(fitB$fit,type='fast',times=times)
exampleSlowB<-exampleCurve(fitB$fit,type='non',times=times)

plotExample<-function(typical,fast,non,ylab='IFNa2 IC50',dat,ifnCol='ic50',addPat=FALSE){
  plotSub<-function(example,col,ylim,time=NULL,ic50=NULL,yAxis=TRUE){
    plot(1,1,type='n',xlim=c(0,85),ylim=ylim,las=1,log='y',yaxt='n',ylab='',xlab='',mgp=c(2.6,.7,0))
    if(yAxis){
      logAxis(las=1)
    }
    #points(time,ic50,bg=sprintf('%s55',col),cex=.5,pch=21,col='#00000099')
    polygon(c(times,rev(times))/7,exp(c(example$mean[,'lower'],rev(example$mean[,'upper']))),col=sprintf('%s33',col),border=NA)#,border=sprintf('%s66',col),lty=2)
    if(addPat)polygon(c(times,rev(times))/7,exp(c(example$pat[,'lower'],rev(example$pat[,'upper']))),col=sprintf('%s33',col),border=NA)#,border=sprintf('%s66',col),lty=3)
    polygon(c(times,rev(times))/7,exp(c(example$iso[,'lower'],rev(example$iso[,'upper']))),col=sprintf('%s33',col),border=NA)#,border=sprintf('%s11',col))
    lines(times/7,exp(example$mean[,'mean']),col=col)
  }
  withAs(xx=dat[dat$pat %in% c('MM14','MM23','MM33','MM34','MM39','MM40'),],plotSub(typical,classCol['typical'],ylim=exp(range(typical$mean)),time=xx$time/7,ic50=xx[,ifnCol]))
  title(main='Typical',line=-1)
  title(ylab=ylab,xpd=NA,mgp=c(2.9,1,0),cex.lab=1.3)
  withAs(xx=dat[dat$pat %in% c('MM55','MM62'),],plotSub(non,classCol['slow'],ylim=exp(range(typical$mean)),yAxis=FALSE,time=xx$time/7,ic50=xx[,ifnCol]))
  title(xlab='Weeks from onset of symptoms',mgp=c(2.,.7,0),xpd=NA,cex.lab=1.3)
  title(main='Non',line=-1)
  withAs(xx=dat[dat$pat %in% c('MM15','WEAU'),],plotSub(fast,classCol['fast'],ylim=exp(range(typical$mean)),yAxis=FALSE,time=xx$time/7,ic50=xx[,ifnCol]))
  title(main='Fast',line=-1)
}
pdf('out/bayesExample.pdf',width=8,height=3)
  layout(matrix(c(0,1:3,0,4:6,0),nrow=1),width=c(.42,1,1,1,.5,1,1,1,.01))
  #drops
  par(mar=c(3,0,0.1,0))
  plotExample(example,exampleFast,exampleSlow,dat=dat)
  text(grconvertX(.001,from='ndc'),grconvertY(.99,from='ndc'),'A',xpd=NA,adj=c(0,1),cex=2)
  plotExample(exampleB,exampleFastB,exampleSlowB,ylab='IFNb IC50',dat=dat,ifnCol='beta')
  text(grconvertX(.52,from='ndc'),grconvertY(.99,from='ndc'),'B',xpd=NA,adj=c(0,1),cex=2)
dev.off()




calcCd4Curve<-function(fit,fakeCd4=seq(-600,500,length.out=200)){
  mat<-as.matrix(fit$fit)
  nadir<-mat[,grep('acute\\[',colnames(mat))]+mat[,grep('nadirChange\\[',colnames(mat))]
  betas<-mat[,grep('cd4Beta\\[',colnames(mat))]
  sigma<-mat[,'sigma']
  preds<-lapply(fit$pats,function(ii){
    preds<-nadir[,ii]+outer(betas[,ii],fakeCd4/100)
    out<-cbind(fakeCd4,exp(t(apply(preds,2,meanCrI))))
    isos<-preds
    isos[,]<-rnorm(prod(dim(isos)),isos,sigma)
    out<-cbind(out,exp(t(apply(isos,2,quantile,c(.025,.975)))))
    colnames(out)<-c('cd4','mean','lower','upper','isoLower','isoUpper')
    return(out)
  })
  return(preds)
}
cd4Curve<-calcCd4Curve(fit)
cd4CurveB<-calcCd4Curve(fitB)
calcDayProbs<-function(fit,days=-50:4000){
  mat<-as.matrix(fit$fit)
  preds<-lapply(fit$pats,function(ii){
    ps<-apply(apply(mat[,grep(sprintf('lp\\[%d,',ii),colnames(mat))],1,softmax),1,mean)
    dayPs<-approx(0:fit$dat$tMax*7,c(0,cumsum(ps)),days,rule=2)$y
  })
  out<-cbind('day'=days,do.call(cbind,preds))
  rownames(out)<-out[,'day']
  return(out)
}
dayProbs<-calcDayProbs(fit)
dayProbsB<-calcDayProbs(fitB)

betas<-dnar::withAs(mat=as.matrix(fit$fit),apply(mat[,grep('cd4Beta\\[',colnames(mat))],2,meanCrI))
betasB<-dnar::withAs(mat=as.matrix(fitB$fit),apply(mat[,grep('cd4Beta\\[',colnames(mat))],2,meanCrI))
nadirs<-(dnar::withAs(mat=as.matrix(fit$fit),apply(mat[,grep('acute\\[',colnames(mat))]+mat[,grep('nadirChange\\[',colnames(mat))],2,meanCrI)))
nadirsB<-(dnar::withAs(mat=as.matrix(fitB$fit),apply(mat[,grep('acute\\[',colnames(mat))]+mat[,grep('nadirChange\\[',colnames(mat))],2,meanCrI)))
colnames(betas)<-colnames(betasB)<-colnames(nadirs)<-colnames(nadirsB)<-names(fit$pats)
ps<-dnar::withAs(mat=as.matrix(fit$fit),do.call(rbind,parallel::mclapply(1:fit$dat$nPatient,function(pat)apply(apply(mat[,grep(sprintf('lp\\[%d,',pat),colnames(mat))],1,softmax),1,mean),mc.cores=10)))
psB<-dnar::withAs(mat=as.matrix(fitB$fit),do.call(rbind,parallel::mclapply(1:fitB$dat$nPatient,function(pat)apply(apply(mat[,grep(sprintf('lp\\[%d,',pat),colnames(mat))],1,softmax),1,mean),mc.cores=10)))
maxNadirTimes<-apply(ps,1,function(xx)sum(1:150*xx))
maxNadirTimesB<-apply(psB,1,function(xx)sum(1:150*xx))
nadirCd4s<-apply(ps*fit$dat$cd4Mat,1,sum)*100
nadirCd4sB<-apply(psB*fit$dat$cd4Mat,1,sum)*100
names(maxNadirTimes)<-names(nadirCd4s)<-names(fit$pats)
names(maxNadirTimesB)<-names(nadirCd4sB)<-names(fitB$pats)
pdf('out/bayesCD4.pdf',width=4,height=7)
plotCD4<-function(dat,nadirCd4s,nadirTimes,cd4Curves,dayProbs,ic50Col='ic50',ylab='IFNa2 IC50',patCols,addArrows=FALSE){
  dat$cd4Diff<-dat$fillCD4-nadirCd4s[dat$pat]
  dat$prob<-apply(dat[,c('pat','time')],1,function(xx)dayProbs[as.character(as.numeric(xx[2])),xx[1]])
  ylim<-range(dat[dat$prob>.1,ic50Col],na.rm=TRUE)
  xlim<-range(dat[dat$prob>.1,'cd4Diff'],na.rm=TRUE)
  par(mar=c(0,0,0,0))
  pats<-patOrder
  for(ii in 1:length(pats)){
    thisDat<-dat[dat$pat==pats[ii]&!dat$qvoa&dat$prob>.05,]
    plot(1,1,log='y',yaxt='n',xaxt='n',xlim=xlim,ylim=ylim,type='n',xlab='',ylab='')
    if(ii >8)axis(1)
    if(ii %% 2==1)dnar::logAxis(las=1,axisMin=10^-3)
    xx<-data.frame('ic50'=tapply(thisDat[,ic50Col],thisDat$time,function(xx)exp(mean(log(xx),na.rm=TRUE))),
      'cd4'=tapply(thisDat$cd4Diff,thisDat$time,unique),'prob'=tapply(thisDat$prob,thisDat$time,unique))
    xx$time<-as.numeric(rownames(xx))
    if(addArrows)arrows(xx$cd4[-nrow(xx)],xx$ic50[-nrow(xx)],xx$cd4[-1],xx$ic50[-1],col=grey(1-(.2+.8*xx$prob[-nrow(xx)])),length=.05)
    thisCurve<-cd4Curves[[pats[ii]]]
    if(!is.null(cd4Curves)){
      lines(thisCurve[,'cd4'],thisCurve[,'mean'],col=patCols[pats[ii]])
      polygon(c(thisCurve[,'cd4'],rev(thisCurve[,'cd4'])),c(thisCurve[,'lower'],rev(thisCurve[,'upper'])),border=NA,col=sprintf('%s33',patCols[pats[ii]]),lty=2)
      polygon(c(thisCurve[,'cd4'],rev(thisCurve[,'cd4'])),c(thisCurve[,'isoLower'],rev(thisCurve[,'isoUpper'])),border=NA,col=sprintf('%s18',patCols[pats[ii]]),lty=2)
    }
    points(thisDat$cd4Diff,thisDat[,ic50Col],pch=21,bg=grey(1-thisDat$prob))#,bg=ifelse(thisDat$afterNadir,patCols[pats[ii]],patCols2[pats[ii]]),col=ifelse(thisDat$afterNadir,'#000000','#00000033'))
    title(main=pats[ii],line=-1,adj=.1)
    box()
    labCex=1.3
    if(ii==5)text(par('usr')[1]-.28*diff(par('usr')[1:2]),10^mean(par('usr')[3:4]),ylab,srt=90,xpd=NA,cex=labCex)
    if(ii==9)text(max(par('usr')[1:2]),10^(par('usr')[3]-.29*diff(par('usr')[3:4])),'Change in CD4 from nadir (cells/ul)',xpd=NA,cex=labCex)
  }
}
layout(rbind(0,cbind(0,matrix(c(1:10),nrow=5,byrow=TRUE),0),0),height=c(.01,rep(1,5),.35),width=c(.33,1,1,.01))
plotCD4(dat,nadirCd4s,maxNadirTimes,cd4Curve,dayProbs,'ic50',patCols=patCols)
text(grconvertX(.001,from='ndc'),grconvertY(.994,from='ndc'),'C',xpd=NA,adj=c(0,1),cex=2)
plotCD4(dat,nadirCd4sB,maxNadirTimesB,cd4CurveB,dayProbsB,'beta',patCols=patCols,ylab='IFNb IC50')
text(grconvertX(.001,from='ndc'),grconvertY(.994,from='ndc'),'D',xpd=NA,adj=c(0,1),cex=2)
dev.off()

pdf('out/bayesCombo.pdf',width=8,height=10)
  lay1<-matrix(c(0,rep(1:3,each=2),0,rep(4:6,each=2),0),nrow=1)
  lay2<-do.call(rbind,c(list(0),lapply(seq(7,15,2),function(xx)matrix(c(0,rep(xx+0:1,each=3),0,rep(xx+10:11,each=3),0),nrow=1)),list(0)))
  lay<-rbind(lay1,lay2)
  layout(lay,height=c(1,.39,rep(1,5),.34),width=c(1,rep(1,6),1.1,rep(1,6),.01))
  #drops
  par(mar=c(0,0,0,0))
  plotExample(example,exampleFast,exampleSlow,dat=dat)
  text(grconvertX(.001,from='ndc'),grconvertY(.999,from='ndc'),'A',xpd=NA,adj=c(0,1),cex=2.2)
  plotExample(exampleB,exampleFastB,exampleSlowB,ylab='IFNb IC50',dat=dat,ifnCol='beta')
  text(grconvertX(.51,from='ndc'),grconvertY(.999,from='ndc'),'B',xpd=NA,adj=c(0,1),cex=2.2)
  plotCD4(dat,nadirCd4s,maxNadirTimes,cd4Curve,dayProbs,'ic50',patCols=patCols)
  text(grconvertX(.51,from='ndc'),grconvertY(.815,from='ndc'),'D',xpd=NA,adj=c(0,1),cex=2.2)
  plotCD4(dat,nadirCd4sB,maxNadirTimesB,cd4CurveB,dayProbsB,'beta',patCols=patCols,ylab='IFNb IC50')
  text(grconvertX(.001,from='ndc'),grconvertY(.815,from='ndc'),'C',xpd=NA,adj=c(0,1),cex=2.2)
dev.off()
file.copy('out/bayesCombo.pdf','out/Fig._S4.pdf',overwrite=TRUE)




pdf('out/bayesCD4_2.pdf',width=4,height=4)
plotCD4Together<-function(dat,nadirCd4s,dayProbs,ic50Col='ic50',ylab='IFNa2 IC50',patCols){
  pats<-patOrder[-9:-10]
  dat$cd4Diff<-dat$fillCD4-nadirCd4s[dat$pat]
  dat<-dat[dat$pat %in% pats,]
  ylim<-range(dat[,ic50Col],na.rm=TRUE)
  xlim<-range(dat[,'cd4Diff'],na.rm=TRUE)
  par(mar=c(2.5,3.5,.1,.1))
  plot(1,1,log='y',yaxt='n',xlim=xlim,ylim=ylim,type='n',ylab='',xlab='',mgp=c(1,.5,0),tcl=-.3)
  title(xlab='Change in CD4 from nadir (cells/ul)',mgp=c(1.6,1,0))
  title(ylab=ylab,mgp=c(2.6,1,0))
  dnar::logAxis(las=1,axisMin=10^-3)
  for(ii in 1:length(pats)){
    thisDat<-dat[dat$pat==pats[ii]&!dat$qvoa,]
    thisDat$prob<-dayProbs[as.character(thisDat$time),pats[ii]]
    xx<-data.frame('ic50'=tapply(thisDat[,ic50Col],thisDat$time,function(xx)exp(mean(log(xx),na.rm=TRUE))),
      'cd4'=tapply(thisDat$cd4Diff,thisDat$time,unique),'prob'=tapply(thisDat$prob,thisDat$time,unique))
    xx$time<-as.numeric(rownames(xx))
    #segments(xx$cd4[-nrow(xx)],xx$ic50[-nrow(xx)],xx$cd4[-1],xx$ic50[-1],col=grey(1-(.2+.8*xx$prob[-nrow(xx)])))
    points(thisDat$cd4Diff,thisDat[,ic50Col],pch=21,bg=grey(1-thisDat$prob))#,bg=ifelse(thisDat$afterNadir,patCols[pats[ii]],patCols2[pats[ii]]),col=ifelse(thisDat$afterNadir,'#000000','#00000033'))
  }
}
plotCD4Together(dat,nadirCd4s,dayProbs,'ic50',patCols=patCols)
plotCD4Together(dat,nadirCd4sB,dayProbsB,'beta',patCols=patCols,ylab='IFNb IC50')
dev.off()


mat<-as.matrix(fit$fit)
matB<-as.matrix(fitB$fit)
message('Acute')
print(exp(meanCrI(mat[,'acuteMean'])))
print(exp(meanCrI(matB[,'acuteMean'])))
message('Fast acute')
print(exp(meanCrI(mat[,'fastAcute'])))
print(exp(meanCrI(matB[,'fastAcute'])))
message('Non acute')
print(exp(meanCrI(mat[,'nonAcute'])))
print(exp(meanCrI(matB[,'nonAcute'])))
message('Nadir change')
print(exp(-meanCrI(mat[,'nadirChangeMean'])))
print(exp(-meanCrI(matB[,'nadirChangeMean'])))
message('CD4 change')
print(exp(-meanCrI(mat[,'cd4BetaMean[1]'])))
print(exp(-meanCrI(matB[,'cd4BetaMean[1]'])))
message('Nadir time')
print((meanCrI(mat[,'nadirTimeMean']))*7)
print((meanCrI(matB[,'nadirTimeMean']))*7)
message('Nadir time fast')
print(meanCrI(mat[,'nadirTimeMean']*exp(mat[,'nadirTimeFast']))*7)
print(meanCrI(matB[,'nadirTimeMean']*exp(matB[,'nadirTimeFast']))*7)
message('Nadir time non')
print(meanCrI(mat[,'nadirTimeMean']*exp(mat[,'nadirTimeNon']))*7)
print(meanCrI(matB[,'nadirTimeMean']*exp(matB[,'nadirTimeNon']))*7)
message('Non nadir change')
print(exp(-meanCrI(mat[,'nadirChangeMean']+mat[,'nonChangeMean'])))
print(exp(-meanCrI(matB[,'nadirChangeMean']+matB[,'nonChangeMean'])))
print(exp(meanCrI(mat[,'nonAcute']+mat[,'nonChangeMean'])))
print(exp(meanCrI(matB[,'nonAcute']+matB[,'nonChangeMean'])))
message('Fast nadir change')
mean(mat[,'nadirChangeMean']+mat[,'fastChangeMean']>0)
mean(matB[,'nadirChangeMean']+matB[,'fastChangeMean']>0)
print(exp(-meanCrI(mat[,'nadirChangeMean']+mat[,'fastChangeMean'])))
print(exp(-meanCrI(matB[,'nadirChangeMean']+matB[,'fastChangeMean'])))
print(exp(meanCrI(mat[,'fastAcute']+mat[,'fastChangeMean'])))
print(exp(meanCrI(matB[,'fastAcute']+matB[,'fastChangeMean'])))
message('Fast vs normal nadir change')
print(exp(meanCrI(mat[,'fastChangeMean'])))
print(exp(meanCrI(matB[,'fastChangeMean'])))
message('CD4 change normal')
print(exp(-meanCrI(mat[,'cd4BetaMean[1]'])))
print(exp(-meanCrI(matB[,'cd4BetaMean[1]'])))
message('CD4 change fast')
print(exp(-meanCrI(mat[,'cd4BetaMean[2]'])))
print(exp(-meanCrI(matB[,'cd4BetaMean[2]'])))
message('CD4 change slow')
print(exp(-meanCrI(mat[,'cd4BetaMean[2]'])))
print(exp(-meanCrI(matB[,'cd4BetaMean[2]'])))

print(fit$fit,pars=c('nadirTimeRaw','nadirChangeRaw','expectedIC50','acuteRaw','lp','vlBetaRaw','cd4BetaRaw','superBetaRaw'),include=FALSE)

if(FALSE){
plotIfn<-function(fit,predIc50,ylab='IFNa2 IC50'){
  for(ii in sort(unique(names(fit$dat$patients)))){
    plot(1,1,type='n',xlim=range(fit$dat$days),ylim=range(exp(fit$dat$ic50))*c(.8,1),main=ii,xlab='DFOSx',ylab='',log='y',yaxt='n',mgp=c(1.8,.5,0),tcl=-.3)
    mtext(ylab,2,line=2.9,cex=.83)
    selector<-names(fit$dat$patients)==ii
    points(fit$dat$days[selector],exp(fit$dat$ic50[selector]))
    dnar::logAxis(las=1)
    selector2<-predIc50$simData$pat==ii
    lines(predIc50$simData[selector2,'day'],exp(apply(predIc50$preds[selector2,],1,mean)))
    quants<-(apply(predIc50$preds[selector2,],1,quantile,c(.025,.975),na.rm=TRUE))
    quants2<-(apply(predIc50$isos[selector2,],1,quantile,c(.025,.975),na.rm=TRUE))
    polygon(c(predIc50$simData[selector2,'day'],rev(predIc50$simData[selector2,'day'])),exp(c(quants[1,],rev(quants[2,]))),col='#00000033',border=NA)
    polygon(c(predIc50$simData[selector2,'day'],rev(predIc50$simData[selector2,'day'])),exp(c(quants2[1,],rev(quants2[2,]))),col='#00000033',border=NA)
    abline(v=superTimes[ii],col='red')
  }
}
pdf('Rplots.pdf',width=17,height=7)
  par(mfrow=c(2,5),mar=c(3.5,4,1,.1))
  plotIfn(fit,predIc50)
  plotIfn(fitB,predIc50B,'IFNb IC50')
  #plotIfn(fitR,predIc50R,'Replicative capacity')
  for(ii in sort(unique(names(fit$dat$patients)))){
    selector2<-predIc50$simData$pat==ii
    plot(predIc50$simData[selector2,'day'],(predIc50$simData[selector2,'cd4']),xlim=range(fit$dat$days),ylim=range((fit$dat$cd4)*100),main=ii,xlab='DFOSx',ylab='CD4',mgp=c(1.8,.5,0),tcl=-.3)
    abline(v=superTimes[ii],col='red')
  }
  for(ii in sort(unique(names(fit$dat$patients)))){
    selector2<-predIc50$simData$pat==ii
    plot(predIc50$simData[selector2,'day'],(predIc50$simData[selector2,'vl']),xlim=range(fit$dat$days),ylim=range(exp(fit$dat$vl)),main=ii,xlab='DFOSx',ylab='VL',mgp=c(1.8,.5,0),tcl=-.3,log='y',yaxt='n')
    abline(v=superTimes[ii],col='red')
    dnar::logAxis(las=1)
  }
dev.off()

apply(mat[,grep('lp\\[',colnames(mat))],2,mean)
pdf('Rplots.pdf',height=20,width=20);plot(apply(mat[,grep('lp\\[',colnames(mat))],2,mean));dev.off()
pdf('Rplots.pdf',height=20,width=20);traceplot(fit$fit,'lp');dev.off()
}

