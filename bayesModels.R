library('rstan')
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source('functions.R')
if(!dir.exists('out'))dir.create('out')

dat<-read.csv('allLongitudinal.csv',stringsAsFactors=FALSE)
dat$week<-as.integer(round(dat$time/7))
meta<-read.csv('allLongitudinalMeta.csv')
maxDates<-tapply(meta$time,meta$mm,max)
artStarts<-tapply(meta$daysBeforeArt+meta$time,meta$mm,unique)

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
      real weeks[nVirus];
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
      vector[nPatient] nadirChange;
      vector[nPatient] cd4Beta;
      matrix[nPatient,tMax] lp;
      acute=acuteMean+acuteRaw*acuteSD+isNon*nonAcute+isFast*fastAcute;
      nadirChange=nadirChangeMean+nadirChangeRaw*nadirChangeSD+fastChangeMean*isFast+nonChangeMean*isNon;
      cd4Beta=(cd4BetaMean[1]*(1-isFast).*(1-isNon)+cd4BetaMean[2]*isFast+cd4BetaMean[3]*isNon)+cd4BetaRaw*cd4BetaSD;
      for (s in 1:tMax){
        for(jj in 1:nPatient)lp[jj,s]=neg_binomial_2_lpmf(s|nadirTimeMean*exp(isFast[jj]*nadirTimeFast)*exp(isNon[jj]*nadirTimeNon),nadirTimePhi);
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
      sigma~gamma(1,.1);
      nadirTimePhi~cauchy(0,10);//gamma(1,.01);
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
    }\n
  '
  ic50Mod <- stan_model(model_code = ic50_bp_mar)

  bayesIC50<-function(mod,ic50,time,timePreArt,patient,chains=50,nIter=30000,fastProgressors=c(),cd4=c(),vl=c(),nonProgressors=FALSE,...){
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
      thisDat<-unique(dat[dat$pat==ii&dat$time<=weekMax[ii]*7,c('time','CD4','vl')])
      newDat[newDat$pat==ii,'vl']<-approx(thisDat$time,thisDat$vl,newDat[newDat$pat==ii,'day'],rule=2)$y
      newDat[newDat$pat==ii,'cd4']<-approx(thisDat$time,thisDat$CD4,newDat[newDat$pat==ii,'day'],rule=2)$y
    }
    cd4Mat<-tapply(newDat$cd4,list(newDat$pat,newDat$week),c)
    vlMat<-tapply(newDat$vl,list(newDat$pat,newDat$week),c)
    # no data beyond last point so presumably shouldn't affect (careful if CD4 data ends far before IC50)
    cd4Mat<-t(apply(cd4Mat,1,fillDown))
    vlMat<-t(apply(vlMat,1,fillDown))
    weekOffset<-1-min(round(time/7))
    tMax<-150 
    dat=list(
      nVirus=length(ic50),
      nArt=max(artId),
      nPatient=max(patientId),
      ic50=log(ic50),
      patients=patientId[patient],
      nSample=max(sampleId),
      sample=sampleId[sample],
      isFast=names(patientId) %in% fastProgressors,
      isNon=names(patientId) %in% nonProgressors,
      cd4=as.numeric(cd4),
      vl=as.numeric(vl),
      weeks=round(time/7),
      vlMat=log(vlMat[names(patientId),weekOffset+1:tMax]),
      cd4Mat=cd4Mat[names(patientId),weekOffset+1:tMax]/100,
      tMax=tMax
    )
    fit <- sampling(mod, data = dat, iter=nIter, chains=chains,thin=2,control=list(adapt_delta=.98,max_treedepth=12),...)
    return(list('fit'=fit,pats=patientId,arts=artId,artStart=artStart,sample=sampleId,dat=dat,simDat=newDat))
  }

  #running 5000 iterations will take a while. 50 chains assumes a computer with large number of cores
  fit<-withAs(xx=dat[!is.na(dat$ic50)&!dat$qvoa,],bayesIC50(ic50Mod,xx$ic50,xx$time,xx$timeBeforeArt,xx$pat,fastProgressors=c('MM15','WEAU'),nonProgressors=c('MM55','MM62'),cd4=(xx$fillCD4)/100,vl=log(xx$fillVl),nIter=5000,chains=50))
  fitB<-withAs(xx=dat[!is.na(dat$beta)&!dat$qvoa,],bayesIC50(ic50Mod,xx$beta,xx$time,xx$timeBeforeArt,xx$pat,fastProgressors=c('MM15','WEAU'),nonProgressors=c('MM55','MM62'),cd4=(xx$fillCD4)/100,vl=log(xx$fillVl),nIter=5000,chains=50))
}
#
print(fit$fit,pars=c('nadirTimeRaw','nadirChangeRaw','expectedIC50','acuteRaw','lp','vlBetaRaw','cd4BetaRaw'),include=FALSE)
print(fitB$fit,pars=c('nadirTimeRaw','nadirChangeRaw','expectedIC50','acuteRaw','lp','vlBetaRaw','cd4BetaRaw'),include=FALSE)



logsumexp<-function(xx)max(xx)+log(sum(exp(xx-max(xx))))
softmax<-function(xx)exp(xx-logsumexp(xx))

meanCrI<-function(xx)c(mean(xx,na.rm=TRUE),quantile(xx,c(.025,.975),na.rm=TRUE))
calcPreds<-function(mat,dat,newDat){
  newDat$scaleCd4<-newDat$cd4/100
  newDat$scaleVl<-log(newDat$vl)
  weekMax<-tapply(dat$week,names(dat$patients),max)
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


plotPointsLine<-function(dat,ic50,ii,ylab,addTitle=TRUE,sims=NULL,addFit=TRUE,filterAfter=TRUE,artStart=NULL,ylim=range(ic50,na.rm=TRUE)){
  plot(dat$time/7,ic50,yaxt='n',log='y',bg=patCols[dat$pat],pch=21,type='n',xlab='',ylab=ylab,xaxt='n',cex=1.4,ylim=ylim)
  if(addTitle)title(ii,line=-1)
  thisDat<-dat[dat$pat==ii,]
  thisIc50<-ic50[dat$pat==ii]
  if(sum(!is.na(thisIc50))==0)return()
  if(addFit){
    if(is.list(sims)&&all(names(sims)==c('preds','isos','simData','crI'))){
      selector2<-sims$simData$pat==ii&sims$simData$day>0
      lines(sims$simData[selector2,'day']/7,exp(sims$crI[selector2,'mean']),col=patCols[ii])
      polygon(c(sims$simData[selector2,'day'],rev(sims$simData[selector2,'day']))/7,exp(c(sims$crI[selector2,'lower'],rev(sims$crI[selector2,'upper']))),,col=patCols2[ii],border=NA)
      polygon(c(sims$simData[selector2,'day'],rev(sims$simData[selector2,'day']))/7,exp(c(sims$crI[selector2,'predL'],rev(sims$crI[selector2,'predU']))),,col=patCols3[ii],border=NA)
      if(!is.null(artStarts)&&!is.na(artStarts[ii])){
        if(ii=='WEAU'){
          nRects<-31
          vertBreaks<-seq(par('usr')[3],par('usr')[4],length.out=nRects+1)
          rect(artStarts[ii]/7,10^vertBreaks[seq(1,nRects,2)],maxDates[ii]/7,10^vertBreaks[seq(2,nRects+1,2)],col='#00000008',border=NA)
          rect(artStarts[ii]/7,10^vertBreaks[1],maxDates[ii]/7,10^vertBreaks[nRects+1],col='#00000004',border=NA)
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
  points(thisDat$time/7,thisIc50,pch=21+thisDat$bulk,bg=c(NA,ifelse(any(!thisDat$bulk&!thisDat$qvoa),'#ffff3399',NA),sprintf('%scc',classCol['qvoa']))[1+thisDat$bulk+2*thisDat$qvoa],cex=1+.2*thisDat$qvoa)
  normCol<-patCols[ii]
  for(isBulk in sort(unique(thisDat$bulk)))withAs(xx=cbind(thisDat,thisIc50)[thisDat$bulk==isBulk,],points(xx$time/7,xx$thisIc50,pch=21+xx$bulk,bg=c(normCol,ifelse(any(!thisDat$bulk&!thisDat$qvoa),'#ffff33',normCol),sprintf('%s',classCol['qvoa']))[1+xx$bulk+2*xx$qvoa],cex=1+.1*xx$qvoa))
}
plotCondenseIfn<-function(dat,ic50,ylab,showLegend=TRUE,sims=NULL,addFit=TRUE,filterAfter=TRUE,subplotLetters=LETTERS[1:3],artStarts=NULL,ylimExpand=c(1,1)){
  par(mar=c(0,0,0,0))
  layout(lay2,width=c(.42,rep(1,2),.01),height=c(.16,c(1,1,1,.2,1,.2,1),ifelse(showLegend,1.3,.32)))
  counter<-1
  for(ii in patOrder){
    plotPointsLine(dat,ic50,ii,ylab,sims=sims,addFit=addFit,filterAfter=filterAfter,ylim=range(ic50,na.rm=TRUE)*ylimExpand)
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
pdf('out/bayesFit.pdf',width=3,height=7)
plotCondenseIfn(dat,dat$ic50,ylab=expression('IFN'*alpha*'2 IC'[50]*' (pg/ml)'),sims=predIc50,filterAfter=TRUE,showLegend=TRUE)
plotCondenseIfn(dat,dat$beta,ylab=expression('IFN'*beta*' IC'[50]*' (pg/ml)'),sims=predIc50B,filterAfter=TRUE,subplotLetters=LETTERS[4:6],showLegend=TRUE,ylimExpand=c(1,2))
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

betas<-withAs(mat=as.matrix(fit$fit),apply(mat[,grep('cd4Beta\\[',colnames(mat))],2,meanCrI))
betasB<-withAs(mat=as.matrix(fitB$fit),apply(mat[,grep('cd4Beta\\[',colnames(mat))],2,meanCrI))
nadirs<-(withAs(mat=as.matrix(fit$fit),apply(mat[,grep('acute\\[',colnames(mat))]+mat[,grep('nadirChange\\[',colnames(mat))],2,meanCrI)))
nadirsB<-(withAs(mat=as.matrix(fitB$fit),apply(mat[,grep('acute\\[',colnames(mat))]+mat[,grep('nadirChange\\[',colnames(mat))],2,meanCrI)))
colnames(betas)<-colnames(betasB)<-colnames(nadirs)<-colnames(nadirsB)<-names(fit$pats)
ps<-withAs(mat=as.matrix(fit$fit),do.call(rbind,parallel::mclapply(1:fit$dat$nPatient,function(pat)apply(apply(mat[,grep(sprintf('lp\\[%d,',pat),colnames(mat))],1,softmax),1,mean),mc.cores=10)))
psB<-withAs(mat=as.matrix(fitB$fit),do.call(rbind,parallel::mclapply(1:fitB$dat$nPatient,function(pat)apply(apply(mat[,grep(sprintf('lp\\[%d,',pat),colnames(mat))],1,softmax),1,mean),mc.cores=10)))
maxNadirTimes<-apply(ps,1,function(xx)sum(1:150*xx))
maxNadirTimesB<-apply(psB,1,function(xx)sum(1:150*xx))
nadirCd4s<-apply(ps*fit$dat$cd4Mat,1,sum)*100
nadirCd4sB<-apply(psB*fit$dat$cd4Mat,1,sum)*100
names(maxNadirTimes)<-names(nadirCd4s)<-names(fit$pats)
names(maxNadirTimesB)<-names(nadirCd4sB)<-names(fitB$pats)
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
    if(ii %% 2==1)logAxis(las=1,axisMin=10^-3)
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


