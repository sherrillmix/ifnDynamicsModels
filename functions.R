library(dnar)

patOrder<-c("MM14","MM23","MM33","MM34","MM39","MM40","MM55","MM62","MM15","WEAU")

classCol<-c('rebound'='#e41a1c','qvoa'='#377eb8','typical'='#ff7f00','slow'='#4daf4a','fast'='#f781bf','postQvoa'='#3fe0d0')
patCols<-structure(rep(classCol[c('typical','slow','fast')],c(6,2,2)),.Names=patOrder)

patCols2<-sprintf('%s44',patCols)
patCols3<-sprintf('%s22',patCols)
names(patCols2)<-names(patCols3)<-names(patCols)


withAs<-function(...,expr=NULL){
  parent<-parent.frame()
  env<-new.env(parent=parent)
  dotVars<-match.call(expand.dots=FALSE)$'...'
  if(missing(expr)){
    expr<-dotVars[[length(dotVars)]]
    dotVars<-dotVars[-length(dotVars)]
  }
  if(is.null(names(dotVars))||any(names(dotVars)==''))stop(simpleError('Unassigned variables passed to withAs'))
  #make sure to get parent above and not within the mapply
  mapply(function(as,val)assign(as,eval(val,envir=parent),env),names(dotVars),dotVars)
  return(eval(substitute(expr), env))
}

fillDown<-function(x,emptyStrings=c(NA,''),errorIfFirstEmpty=TRUE){
  #depending on %in% to catch NAs if necessary
  isEmpty<-x %in% emptyStrings
  if(isEmpty[1]&errorIfFirstEmpty)stop(simpleError('First value empty'))
  #if first is empty and we don't want errors then have to just fill down from it anyway
  isEmpty[1]<-FALSE
  ids<-1:length(x)
  ids[isEmpty]<-0
  ids<-cummax(ids)
  return(x[ids])
}



logAxis<-function(side=2,exponent=TRUE,addExtra=!exponent,minorTcl=-.2,axisMin=-Inf,axisMax=Inf,offset=0,col.ticks='black',axisVals=NULL,...){
  if(side %in% c(2,4)) parX<-sort(graphics::par('usr')[3:4])
  else parX<-sort(graphics::par('usr')[1:2])
  minX<-max(10^parX[1],axisMin)
  maxX<-min(10^parX[2],axisMax)
  if(log10(maxX)-log10(minX)>400)stop(simpleError('Huge range in logged axis'))
  allTicks<-unlist(lapply(floor(log10(minX)):ceiling(log10(maxX)),function(x)1:9*10^x))
  allTicks<-allTicks[allTicks<=maxX & allTicks>=minX]
  graphics::axis(side,allTicks+offset,rep('',length(allTicks)),tcl=minorTcl,col.ticks=col.ticks)
  if(ceiling(log10(minX))<=floor(log10(maxX)))prettyY<-seq(ceiling(log10(minX)),floor(log10(maxX)),1)
  else prettyY<-c()
  graphics::axis(side,10^prettyY+offset,rep('',length(prettyY)),tcl=minorTcl*2,col.ticks=col.ticks)
  if(length(prettyY)>7)prettyY<-pretty(prettyY)
  if(length(prettyY)==0)prettyY<-c(ceiling(log10(minX)),floor(log10(maxX)))
  if(addExtra){
    origPretty<-prettyY
    if(sum(prettyY>=log10(minX)&prettyY<=log10(maxX))<4)prettyY<-unique(c(prettyY,origPretty+log10(5),origPretty-log10(10/5)))
    if(sum(prettyY>=log10(minX)&prettyY<=log10(maxX))<4)prettyY<-unique(c(prettyY,origPretty+log10(2),origPretty-log10(10/2)))
    if(sum(prettyY>=log10(minX)&prettyY<=log10(maxX))<4)prettyY<-unique(c(prettyY,origPretty+log10(3),origPretty-log10(10/3)))
    if(sum(prettyY>=log10(minX)&prettyY<=log10(maxX))<4)prettyY<-unique(c(prettyY,origPretty+log10(7),origPretty-log10(10/7)))
  }
  if(exponent){
    if(any(prettyY%%1!=0))labs<-sapply(prettyY,function(x)as.expression(bquote(.(10^(x%%1))%*%10^.(floor(x)))))
    else labs<-ifelse(prettyY==0,1,sapply(prettyY,function(x)as.expression(bquote(10^.(floor(x))))))
  }
  else labs<-10^prettyY
  if(is.null(axisVals))axisVals<-prettyY
  graphics::axis(side,10^prettyY+offset,ifelse(prettyY %in% axisVals,labs,''),col.ticks=col.ticks,...)
  return(invisible(list('minor'=allTicks,'major'=10^prettyY)))
}


plotQvoa2<-function(ic50,label,pos,class,study,speed,ylab='IFNa2 IC50 (pg/ml)',mar=c(6.9,4,.1,3.9),cex.axis=1.25,startDown=FALSE,pats=NULL,classCols,labelXAxis=TRUE,log='y',ylim=range(ic50),addStudyLines=TRUE,mgp=c(2.4,1,0)){
  speedPch<-c('Fast'=23,'Slow'=21,'Standard'=21,'Other'=21)
  speedCex<-c('Fast'=1,'Slow'=1,'Standard'=1,'Other'=1)
  deemphasize<-c("Acute","6 Month","Donor","Nadir","Last","Acute Recipients","Chronic Donors","Acute Recipient","Chronic Donor",'1 Year','Acute','Chronic','Early','Late')
  spread<--vipor::offsetX(log10(ic50),label,width=.4)
  spread[ic50>.9&class=='Post-ATI']<-0 #manual tweak to improve figure
  ii<-'all'
  marSpace<-0
  if(ii!='all')selector<-label %in% subsets[[ii]]
  else selector<-rep(TRUE,length(label))
  par(mar=mar)
  plot(pos[label[selector]]+spread[selector],ic50[selector],log=log,yaxt='n',ylab=ylab,xlab='',xaxt='n',type='n',cex.lab=1.2,ylim=ylim,mgp=mgp,xlim=range(pos)+c(-1,1),xaxs='i')
  spread2<-ave(ic50,label,FUN=function(xx)beeswarm::swarmx(rep(0,length(xx)),xx,cex=.4,priority=ifelse(length(xx)>2,'density','ascending'))$x)
  spread<-ifelse(ave(ic50,label,FUN=length)>5,spread,spread2)
  if(ii=='all'&exists('subsets'))marSpaces<-sapply(subsets,function(xx)diff(convertUserToLine(1:0,2))*sum(!names(pos) %in% xx))
  if(!is.null(pats)){
    ranges<-tapply(pos[label],pats,range)
    #print(ranges[sapply(ranges,diff)>0.1])
    axFunc<-if(log=='y')function(xx)10^xx else function(xx)xx
    lapply(ranges[!is.na(sapply(ranges,diff))&sapply(ranges,diff)>0.1&!grepl('MM|^$',names(ranges))],function(xx)rect(min(xx)-.4,axFunc(par('usr')[3]),max(xx)+.4,axFunc(par('usr')[4]),col='#00000014',border=NA))
  }
  slantSelect<-!grepl('Outgrowth|VOA',names(pos)) &grepl('VOA MM|Outgrowth MM|Acute|6 Month|1 Year|Nadir|Last|Chronic',names(pos))
  #textOffsets<-c('Acute'=-.2,'6 Month'=-.2,'Nadir'=0,'Last'=.2,'Acute Recipients'=-.3,'Chronic Donors'=.3)[names(pos)[slantSelect]]
  if(any(slantSelect)){
    textOffsets<-c('Acute'=0.3,'6 Month'=0,'Nadir'=0,'Last'=0,'Chronic'=.3,'Acute Recipients'=-.1,'Chronic Donors'=.1)[names(pos)[slantSelect]]
    textOffsets[is.na(textOffsets)]<-0
    slantAxis(1,pos[slantSelect],names(pos)[slantSelect],srt=-45,cex=cex.axis,location=.8,xpd=NA,textOffsets=textOffsets)
  }
  if(log=='y'){
    if(diff(par('usr')[3:4])>1)logAxis(las=1,mgp=c(1,.7,0))
    else axis(2,las=1,mgp=c(3,.5,0),tcl=-.3)
  }else{
    axis(2,las=1,mgp=c(3,.5,0),tcl=-.3)
  }
  abline(v=pos,col='#00000055',lty=3)
  points(pos[label[selector]]+spread[selector],ic50[selector],pch=speedPch[speed[selector]],bg=classCols[class[selector]],col=ifelse(speedPch[speed[selector]]>20,ifelse(class[selector] %in% deemphasize,'#00000066','#00000088'),classCols[class[selector]]),lwd=ifelse(class[selector]=='Rebound',1,1),cex=ifelse(class[selector]=='QVOA',max(speedCex),speedCex[speed[selector]]))
  abline(v=pos['Acute Recipients']-(pos['Acute Recipients']-pos[which(names(pos)=='Acute Recipients')-1])/2,lty=1,col='#000000')
  #abline(v=pos['Acute']-(pos['Acute']-pos[which(names(pos)=='Acute')-1])/2,lty=1,col='#000000')
  #abline(v=pos['Outgrowth MM14']-(pos['Outgrowth MM14']-pos[which(names(pos)=='Outgrowth MM14')-1])/2,lty=1,col='#00000099')
  #abline(v=pos['Outgrowth MM14']-(pos['Outgrowth MM14']-pos[which(names(pos)=='Outgrowth MM14')-1])/2,lty=1,col='#00000099')
  #abline(v=pos['VOA MM14']-(pos['VOA MM14']-pos[which(names(pos)=='VOA MM14')-1])/2,lty=1,col='#000000')
  cols<-c('#00000033','#00000000')
  studyOrder<-unique(sapply(names(pos),function(xx)study[label==xx][1]))
  studyOrder<-studyOrder[studyOrder!='Transmission'&!is.na(studyOrder)]
  counter<-1
  for(ii in studyOrder){
    if(ii=='MM'){
      minPos<-pos[min(which(names(pos) %in% label[study==ii&class %in% c('VOA','QVOA','Outgrowth')]))]
      maxPos<-pos[max(which(names(pos) %in% label[study==ii&class %in% c('VOA','QVOA','Outgrowth')]))]
    }else{
      minPos<-pos[min(which(names(pos) %in% label[study==ii]))]
      maxPos<-pos[max(which(names(pos) %in% label[study==ii]))]
      nextPos<-pos[max(which(names(pos) %in% label[study==ii]))+1]
      if(ii!='Reservoir'&addStudyLines)abline(v=(maxPos+nextPos)/2,lty=2)
    }
    #rect(minPos-.5,10^par('usr')[3],maxPos+.5,10^par('usr')[4],col=cols[counter%%2+1],border=NA)
    counter<-counter+1
    if(is.na(minPos))browser()
    if(labelXAxis)axis(1,mean(c(minPos,maxPos)),sub('MM','Outgrowth',sub('ATI','Interrupt',sub('BEAT','IFNa2',sub('/','/\n',ii)))),padj=1,mgp=c(3,.2+2*(startDown+counter)%%2,0),tcl=-.7+-2*(startDown+counter)%%2,cex.axis=cex.axis)
  }
  return(list(ylim=ylim,pos=pos))
}

meanCrI<-function(xx)c(mean(xx,na.rm=TRUE),quantile(xx,c(.025,.975),na.rm=TRUE))
logsumexp<-function(xx)max(xx)+log(sum(exp(xx-max(xx))))
softmax<-function(xx)exp(xx-logsumexp(xx))

insetScale<-function(breaks,col,insetPos=c(.025,.015,.04,.25),main='',offset=1e-3,at=NULL,labels=NULL,cex=1,labXOffset=0,labYOffset=0){
  if(length(breaks)!=length(col)+1)stop('Number of breaks must be one more than colors')
  insetPos<-c(graphics::grconvertY(insetPos[1],'nfc','user'),graphics::grconvertX(insetPos[2],'nfc','user'),graphics::grconvertY(insetPos[3],'nfc','user'),graphics::grconvertX(insetPos[4],'nfc','user'))
  breakPos<-((breaks)-(min(breaks)))/max((breaks)-(min(breaks)))*(insetPos[4]-insetPos[2])+insetPos[2]
  #add a bit of offset to avoid pdf viewers displaying breaks between exact rectangle border meeting
  offsetPos<-breakPos[-1]+c(rep(offset*diff(range(breakPos)),length(breakPos)-2),0)
  graphics::rect(breakPos[-length(breakPos)],insetPos[1],offsetPos,insetPos[3],col=col,xpd=NA,border=NA)
  graphics::rect(insetPos[2],insetPos[1],insetPos[4],insetPos[3],xpd=NA)
  if(is.null(at)){
    at<-pretty(breaks)
    at<-at[at<=max(breaks)&at>=min(breaks)]
  }
  if(is.null(labels))labels<-at
  convertPos<-(at-(min(breaks)))/((max(breaks))-(min(breaks)))*(insetPos[4]-insetPos[2])+insetPos[2]
  graphics::segments(convertPos,insetPos[1],convertPos,insetPos[1]-diff(insetPos[c(1,3)])*.1,xpd=NA)
  graphics::text(convertPos+labXOffset*diff(insetPos[c(2,4)]),insetPos[1]-diff(insetPos[c(1,3)])*.175+labYOffset*diff(insetPos[c(1,3)]),labels,xpd=NA,adj=c(.5,1),cex=.85*cex)
  graphics::text(mean(insetPos[c(2,4)]),insetPos[3]+diff(insetPos[c(1,3)])*.45,main,xpd=NA,adj=c(.5,0),cex=cex)
  invisible(NULL)
}
