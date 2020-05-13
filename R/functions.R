
library(devtools)
library(roxygen2)
library(reshape2) 
library(ggplot2) 
library(DescTools) 
library(ggpubr) 

#'  Align the curves 
#' 	@param data Do you love cats? Defaults to TRUE.
#' 	@param days Do you love cats? Defaults to TRUE.
#' 	@param k Do you love cats? Defaults to TRUE.
#' 	@keywords align.curves
#' 	@export
#' 	@examples
#' align.curves()
align.curves<-function(data, days = 1:263, k=3, n=5){
	title='Mean clinical score'

	m=melt(data[,c('condition','sex','souris',days)],id.vars=c('condition','sex','souris'))
	m=na.omit(m)
	m$key = paste(m$souris,m$sex,m$condition,m$variable,sep='.')
	m$souris.id = paste(m$souris,m$sex,m$condition,sep='.')
	
	vec=c()
	ages = c()
	vars = c()
	sums=c()
	m=m[!duplicated(m$key),]
	m=m[order(m$souris.id,m$variable),]
	for(mouse in unique(m$souris.id)){
		sub = subset(m, souris.id==mouse)	
		sick = F
		check = T
		for(day in unique(as.character(sub$variable))){
			if(check){
				if(is.na(sub[which(sub$variable==day),'value'])){
					if(!is.na(sub[which(sub$variable==as.character(as.numeric(day)+1)),'value'])){
						sub[which(sub$variable==day),'value']=0	
					}
				}
				if(sub[which(sub$variable==day),'value']>0){
					ages=c(ages,day)
					sick=T
					check=F
					vec[(length(vec)-n):(length(vec))]=1
				}
			}
			if(sick){
				vec=c(vec,1)
			}else{
				vec=c(vec,0)
			}
		}
	}
	m$sick = vec
	m.sick=subset(m, sick==1)
	m.sick$variable=as.numeric(as.character(m.sick$variable))
	df=data.frame(matrix(ncol=264))
	m.sick$souris.id=paste(m.sick$souris,m.sick$sex,m.sick$condition,sep='.')
	for(mouse in unique(m.sick$souris.id)){
		vec=rep(NA,264)
		names(vec)=1:263
		colnames(df)=names(vec)
		sub = subset(m.sick, souris.id==mouse)
		sums=c(sums, sum(as.numeric(sub$value)))
		vars=c(vars, var(as.numeric(sub$value)))
		sub$variable=(sub$variable-min(sub$variable))+1
		vec[sub$variable]=sub$value
		df=rbind(df,vec)
	
	}
	df=df[2:nrow(df),]
	df=data.frame(cbind(unique(m.sick$souris.id),df))
	colnames(df)=gsub('X','',colnames(df))
	df$souris = unlist(strsplit(as.character(df[,1]),'\\.'))[seq(1,nrow(df)*3,3)]
	df$sex = unlist(strsplit(as.character(df[,1]),'\\.'))[seq(2,nrow(df)*3,3)]
	df$condition = unlist(strsplit(as.character(df[,1]),'\\.'))[seq(3,nrow(df)*3,3)]
	df=df[,c('souris','sex','condition',1:263)]
	df$max_score = as.numeric(apply(df[,4:266],1,FUN=max,na.rm=T))
	df$age.at.onset = ages
	rownames(df)=paste(df[,1],df[,2],df[,3],sep='.')
	hc=hclust(dist(df[,as.character(1:263)]))
	df$cluster = cutree(hc,k=k)
	df$var = vars
	df$score_sum=sums
	#df[is.na(df)]=5
	return(df)
}

#'  Plot the means along a period of time 
#' 	@param data Do you love cats? Defaults to TRUE.
#' 	@param days Do you love cats? Defaults to TRUE.
#' 	@param min.score Do you love cats? Defaults to TRUE.
#' 	@keywords plot.mean
#' 	@export
#' 	@examples
#' plot.mean()
plot.mean<-function(data,days=1:263,min.score=0){
	data=subset(data, max_score>min.score)
	title='Mean clinical score'
	m=melt(data[,c('condition','sex','souris',days)],id.vars=c('condition','sex','souris'))
	m$value=as.numeric(as.character(m$value))
	if(min.score>0){
		title = paste0(title,' for mice that reach at least ',min.score)
	}
	a=aggregate(m$value,by=list(m$souris,m$condition,m$sex),FUN=mean,na.rm=T)	 
	g=ggplot(a,aes(Group.3, x,fill=Group.2))+geom_point(stroke=1,size=3,shape=21,position = position_jitterdodge(jitter.width=.1))+geom_boxplot(size=1,alpha=.6,outlier.shape=NA)+theme_bw()+scale_fill_manual(values=c('firebrick3','grey80'))+ylab(paste0('Mean clinical score: days ',days[1],'-',days[length(days)]))+ggtitle(title)+xlab('')
	 return(g)
}

#'  Plot the curves for the groups by mice 
#' 	@param data Do you love cats? Defaults to TRUE.
#' 	@param group Do you love cats? Defaults to TRUE.
#' 	@param min.score Do you love cats? Defaults to TRUE.
#' 	@param days Do you love cats? Defaults to TRUE.
#' 	@keywords plot.curves
#' 	@export
#' 	@examples
#' plot.curves()
plot.curves<-function(data,group,min.score=0,days=1:263){
	data=subset(data, max_score>min.score)
	data = data[,c('condition','sex','souris',as.character(days))]
	m=melt(data,id.vars=c('condition','sex','souris'))
	m$value=as.numeric(as.character(m$value))
	a=aggregate(m$value,by=list(m$souris,m$condition,m$sex,m$variable),FUN=mean,na.rm=T)
	m$variable=as.numeric(as.character(m$variable))
	m$value=as.numeric(as.character(m$value))
	a$var=paste0(a[,3],a[,2])
	a$Group.4=as.numeric(as.character(a$Group.4))	
	g=ggplot(subset(a,var %in% group),aes(Group.4,x))+geom_path(aes(group=Group.1,col=Group.1))+facet_wrap(Group.1~.,ncol=4)+theme_bw()+geom_vline(xintercept=(1:(263/50)*50),linetype='dashed',size=.2)+theme(strip.text=element_text(size=5))+guides(col=F)+xlab('Days')+ylab('Clinical Score')
	return(g)
}

#'  Plot the mean curves for the groups 
#' 	@param dataf Do you love cats? Defaults to TRUE.
#' 	@param days Do you love cats? Defaults to TRUE.
#' 	@param show Do you love cats? Defaults to TRUE.
#' 	@param min.score Do you love cats? Defaults to TRUE.
#' 	@param by Do you love cats? Defaults to TRUE.
#' 	@keywords plot.curve
#' 	@export
#' 	@examples
#' plot.curve()
plot.curve<-function(dataf,days=1:263,show='None',min.score=0,by='sex'){
	if(by=='sex'){
		col='condition'
	}else{
		col='sex'
	}
	dataf=subset(dataf, max_score>min.score)
	dataf$group = paste0(dataf$sex, dataf$condition)
	first=T
	survival = list()
	for(g in rev(unique(dataf$group))){
		sub = subset(dataf, group==g)
		for(i in 1:263){
			t = as.numeric(as.character(sum(na.omit(as.numeric(!is.na(sub[,as.character(i)]))))/nrow(sub)))
			vec=c(t, i,g)
			if(first){
				first=F 
				df = vec
			}else{
				df=rbind(df,  vec)
			}
		}
	}
	m=melt(dataf[,c('condition','sex','souris',days)],id.vars=c('condition','sex','souris'))
	m$value=as.numeric(as.character(m$value))
	m$variable=as.numeric(as.character(m$variable))
	m$value=as.numeric(as.character(m$value))
	m$var=paste0(m[,2],m[,1])
	m$souris.id = paste(m$souris,m$sex,m$condition,sep='.')

	a=aggregate(m$value,by=list(m$condition,m$sex,m$variable),FUN=mean,na.rm=T)
	a$key = paste0(a[,2],a[,1],a[,3])
	rownames(a)=a$key
	df=data.frame(df)
	rownames(df)=paste0(df[,3],df[,2])
	a$group=paste0(a[,1],a[,2])

	a$survival = as.numeric(as.character(df[rownames(a),1]))

	a$Group.3=as.numeric(as.character(a$Group.3))
	colnames(a)=c('condition','sex','day','score','key','group','survival')
	a$score = as.numeric(as.character(a$score))
	sums=c()
	a=a[order(a$group,a$day),]
	a=na.omit(a)
	for(souris in unique(a$group)){
		sub = a[which(a$group==souris),]
		sum.cum = 0
		for(day in 1:nrow(sub)){
				sum.cum=sum.cum+sub[day,'score']
				sums=c(sums, sum.cum)
			}
	}
	a$sums = sums
	if(col == 'sex'){
		cols=c('lightpink2','deepskyblue')	 
	}else{
		cols=c('darkorange','olivedrab2')
	}
	if(show=='None'){
	g=ggplot(a,aes(day))+geom_path(aes(group=group,y=score,col=a[,col]),size=1)+facet_grid(a[,by]~.)+theme_bw()+scale_color_manual(values=cols)+xlab('Days')+ylab('Clinical Score')
	}else if(show=='survival'){
	a[is.nan(a$score),'score']=NA
	scaleFactor <- max(na.omit(a$score)) / max(a$survival)
	g=ggplot(a,aes(day))+geom_path(aes(group=group,y=score,col=a[,col]),size=1)+facet_grid(a[,by]~.)+theme_bw()+scale_color_manual(values=cols)+xlab('Days')+ylab('Clinical Score')+geom_path(aes(y=survival*scaleFactor,linetype=a[,col]))+scale_y_continuous(name="score", sec.axis=sec_axis(~./scaleFactor, name="survival"))
	}else if(show=='sums'){
		a[is.nan(a$score),'score']=NA
		scaleFactor <- max(na.omit(a$score)) / max(a$sums)
		g=ggplot(a,aes(day))+geom_path(aes(group=group,y=score,col=a[,col]),size=1)+facet_grid(a[,by]~.)+theme_bw()+scale_color_manual(values=cols)+xlab('Days')+ylab('Clinical Score')+geom_path(aes(y=sums*scaleFactor,linetype=a[,col]))+scale_y_continuous(name="score", sec.axis=sec_axis(~./scaleFactor, name="cumsum"))
	}else if(show=='log.sums'){
		a[is.nan(a$score),'score']=NA
		scaleFactor <- max(na.omit(a$score)) / max(log1p(a$sums))
		g=ggplot(a,aes(day))+geom_path(aes(group=group,y=score,col=a[,col]),size=1)+facet_grid(a[,by]~.)+theme_bw()+scale_color_manual(values=cols)+xlab('Days')+ylab('Clinical Score')+geom_path(aes(y=log1p(sums)*scaleFactor,linetype=a[,col]))+scale_y_continuous(name="score", sec.axis=sec_axis(~./scaleFactor, name="cumsum"))
	}
	return(g)
}

#'  Compare the mean scores between 2 days 
#' 	@param data Do you love cats? Defaults to TRUE.
#' 	@param days Do you love cats? Defaults to TRUE.
#' 	@param min.score Do you love cats? Defaults to TRUE.
#' 	@param test Do you love cats? Defaults to TRUE.
#' 	@keywords stat.comp
#' 	@export
#' 	@examples
#' stat.comp()
stat.comp <- function(data, days=c(1,263),min.score=0,test='t.test'){
	data=subset(data, max_score>min.score)
	sub = data[,c('condition','sex','souris',days)]
	m = melt(data[,c('condition','sex','souris',days)],id.vars=c('condition','sex','souris'))
	m$variable = as.numeric(as.character(m$variable))
	m$value = as.numeric(as.character(m$value))
	g=ggplot(m,aes(factor(variable),value,fill=factor(variable)))+geom_point(position=position_jitter(width=.1))+geom_boxplot(alpha=.6,outlier.shape=NA,size=1)+stat_compare_means(method=test)+facet_grid(sex~condition)+theme_bw()+scale_fill_manual(values=c('grey80','grey50'))
	return(g)
}

#'  Compute AUC 
#' 	@param data Do you love cats? Defaults to TRUE.
#' 	@param days Do you love cats? Defaults to TRUE.
#' 	@param min.score Do you love cats? Defaults to TRUE.
#' 	@keywords find.auc
#' 	@export
#' 	@examples
#' find.auc()
find.auc <-function(data, days=1:263, min.score=0){

	data=subset(data, max_score>min.score)
	first=T
	for(i in as.character(unique(data$souris))){
		x=as.numeric(days[!is.na(data[which(data$souris==i),days+3])])
		y=as.numeric(data[which(data$souris==i),days+3][!is.na(data[which(data$souris==i),days+3])])
		auc=AUC(x,y)
		if(first){
			df=c(as.numeric(auc),i)
			first=F
		}else{
			df=rbind(df,c(as.numeric(auc),i))
		}
	}
	df=data.frame(df)
	return(df)
}

#'  Compare AUCs
#' 	@param data Do you love cats? Defaults to TRUE.
#' 	@param meta Do you love cats? Defaults to TRUE.
#' 	@param by Do you love cats? Defaults to TRUE.
#' 	@param test Do you love cats? Defaults to TRUE.
#' 	@keywords compare.auc
#' 	@export
#' 	@examples
#' compare.auc()
compare.auc<-function(dataf,meta,by='sex',test='t.test'){
	dataf$X1 = as.numeric(as.character(dataf$X1))
	rownames(dataf)=dataf[,2]
	meta = meta[rownames(dataf),]
	dataf$sex = meta[rownames(meta),'sex']
	dataf$condition = meta[rownames(meta),'condition']
	dataf$condition = factor(dataf$condition,levels=c('S','O'))
	if(by=='sex'){
		g=ggplot(dataf,aes(condition,X1,fill=condition))+geom_point(shape=21,size=3,stroke=1,position=position_jitter(width=.1))+geom_boxplot(alpha=.6,size=1)+facet_wrap(~sex)+stat_compare_means(method=test)+scale_fill_manual(values=c('firebrick3','grey80'))+theme_bw()+xlab('')+ylab('AUC')+guides(fill=F)
	}else{
		g=ggplot(dataf,aes(sex,X1,fill=sex))+geom_point(shape=21,size=3,stroke=1,position=position_jitter(width=.1))+geom_boxplot(alpha=.6,size=1)+facet_wrap(~condition)+stat_compare_means(method=test)+scale_fill_manual(values=c('lightpink2','deepskyblue'))+theme_bw()+xlab('')+ylab('AUC')+guides(fill=F)
	}
	return(g)
}

#'  Compute cumulative sums
#' 	@param data Do you love cats? Defaults to TRUE.
#' 	@param days Do you love cats? Defaults to TRUE.
#' 	@keywords cumsum
#' 	@export
#' 	@examples
#' cumsum()
cumsum <- function(data,days=1:263){
	m=melt(data[,c('condition','sex','souris',days)],id.vars=c('condition','sex','souris'))
	m=na.omit(m)
	m$key = paste(m$souris,m$sex,m$condition,m$variable,sep='.')
	m$souris.id = paste(m$souris,m$sex,m$condition,sep='.')
	
	vec=c()
	ages = c()
	vars = c()
	sums=c()
	m=m[!duplicated(m$key),]
	m=m[order(m$souris.id,m$variable),]
	sums = c()
	for(souris in unique(m$souris.id)){
		sub = m[which(m$souris.id==souris),]
		sum.cum = 0
		for(day in 1:nrow(sub)){
			sum.cum=sum.cum+sub[day,'value']
			sums=c(sums, sum.cum)
			}
	}
	m$sums=sums
	return(m)

}

#'  Return proportion of sick mice per day
#' 	@param data Do you love cats? Defaults to TRUE.
#' 	@param days Do you love cats? Defaults to TRUE.
#' 	@keywords sickMice.proportion
#' 	@export
#' 	@examples
#' sickMice.proportion()
sickMice.proportion <- function(data,days=1:263){
	first=T
	for(i in days){
		t = as.numeric(as.character(sum(na.omit(as.numeric(data[,as.character(i)]))>0|is.na(data[,as.character(i)]))/nrow(data)))
		vec=c(t, i)
		if(first){
			first=F 
			df = vec
		}else{
			df=rbind(df,  vec)
		}
	}
	return(df)
}

#'  Cluster mice based on cumulative sum of scores 
#' 	@param c Do you love cats? Defaults to TRUE.
#' 	@keywords cluster.on.sums
#' 	@export
#' 	@examples
#' cluster.on.sums()
cluster.on.sums <- function(c){
	first=T
	for(s in unique(c$souris.id)){
		vals=c()
		for(i in unique(c$souris.id)){
			k=ks.test(subset(c,souris.id==s)$sums,subset(c,souris.id==i)$sums)
			vals = c(vals, -log10(k$p.value))
			}
				names(vals) = unique(c$souris.id)
	
		if(first){
			first=F
			df = vals
		}else{
			df=rbind(df,vals)
		}
	}
	pheatmap(df)
}