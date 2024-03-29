#enhancerPromoter.R

#plotting functions for enhancerPromoter code

#=========================================================
#===========================JOB PARAMETERS================
#=========================================================


args = commandArgs()
print(args)
geneTablePath = args[6]
outputFolder = args[7]
analysisName = args[8]
top = as.numeric(args[9])

#=========================================================
#===========================DEBUG SECTION=================
#=========================================================


# setwd('~/Dropbox/mycn_cyl')

# #needed arguments

# #gene table
# geneTablePath = './enhancerPromoter/HG19_P493-6_MYC_T24_-0_+0/HG19_P493-6_MYC_T24_-0_+0_GENE_TABLE.txt'

# #output folder
# outputFolder = './enhancerPromoter/HG19_P493-6_MYC_T24_-0_+0/output/'

# top = 5000


# analysisName = 'P493-6_MYC_T24'

#=========================================================
#========================HELPER FUNCTIONS=================
#=========================================================


plotContribution <- function(geneTable,analysisName,outputFolder,top=0,nBins=100){
	
	if(top == 0){
		top = nrow(geneTable)
		topString = 'all'
	}
	else if(top > nrow(geneTable)){
		top = nrow(geneTable)
		topString = 'all'			
	}else{
		topString = as.character(top)
	}

	#get the signal for tss, distal and total
	promoterSignal =geneTable[,2]
	enhancerSignal = geneTable[,3]
	totalSignal = promoterSignal + enhancerSignal

	#order by total signal
	totalOrder = order(totalSignal,decreasing=FALSE)[(length(totalSignal)-top+1):length(totalSignal)]


	#do a simpile moving average w/ a step size set by number of bins
	enhancerVector = c()
	promoterVector = c()
	
	i = 1
	stepSize = length(totalOrder)/nBins
	
	while(i < (length(totalOrder) - stepSize)){
	enhancerVector = c(enhancerVector,mean(enhancerSignal[totalOrder][i:(i+stepSize)]))
	promoterVector = c(promoterVector,mean(promoterSignal[totalOrder][i:(i+stepSize)]))
	i = i + stepSize/2
	}
	
	#check for crazy outliers
	if(max(totalSignal)/quantile(totalSignal,.99) > 2){
		yMax = as.numeric(quantile(totalSignal,.99)*1.5)
	}else{
		yMax = max(totalSignal)
	}
	
	#set the linewidth 
	linewidth = max(0.001,100/length(totalOrder))
	print(linewidth)
	
	plotPath = paste(outputFolder,analysisName,'_TOP_', topString,'_ORDERED_CONTRIB.pdf',sep='')
	pdf(file=plotPath,width = 11,height=8.5)
	par(mfrow=c(2,1))
	#we want to plot the bigger factor as the background
	if(max(enhancerSignal) > max(promoterVector)){
		plot(totalSignal[totalOrder],type='h',ylim =c(0,yMax),col='blue',xlab='All genes ranked by total signal',main=analysisName,lwd=0.25)
		lines(promoterSignal[totalOrder],type='h',col=rgb(1,0,0,0.3),lwd= linewidth)
	}else{
		plot(totalSignal[totalOrder],type='h',ylim =c(0,yMax),col='red',xlab='All genes ranked by total signal',main=analysisName,lwd=0.25)
		lines(enhancerSignal[totalOrder],type='h',col=rgb(0,0,1,0.3),lwd= linewidth)		
	}
	legend(0,0.75* yMax,c('promoter contribution','enhancer contribution'),fill=c('red','blue'))			

	enhancerContribVector = enhancerVector/(promoterVector+enhancerVector)*100
	enhancerContribVector[is.na(enhancerContribVector)] <- 0
	plot(1:length(enhancerContribVector),enhancerContribVector,type='p',col='blue',pch=16,ylab='% enhancer contribution to total signal',xaxt='n',xlab='',ylim =c(0,max(enhancerContribVector*1.15)))
	x=1:length(enhancerContribVector)
	lw1 = loess(enhancerContribVector ~x)
	lines(x,lw1$fitted,col='blue',lwd=2)
	dev.off()
	
	
}

addTicks <- function(plotTable,geneList){
	geneTicks = c()
	geneNames = c()
	for(gene in geneList){
		xIndex = which(plotTable[,1]==gene)
		if(length(xIndex) == 1){
			geneTicks = c(geneTicks,xIndex)
			geneNames = c(geneNames,gene)
		}
	}
	axis(3,geneTicks,geneNames,las=2)	
	
}

runWaterfall <- function(geneTable,analysisName,outputFolder,top=0,geneList = c()){
	if(top == 0){
		top = nrow(geneTable)
		topString = 'all'
	}else if(top > nrow(geneTable)){
		top = nrow(geneTable)
		topString = 'all'			
	}else{
		topString = as.character(top)
	}
	#get the signal for tss, distal and total
	promoterSignal =geneTable[,2]
	enhancerSignal = geneTable[,3]
	totalSignal = promoterSignal + enhancerSignal

	#order by total signal and get top N
	totalOrder = order(totalSignal,decreasing=FALSE)[(length(totalSignal)-top+1):length(totalSignal)]

	topTable = geneTable[totalOrder,]
	
	#get the % enhancer contribution
	topEnhancerContrib = topTable[,3]/(topTable[,2]+topTable[,3])

	#now in log2 form w/ limits at +/-9
	topEnhancerContrib_log = log2(topTable[,3]/topTable[,2])
	topEnhancerContrib_log[which(topEnhancerContrib_log== -Inf)] <- -8
	topEnhancerContrib_log[which(topEnhancerContrib_log== +Inf)] <- +8
	topEnhancerContrib_log[which(topEnhancerContrib_log < -8)] <- -8
	topEnhancerContrib_log[which(topEnhancerContrib_log > 8)] <- 8	
	topEnhancerContrib_log[is.na(topEnhancerContrib_log)] <- 1	

	#now get the difference
	topDiffVector = topTable[,3]-topTable[,2]
	topDiffOrder = order(topDiffVector,decreasing=FALSE)
	
	#enhancerContribOrder = order(topEnhancerContrib)
	topContribOrderedTable = cbind(topTable[topDiffOrder,],topDiffVector[topDiffOrder],topEnhancerContrib[topDiffOrder], topEnhancerContrib_log[topDiffOrder])
	colnames(topContribOrderedTable)[4] = 'ENHANCER_PROMOTER_DIFFERENCE'
	colnames(topContribOrderedTable)[5] = 'ENHANCER_CONTRIBUTION'
	colnames(topContribOrderedTable)[6] = 'ENHANCER_PROMOTER_RATIO'

	#first plot by contribution linear
	plotPath = paste(outputFolder,analysisName,'_TOP_', topString,'_WATERFALL_CONTRIBUTION.pdf',sep='')
	pdf(file=plotPath,width = 8,height =6)
	
	enhancerContribOrder = order(topContribOrderedTable[,5])
	plot(1:top,1-topContribOrderedTable[enhancerContribOrder,5],type='h',col='red',ylim=c(-1,1),ylab='Relative enhancer/promoter contribution',xlab=paste('Top',top,'genes as ranked by total signal',sep=' '))
	lines(-1*topContribOrderedTable[enhancerContribOrder,5],type='h',col='blue')
	
	if(length(geneList)>0){
		addTicks(topContribOrderedTable[enhancerContribOrder,],geneList)
	}
	
	dev.off()
	
	#plotting the log2 waterfall
	#setting the color for the log2 waterfall
	colorSpectrum <- colorRampPalette(c("red","grey","grey","blue"))(100)

	#setting a color data range
	minValue <- -8
	maxValue <- 8
	color_cuts <- seq(minValue,maxValue,length=100)
	color_cuts <- c(min(topEnhancerContrib_log,na.rm=TRUE), color_cuts,max(topEnhancerContrib_log,na.rm=TRUE))


	#add one extra min color to even out sampling
	colorSpectrum <- c(colorSpectrum[1],colorSpectrum[1],colorSpectrum)

	colorVector = c()
	for(i in enhancerContribOrder){
		delta = topContribOrderedTable[i,6]
		color = colorSpectrum[max(which(color_cuts <= delta))]
		colorVector =c(colorVector,color)	
	}	
	plotPath = paste(outputFolder,analysisName,'_TOP_', topString,'_WATERFALL_LOG.pdf',sep='')
	pdf(file=plotPath,width = 8,height =6)
	plot(1:top,topContribOrderedTable[enhancerContribOrder,6],type='h',col=colorVector,ylim=c(-8,8),ylab='log2 enhancer/promoter ratio',xlab=paste('Top',top,'genes as ranked by total signal',sep=' '))
	
	if(length(geneList)>0){
		addTicks(topContribOrderedTable[enhancerContribOrder,],geneList)
	}
	dev.off()	
	
	#now plot by difference
	#flipping this plot for aesthetic purposes
	colorVector = rep('red',nrow(topContribOrderedTable))
	colorVector[which(topContribOrderedTable[,4]>0)] <- 'blue'
	yMin = -1*quantile(topContribOrderedTable[,4],probs=(0.99))
	yMax = -1*quantile(topContribOrderedTable[,4],probs=(0.01))
	
	plotPath = paste(outputFolder,analysisName,'_TOP_', topString,'_WATERFALL_DIFFERENCE.pdf',sep='')
	pdf(file=plotPath,width = 8,height =6)
	plot(-1*topContribOrderedTable[,4],type='h',col=colorVector,ylim = c(yMin,yMax),ylab='Net promoter - enhancer signal (total rpm)',xlab=paste('Top',top,'genes as ranked by promoter enhancer difference'))
	if(length(geneList)>0){
		addTicks(topContribOrderedTable[enhancerContribOrder,],geneList)
	}
	dev.off()	

	#writing out the rank ordered table
	tablePath = paste(outputFolder,analysisName,'_TOP_', topString,'_ORDERED.txt',sep='')
	write.table(topContribOrderedTable,file= tablePath,quote=FALSE,row.names=FALSE,sep='\t')
	
	#making the gct
	filename_gct= paste(outputFolder,analysisName,'_top_', topString,'.gct',sep='')
	gctMatrix =matrix(ncol=4,nrow=nrow(topContribOrderedTable))
	colnames(gctMatrix) = c('NAME','DESCRIPTION','PROMOTER','DISTAL')
	gctMatrix[,1]= as.character(topContribOrderedTable[,1])
	gctMatrix[,3]= topContribOrderedTable[,2]
	gctMatrix[,4]= topContribOrderedTable[,3]
	
	gctHeader = matrix(data='',ncol=4,nrow=3)
	gctHeader[1,1]='#1.2'
	gctHeader[2,1]=nrow(topContribOrderedTable)
	gctHeader[2,2]='2'
	gctHeader[3,]=c('NAME','DESCRIPTION','PROMOTER','DISTAL')
	gctCombined = rbind(gctHeader,gctMatrix)
	write.table(gctCombined,file=filename_gct,quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)	

	
	#making the cls
	filename_cls= paste(outputFolder,analysisName,'_top_', topString,'.cls',sep='')
	clsTable = matrix(data='',ncol=3,nrow=3)
	clsTable[1,] =c(2,2,1)
	clsTable[2,1]=paste('#','PROMOTER',sep='')
	clsTable[2,2]='DISTAL'
	clsTable[3,1]='PROMOTER'
	clsTable[3,2]='DISTAL'
	print('WRITING .cls OUTPUT TO:')
	print(filename_cls)
	write.table(clsTable,file=filename_cls,quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)		


	#making the gct and cls for total contribution
	filename_gct= paste(outputFolder,analysisName,'_top_', topString,'_total_contrib.gct',sep='')
	gctMatrix =matrix(ncol=4,nrow=nrow(topContribOrderedTable))
	colnames(gctMatrix) = c('NAME','DESCRIPTION','BACKGROUND','SIGNAL')
	gctMatrix[,1]= as.character(topContribOrderedTable[,1])
	gctMatrix[,3]= as.numeric(rep(1,nrow(gctMatrix)))
	gctMatrix[,4]= topContribOrderedTable[,3] + topContribOrderedTable[,2]
	
	gctHeader = matrix(data='',ncol=4,nrow=3)
	gctHeader[1,1]='#1.2'
	gctHeader[2,1]=nrow(topContribOrderedTable)
	gctHeader[2,2]='2'
	gctHeader[3,]=c('NAME','DESCRIPTION','BACKGROUND','SIGNAL')
	gctCombined = rbind(gctHeader,gctMatrix)
	write.table(gctCombined,file=filename_gct,quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)	

	# making the cls
	filename_cls= paste(outputFolder,analysisName,'_top_', topString,'_total_contrib.cls',sep='')
	clsTable = matrix(data='',ncol=3,nrow=3)
	clsTable[1,] =c(2,2,1)
	clsTable[2,1]=paste('#','BACKGROUND',sep='')
	clsTable[2,2]='SIGNAL'
	clsTable[3,1]='BACKGROUND'
	clsTable[3,2]='SIGNAL'	
	write.table(clsTable,file=filename_cls,quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)


	# #making the gct and cls in log space <- no longer used in the analysis
	# #making the gct
	# filename_gct= paste(outputFolder,analysisName,'_top_', topString,'_log2.gct',sep='')
	# gctMatrix =matrix(ncol=4,nrow=nrow(topContribOrderedTable))
	# colnames(gctMatrix) = c('NAME','DESCRIPTION','PROMOTER','DISTAL')
	# for(i in 1:length(enhancerContribOrder)){
	# 	row = enhancerContribOrder[i]
	# 	gctMatrix[i,1] = as.character(topContribOrderedTable[row,1])
	# 	if(topContribOrderedTable[row,6] <=0){
	# 		gctMatrix[i,3] = 1/2^topContribOrderedTable[row,6]
	# 		gctMatrix[i,4] = 1
	# 	}else{
	# 		gctMatrix[i,3] = 1
	# 		gctMatrix[i,4] = 2^topContribOrderedTable[row,6]		
	# 	}
	# }
	# gctMatrix[,1]= as.character(topContribOrderedTable[enhancerContribOrder,1])
	# gctMatrix[,3]= topContribOrderedTable[,2]
	# gctMatrix[,4]= topContribOrderedTable[,3]
	
	# gctHeader = matrix(data='',ncol=4,nrow=3)
	# gctHeader[1,1]='#1.2'
	# gctHeader[2,1]=nrow(topContribOrderedTable)
	# gctHeader[2,2]='2'
	# gctHeader[3,]=c('NAME','DESCRIPTION','PROMOTER','DISTAL')
	# gctCombined = rbind(gctHeader,gctMatrix)
	# write.table(gctCombined,file=filename_gct,quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)	
	
	
	# #making the cls
	# filename_cls= paste(outputFolder,analysisName,'_top_', topString,'_log2.cls',sep='')
	# clsTable = matrix(data='',ncol=3,nrow=3)
	# clsTable[1,] =c(2,2,1)
	# clsTable[2,1]=paste('#','PROMOTER',sep='')
	# clsTable[2,2]='DISTAL'
	# clsTable[3,1]='PROMOTER'
	# clsTable[3,2]='DISTAL'
	# write.table(clsTable,file=filename_cls,quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)
	
	
	
}
#========================================================
#========================DATA INPUT======================
#========================================================

geneTable = read.delim(geneTablePath,sep='\t')




#========================================================
#=================PLOTTING CONTRIBUTION==================
#========================================================

#for all
plotContribution(geneTable,analysisName,outputFolder)
runWaterfall(geneTable,analysisName,outputFolder)


#top N
print('working on top genes')
print(top)
plotContribution(geneTable,analysisName,outputFolder,as.numeric(top)) 
runWaterfall(geneTable,analysisName,outputFolder,as.numeric(top))
