library(ggfortify)
library(sva)
meta = read.table("meta.txt",head=T,row.names=1)
psi = read.table("vHIP_CSDS_raw_psi.txt.gz",head=T,row.names=1)
pdf("vHIP_CSDS_raw_pca.pdf")
pc<-prcomp(t(psi))
pca = autoplot(pc,data=meta,colour="Treatment")
plot(pca)
dev.off()
eigs<-pc$sdev^2 # calculate how much variance explained by each PC.
percVar<-eigs/sum(eigs)

##pvca
library(lme4)
    variance_1perc<-sum(percVar[percVar>0.01])
    pct_threshold = variance_1perc # Amount of variability desired to be explained by the principal components.  Set to match the results in book chapter and SAS code.  User can adjust this to a higher (>= 0.8) number but < 1.0

    # Load matrix data
    theDataMatrix<-psi
    dataRowN <- nrow(theDataMatrix)
    dataColN <- ncol(theDataMatrix)

    # Center the data (center rows)
    theDataMatrixCentered <- matrix(data = 0, nrow = dataRowN, ncol = dataColN)
    theDataMatrixCentered_transposed = apply(theDataMatrix, 1, scale, center = TRUE, scale = FALSE)
    theDataMatrixCentered = t(theDataMatrixCentered_transposed)

    # Compute correlation matrix
    theDataCor <- cor(theDataMatrixCentered)

    # meta info
    expDesignRowN <- nrow(meta)
    expDesignColN <- ncol(meta)
    myColNames <- names(meta)

    # Obtain eigenvalues
    eigenData <- eigen(theDataCor)
    eigenValues = eigenData$values
    ev_n <- length(eigenValues)
    eigenVectorsMatrix = eigenData$vectors
    eigenValuesSum = sum(eigenValues)
    percents_PCs = eigenValues /eigenValuesSum 

    # Merge experimental file and eigenvectors for n components
    my_counter_2 = 0
    my_sum_2 = 1
    for (i in ev_n:1){
        my_sum_2  = my_sum_2 - percents_PCs[i]
        if ((my_sum_2) <= pct_threshold ){
            my_counter_2 = my_counter_2 + 1
        }
    }
    if (my_counter_2 < 3){
        pc_n  = 3 # pc_n is the number of principal components to model
    }else {
        pc_n = my_counter_2 
    }
    pc_data_matrix <- matrix(data = 0, nrow = (expDesignRowN*pc_n), ncol = 1)
    mycounter = 0
    for (i in 1:pc_n){
        for (j in 1:expDesignRowN){
        mycounter <- mycounter + 1
            pc_data_matrix[mycounter,1] = eigenVectorsMatrix[j,i]
        }
    }
    AAA <- meta[rep(1:expDesignRowN,pc_n),]
    Data <- cbind(AAA,pc_data_matrix)

    # Edit these variables according to your factors
    variables <- c(colnames(meta))
    # Mixed linear model
    op <- options(warn = (-1))
    effects_n = expDesignColN + 1
    randomEffectsMatrix <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
    model.func <- c()
    index <- 1
    for (i in 1:length(variables)) {
        mod = paste("(1|", variables[i], ")", sep = "")
        model.func[index] = mod
        index = index + 1
    }
    function.mods <- paste(model.func, collapse = " + ")
    for (i in 1:pc_n) {
        y = (((i - 1) * expDesignRowN) + 1)
        funct <- paste("pc_data_matrix", function.mods, sep = " ~ ")
        Rm1ML <- lmer(funct, Data[y:(((i - 1) * expDesignRowN) + 
            expDesignRowN), ], REML = TRUE, control=lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"),verbose = FALSE, na.action = na.omit)
        randomEffects <- Rm1ML
        randomEffectsMatrix[i, ] <- c(unlist(VarCorr(Rm1ML)), 
            resid = sigma(Rm1ML)^2)
    }
    effectsNames <- c(names(getME(Rm1ML, "cnms")), "resid")

    # Standardize Variance
    randomEffectsMatrixStdze <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
    for (i in 1:pc_n){
        mySum = sum(randomEffectsMatrix[i,])
        for (j in 1:effects_n){
            randomEffectsMatrixStdze[i,j] = randomEffectsMatrix[i,j]/mySum	
        }
    }

    # Compute Weighted Proportions
    randomEffectsMatrixWtProp <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
    for (i in 1:pc_n){
        weight = eigenValues[i]/eigenValuesSum
        for (j in 1:effects_n){
            randomEffectsMatrixWtProp[i,j] = randomEffectsMatrixStdze[i,j]*weight
        }
    }

    # Compute Weighted Ave Proportions
    randomEffectsSums <- matrix(data = 0, nrow = 1, ncol = effects_n)
    randomEffectsSums <-colSums(randomEffectsMatrixWtProp)
    totalSum = sum(randomEffectsSums)
    randomEffectsMatrixWtAveProp <- matrix(data = 0, nrow = 1, ncol = effects_n)
    for (j in 1:effects_n){
        randomEffectsMatrixWtAveProp[j] = randomEffectsSums[j]/totalSum 	
        
    }

    # Draw
    #pdf("pvca.qn1.batch1.degradation1.sva.pdf")
    pdf("psi.pvca.raw.pdf",width=7,height=3)
    bp <- barplot(randomEffectsMatrixWtAveProp,  xlab = "Covariates", ylab = "Weighted average proportion variance", ylim= c(0,0.5),col = c("blue"), las=2, cex.axis=0.5, cex.lab=0.75)
    #axis(1, at = bp, labels = effectsNames, xlab = "Covariates", cex.axis = 0.5, las=2)  # vertical x-axis
    text(bp, par("usr")[3]-0.02, srt = 45, adj = 1,labels = effectsNames, xpd = TRUE,cex=0.5)  # rotate 45 x-axis

    values = randomEffectsMatrixWtAveProp
    new_values = round(values , 3)
    text(bp,randomEffectsMatrixWtAveProp,labels = new_values, pos=3, cex = 0.5) # place numbers on top of bars 
    dev.off()
	
##sva hidden factors
Lmbeta<- as.matrix(psi)
    #mod = model.matrix( ~ meta$BrainWeight + meta$RIN + meta$YearAutopsy + meta$Ethnicity + meta$SequencingPlatform + meta$TissueState + meta$Sex + meta$Hemisphere)
####select top 5 covariates or the proportion value of covariates > 0.1
mod = model.matrix( ~ meta$ + meta$ + meta$ + meta$ + meta$)
n.sv = num.sv(Lmbeta, mod)   # n.sv=30
svobj = sva(Lmbeta, mod, n.sv=n.sv)
modSv = cbind(mod, svobj$sv)
rownames(modSv)=rownames(meta)
modSv = modSv[,-1]
write.table(modSv,"group.txt",col.names=T,row.names=T,quote=F,sep="\t")
