#### divfield_pams ####
# Function to generate species' diversity field pams. That is, get a presence-absence matrix of the species co-occurring within a focal sp's range
# diversity field pams are saved in the WD without keeping them in the R environment

divfield_pams <- function(pam_obs){
	
  sppnames<-colnames(pam_obs)
  spp.divfields <- matrix(NA,nrow=ncol(pam_obs),ncol=3)
  rownames(spp.divfields) <- sppnames
  colnames(spp.divfields) <- c("RangeSize","TotalRich", "MeanRich")
  
  RangeSize<- as.matrix(colSums(pam_obs))
  SpeciesRichness<- as.matrix(rowSums(pam_obs))
  D_Volume<- t(pam_obs)%*%SpeciesRichness
  RangeRichness<- D_Volume/RangeSize
  
	for (i in 1:ncol(pam_obs)){
		
		sp_sites <- which(pam_obs[,i]==1)
		
			if(length(sp_sites)>1){
				
					sp_pam <- pam_obs[sp_sites,which(colSums(pam_obs[sp_sites,])!=0)]
				
			} else {sp_pam <- pam_obs[sp_sites,which(pam_obs[sp_sites,]!=0)]
			  }
		
        if(is.null(colnames(sp_pam)) & is.null(names(sp_pam))){
          
          spp.divfields[i,1] <- RangeSize[i]
          spp.divfields[i,2] <- 1
          spp.divfields[i,3] <- RangeRichness[i]
          i<-i+1
          
        }	else {
          
          spp.divfields[i,1] <- RangeSize[i]
          
            if(is.null(dim(sp_pam))){
              spp.divfields[i,2] <- length(sp_pam)
            }else{spp.divfields[i,2] <- ncol(sp_pam)}
                    
          spp.divfields[i,3] <- RangeRichness[i]
          
          write.table(sp_pam, paste(sppnames[i],'.txt', sep=""), sep="\t")
          }
	}
	return(spp.divfields)
}


### sp_assemb function###
# Function to create and save a focal species assemblage (based on the sp_pams from the function 'divfield_pams'). objects are saved to the WD instead of keeping them in the R environment
# only saves assemblages for those species with total within-range richness > 1 (i.e. if a focal sp cooccurs with at least other sp)

sp_assemb <- function(spp.divfields){
	
  sppnames<-rownames(spp.divfields[which(spp.divfields[,2]>1),])
  
		for (i in 1:length(sppnames)){

			sp_pam <- read.table(paste(sppnames[i],'.txt', sep=""), h=T)
			
				if(ncol(sp_pam)==1){
					
					sp_pam <- t(sp_pam)
					write.table(sp_pam, paste(sppnames[i],'_assemblage','.txt', sep=""), sep="\t")
					
				} else	{ sp_assemblage <- sp_pam[1,]	
			
				sp_assemblage[1,] <- 1
			
			write.table(sp_assemblage, paste(sppnames[i],'_assemblage','.txt', sep=""), sep="\t")

			}	
		
		}
}

#### modified function "psv" from PICANTE to give a single value ####
# this function will be used within the phylofields_psv function below
psv_modified <- function (samp, tree, compute.var = TRUE) {

    psv_values <- as.matrix(0,nrow=1,ncol=3)

	samp[samp > 0] <- 1
    flag = 0
    if (is.null(dim(samp))) {
        samp <- rbind(samp, samp)
        flag = 2
    }
    if (is(tree)[1] == "phylo") {
        if (is.null(tree$edge.length)) {
            tree <- compute.brlen(tree, 1)
        }
        tree <- prune.sample(samp, tree)
        samp <- samp[, tree$tip.label]
        Cmatrix <- vcv.phylo(tree, cor = TRUE)
    }
    else {
        Cmatrix <- tree
        species <- colnames(samp)
        preval <- colSums(samp)/sum(samp)
        species <- species[preval > 0]
        Cmatrix <- Cmatrix[species, species]
        samp <- samp[, colnames(Cmatrix)]
    }
    SR <- rowSums(samp)
    nlocations <- dim(samp)[1]
    nspecies <- dim(samp)[2]
    PSVs <- NULL
    for (i in 1:nlocations) {
        index <- seq(1:nrow(Cmatrix))[samp[i, ] == 1]
        n <- length(index)
        if (n > 1) {
            C <- Cmatrix[index, index]
            PSV <- (n * sum(diag(as.matrix(C))) - sum(C))/(n * 
                (n - 1))
        }
        else {
            PSV <- NA
        }
        PSVs <- c(PSVs, PSV)
    }
    PSVout <- cbind(PSVs, SR)
    if (flag == 2) {
        PSVout <- PSVout[-2, ]
        return(PSVout)
    }
    else {
        if (compute.var == FALSE) {
            return(data.frame(PSVout))
        }
        else {
            PSVvar <- NULL
            X <- Cmatrix - (sum(sum(Cmatrix - diag(nspecies))))/(nspecies * 
                (nspecies - 1))
            X <- X - diag(diag(X))
            SS1 <- sum(X * X)/2
            SS2 <- 0
            for (i in 1:(nspecies - 1)) {
                for (j in (i + 1):nspecies) {
                  SS2 <- SS2 + X[i, j] * (sum(X[i, ]) - X[i, 
                    j])
                }
            }
            SS3 <- -SS1 - SS2
            S1 <- SS1 * 2/(nspecies * (nspecies - 1))
            S2 <- SS2 * 2/(nspecies * (nspecies - 1) * (nspecies - 
                2))
            if (nspecies == 3) {
                S3 <- 0
            }
            else {
                S3 <- SS3 * 2/(nspecies * (nspecies - 1) * (nspecies - 
                  2) * (nspecies - 3))
            }
            for (n in 2:nspecies) {
                approxi <- 2/(n * (n - 1)) * (S1 + (n - 2) * 
                  S2 + (n - 2) * (n - 3) * S3)
                PSVvar <- rbind(PSVvar, c(n, approxi))
            }
            vars <- rep(0, nlocations)
            PSVout <- cbind(PSVout, vars)
            for (g in 1:nlocations) {
                if (PSVout[g, 2] > 1) {
                  PSVout[g, 3] <- PSVvar[PSVout[g, 2] - 1, 2]
                }
                else {
                  PSVout[g, 3] <- NA
                }
            }
            #return(data.frame(PSVout))
			psv_values[1,1] <- PSVs # valor original es [1,1]
			#psv_values[1,2] <- ncol(samp)
			#psv_values[1,3] <- vars
			
			return(psv_values)
        }
    }
}



#### phylofields_psv function####
# Function to iterate the calculation of the phylogenetic field (single-value) of a species; using PSV (phylogenetic species variability) --function from PICANTE (needs to be loaded)

phylofields_psv <- function(spp.divfields, tree){

  sppnames<-rownames(spp.divfields[which(spp.divfields[,2]>1),])
  
  phylo_fields <- matrix(0,nrow=length(sppnames), ncol=2)
  colnames(phylo_fields) <- c("PSV","SppRichness")
  rownames(phylo_fields) <- sppnames
	
	for (i in 1:length(sppnames)){
		
		sp_assemblage <- read.table(paste(sppnames[i],'_assemblage', '.txt', sep=""), h=T)
		
		if (ncol(sp_assemblage)==1){
			
			phylo_fields[i,1] <- NA
			phylo_fields[i,2] <- ncol(sp_assemblage)
			
		} else {
				
		#sp_assemblage <- sp_pam[1,]
		#sp_assemblage[1,] <- 1
		
		phylo_fields[i,1] <- psv_modified(sp_assemblage,tree)
		phylo_fields[i,2] <- ncol(sp_assemblage)
				
		}
	
	write.table(phylo_fields, paste("phylo_fields.txt", sep=""), sep="\t")
	
	}
	return(phylo_fields)
}

