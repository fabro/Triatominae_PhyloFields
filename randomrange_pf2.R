### Function to create random phylogenetic fields. Sample sp's assemblages from the global PAM (rows are sites, columns are species), keeping the observed sp's assemblage richness (# of spp coexisting within its range)
#### and calculate the phylogenetic field (PSV) for each. Save all this on a table to compare with observed values.
# PAM: presence-absence matrix, sites are rows and species are columns
# spp.divfields: object with three columns derived with the function divfield_pams()
# patch: file path to the folder where the individual species assemblages are stored (results from 'sp_assem' function)

## In this restricted version, considered the range size of spp as its probability of being chosen (larger range, more likely to be chosen)
### It's necessary to create a probability vector (according to range sizes) and "feed" that into the function. Important: the spp have to be in the same order as the columns in the global PAM

randomrange_phylofields <- function(PAM, spp.divfields, path, tree, sims){


# create a matrix to save results
sppnames<-rownames(spp.divfields)

random_phylfields <- matrix(0, nrow=sims, ncol=length(sppnames))
colnames(random_phylfields) <- sppnames  

	for (j in 1:length(sppnames)){
	
	  rich <- spp.divfields[j,2]
	  
	  if(rich<3){
	    
	    #random_phylfields <- 0
      j=j+1
	    
	  } else {
	    
	    PAM.nofocal <- as.matrix(PAM[,-c(j)])
	    sp.focal <- as.matrix(PAM[1,j])
      
	    # create the probability vector
	    ranges <- as.matrix(colSums(PAM.nofocal))
	    sites <- nrow(PAM.nofocal)
	    prob_spp <- as.vector(ranges/sites)
    
		for (i in 1:sims){
			# read the observed assemblage
		  sp_assemb <- read.table(paste(path,sppnames[j],'_assemblage','.txt', sep=""), h=T)
			# create a random assemblage, keeping the observed spp richness
			sample_assemb <- as.matrix(sample(PAM.nofocal[1,],(rich-1),rep=F,prob=prob_spp))
			sample_assemb <- t(sample_assemb)
			sample_assembfocal <- as.data.frame(cbind(sample_assemb,t(sp.focal)))
			sample_assembfocal[1,] <- 1
			# calculate the PSV of the random assemblage
			random_phylfields[i,j] <- psv_modified(sample_assembfocal,tree)  
	

		}
		
	
	}	
  
	}
	
	write.table(random_phylfields, paste("random_phylofields_constrained.txt", sep=""), sep="\t")
	return(random_phylfields)

}