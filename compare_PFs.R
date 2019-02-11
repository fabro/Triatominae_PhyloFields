##### compare phylofields observed vs null ####
#### Compare observed values against nulls (created with the above functions) ####

compare_phylofields <- function(observed,nulls){
	#read the table with the observed values (remember that this table has two columns: PSV and SppRich)
	phylofields <- read.table(observed)
	#read the table with the null values
	null_phylofields <- read.table(nulls)
# 	#order (ascending) null values, to ease comparison later
# 	ordered_nulls <- apply(null_phylofields,2,sort)

	compare_table <- matrix(NA,nrow=nrow(phylofields),ncol=1)
	
	for (i in 1:nrow(phylofields)){
		#read first observed value (for the first species. Order corresponds to the global PAM)
		#check if sp's PSV == 1 and SppRich < 3. If so, no comparions can be made and leave NA
    if(phylofields[i,1]==1 & phylofields[i,2]<3){
      i <- i+1   
    } else {
      obs <- phylofields[i,1]
      #get the quartiles for comparison
      lower.quartile <- quantile(x = null_phylofields[,i],probs = 0.025)
      higher.quartile <- quantile(x = null_phylofields[,i],probs = 0.975)
      #check if observed value is < 0.25% or > 0.75%  of the null values (clustered or overdispersed, respectively)
      if (obs <= lower.quartile){
        
        compare_table[i,] <- "clustered"
      }
      else {
        if (obs >= higher.quartile){
          compare_table[i,] <- "overdispersed"
        }
        else {
          compare_table[i,] <- "random"
        }
      }
    }
    
		
	}
	rownames(compare_table) <- rownames(phylofields)
	write.table(compare_table, paste("compare_phylofields.txt", sep=""),sep="\t")
	
}