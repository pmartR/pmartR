rrollup_r <- function(d){
#### Perform R_Rollup ####
## store number of peptides ##
num_peps = nrow(d)

res = matrix(NA, nrow = 1, ncol =  ncol(d))
## if only 1 peptide, set the protein value to the peptide ##
if(num_peps==1){
	protein_val = unlist(d)
}else{
	## Step 1: Select Reference Peptide -- peptide with least amount of missing data ##
	na.cnt = apply(is.na(d),1,sum)
	least.na = which(na.cnt == min(na.cnt))
	
	## If tied, select one with highest median abundance##
	if(length(least.na)>1){
    	mds = apply(d,1,median,na.rm=T)[least.na]
		least.na = least.na[which(mds==max(mds))]		
	}
	prot_val = unlist(d[least.na,])

	## Step 2: Ratio all peptides to the reference.  Since the data is on the log scale, this is the difference ##
	scaling_factor = apply(matrix(prot_val, nrow = num_peps, ncol = ncol(d), byrow=T) - d,1,median,na.rm=T)

	## Step 3: Use the median of the ratio as a scaling factor for each peptide ##
	x_scaled = d + matrix(scaling_factor, nrow = num_peps, ncol = ncol(d))
	
	## Step 4: Set Abundance as Median Peptide Abundance ##
	protein_val = apply(x_scaled, 2, median, na.rm=T)
}

res[1,] = protein_val
res2 = data.frame(cbind(num_peps, res))
names(res2) = c("N_peptides", names(d))
return(res2)
}