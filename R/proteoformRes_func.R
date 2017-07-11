# function to use with lapply on proteoformRes

proteoformRes_func<- function(df){
  temp<- vector("list", max(df$proteoformID))
  
  for(i in 1:max(df$proteoformID)){
    cur_subset<- df[which(df$proteoformID == i), ]
    new_df<- data.frame(cur_subset$Protein, paste(cur_subset$Protein, cur_subset$proteoformID, sep = ';'), as.numeric(cur_subset$Mass_Tag_ID))
    names(new_df)<- c(names(df)[1], "Protein_Isoform", names(df)[2])
    temp[[i]]<- new_df
  }
  
  do.call(rbind, temp)
}