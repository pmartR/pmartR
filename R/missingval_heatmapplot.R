
#here we will make a temporary function called missingval_heatmapplot, it takes in omicsData and returns a heatmap of omicsData$e_data

missingval_heatmapplot<- function(omicsData){

#here we call 'missingval_result' function that will give us an object of type naRes containing informaiton on the number of missing values per molecule   
na_Res<- missingval_result(omicsData)   
by_molecule<- na_Res$na.by.molecule
by_molecule<- by_molecule[order(by_molecule$num_NA),]

edata_cname<- attr(omicsData, "cnames")$edata_cname 

e_data<- omicsData$e_data

edata_melt<- melt(e_data, id.vars = edata_cname)
edata_melt[[edata_cname]]<- factor(edata_melt[[edata_cname]], levels = rev(by_molecule[[edata_cname]]))  
names(edata_melt)[1]<- "edata_cname"

p <- ggplot(edata_melt, aes(x=variable, y = edata_cname)) + geom_tile(aes(fill = value), colour = "white") + 
     theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + ylab("Molecule")
  
return(p)

}
