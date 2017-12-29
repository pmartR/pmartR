#msnset to pepData object function 

library("MSnbase")
data("msnset")

to_pepData<- function(msnset_object, edata_cname = "UniqueID", fdata_cname = "SampleID", emeta_cname = "UniqueID", data_scale = "abundance"){
  
  #check that msnset_object is of correct class
  if(class(msnset_object)!= "MSnSet") stop("msnset_object must be of class 'MSnSet'")
  
  msnset_edata<- msnset_object@assayData$exprs
  msnset_edata<- as.data.frame(msnset_edata)
  msnset_edata<- cbind(row.names(msnset_edata), msnset_edata)
  row.names(msnset_edata)<- NULL
  names(msnset_edata)[1]<- "UniqueID"
  
  msnset_fdata<- msnset_object@phenoData@data
  msnset_fdata<- as.data.frame(msnset_fdata)
  msnset_fdata<- cbind(row.names(msnset_fdata), msnset_fdata)
  row.names(msnset_fdata)<- NULL
  names(msnset_fdata)[1]<- "SampleID"
  
  
  msnset_emeta<- msnset_object@featureData@data
  msnset_emeta<- as.data.frame(msnset_emeta)
  msnset_emeta<- cbind(row.names(msnset_emeta), msnset_emeta)
  row.names(msnset_emeta)<- NULL
  names(msnset_emeta)[1]<- "UniqueID"
  
  res<- as.pepData(e_data = msnset_edata, f_data = msnset_fdata, e_meta = msnset_emeta, edata_cname = edata_cname, fdata_cname = fdata_cname, emeta_cname = emeta_cname, data_scale = data_scale)
  return(res)
  
}















