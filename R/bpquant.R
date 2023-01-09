#' Runs BP-Quant
#'
#' Applies BP-Quant to a pepData object 
#'
#' @param statRes an object of the class 'statRes'
#' @param pepData an omicsData object of the class 'pepData' that includes the
#'   e_meta component
#' @param pi_not numeric value between 0 and 1 indicating the background
#'   probability/frequency of a zero signature
#' @param max_proteoforms a numeric value corresponding to the maximum threshold
#'   for the number of possible proteoforms
#'
#' @return a list of data frames, one for each unique protein. The data frames
#'   have three columns, a protein identifier, a peptide identifier, and a
#'   "ProteoformID". The class of this list is 'isoformRes'.
#'
#' @details The result of this function can be used as one the \code{isoformRes}
#'   input argument to \code{\link{protein_quant}}. The \code{bpquant} function
#'   itself operates as follows: The statRes object contains the signatures data
#'   frame, the pepData object is used for its e_meta data frame. Next the
#'   signatures data frame and e_meta are merged by their edata_cname (e.g.
#'   peptide identifier) columns, this new data frame called protein_sig_data
#'   will be input to bpquant_mod in a “foreach” statement. “Foreach” will
#'   subset protein_sig_data for each unique protein and apply bpquant_mod to
#'   each subset and store the results.
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#'
#' mypepData <- group_designation(omicsData = pep_object,
#'                                main_effects = c("Phenotype"))
#' mypepData = edata_transform(omicsData = mypepData, data_scale = "log2")
#'
#' imdanova_Filt <- imdanova_filter(omicsData = mypepData)
#' mypepData <- applyFilt(filter_object = imdanova_Filt,
#'                        omicsData = mypepData,
#'                        min_nonmiss_anova=2)
#'
#' imd_anova_res <- imd_anova(omicsData = mypepData,
#'                            test_method = 'combined',
#'                            pval_adjust_a ='bon',
#'                            pval_adjust_g = 'bon')
#'
#' result = bpquant(statRes = imd_anova_res, pepData = mypepData)
#'
#' }
#'
#' @export
#'
bpquant<- function (statRes, pepData, pi_not = .9,
                    max_proteoforms = 5, parallel = TRUE) {

  # some checks
  
  if(is.null(pepData$e_meta)) {
    stop("pepData object must contain e_meta element")
  }

  if (!is.numeric(pi_not)) {
    # Again! When will you look at the examples and follow them?
    stop ("pi_not must be numeric.")
  }
  
  if(pi_not < 0 || pi_not > 1) stop("pi_not must be between 0 and 1")

  if (!is.numeric(max_proteoforms)) {
    # Hmmmmm. I can totally see that a number should be input as a character.
    stop ("max_proteoforms must be numeric.")
  }
  
  if(max_proteoforms <= 0) stop("max_proteoforms must be 1 or greater")

  if(!inherits(statRes, "statRes"))
    stop("statRes must be an object of class statRes")
  if (!inherits(pepData, "pepData"))
    stop("pepData must be an object of class pepData")

  # pulling e_meta from pepData, aswell as cnames
  e_meta <- pepData$e_meta
  edata_cname <- attr(pepData, "cnames")$edata_cname
  emeta_cname <- attr(pepData, "cnames")$emeta_cname

  # Check if the bpFlags attribute is NULL. If it is then only the G-test was
  # run and the bpquant function should not be run with just G-test flags.
  if (is.null(attr(statRes, "bpFlags"))) {

    stop (paste("The bpquant function cannot be run with just G-test flags.",
                "Run the imd_anova function with either 'anova' or 'combined'",
                "as the input to test_method.",
                sep = " "))

  }

  # Pull signatures (which are the flags) from the statRes object.
  signatures <- attr(statRes, "bpFlags")

  #check that rows in signatures equals rows in e_data
  if(nrow(signatures) != nrow(pepData$e_data)) stop("rows in signatures must equal rows in e_data")

  #merging e_meta with signatures from statRes
  protein_sig_data<- merge(e_meta, signatures, by = edata_cname, all.x = TRUE)

  protein_sig_data<- as.data.frame(protein_sig_data)

  #removing any signatures with NA's
  protein_sig_data<- na.omit(protein_sig_data)

 #pull protein column from protein_sig_data and apply unique function
  unique_proteins<- unique(protein_sig_data[[emeta_cname]])

  # set up parallel backend
  if(parallel == TRUE){
    cores<- parallel::detectCores()
    cl<- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl))
    doParallel::registerDoParallel(cl)
  } else {
    foreach::registerDoSEQ()
  }

  isoformRes <- foreach::foreach(i=1:length(unique_proteins),
                                 .export = "bpquant_mod")%dopar%{

    row_ind<- which(protein_sig_data[, emeta_cname] == unique_proteins[i])
    cur_protein<- protein_sig_data[row_ind, ]
    cur_protein_sigs = cur_protein[which(names(cur_protein) %in% colnames(signatures[-which(colnames(signatures) == edata_cname)]))]
    colnames(cur_protein_sigs) = paste(rep("flags", ncol(cur_protein_sigs)), 1:ncol(cur_protein_sigs), sep = "")

    result<- bpquant_mod(protein_sig = cur_protein_sigs,
                                  pi_not = pi_not,
                                  max_proteoforms = max_proteoforms)

    peptide_id<- result$peptide_idx

    temp<- data.frame(Protein = as.character(cur_protein[, emeta_cname]),
                      Mass_Tag_ID = as.character(cur_protein[, edata_cname]),
                      proteoformID = peptide_id,
                      stringsAsFactors = FALSE)
    names(temp)<- c(emeta_cname, edata_cname, "proteoformID")

    temp

 }

 #filter out isoformRes with max(proteoformID) of zero
 list_inds = lapply(isoformRes,function(x){max(x$proteoformID)})
 inds = which(list_inds == 0)

 # only get rid of isoform results if we find something with no peptides #
 if(length(inds)>0){
 isoformRes = isoformRes[-inds]
 }

 isoformRes2<- lapply(isoformRes,
                      isoformRes_func,
                      emeta_cname = emeta_cname,
                      edata_cname = edata_cname)
 isoformRes2<- data.frame(data.table::rbindlist(isoformRes2))

 attr(isoformRes, "isoformRes_subset")<- isoformRes2

 class(isoformRes)<- "isoformRes"
 return(isoformRes)
}

#' bpquant_mod function
#'
#' The function is written to take input from one protein at a time and requires
#' three inputs: protein_sig, pi_not and max_proteforms
#'
#' @param protein_sig  is a matrix or data.frame with p rows and n columns,
#'   where p is the number of peptides mapped to the protein of interest and n
#'   is the number of tests conducted to generate signatures made up of values
#'   0, 1, and -1.
#' @param pi_not is a numeric value between 0 and 1 indicating the background
#'   probability/frequency of a zero signature.
#' @param max_proteoforms a numeric value, a maximum threshold for the number of
#'   possible proteoforms.
#'
#' @return  a list of five items: num_proteoforms, unique_sigs,
#'   proteoform_configs, post_prob and peptide_idx
#'
#' @details \tabular{ll}{ num_proteoforms \tab the number of proteoforms as
#' identified by bpquant \cr unique_sigs \tab matrix of unique signatures
#' observed  \cr proteoform_configs \tab matrix of 0/1 values indicating
#' scenarios of proteoform absence/presence scenarios \cr post_prob \tab vector
#' of posterior probabilities corresponding to each proteoform configuration in
#' "proteoform_configs" \cr peptide_idx \tab vector of 0, 1, 2, . . . values
#' indicating which proteoform each peptide belongs to \cr }
#'
#' @rdname bpquant_mod
#'
bpquant_mod <- function (protein_sig, pi_not, max_proteoforms) {

  ## implement some checks ##

  ## signatures can only contain values equal to -1, 1, or 0
  if(sum(apply(protein_sig, 1, function(x)
    sum(!(x %in% c(
      1, -1, 0
    )))))> 0)
    stop("Entries in the signatures matrix may only take values of -1, 1, or 0")
  

  ### generate parameters associated with signature counts and probabilities ###

  ## generate a list of unique signatures ##
  sigs = unique(protein_sig)

  ## store number of unique signatures ##
  uniq_sigs = nrow(sigs)

  #### count occurences of each signature ####
  ## save current signatures as character strings ##
  sig_str = apply(protein_sig,1,paste,collapse="")

  counts = NULL
  for(j in 1:uniq_sigs){
    ## count occurences of each unique signatures ##
    counts[j] = sum(sig_str==paste(sigs[j,],collapse=""))
  }

  ## order signatures by count##
  cnt.ord = counts[order(counts, decreasing=TRUE)]
  if(nrow(sigs) > 1){
    sig.ord = as.data.frame(sigs[order(counts, decreasing=TRUE),])
  }else{
    sig.ord = sigs
  }

  ## set zero signature first ##
  z.sig = which(apply(abs(sig.ord),1,sum)==0)
  if(length(z.sig) == 0){
    pi_probs = rep((1 - pi_not)/uniq_sigs, uniq_sigs)
  }else{
    sig.ord = rbind(sig.ord[z.sig,],sig.ord[-z.sig,])
    cnt.ord = c(cnt.ord[z.sig],cnt.ord[-z.sig])

    ## calculate prior probabilities associated with non-zero signatures ##
    pi_probs = rep((1 - pi_not)/(uniq_sigs - 1), uniq_sigs)
    pi_probs[1] = pi_not
  }

  ## identify which unique signature each row matches ##
  counts_ids = NULL
  for(j in 1:nrow(sig.ord)){
    counts_ids[sig_str==paste(sig.ord[j,],collapse="")] = j
  }

  ############## generate possible proteoforms ####################
  ## specify the number of unique signatures observed ##
  nu = length(pi_probs)

  #### build matrix of possible proteoforms ####
  ## number of possible proteoform presence/absence combinations ##

  ## if the number of possible proteoforms exceeds 5, then cap at 5 ##
  if(nu > max_proteoforms){nu = max_proteoforms}

  n_combos = 2^nu

  ## number of reps for each iteration ##
  n_reps = rep(1,nu)
  if(nu > 1){
    for(j in 2:nu){
      n_reps[j] = n_reps[j-1]*2
    }
  }

  ## construct matrix of all possible proteoform configurations ##
  mat = matrix(0, nrow = n_combos, ncol = nu)
  mat[,1] = rep(c(0,1),each = n_reps[nu])
  if(nu > 1){
    for(j in 2:nu){
      mat[,j] = rep(c(0,1),each = n_reps[nu-j + 1],length=n_combos)
    }
  }else{ mat = matrix(mat,nrow=2,ncol=1)}

  ## order configurations for readability ##
  ## rows are first ordered by sum (ascending)
  ## then first column value (descending), 
  ## second column value (descending), etc. ##

  if(nu == 1){p_configs = mat}else{
    mat.list = list()
    mat.list[[1]] = mat[apply(mat,1,sum)==0,]
    mat.list[[nu+1]] = mat[apply(mat,1,sum)==nu,]
    for(j in 1:(nu-1)){
      tmp.rws = mat[apply(mat,1,sum)==j,]
      ord.var = apply(tmp.rws,1,function(x) as.numeric(paste(x,collapse="")))
      mat.list[[j+1]] = tmp.rws[order(ord.var, decreasing=TRUE),]
    }
    p_configs = do.call(rbind,mat.list)
  }

  ############## calculate posterior probability vector ##################
  ## specify number of peptides observed ##
  n_peps = nrow(protein_sig)

  ## define matrix of Binomial CDF values. first column gives probability of
  ## observing proteoform at least as many times as it has been for this data.
  ## second column gives probability of observing proteoform fewer times than it
  ## has been for this data.
  x_mat = matrix(0, nrow=nu, ncol=2)
  for (j in 1:nu){
    x_mat[j,2] = pbinom(cnt.ord[j]-1, n_peps, pi_probs[j])
  }
  x_mat[,1] = 1 - x_mat[,2]

  ## calculate posterior distribution for each proteoform configuration##
  post_prob = rep(1, nrow(p_configs))
  for (j in 1:nrow(p_configs)){
    x = p_configs[j,]
    for(k in 1:nu){
      prior = (pi_probs[k]^x[k])*(1-pi_probs[k])^(1-x[k])
      post_prob[j] = post_prob[j]*prior*x_mat[k,x[k]+1]
    }
  }

  post_prob = post_prob/sum(post_prob)

  ## store which proteoform configuration has highest posterior probability ##
  id.max = which.max(post_prob)

  ## pull maximum configuration ##
  tmp = p_configs[id.max,]

  ## determine number of proteoforms and map to protein signatures in observed
  ## peptides if an all zero signature pull signatures of counts_ids==1
  if(sum(tmp)==0){
    peptide_ids = rep(0, length(counts_ids))
    peptide_ids[which(counts_ids==1)] = 1
    num_proteoforms = 1
  }else{
    num_proteoforms = sum(tmp)
    tid = which(tmp > 0)

    peptide_ids = rep(0,length(counts_ids))
    for(q in 1:num_proteoforms){
      peptide_ids[which(counts_ids == tid[q])] = q
    }
  }

  return(list(post_prob = post_prob,
              peptide_idx = peptide_ids,
              unique_sigs = sigs,
              num_proteoforms = num_proteoforms,
              proteoform_configs = p_configs))

}


# function to use with lapply on isoformRes
isoformRes_func<- function(df, emeta_cname, edata_cname){
  temp<- vector("list", max(df$proteoformID))

  for(i in 1:max(df$proteoformID)){
    cur_subset<- df[which(df$proteoformID == i), ]

    # append proteoform number if there's more than one #
    if(max(df$proteoformID) > 1){
    new_df<- data.frame(cur_subset[[emeta_cname]],
                        paste(cur_subset[[emeta_cname]],
                              cur_subset$proteoformID,
                              sep = ';'),
                        cur_subset[[edata_cname]])
    }else{
    new_df = data.frame(cur_subset[[emeta_cname]],
                        cur_subset[[emeta_cname]],
                        cur_subset[[edata_cname]])
    }
    names(new_df)<- c(names(df)[1], "Protein_Isoform", names(df)[2])
    temp[[i]]<- new_df
  }

  data.frame(data.table::rbindlist(temp))
}
