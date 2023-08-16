#' Creates a list of six Data Frames, one for each summarizing metric
#'
#' This function takes in an omicsData object and returns a summary of the
#' e_data component. The six summarizing metrics include the mean, standard
#' deviation, median, percent observed, minimum, and maximum.
#'
#' @param omicsData object of the class 'lipidData', 'metabData', 'pepData',
#'   'proData', or 'nmrData' created by \code{\link{as.lipidData}},
#'   \code{\link{as.metabData}}, \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.nmrData}}, respectively.
#' @param by character string indicating whether summarizing metrics will be
#'   applied by 'sample' or by 'molecule'. Defaults to 'sample'.
#' @param groupvar a character vector with no more than two variable names that
#'   should be used to determine group membership of samples. The variable name
#'   must match a column name from \code{f_data}. Defaults to NULL, in which
#'   case group_DF attribute will be used.
#'
#' @details If groupvar is NULL and group_designation has not been applied to
#'   omicsData, then the metrics will be applied to each column of e_data (when
#'   by = 'sample) or to each row of e_data (when by = 'molecule'). When
#'   groupvar is provided, it must match a column name from \code{f_data}, this
#'   column of f_data is used to group e_data in order to apply the metrics.
#'
#' @return A list of six data frames, of class 'dataRes' (data Result), which
#'   are the results of applying the metrics (mean, standard deviation, median,
#'   percent observed, minimum and maximum) to omicsData$e_data.
#'
#' @export
#'
#' @examples
#' library(pmartRdata)
#'
#' mylipid <- edata_transform(omicsData = lipid_pos_object, data_scale = "log2")
#' mylipid <- group_designation(omicsData = mylipid, main_effects = "Virus")
#' result <- edata_summary(omicsData = mylipid, by = "sample", groupvar = NULL)
#'
edata_summary <- function(omicsData, by = 'sample', groupvar = NULL) {
  # some checks
  if (!inherits(omicsData, c(
    'pepData', 'proData', 'lipidData',
    'metabData', 'nmrData', 'seqData'
  )))
    stop("omicsData must be an object of class pepData, proData, lipidData, metabData, nmrData, or seqData")
  if (!isTRUE(by %in% c('sample', 'molecule')))
    stop("by must be either sample or molecule")
  if (isTRUE(groupvar == attr(omicsData, "cnames")$fdata_cname))
    stop("The sample ID column in f_data cannot be used as a grouping column. Specify by = 'sample' to see a by-sample summary of the data")
  if (!all(groupvar %in% names(omicsData$f_data))) {
    # Why do you think you can group by variables that do not exist?!
    stop("The variable(s) in groupvar are not present in f_data.")
  }

  # pull cnames attr from omicsData
  edata = omicsData$e_data
  fdata = omicsData$f_data
  edata_cname = attr(omicsData, "cnames")$edata_cname
  fdata_cname = attr(omicsData, "cnames")$fdata_cname
  edata_cname_id = which(names(edata) == edata_cname)
  fdata_cname_id = which(names(fdata) == fdata_cname)
  groupDF = attr(omicsData, "group_DF")

  # all groupvars must be present in fdata
  if (any((groupvar %in% names(fdata)) == F)) {
    stop(paste0(
      "The following variables were not present in f_data: '",
      groupvar[which((groupvar %in% names(fdata)) == F)],
      "'"
    ))
  }
  if (length(groupvar) > 2) stop("No more than two groupvar can be provided")

  if (by == 'sample') {
    # check that groupvar is NULL, groupvar is only used when by == 'molecule'
    if (!is.null(groupvar)) stop("groupvar is only used when by == 'molecule'")

    avg = as.data.frame(apply(
      edata[, -edata_cname_id], 2,
      function(x) {
        if (all(is.na(x))) {
          mean(x)
        } else {
          mean(x, na.rm = TRUE)
        }
      }
    ))
    avg = cbind(names(edata[, -edata_cname_id]), avg)
    names(avg) <- c("sample", "mean")
    rownames(avg) <- NULL

    sd = as.data.frame(apply(edata[, -edata_cname_id], 2, sd, na.rm = TRUE))
    sd = cbind(names(edata[, -edata_cname_id]), sd)
    names(sd) <- c("sample", "sd")
    rownames(sd) <- NULL

    mds = as.data.frame(apply(edata[, -edata_cname_id], 2, median, na.rm = TRUE))
    mds = cbind(names(edata[, -edata_cname_id]), mds)
    names(mds) <- c("sample", "median")
    rownames(mds) <- NULL

    if (inherits(omicsData, "seqData")) {
      pct_obs = as.data.frame(apply(
        edata[, -edata_cname_id], 2,
        function(x) {
          sum(x != 0) / length(x)
        }
      ))
      pct_obs = cbind(names(edata[, -edata_cname_id]), pct_obs)

      names(pct_obs) <- c("sample", "pct_nonzero_obs")
      rownames(pct_obs) <- NULL
    } else {
      pct_obs = as.data.frame(apply(
        edata[, -edata_cname_id], 2,
        function(x) {
          sum(!is.na(x)) / length(x)
        }
      ))
      pct_obs = cbind(names(edata[, -edata_cname_id]), pct_obs)

      names(pct_obs) <- c("sample", "pct_obs")
      rownames(pct_obs) <- NULL
    }

    min = as.data.frame(apply(edata[, -edata_cname_id], 2, min, na.rm = TRUE))
    min = cbind(names(edata[, -edata_cname_id]), min)
    names(min) <- c("sample", "min")
    rownames(min) <- NULL

    max = as.data.frame(apply(edata[, -edata_cname_id], 2, max, na.rm = TRUE))
    max = cbind(names(edata[, -edata_cname_id]), max)
    names(max) <- c("sample", "max")
    rownames(max) <- NULL

    res_list = list(
      mean = avg,
      sd = sd,
      median = mds,
      pct_obs = pct_obs,
      min = min,
      max = max
    )
    class(res_list) <- "dataRes"
    attr(res_list, "by") <- by
    attr(res_list, "groupvar") <- groupvar
    attr(res_list, "cnames") <- list(
      "edata_cname" = edata_cname,
      "fdata_cname" = fdata_cname
    )
    attr(res_list, "data_scale") <- get_data_scale(omicsData)
  }

  if (by == "molecule") {
    if (is.null(groupvar)) {
      if (is.null(groupDF)) {
        # when groupvar is NULL and group_designation is NULL, we calculate
        # metric for each row
        avg = apply(
          edata[, -edata_cname_id], 1,
          function(x) {
            if (all(is.na(x))) {
              mean(x)
            } else {
              mean(x, na.rm = TRUE)
            }
          }
        )
        avg = data.frame(
          molecule = edata[, edata_cname_id],
          mean = avg,
          stringsAsFactors = F
        )
        names(avg)[1] <- edata_cname
        sd = apply(edata[, -edata_cname_id], 1, sd, na.rm = TRUE)
        sd = data.frame(
          molecule = edata[, edata_cname_id],
          sd = sd,
          stringsAsFactors = F
        )
        names(sd)[1] <- edata_cname
        mds = apply(edata[, -edata_cname_id], 1, median, na.rm = TRUE)
        mds = data.frame(
          molecule = edata[, edata_cname_id],
          median = mds,
          stringsAsFactors = F
        )
        names(mds)[1] <- edata_cname
        pct_obs = apply(
          edata[, -edata_cname_id], 1,
          function(x) {
            sum(!is.na(x)) / length(x)
          }
        )
        pct_obs = data.frame(
          molecule = edata[, edata_cname_id],
          pct_obs = pct_obs,
          stringsAsFactors = F
        )
        names(pct_obs)[1] <- edata_cname
        min = apply(edata[, -edata_cname_id], 1, min, na.rm = TRUE)
        min = data.frame(
          molecule = edata[, edata_cname_id],
          min = min,
          stringsAsFactors = F
        )
        names(min)[1] <- edata_cname
        max = apply(edata[, -edata_cname_id], 1, max, na.rm = TRUE)
        max = data.frame(
          molecule = edata[, edata_cname_id],
          max = max,
          stringsAsFactors = F
        )
        names(max)[1] <- edata_cname

        res_list = list(
          mean = avg,
          sd = sd,
          median = mds,
          pct_obs = pct_obs,
          min = min,
          max = max
        )
        class(res_list) <- "dataRes"
        attr(res_list, "by") <- by
        attr(res_list, "groupvar") <- groupvar
        attr(res_list, "cnames") <- list(
          "edata_cname" = edata_cname,
          "fdata_cname" = fdata_cname
        )
        attr(res_list, "group_DF") <- groupDF
        attr(res_list, "data_scale") <- get_data_scale(omicsData)
      } else if (!is.null(groupDF)) {
        # when groupvar is NULL but group_designation has been run

        # check that there are atleast 2 samples in each group and remove groups
        # that have less than two samples per group
        n_per_grp = as.data.frame(groupDF %>%
          dplyr::group_by(Group) %>%
          dplyr::summarise(count = dplyr::n()))
        remove_group = as.character(
          n_per_grp[which(n_per_grp$count < 2), "Group"]
        )
        if (length(remove_group) > 0) {
          n_per_grp = n_per_grp[-which(n_per_grp$count < 2), ]
          groupDF = groupDF[-which(groupDF$Group %in% remove_group), ]
        }

        # rearranging edata
        edata_melt <- tidyr::pivot_longer(
          edata, -!!edata_cname,
          cols_vary = "slowest",
          names_to = fdata_cname,
          values_to = "value"
        )
        edata_melt <- merge.data.frame(edata_melt, groupDF, by = fdata_cname)
        edata_melt <- edata_melt[, -which(names(edata_melt) == fdata_cname)]

        # checking that n_per_grp group order matches that of edata_melt
        n_per_grp = n_per_grp[match(
          unique(edata_melt$Group),
          n_per_grp$Group
        ), ]

        dcast_fn <- function(fun) {
          edata_melt %>%
            dplyr::group_by(
              !!dplyr::sym(edata_cname), Group
            ) %>%
            dplyr::summarise(
              value = fun(value),
              .groups = "drop"
            ) %>%
            tidyr::pivot_wider(
              id_cols = !!edata_cname,
              names_from = Group,
              names_vary = "slowest",
              values_from = value
            ) %>%
            data.frame
        }
        
        avg <- dcast_fn(\(x) mean(x, na.rm = !all(is.na(x))))
        std_dev <- dcast_fn(\(x) sd(x, na.rm = TRUE))
        mds <- dcast_fn(\(x) median(x, na.rm = !all(is.na(x))))
        pct_obs <- dcast_fn(\(x) sum(!is.na(x)) / length(x))
        mins <- dcast_fn(\(x) min(x, na.rm = !all(is.na(x))))
        maxs <- dcast_fn(\(x) max(x, na.rm = !all(is.na(x))))

        res_list = list(
          n_per_grp = n_per_grp,
          mean = avg,
          sd = std_dev,
          median = mds,
          pct_obs = pct_obs,
          min = mins,
          max = maxs
        )
        class(res_list) <- "dataRes"
        attr(res_list, "by") <- by
        attr(res_list, "groupvar") <- groupvar
        attr(res_list, "cnames") <- list(
          "edata_cname" = edata_cname,
          "fdata_cname" = fdata_cname
        )
        attr(res_list, "group_DF") <- groupDF
        attr(res_list, "data_scale") <- get_data_scale(omicsData)
      }
    } else if (length(groupvar) == 1) {
      #### case where groupvar is provided and has length 1####
      temp_fdata = fdata[, c(fdata_cname_id, which(names(fdata) == groupvar))]
      names(temp_fdata)[2] <- "Group"

      # check that there are atleast 2 samples in each group and remove groups
      # that have less than two samples per group
      n_per_grp = as.data.frame(temp_fdata %>%
        dplyr::group_by(Group) %>%
        dplyr::summarise(count = dplyr::n()))
      remove_group = as.character(n_per_grp[which(n_per_grp$count < 2), "Group"])
      if (length(remove_group) > 0) {
        n_per_grp = n_per_grp[-which(n_per_grp$count < 2), ]
        temp_fdata = temp_fdata[-which(temp_fdata$Group %in% remove_group), ]
      }

      # check to see if grouping structure was 1 sample per group
      if (nrow(temp_fdata) == 0) stop("The grouping variable must assign more than 1 sample to at least one group level")

      # rearranging edata
      edata_melt <- tidyr::pivot_longer(
        omicsData$e_data, -!!edata_cname,
        cols_vary = "slowest",
        names_to = fdata_cname,
        values_to = "value"
      )
      edata_melt = merge.data.frame(edata_melt, temp_fdata, by = fdata_cname)
      edata_melt = edata_melt[, -which(names(edata_melt) == fdata_cname)]

      # checking that n_per_grp group order matches that of edata_melt
      n_per_grp = n_per_grp[match(unique(edata_melt$Group), n_per_grp$Group), ]

      dcast_fn <- function(fun) {
        edata_melt %>%
          dplyr::group_by(
            !!dplyr::sym(edata_cname), Group
          ) %>%
          dplyr::summarise(
            value = fun(value),
            .groups = "drop"
          ) %>%
          tidyr::pivot_wider(
            id_cols = !!edata_cname,
            names_from = Group,
            names_vary = "slowest",
            values_from = value
          ) %>%
          data.frame
      }
      
      avg <- dcast_fn(\(x) mean(x, na.rm = !all(is.na(x))))
      std_dev <- dcast_fn(\(x) sd(x, na.rm = TRUE))
      mds <- dcast_fn(\(x) median(x, na.rm = !all(is.na(x))))
      pct_obs <- dcast_fn(\(x) sum(!is.na(x)) / length(x))
      mins <- dcast_fn(\(x) min(x, na.rm = !all(is.na(x))))
      maxs <- dcast_fn(\(x) max(x, na.rm = !all(is.na(x))))

      res_list = list(
        n_per_grp = n_per_grp,
        mean = avg,
        sd = std_dev,
        median = mds,
        pct_obs = pct_obs,
        min = mins,
        max = maxs
      )
      class(res_list) <- "dataRes"
      attr(res_list, "by") <- by
      attr(res_list, "groupvar") <- groupvar
      attr(res_list, "cnames") <- list(
        "edata_cname" = edata_cname,
        "fdata_cname" = fdata_cname
      )
      attr(res_list, "data_scale") <- get_data_scale(omicsData)
    } else if (length(groupvar) == 2) {
      #### case where length of groupvar is 2####

      group_vars = fdata[names(fdata) %in% groupvar]
      group_vars <- apply(group_vars, 2, as.character)

      # create a group variable and paste grouvar levels together for samples #
      # samples with a value of NA for either groupvar will have a Group value
      # of NA #
      Group = rep(NA, nrow(fdata))

      # identify samples that will have a Group membership that is not missing #
      nonna.group = (!is.na(group_vars[, 1]) & !is.na(group_vars[, 2]))
      Group[nonna.group] = paste(as.character(group_vars[nonna.group, 1]),
        as.character(group_vars[nonna.group, 2]),
        sep = "_"
      )

      # create output formatted with first column being fdata_cname and second
      # column group id #
      output = data.frame(Sample.ID = fdata[, fdata_cname], Group = Group)
      names(output)[1] = fdata_cname

      # check that there are atleast 2 samples in each group and remove groups
      # that have less than two samples per group
      n_per_grp = as.data.frame(output %>%
        dplyr::group_by(Group) %>%
        dplyr::summarise(count = dplyr::n()))
      remove_group = as.character(n_per_grp[which(n_per_grp$count < 2), "Group"])
      if (length(remove_group) > 0) {
        n_per_grp = n_per_grp[-which(n_per_grp$count < 2), ]
        output = output[-which(output$Group %in% remove_group), ]
      }

      if (nrow(output) == 0) stop("The grouping structure must assign more than 1 sample to at least one group level")

      # groupvar was provided, rearranging edata
      edata_melt <- tidyr::pivot_longer(
        omicsData$e_data, -!!edata_cname,
        cols_vary = "slowest",
        names_to = fdata_cname,
        values_to = "value"
      )
      edata_melt = merge.data.frame(edata_melt, output, by = fdata_cname)
      edata_melt = edata_melt[, -which(names(edata_melt) == fdata_cname)]

      # checking that n_per_grp group order matches that of edata_melt
      n_per_grp = n_per_grp[match(unique(edata_melt$Group), n_per_grp$Group), ]

      dcast_fn <- function(fun) {
        edata_melt %>%
          dplyr::group_by(
            !!dplyr::sym(edata_cname), Group
          ) %>%
          dplyr::summarise(
            value = fun(value),
            .groups = "drop"
          ) %>%
          tidyr::pivot_wider(
            id_cols = !!edata_cname,
            names_from = Group,
            names_vary = "slowest",
            values_from = value
          ) %>%
          data.frame
      }
      
      avg <- dcast_fn(\(x) mean(x, na.rm = !all(is.na(x))))
      std_dev <- dcast_fn(\(x) sd(x, na.rm = TRUE))
      mds <- dcast_fn(\(x) median(x, na.rm = !all(is.na(x))))
      pct_obs <- dcast_fn(\(x) sum(!is.na(x)) / length(x))
      mins <- dcast_fn(\(x) min(x, na.rm = !all(is.na(x))))
      maxs <- dcast_fn(\(x) max(x, na.rm = !all(is.na(x))))

      res_list = list(
        n_per_grp = n_per_grp,
        mean = avg,
        sd = std_dev,
        median = mds,
        pct_obs = pct_obs,
        min = mins,
        max = maxs
      )
      class(res_list) <- "dataRes"
      attr(res_list, "by") <- by
      attr(res_list, "groupvar") <- groupvar
      attr(res_list, "cnames") <- list(
        "edata_cname" = edata_cname,
        "fdata_cname" = fdata_cname
      )
      attr(res_list, "data_scale") <- get_data_scale(omicsData)
    }
  }

  return(res_list)
}
