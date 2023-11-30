#' print.moleculeFilt
#'
#' For printing an S3 object of type 'moleculeFilt'
#' @param x An object of type 'moleculeFilt'
#' @param ...  further arguments passed to or from other methods
#'
#' @return No return value, prints details about x
#'
#' @rdname print-moleculeFilt
#' @export
#'
print.moleculeFilt <- function(x, ...) {
  filter_object <- x

  if (!inherits(filter_object, "moleculeFilt")) stop("filter object must be of the class 'moleculeFilt'")
  filter_object <- as.data.frame(lapply(filter_object, as.character), stringsAsFactors = FALSE)

  filter_object_head = head(filter_object, 4)[, 1:ncol(filter_object)]
  filter_object_tail = tail(filter_object, 4)[, 1:ncol(filter_object)]
  blank_row = rep("---", ncol(filter_object))

  result <- rbind(filter_object_head, blank_row, filter_object_tail)

  cat("moleculeFilt object\n")
  cat(capture.output(result), sep = "\n")
  cat("\n")
}

#' print.totalCountFilt
#'
#' For printing an S3 object of type 'totalCountFilt'
#' @param x An object of type 'totalCountFilt'
#' @param ...  further arguments passed to or from other methods
#'
#' @return No return value, prints details about x
#'
#' @rdname print-totalCountFilt
#' @export
#'
print.totalCountFilt <- function(x, ...) {
  filter_object <- x

  if (!inherits(filter_object, "totalCountFilt")) stop("filter object must be of the class 'totalCountFilt'")
  filter_object <- as.data.frame(lapply(filter_object, as.character), stringsAsFactors = FALSE)

  filter_object_head = head(filter_object, 4)[, 1:ncol(filter_object)]
  filter_object_tail = tail(filter_object, 4)[, 1:ncol(filter_object)]
  blank_row = rep("---", ncol(filter_object))

  result <- rbind(filter_object_head, blank_row, filter_object_tail)

  cat("totalCountFilt object\n")
  cat(capture.output(result), sep = "\n")
  cat("\n")
}

#' print.RNAFilt
#'
#' For printing an S3 object of type 'RNAFilt'
#' @param x An object of type 'RNAFilt'
#' @param ...  further arguments passed to or from other methods
#'
#' @return No return value, prints details about x
#'
#' @rdname print-RNAFilt
#' @export
#'
print.RNAFilt <- function(x, ...) {
  filter_object <- x

  if (!inherits(filter_object, "RNAFilt")) stop("filter object must be of the class 'RNAFilt'")
  filter_object <- as.data.frame(lapply(filter_object, as.character), stringsAsFactors = FALSE)

  filter_object_head = head(filter_object, 4)[, 1:ncol(filter_object)]
  filter_object_tail = tail(filter_object, 4)[, 1:ncol(filter_object)]
  blank_row = rep("---", ncol(filter_object))

  result <- rbind(filter_object_head, blank_row, filter_object_tail)

  cat("RNAFilt object\n")
  cat(capture.output(result), sep = "\n")
  cat("\n")
}

#' print.proteomicsFilt
#'
#' For printing an S3 object of type 'proteomicsFilt'
#' @param x An object of type 'proteomicsFilt'
#' @param ...  further arguments passed to or from other methods
#'
#' @return No return value, prints details about x
#' 
#' @rdname print-proteomicsFilt
#' @export
#'
print.proteomicsFilt <- function(x, ...) {
  filter_object <- x

  if (!inherits(filter_object, "proteomicsFilt")) stop("filter object must be of the class 'proteomicsFilt'")

  counts_by_pep <- as.data.frame(lapply(filter_object$counts_by_pep, as.character), stringsAsFactors = FALSE)
  counts_by_pro <- as.data.frame(lapply(filter_object$counts_by_pro, as.character), stringsAsFactors = FALSE)

  counts_by_pep_head <- head(counts_by_pep, 4)[1:ncol(counts_by_pep)]
  counts_by_pep_tail <- tail(counts_by_pep, 4)[1:ncol(counts_by_pep)]
  blank_row = rep("---", ncol(counts_by_pep))

  counts_by_pro_head <- head(counts_by_pro, 4)[1:ncol(counts_by_pro)]
  counts_by_pro_tail <- tail(counts_by_pro, 4)[1:ncol(counts_by_pro)]
  blank_row2 = rep("---", ncol(counts_by_pro))

  bypep <- rbind(counts_by_pep_head, blank_row, counts_by_pep_tail)
  bypro <- rbind(counts_by_pro_head, blank_row2, counts_by_pro_tail)

  cat("proteomicsFilt object\n")
  cat("counts_by_pep\n")
  cat(capture.output(bypep), sep = "\n")
  cat("\n")

  cat("counts_by_pro\n")
  cat(capture.output(bypro), sep = "\n")
  cat("\n")
}

#' print.imdanovaFilt
#'
#' For printing an S3 object of type 'imdanovaFilt'
#' @param x An object of type 'imdanovaFilt'
#' @param ...  further arguments passed to or from other methods
#'
#' @return No return value, prints details about x
#' 
#' @rdname print-imdanovaFilt
#' @export
#'
print.imdanovaFilt <- function(x, ...) {
  filter_object <- x

  if (!inherits(filter_object, "imdanovaFilt")) stop("filter object must be of the class 'imdanovaFilt'")
  filter_object <- as.data.frame(lapply(filter_object, as.character), stringsAsFactors = FALSE)

  filter_object_head = head(filter_object, 4)[, 1:ncol(filter_object)]
  filter_object_tail = tail(filter_object, 4)[, 1:ncol(filter_object)]
  blank_row = rep("---", ncol(filter_object))

  result <- rbind(filter_object_head, blank_row, filter_object_tail)

  cat("imdanovaFilt object\n")
  cat(capture.output(result), sep = "\n")
  cat("\n")
}

#' print.rmdFilt
#'
#' For printing an S3 object of type 'rmdFilt'
#' @param x An object of type 'rmdFilt'
#' @param ...  further arguments passed to or from other methods
#'
#' @return No return value, prints details about x
#'
#' @rdname print-rmdFilt
#' @export
#'
print.rmdFilt <- function(x, ...) {
  filter_object <- x

  if (!inherits(filter_object, "rmdFilt")) stop("filter object must be of the class 'rmdFilt'")
  filter_object <- as.data.frame(lapply(filter_object, as.character), stringsAsFactors = FALSE)
  num_cols <- ncol(filter_object)

  filter_object_head = head(filter_object, 4)[, 1:min(num_cols, 5)]
  filter_object_tail = tail(filter_object, 4)[, 1:min(num_cols, 5)]
  blank_row = rep("---", ncol(filter_object))

  result <- rbind(filter_object_head, blank_row, filter_object_tail)

  if (num_cols > 5) message("only first 5 columns are shown")
  cat("rmdFilt object\n")
  cat(capture.output(result), sep = "\n")
  cat("\n")
}

#' print.cvFilt
#'
#' For printing an S3 object of type 'cvFilt'
#' @param x An object of type 'cvFilt'
#' @param ...  further arguments passed to or from other methods
#'
#' @return No return value, prints details about x
#'
#' @rdname print-cvFilt
#' @export
#'
print.cvFilt <- function(x, ...) {
  filter_object <- x

  if (!inherits(filter_object, "cvFilt")) stop("filter object must be of the class 'cvFilt'")
  filter_object <- as.data.frame(lapply(filter_object, as.character), stringsAsFactors = FALSE)

  filter_object_head = head(filter_object, 4)[, 1:ncol(filter_object)]
  filter_object_tail = tail(filter_object, 4)[, 1:ncol(filter_object)]
  blank_row = rep("---", ncol(filter_object))

  result <- rbind(filter_object_head, blank_row, filter_object_tail)

  cat("cvFilt object\n")
  cat(capture.output(result), sep = "\n")
  cat("\n")
}

#' print.customFilt
#'
#' For printing an S3 object of type 'customFilt'
#' @param x An object of type 'customFilt'
#' @param ...  further arguments passed to or from other methods
#'
#' @return No return value, prints details about x
#'
#' @rdname print-customFilt
#' @export
#'
print.customFilt <- function(x, ...) {
  filter_object <- x

  if (!inherits(filter_object, "customFilt")) stop("filter object must be of the class 'customFilt'")

  if (!is.null(filter_object$e_data_remove) | !is.null(filter_object$f_data_remove) | !is.null(filter_object$e_meta_remove)) {
    if (!is.null(filter_object$e_data_remove)) {
      edata_rm <- filter_object$e_data_remove

      if (length(edata_rm) < 8) {
        out <- cbind(edata_rm)
        colnames(out) <- 'e_data_remove'
        cat(capture.output(out), sep = '\n')
        cat('\n')
      } else {
        edata_rm_hd <- cbind(head(edata_rm, 4))
        edata_rm_tl <- cbind(tail(edata_rm, 4))
        blank_row <- '---'

        out <- rbind(edata_rm_hd, blank_row, edata_rm_tl)
        colnames(out) <- 'e_data_remove'
        rownames(out) <- NULL
        cat(capture.output(out), sep = '\n')
        cat('\n')
      }
    }

    if (!is.null(filter_object$f_data_remove)) {
      fdata_rm <- filter_object$f_data_remove

      if (length(fdata_rm) < 8) {
        out <- cbind(fdata_rm)
        colnames(out) <- 'f_data_remove'
        cat(capture.output(out), sep = '\n')
        cat('\n')
      } else {
        fdata_rm_hd <- cbind(head(fdata_rm, 4))
        fdata_rm_tl <- cbind(tail(fdata_rm, 4))
        blank_row <- '---'

        out <- rbind(fdata_rm_hd, blank_row, fdata_rm_tl)
        colnames(out) <- 'f_data_remove'
        rownames(out) <- NULL
        cat(capture.output(out), sep = '\n')
        cat('\n')
      }
    }

    if (!is.null(filter_object$e_meta_remove)) {
      emeta_rm <- filter_object$e_meta_remove

      if (length(emeta_rm) < 8) {
        out <- cbind(emeta_rm)
        colnames(out) <- 'e_meta_remove'
        cat(capture.output(out), sep = '\n')
        cat('\n')
      } else {
        emeta_rm_hd <- cbind(head(emeta_rm, 4))
        emeta_rm_tl <- cbind(tail(emeta_rm, 4))
        blank_row <- '---'

        out <- rbind(emeta_rm_hd, blank_row, emeta_rm_tl)
        colnames(out) <- 'e_meta_remove'
        rownames(out) <- NULL
        cat(capture.output(out), sep = '\n')
        cat('\n')
      }
    }
  }
  if (!is.null(filter_object$e_data_keep) | !is.null(filter_object$f_data_keep) | !is.null(filter_object$e_meta_keep)) {
    if (!is.null(filter_object$e_data_keep)) {
      edata_kp <- filter_object$e_data_keep

      if (length(edata_kp) < 8) {
        out <- cbind(edata_kp)
        colnames(out) <- 'e_data_keep'
        cat(capture.output(out), sep = '\n')
        cat('\n')
      } else {
        edata_kp_hd <- cbind(head(edata_kp, 4))
        edata_kp_tl <- cbind(tail(edata_kp, 4))
        blank_row <- '---'

        out <- rbind(edata_kp_hd, blank_row, edata_kp_tl)
        colnames(out) <- 'e_data_keep'
        rownames(out) <- NULL
        cat(capture.output(out), sep = '\n')
        cat('\n')
      }
    }

    if (!is.null(filter_object$f_data_keep)) {
      fdata_kp <- filter_object$f_data_keep

      if (length(fdata_kp) < 8) {
        out <- cbind(fdata_kp)
        colnames(out) <- 'f_data_keep'
        cat(capture.output(out), sep = '\n')
        cat('\n')
      } else {
        fdata_kp_hd <- cbind(head(fdata_kp, 4))
        fdata_kp_tl <- cbind(tail(fdata_kp, 4))
        blank_row <- '---'

        out <- rbind(fdata_kp_hd, blank_row, fdata_kp_tl)
        colnames(out) <- 'f_data_keep'
        rownames(out) <- NULL
        cat(capture.output(out), sep = '\n')
        cat('\n')
      }
    }

    if (!is.null(filter_object$e_meta_keep)) {
      emeta_kp <- filter_object$e_meta_keep

      if (length(emeta_kp) < 8) {
        out <- cbind(emeta_kp)
        colnames(out) <- 'e_meta_keep'
        cat(capture.output(out), sep = '\n')
        cat('\n')
      } else {
        emeta_kp_hd <- cbind(head(emeta_kp, 4))
        emeta_kp_tl <- cbind(tail(emeta_kp, 4))
        blank_row <- '---'

        out <- rbind(emeta_kp_hd, blank_row, emeta_kp_tl)
        colnames(out) <- 'e_meta_keep'
        rownames(out) <- NULL
        cat(capture.output(out), sep = '\n')
        cat('\n')
      }
    }
  }
}
