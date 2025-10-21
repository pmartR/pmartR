# Version 2.5.0 April, 2025

Fixed NOTE's about updating R-version to > 4.1.0 to accommodate base pipe operator and new lambda function syntax.

# Version 2.5.1 October, 2025

Fixed package deprecations that arose from an update to RcppArmadillo. Includes a new package flag during installation (-DARMA_USE_CURRENT) in addition to replacing functions arma::conv_to::from and arma::is_finite to arma::as_scalar and std::isfinite respectively.

