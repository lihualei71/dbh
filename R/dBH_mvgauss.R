#' The dependence-adjusted Benjamini-Hochberg (BH) procedure and step-up procedures for multivariate z-statistics
#'
#' \code{dBH_mvgauss} computes the rejection set of \eqn{dBH_\gamma(\alpha)} and \eqn{dBH^2_\gamma(\alpha)}, as well as
#' \eqn{dSU_{\gamma, \Delta}(\alpha)} and \eqn{dSU^2_{\gamma, \Delta}(\alpha)} a broad class of step-up procedures for
#' multivariate z-statistics \eqn{z\sim N(\mu, \Sigma)}. \code{dBH_mvgauss} can handle both one-sided tests
#' \eqn{H_i: \mu_i \le 0} or \eqn{H_i: \mu_i \ge 0}, and two-sided tests \eqn{H_i: \mu_i = 0}.
#'
#' @details \code{dBH_mvgauss} supports two types of inputs for the covariance matrix.
#' \itemize{
#' \item When the covariance matrix fits into the memory, it can be inputted through \code{Sigma}. The diagonals do not
#' have to be 1. In this case, \code{Sigmafun} and \code{vars} should be left as their default;
#' \item When the covariance matrix does not fit into the memory, it can be inputted through \code{Sigmafun}, for the
#' correlation matrix, and \code{vars}, for marginal variances. \code{Sigmafun} should be a function with a single input
#' \code{i} that gives the row index and a vector output \eqn{\Sigma_{i, }/\sqrt{\Sigma_{i, i}\Sigma_{j, j}}} where
#' \eqn{\Sigma_{i, }} is the i-th row of \eqn{\Sigma}. The marginal variances \eqn{(\Sigma_{1, 1}, \ldots, \Sigma_{m, m})}
#' are inputted through \code{vars}. If all marginal variances are 1, \code{vars} can be left as its default.
#' }
#'
#' \code{dBH_mvgauss} can handle all dSU procedures that are defined in Appendix C.1. with
#' \deqn{\Delta_{\alpha}(r) = \frac{\alpha a_{\ell}}{m}, r\in [a_{\ell}, a_{\ell + 1}), \ell = 0, 1, \ldots, L}
#' for any set of integer a-values \eqn{1\le a_1 < \ldots < a_L\le m} (with \eqn{a_0 = 0, a_{L+1} = m+1}). There are two-ways to input the a-values.
#' \itemize{
#' \item (Recommended) use \code{avals_type} while leave \code{avals} as its default. Three types of \code{avals_type}
#' are supported.
#'     \itemize{
#'     \item When \code{avals_type = "BH"}, the BH procedure is used (i.e. \eqn{a_\ell = \ell, \ell = 1, \ldots, m});
#'     \item When \code{avals_type = "geom"}, the geometrically increasing a-values that are defined in Appendix C.1.
#' are used with the growing rate specified by \code{geom_fac}, whose default is 2;
#'     \item When \code{avals_type = "bonf"}, the Bonferroni procedure is used (i.e. \eqn{a_1 = 1, L = 1}).
#'     }
#' \item (Not recommended) use \code{avals} while set \code{avals_type = "manual"}. This option allows any set of
#' increasing a-values. However, it can be much slower than the recommended approach.
#' }
#'
#' @param zvals a vector of z-values. The marginal variances do not have to be 1.
#' @param Sigma the covariance matrix. For very large m, it is preferable to use \code{Sigmafun} while leave
#' \code{Sigma = NULL} as its default. See Details.
#' @param Sigmafun the function that takes the row index as the input and outputs the i-th row of the covariance matrix
#' divided by \eqn{\Sigma_{i, i}}. This is an alternative for \code{Sigma} when m is very large. See Details.
#' @param vars the vector of marginal variances \eqn{(\sigma_1^2, \ldots, \sigma_m^2)}. Used only with \code{Sigmafun}.
#' See Details.
#' @param side a string that takes values in \{"right", "left", "two"\}, with "right" for one-sided tests
#' \eqn{H_i: \mu_i \le 0}, "left" for one-sided tests, \eqn{H_i: \mu_i \ge 0}, and "two" for two-sided tests
#' \eqn{H_i: \mu_i = 0}.
#' @param alpha the target FDR level.
#' @param gamma the parameter for the dBH and dSU procedures. The default is \code{NULL}, which gives dBY or the safe
#' dSU with gamma = 1 / Lm defined in Appendix C.1.
#' @param niter the number of iterations. In the current version it can only be 1, for dBH/dSU, or 2, for
#' \eqn{dBH^2}/\eqn{dSU^2}.
#' @param tautype the type of tau function. In the current version, only "QC" (q-value-calibration) is supported with
#' \eqn{\tau(c; X) = cR_{BH}(c) / m}.
#' @param avals the a-values in the step-up procedures defined in Appendix C.1. The default is NULL, in which case
#' \code{avals} is determined by \code{avals_type}. See Details.
#' @param avals_type a string that takes values in \{"BH", "geom", "bonf", "manual"\}, which determines the type of
#' thresholds in the step-up procedures. See Details.
#' @param geom_fac a real number that is larger than 1. This is the growing rate of thresholds when
#' \code{avals_type = "geom"}. See Details.
#' @param eps a real number in [0, 1], which is used to determine \eqn{t_{hi}} discussed in Section 4.2.
#' @param qcap a real number that is larger than 1. It is used to filter out hypotheses with q-values above
#' \code{qcap} X \code{alpha}, as discussed in Appendix C.2.2.
#' @param gridsize an integer for the size of the grid used when \code{niter = 2}. For two-sided tests, \code{gridsize}
#' knots will be used for both the positive and negative sides.
#' @param exptcap a real number in [0, 1]. It is used by \eqn{dBH^2}/\eqn{dSU^2} to filter out hypotheses with
#' \eqn{g_i*(q_i | S_i)} below \code{exptcap} X \code{alpha} / m in their initializations, as discussed in Appendix C.2.5.
#' @param is_safe a logical or NULL indicating whether the procedure is taken as safe. DON'T set \code{is_safe = TRUE}
#' unless the covariance structure is known to be CPRDS. The default is \code{NULL}, which sets \code{is_safe = TRUE}
#' if \code{gamma = NULL} or \code{gamma} is below 1 / Lm, and \code{is_safe = FALSE} otherwise.
#' @param verbose a logical indicating whether a progress bar is shown.
#'
#' @return a list with the following attributes
#' \itemize{
#' \item \code{rejs}: the indices of rejected hypotheses (after the randomized pruning step if any);
#' \item \code{initrejs}: the indices of rejected hypotheses (before the randomized pruning step if any);
#' \item \code{cand}: the set of candidate hypotheses for which \eqn{g_i^{*}(q_i | S_i)} is evaluated;
#' \item \code{expt}: \eqn{g_i^{*}(q_i | S_i)} for each hypothesis in \code{cand};
#' \item \code{safe}: TRUE iff the procedure is safe;
#' \item \code{secBH}: TRUE iff the randomzied pruning step (a.k.a. the secondary BH procedure) is invoked;
#' \item \code{secBH_fac}: a vector that gives \eqn{\hat{R}_i / R_{+}} that is defined in Section 2.2. It only shows up
#' in the output if \code{secBH = TRUE}.
#' }
#'
#' @seealso
#' \code{\link{dBH_mvt}}, \code{\link{dBH_lm}}
#'
#' @examples
#' \donttest{# Generate mu and Sigma for an AR process
#' n <- 100
#' rho <- 0.8
#' Sigma <- rho^(abs(outer(1:n, 1:n, "-")))
#' mu1 <- 2.5
#' nalt <- 10
#' mu <- c(rep(mu1, nalt), rep(0, n - nalt))
#'
#' # Generate the z-values
#' set.seed(1)
#' zvals <- rep(NA, n)
#' zvals[1] <- rnorm(1)
#' for (i in 2:n){
#'     zvals[i] <- zvals[i - 1] * rho + rnorm(1) * sqrt(1 - rho^2)
#' }
#' zvals <- zvals + mu
#'
#' # Run dBH_1(\alpha) for one-sided tests
#' alpha <- 0.05
#' res <- dBH_mvgauss(zvals = zvals, Sigma = Sigma, side = "right", alpha = alpha,
#'                    gamma = 1, niter = 1, avals_type = "BH")
#'
#' # Run dBH_1(\alpha) for two-sided tests
#' res <- dBH_mvgauss(zvals = zvals, Sigma = Sigma, side = "right", alpha = alpha,
#'                    gamma = 1, niter = 1, avals_type = "BH")
#'
#' # Run dBH^2_1(\alpha) for one-sided tests
#' alpha <- 0.05
#' res <- dBH_mvgauss(zvals = zvals, Sigma = Sigma, side = "right", alpha = alpha,
#'                    gamma = 1, niter = 2, avals_type = "BH")
#'
#' # Run dBH^2_1(\alpha) for one-sided tests
#' res <- dBH_mvgauss(zvals = zvals, Sigma = Sigma, side = "right", alpha = alpha,
#'                    gamma = 1, niter = 2, avals_type = "BH")
#'
#' # Run dSU_1(\alpha) with the geometrically increasing a-values for one-sided tests
#' res <- dBH_mvgauss(zvals = zvals, Sigma = Sigma, side = "right", alpha = alpha,
#'                    gamma = 1, niter = 1, avals_type = "geom", geom_fac = 2)
#'
#' # Run dBY(\alpha) for one-sided tests
#' res <- dBH_mvgauss(zvals = zvals, Sigma = Sigma, side = "right", alpha = alpha,
#'                    gamma = NULL, niter = 1, avals_type = "BH")
#'
#' # Run dBY^2(\alpha) for one-sided tests
#' res <- dBH_mvgauss(zvals = zvals, Sigma = Sigma, side = "right", alpha = alpha,
#'                    gamma = NULL, niter = 2, avals_type = "BH")
#'
#' # Input the covariance matrix through \code{Sigmafun}
#' Sigmafun <- function(i){
#'    rho^(abs(1:n - i))
#'}
#' vars <- rep(1, n)
#' res1 <- dBH_mvgauss(zvals = zvals, Sigmafun = Sigmafun, vars = vars, side = "right",
#'                     alpha = alpha, gamma = NULL, niter = 1, avals_type = "BH")
#' res2 <- dBH_mvgauss(zvals = zvals, Sigma = Sigma, side = "right",
#'                     alpha = alpha, gamma = NULL, niter = 1, avals_type = "BH")
#' identical(res1, res2)
#'
#' }
#'
#' @export
dBH_mvgauss <- function(zvals,
                        Sigma = NULL,
                        Sigmafun = NULL,
                        vars = NULL,
                        side = c("right", "left", "two"),
                        alpha = 0.05, gamma = NULL,
                        niter = 1,
                        tautype = "QC",
                        avals = NULL,
                        avals_type = c("BH", "geom", "bonf", "manual"),
                        geom_fac = 2,
                        eps = 0.05,
                        qcap = 2,
                        gridsize = 20,
                        exptcap = 0.9,
                        is_safe = NULL,
                        verbose = FALSE){
    if (niter > 2){
        stop("\'niter\' can only be 1 or 2.")
    }

    side <- side[1]
    avals_type <- avals_type[1]
    tautype <- tautype[1]
    n <- length(zvals)
    if (is.null(avals)){
        if (avals_type == "manual"){
            stop("avals must be inputted when avals_type = \"manual\"")
        } else if (avals_type == "geom" && geom_fac <= 1){
            stop("geom_fac must be larger than 1 when avals_type = \"geom\"")
        }
        avals <- switch(avals_type,
                        BH = 1:n,
                        geom = geom_avals(geom_fac, n),
                        bonf = 1)
    } else if (is.null(avals_type)) {
        if (avals[1] != 1){
            stop("The first element of avals must be 1.")
        }
        avals_type <- "manual"
        warning("avals is inputted and avals_type is set to be \"manual\" by default. This may slow down the code. Use the built-in avals_type (\"BH\", \"geom\" or \"bonf\") instead unless there is a good reason to use the inputted avals.")
    } else {
        stop("Set avals = NULL when avals_type is specified")
    }

    if (is.null(gamma)){
        gamma <- 1 / normalize(avals)
        is_safe <- TRUE
    } else if (is.null(is_safe)){
        is_safe <- (gamma <= 1 / normalize(avals))
    } else if (!is_safe){
        warning("Set is_safe = TRUE only if you can prove the procedure is safe (e.g. CPRDS case with gamma = 1)")
    }

    if (!is.null(Sigma) && any(diag(Sigma) != 1)){
        vars <- diag(Sigma)
        Sigma <- cov2cor(Sigma)
    } else {
        if (is.null(vars)){
            vars <- rep(1, n)
        }
    }
    zvals <- zvals / sqrt(vars)
    if (side == "left"){
        zvals <- -zvals
        side <- "one"
    } else if (side == "right"){
        side <- "one"
    }

    if (niter == 1){
        if (tautype == "QC"){
            dBH_mvgauss_qc(zvals = zvals,
                           Sigma = Sigma,
                           Sigmafun = Sigmafun,
                           side = side,
                           alpha = alpha,
                           gamma = gamma,
                           is_safe = is_safe,
                           avals = avals,
                           avals_type = avals_type,
                           geom_fac = geom_fac,
                           eps = eps,
                           qcap = qcap,
                           verbose = verbose)
        }
    } else if (niter == 2){
        if (tautype == "QC"){
            dBH_mvgauss_qc_grid(zvals = zvals,
                                Sigma = Sigma,
                                Sigmafun = Sigmafun,
                                side = side,
                                alpha = alpha,
                                gamma = gamma,
                                is_safe = is_safe,
                                avals = avals,
                                avals_type = avals_type,
                                geom_fac = geom_fac,
                                eps = eps,
                                qcap = qcap,
                                gridsize = gridsize,
                                exptcap = exptcap,
                                verbose = verbose)
        }
    }
}
