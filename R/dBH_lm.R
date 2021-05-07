#' The dependence-adjusted Benjamini-Hochberg (BH) procedure and step-up procedures for fixed-design
#' homoscedastic Gaussian linear models
#'
#' \code{dBH_lm} computes the rejection set of \eqn{dBH_\gamma(\alpha)} and \eqn{dBH^2_\gamma(\alpha)}, as well as
#' \eqn{dSU_{\gamma, \Delta}(\alpha)} and \eqn{dSU^2_{\gamma, \Delta}(\alpha)} a broad class of step-up procedures for
#' linear models \eqn{y = X\beta + \epsilon} where \eqn{X} is a fixed-design matrix and \eqn{\epsilon} has i.i.d.
#' components drawn from \eqn{N(0, \sigma^2)} for some unknonw variance \eqn{\sigma^2}. \code{dBH_lm} can handle
#' both one-sided tests \eqn{H_i: \beta_i \le 0} or \eqn{H_i: \beta_i \ge 0}, and two-sided tests \eqn{H_i: \beta_i = 0}.
#'
#' @details \code{dBH_lm} can handle all dSU procedures that are defined in Appendix C.1. with
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
#' @param y the outcome vector.
#' @param X the design matrix. The intercept term should be added manually into \code{X} to be included. See Details.
#' @param subset a subset of \code{1:ncol(X)}, which specifies the subset of coefficients to be tested. The default is
#' \code{1:ncol(X)}, which tests all coefficients.
#' @param intercept a logical indicating whether an intercept is included in the linear regression. By default,
#' the intercept term is not tested. If it needs to be tested together with other variables, manually add it into
#' \code{X} and set \code{intercept = FALSE}.
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
#' \code{\link{dBH_mvgauss}}, \code{\link{dBH_mvt}}
#' 
#' @examples
#' \donttest{# Generate beta
#' n <- 150
#' p <- 100
#' beta1 <- 0.5
#' nalt <- 10
#' beta <- c(rep(beta1, nalt), rep(0, p - nalt))
#'
#' # Generate X and y
#' set.seed(1)
#' X <- matrix(rnorm(n * p), nrow = n)
#' y <- X %*% beta + rnorm(n)
#'
#' # Run dBH_1(\alpha) for one-sided tests
#' alpha <- 0.05
#' res <- dBH_lm(y = y, X = X, side = "right", alpha = alpha,
#'               gamma = 1, niter = 1, avals_type = "BH")
#' 
#' # Run dBH_1(\alpha) for two-sided tests
#' res <- dBH_lm(y = y, X = X, side = "right", alpha = alpha,
#'               gamma = 1, niter = 1, avals_type = "BH") 
#'
#' # Run dBH^2_1(\alpha) for one-sided tests
#' alpha <- 0.05
#' res <- dBH_lm(y = y, X = X, side = "right", alpha = alpha,
#'               gamma = 1, niter = 2, avals_type = "BH")
#' 
#' # Run dBH^2_1(\alpha) for one-sided tests
#' res <- dBH_lm(y = y, X = X, side = "right", alpha = alpha,
#'               gamma = 1, niter = 2, avals_type = "BH")
#' 
#' # Run dSU_1(\alpha) with the geometrically increasing a-values for one-sided tests
#' res <- dBH_lm(y = y, X = X, side = "right", alpha = alpha,
#'               gamma = 1, niter = 1, avals_type = "geom",
#'               geom_fac = 2)
#'
#' # Run dBY(\alpha) for one-sided tests
#' res <- dBH_lm(y = y, X = X, side = "right", alpha = alpha,
#'               gamma = NULL, niter = 1, avals_type = "BH")
#'
#' # Run dBY^2(\alpha) for one-sided tests
#' res <- dBH_lm(y = y, X = X, side = "right", alpha = alpha,
#'               gamma = NULL, niter = 2, avals_type = "BH")
#' 
#' }
#'
#' @export
dBH_lm <- function(y, X,
                   subset = 1:ncol(X),
                   intercept = TRUE,
                   side = c("right", "left", "two"),
                   alpha = 0.05, gamma = NULL,
                   tautype = "QC",
                   niter = 1,                   
                   avals = NULL,
                   avals_type = c("BH", "geom", "bonf", "manual"),
                   geom_fac = 2,
                   eps = 0.05,
                   qcap = 2,
                   gridsize = 20,
                   exptcap = 0.9,
                   is_safe = NULL,
                   verbose = FALSE){
    stats <- lm_mvt(y, X, subset, intercept)
    dBH_mvt(tvals = stats$tvals, 
            df = stats$df,
            Sigma = stats$Sigma,
            Sigmafun = NULL,
            vars = NULL,
            side = side,
            alpha = alpha,
            gamma = gamma,
            tautype = tautype,
            niter = niter,
            avals = avals,
            avals_type = avals_type,
            geom_fac = geom_fac,
            eps = eps,
            qcap = qcap,
            gridsize = gridsize,
            exptcap = exptcap,
            is_safe = is_safe,
            verbose = verbose)
}
