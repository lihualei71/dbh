dBH_mvgauss_qc_grid <- function(zvals,
                                Sigma = NULL,
                                Sigmafun = NULL,
                                side = c("right", "left", "two"),
                                alpha = 0.05, gamma = NULL,
                                is_safe = FALSE,
                                avals = NULL,
                                avals_type = c("BH", "geom", "bonf", "manual"),
                                geom_fac = 2,
                                eps = 0.05,
                                qcap = 2,
                                gridsize = 20,
                                exptcap = 0.9,
                                verbose = FALSE){
    n <- length(zvals)
    ntails <- ifelse(side == "two", 2, 1)
    high <- qnorm(alpha * eps / n / ntails, lower.tail = FALSE)
    high <- abs(high) # just in case high is negative
    pvals <- zvals_pvals(zvals, side)
    qvals <- qvals_BH_reshape(pvals, avals)

    params_root <- list(Sigma = Sigma,
                        Sigmafun = Sigmafun,
                        side = side,
                        alpha = alpha, gamma = gamma,
                        avals = avals,
                        avals_type = avals_type,
                        geom_fac = geom_fac,
                        eps = eps,
                        qcap = qcap,
                        verbose = FALSE)
    params <- c(params_root, list(zvals = zvals))
    res_init <- do.call(dBH_mvgauss_qc, params)
    Rinit <- rep(length(res_init$initrejs) + 1, n)
    Rinit[res_init$initrejs] <- Rinit[res_init$initrejs] - 1

    cand <- res_init$cand
    if (res_init$safe){
        init_rejlist <- res_init$initrejs
    } else {
        init_rejlist <- which(qvals <= alpha / max(avals))
        if (!is.null(exptcap)){
            init_rejlist <- union(cand[res_init$expt <= exptcap * alpha], init_rejlist)
        }
    }
    cand <- setdiff(cand, init_rejlist)

    if (length(cand) == 0){
        return(list(rejs = init_rejlist,
                    initrejs = init_rejlist,
                    cand = numeric(0),
                    expt = numeric(0),
                    safe = res_init$safe,
                    secBH = FALSE))
    }

    if (verbose){
        pb <- txtProgressBar(style=3)
    }
    ncands <- length(cand)
    cand_info <- sapply(1:ncands, function(id){
        i <- cand[id]
        low <- qnorm(qvals[i] * max(avals) / n / ntails, lower.tail = FALSE)
        if (!is.null(Sigma)){
            cor <- Sigma[-i, i]
        } else {
            cor <- Sigmafun(i)[-i]
        }
        s <- zvals[-i] - cor * zvals[i]

        ## RBH function with alpha = qi
        res_q <- compute_knots_mvgauss(
            zstat = zvals[i],
            zminus = zvals[-i],
            cor = cor,
            alpha = qvals[i],
            side = side,
            low = low,
            high = high,
            avals = avals,
            avals_type = avals_type,
            geom_fac = geom_fac)
        res_q <- lapply(res_q, function(re){
            RBH <- RejsBH(re$posit, re$sgn, re$RCV, avals)
            knots <- c(re$low, re$knots)
            RBH <- rle(RBH)
            nrejs <- RBH$values
            cutinds <- c(1, cumsum(RBH$lengths) + 1)
            knots <- c(knots, re$high)
            knots <- knots[cutinds]
            if (knots[1] < 0 && side == "two"){
                knots <- rev(abs(knots))
                ## This requires the null distribution to be symmetric
                nrejs <- rev(nrejs)
            }
            if (avals_type == "BH"){
                thra <- nrejs
            } else if (avals_type == "geom"){
                thra <- find_ind_geom_avals(geom_fac, nrejs, "max")
                ## 0 rejection should return aval = 0
                thra[thra == 0] <- NA
                thra <- avals[thra]
                thra[is.na(thra)] <- 0
            } else if (avals_type == "manual"){
                thra <- fill_int_general(nrejs, avals)
                ## 0 rejection should return aval = 0
                thra[thra == 0] <- NA
                thra <- avals[thra]
                thra[is.na(thra)] <- 0
            } else if (avals_type == "bonf"){
                thra <- rep(1, length(nrejs))
            }
            thr <- qnorm(thra * qvals[i] / n / ntails, lower.tail = FALSE)
            list(knots = knots, thr = thr)
        })

        ## Get the intervals of knots at which the numerator
        ## is 1
        grids <- lapply(1:ntails, function(k){
            ints <- find_int_above_thr(res_q[[k]]$knots, res_q[[k]]$thr)
            ## Filter out very small intervals
            if (k == 2){
                lapply(ints, function(int){
                    -rev(int)
                })
            } else {
                ints
            }
        })
        grids <- do.call(c, grids)
        prob <- sapply(grids, function(int){
            diff(pnorm(int))
        })
	      if (length(prob) < 1 || !is.numeric(prob) || sum(prob) * n <= alpha){
            return(c(1, NA))
        }

        ## Add extra knots
        nknots <- ceiling(prob / sum(prob) * gridsize * ntails)
        grids <- lapply(1:length(grids), function(k){
            seq(grids[[k]][1], grids[[k]][2],
                length.out = nknots[k] + 1)
        })

        if (verbose){
            cumsum_nknots <- cumsum(nknots)
            sum_nknots <- tail(cumsum_nknots, 1)
            cumsum_nknots <- c(0, head(cumsum_nknots, -1))
        }

        ## Create grid for the denominator
        expt <- sapply(1:length(grids), function(k){
            grid <- grids[[k]]
            prob <- diff(pnorm(grid))
            ex <- sapply(1:(length(grid) - 1), function(j){
                pr <- prob[j]
                if (any(grid > 0)){
                    j <- j + 1
                }
                tmp <- recover_stats_mvgauss(zvals[i], grid[j], s, cor)
                zvals_tmp <- rep(0, n)
                zvals_tmp[i] <- tmp[1]
                zvals_tmp[-i] <- tmp[-1]
                params <- c(params_root, list(zvals = zvals_tmp))
                res <- do.call(dBH_mvgauss_qc, params)
                nrejs <- length(res$initrejs) + !(i %in% res$initrejs)
                if (verbose){
                    tmp <- id - 1 + (j + cumsum_nknots[k]) / sum_nknots
                    setTxtProgressBar(pb, tmp / ncands)
                }
                return(pr / nrejs)
            })
            sum(ex)
        })
        expt <- sum(expt) * n
        ifrej <- expt <= alpha
        return(c(ifrej, expt))
    })

    if (verbose){
        cat("\n")
    }

    ifrej <- as.logical(cand_info[1, ])
    rejlist <- which(ifrej)
    rejlist <- c(init_rejlist, cand[rejlist])
    expt <- cand_info[2, ]

    if (length(rejlist) == 0){
        return(list(rejs = numeric(0),
                    initrejs = numeric(0),
                    cand = cand,
                    expt = expt,
                    safe = res_init$safe,
                    secBH = FALSE))
    }

    rejlist <- sort(rejlist)
    Rplus <- length(rejlist)
    if (Rplus >= max(Rinit[rejlist])){
        return(list(rejs = rejlist,
                    initrejs = rejlist,
                    cand = cand,
                    expt = expt,
                    safe = res_init$safe,
                    secBH = FALSE))
    }

    uvec <- runif(Rplus)
    secBH_fac <- Rinit[rejlist] / Rplus
    tdp <- uvec * secBH_fac
    nrejs <- nrejs_BH(tdp, 1)
    thr <- max(nrejs, 1) / Rplus
    secrejs <- which(tdp <= thr)
    rejs <- rejlist[secrejs]
    return(list(rejs = rejs,
                initrejs = rejlist,
                cand = cand,
                expt = expt,
                safe = FALSE,
                secBH = TRUE,
                secBH_fac = secBH_fac))
}
