dBH_mvgauss_qc_optimal <- function(zvals,
                           Sigma = NULL,
                           Sigmafun = NULL,
                           side = c("one", "two"),
                           MC = NULL,
                           groups = NULL,
                           lfdrinv_type = "max",
                           alpha = 0.05, gamma = NULL,
                           is_safe = FALSE,
                           avals = NULL, 
                           avals_type = c("BH", "geom", "bonf", "manual"),
                           geom_fac = 2,
                           eps = 0.05,
                           qcap = 2,
                           kappa = 0.5,
                           verbose = FALSE) {
    n <- length(zvals)
    alpha0 <- gamma * alpha
    ntails <- ifelse(side == "two", 2, 1)    
    pvals <- zvals_pvals(zvals, side)
    pvals[which(pvals >= kappa)] = Inf

    init_weights <- adaptive_optimal.weights(
                    groups = groups, 
                    zvals = zvals, 
                    alpha = alpha, 
                    side = side,
                    type = lfdrinv_type)
    init_qvals <- qvals_BH_reshape(pvals/init_weights, avals)
    init_acclist <- which(init_qvals >= qcap * alpha | pvals >= kappa)
    if (length(init_acclist) > 0){
        cand <- (1:n)[-init_acclist]
    }

    if (length(cand) == 0){
        return(list(rejs = numeric(0),
                    initrejs = numeric(0),
                    cand = numeric(0),
                    expt = numeric(0),
                    safe = is_safe,
                    secBH = FALSE))
    }

    if (verbose){
        pb <- txtProgressBar(style=3)
    }

    ncands <- length(cand)
    
    cand_info <- sapply(1:ncands, function(id){
        i <- cand[id]
        if (!is.null(Sigma)){
            cor <- Sigma[-i, i]
        } else {
            cor <- Sigmafun(i)[-i]
        }
        zstat = zvals[i]
        zminus = zvals[-i]
        s <- zminus - cor * zstat
        group <- groups[i]
        groupminus <- groups[-i]

        mcz <- rnorm(MC)
        w <- sapply(mcz, function(z){
          adaptive_optimal.weights(groups = c(group, groupminus), 
                                   zvals = c(z, s + cor * z), 
                                   alpha = alpha, 
                                   side = side,
                                   type = lfdrinv_type)
        })

        weights <- rowMeans(w)
        weight <- weights[1]
        reconsp <- c(pvals[i], pvals[-i])
        qvals <- qvals_BH_reshape(reconsp/weights, avals)
        qval <- qvals[1]
        dBH_rej0 <- which(qvals <= alpha0)
        Rinit <- length(dBH_rej0)
        #Rinit <- ifelse(qval <= alpha0, length(dBH_rej0), length(dBH_rej0)+1)

        if ((qval <= alpha0 & is_safe) | qval <= alpha / max(avals)) {
            initrej = T
            ifrej = T
            expt = NA
            weight = weight 
            return(c(ifrej, initrej, expt, Rinit, weight))
        }
        
        if (qval > qcap * alpha) {
            initrej = F
            ifrej = F
            expt = NA
            weight = weight
            return(c(ifrej, initrej, expt, Rinit, weight))
        }

        initrej = F
        low <- qnorm(min(qval * weight * max(avals) / n / ntails, kappa), lower.tail = FALSE)
        high <- qnorm(weight * alpha * eps / n / ntails, lower.tail = FALSE)

        
        ## RBH function with alpha = qi
        res_q <- compute_knots_mvgauss(
            zstat = zstat,
            zminus = zminus,
            cor = cor,
            alpha = qval,
            side = side,            
            low = low,
            high = high,
            avals = avals,
            avals_type = avals_type,
            geom_fac = geom_fac,
            kappa = kappa,
            weight = weights[1],
            weightminus = weights[-1])
        
        res_q <- lapply(res_q, function(re){
            RBH <- RejsBH(re$posit, re$sgn, re$RCV, avals)
            knots <- c(re$low, re$knots)
            RBH <- rle(RBH)
            nrejs <- RBH$values
            cutinds <- c(1, cumsum(RBH$lengths) + 1)
            knots <- c(knots, re$high)        
            knots <- knots[cutinds]
            if (knots[1] < 0){
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
            thr <- qnorm(thra * qval * weight / n / ntails, lower.tail = FALSE)
            list(knots = knots, thr = thr)
        })

        ## RBH function with alpha = alpha0
        res_alpha0 <- compute_knots_mvgauss(
            zstat = zstat,
            zminus = zminus,
            cor = cor,
            alpha = alpha0,
            side = side,            
            low = low, 
            high = high,
            avals = avals,
            avals_type = avals_type,
            geom_fac = geom_fac,
            kappa = kappa,
            weight = weights[1],
            weightminus = weights[-1])
        
        res_alpha0 <- lapply(res_alpha0, function(re){
            RBH <- RejsBH(re$posit, re$sgn, re$RCV, avals)
            knots <- c(re$low, re$knots)
            RBH <- rle(RBH)
            nrejs <- RBH$values
            cutinds <- c(1, cumsum(RBH$lengths) + 1)
            knots <- c(knots, re$high)        
            knots <- knots[cutinds]
            if (knots[1] < 0){
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
            thr <- qnorm(thra * alpha0 * weight / n / ntails, lower.tail = FALSE)
            knots_lo <- head(knots, -1)
            knots_hi <- tail(knots, -1)
            nrejs <- nrejs + ((knots_lo + knots_hi) / 2 < thr)
            list(knots = knots, nrejs = nrejs)
        })

        ## Combine two RBH functions
        res <- RBHfun_combine(res_q, res_alpha0)

        ## Compute conditional expectation
        expt <- sapply(res, function(re){
            compute_cond_exp(abs(zstat), re$knots, re$nrejs, re$thr, dist = pnorm)
        })

        expt <- sum(expt)
        ifrej <- (expt <= alpha * weight / n)

        if (verbose){
            setTxtProgressBar(pb, id / ncands)
        }
        
        #return(c(ifrej, expt))
        return(c(ifrej, initrej, expt, Rinit, weight))
    })
    rownames(cand_info) <- c("ifrej", "ifrej_init", "exp", "Rinit", "weight")

    if (verbose){
        cat("\n")
    }
    
    ifrej <- as.logical(cand_info[1, ])
    rejindex <- which(ifrej)
    rejlist <- cand[rejindex]
    rejthoCal <- cand_info[2, ]
    expt <- cand_info[3, ]
    Rinit <- cand_info[4, ]
    weights <- cand_info[5, ]

    if (length(rejlist) == 0){
        return(list(rejs = numeric(0),
                    initrejs = numeric(0),
                    rejthoCal = numeric(0),
                    cand = cand,
                    expt = expt,
                    safe = is_safe,
                    secBH = FALSE))
    }

    rejlist <- sort(rejlist)
    Rplus <- length(rejlist)
    if (Rplus >= max(Rinit[rejindex])){
        return(list(rejs = rejlist,
                    initrejs = rejlist,
                    rejthoCal = rejthoCal,
                    cand = cand,
                    expt = expt,
                    weights = weights,
                    safe = is_safe,
                    secBH = FALSE))
    }

    uvec <- runif(Rplus)
    secBH_fac <- Rinit[rejindex] / Rplus
    tdp <- uvec * secBH_fac
    nrejs <- nrejs_BH(tdp, 1)
    thr <- max(nrejs, 1) / Rplus
    secrejs <- which(tdp <= thr)
    rejs <- rejlist[secrejs]
    return(list(rejs = rejs,
                initrejs = rejlist,
                rejthoCal = rejthoCal,
                cand = cand,
                expt = expt,
                weights = weights,
                safe = FALSE,
                secBH = TRUE,
                secBH_fac = secBH_fac))
}




