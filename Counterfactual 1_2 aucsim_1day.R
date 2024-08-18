########################################################################################
# This file aims to do Counterfactual 1
# using subcategories clicked by users in the past X days (x=1) in the matching


########################################################################################
library(bayesm)
library(msm)
library(LearnBayes)
library(MCMCpack)
library(numDeriv)
library(fastGHQuad)
library(doParallel)
library(apollo)



########################################################################################
rm(list = ls())


########################################################################################
library(apollo)
library(stargazer)
########################################################################################
### Initialize code


apollo_initialise()

apollo_control = list(
  modelName  ="Demand eMDC model with latent class",
  modelDescr ="Demand eMDC model with latent class",
  indivID    ="user_id3",
  nCores     = 7,
  workInLogs=F
)




#####This function calculates the ad click probability using subcategory information, predicted clicks from MDCEV and choice probability
add_restrict_imputevalue = function(x1) {
  x = unlist(x1[1:8]) 
  y = x1[9:16] 
  z = x1[17:24] 
  
  for (i in 1:length(unique(x))) {
    temp = max(y[x == unique(x)[i]]) 
    y[x == unique(x)[i]] = min(temp, sum( x == unique(x)[i] ) )
    z[x == unique(x)[i] ] = z[x == unique(x)[i] ] / sum(z[x == unique(x)[i] ]) * y[x == unique(x)[i]]
  }
  z[z>1] = 1
  return(z)
}

####This function makes prediction
choiceprobpred = function(database) {
  apollo_mdcev = function (mdcev_settings, functionality)
  {
    modelType = "MDCEV"
    if (is.null(mdcev_settings[["componentName"]])) {
      mdcev_settings[["componentName"]] = ifelse(!is.null(mdcev_settings[["componentName2"]]),
                                                 mdcev_settings[["componentName2"]], modelType)
      test <- functionality == "validate" && mdcev_settings[["componentName"]] !=
        "model" && !apollo_inputs$silent
      if (test)
        apollo_print(paste0("Apollo found a model component of type ",
                            modelType, " without a componentName. The name was set to \"",
                            mdcev_settings[["componentName"]], "\" by default."))
    }
  
    apollo_inputs = tryCatch(get("apollo_inputs", parent.frame(),
                                 inherits = FALSE), error = function(e) return(list(apollo_control = list(cpp = FALSE))))
    if (!is.null(apollo_inputs[[paste0(mdcev_settings$componentName,
                                       "_settings")]]) && (functionality != "preprocess")) {
      tmp <- apollo_inputs[[paste0(mdcev_settings$componentName,
                                   "_settings")]]
      if (is.null(tmp$V))
        tmp$V <- mdcev_settings$V
      if (is.null(tmp$alpha))
        tmp$alpha <- mdcev_settings$alpha
      if (is.null(tmp$gamma))
        tmp$gamma <- mdcev_settings$gamma
      if (is.null(tmp$sigma))
        tmp$sigma <- mdcev_settings$sigma
      mdcev_settings <- tmp
      rm(tmp)
    }
    else {
      mdcev_settings <- apollo_preprocess(mdcev_settings, modelType,
                                          functionality, apollo_inputs)
      if (apollo_inputs$apollo_control$cpp & !apollo_inputs$silent)
        apollo_print("No C++ optimisation available for OL components.")
      if (mdcev_settings$hasOutside)
        mdcev_settings$probs_MDCEV <- function(inputs) {
          inputs$V <- mapply(function(v, a) apollo_setRows(v,
                                                           !a, 0), inputs$V, inputs$avail, SIMPLIFY = FALSE)
          inputs$alpha <- mapply(function(l, a) apollo_setRows(l,
                                                               !a, 0), inputs$alpha, inputs$avail, SIMPLIFY = FALSE)
          inputs$V[[1]] = (inputs$alpha[[1]] - 1) * log(inputs$continuousChoice[[1]])
          for (j in 2:inputs$nAlt) {
            if (inputs$minX) {
              tmp <- inputs$continuousChoice[[j]] - (inputs$continuousChoice[[j]] >=
                                                       inputs$minConsumption[[j]]) * inputs$minConsumption[[j]]
              inputs$V[[j]] = inputs$V[[j]] + inputs$avail[[j]] *
                ((inputs$alpha[[j]] - 1) * log((tmp/inputs$gamma[[j]]) +
                                                 1) - log(inputs$cost[[j]]))
            }
            else {
              inputs$V[[j]] = inputs$V[[j]] + inputs$avail[[j]] *
                ((inputs$alpha[[j]] - 1) * log((inputs$continuousChoice[[j]]/inputs$gamma[[j]]) +
                                                 1) - log(inputs$cost[[j]]))
            }
          }
          term1 = (1 - inputs$totalChosen) * log(inputs$sigma)
          logfi = list()
          for (j in 2:inputs$nAlt) {
            if (inputs$minX) {
              tmp <- inputs$continuousChoice[[j]] - (inputs$continuousChoice[[j]] >=
                                                       inputs$minConsumption[[j]]) * inputs$minConsumption[[j]]
              logfi[[j - 1]] = inputs$avail[[j]] * (log(1 -
                                                          inputs$alpha[[j]]) - log(tmp + inputs$gamma[[j]]))
            }
            else {
              logfi[[j - 1]] = inputs$avail[[j]] * (log(1 -
                                                          inputs$alpha[[j]]) - log(inputs$continuousChoice[[j]] +
                                                                                     inputs$gamma[[j]]))
            }
          }
          term2 = log(1 - inputs$alpha[[1]]) - log(inputs$continuousChoice[[1]])
          for (j in 2:inputs$nAlt) term2 = term2 + inputs$avail[[j]] *
            (logfi[[j - 1]] * inputs$discrete_choice[[j]])
          term3 = inputs$continuousChoice[[1]]/(1 - inputs$alpha[[1]])
          for (j in 2:inputs$nAlt) term3 = term3 + inputs$avail[[j]] *
            (inputs$cost[[j]]/exp(logfi[[j - 1]]) * inputs$discrete_choice[[j]])
          term3 = log(term3)
          term4_1 = inputs$V[[1]]/inputs$sigma
          term4_2 = exp(inputs$V[[1]]/inputs$sigma)
          for (j in 2:inputs$nAlt) {
            term4_1 = term4_1 + inputs$avail[[j]] * (inputs$V[[j]]/inputs$sigma *
                                                       inputs$discrete_choice[[j]])
            term4_2 = term4_2 + inputs$avail[[j]] * exp(inputs$V[[j]]/inputs$sigma)
          }
          term4_2 = inputs$totalChosen * log(term4_2)
          term4 = term4_1 - term4_2
          rm(term4_1, term4_2)
          term5 = lfactorial(inputs$totalChosen - 1)
          P = exp(term1 + term2 + term3 + term4 + term5)
          rm(term1, term2, term3, term4, term5)
          if (any(inputs$chosenUnavail))
            P <- apollo_setRows(P, inputs$chosenUnavail,
                                0)
          return(P)
        }
      if (!mdcev_settings$hasOutside)
        mdcev_settings$probs_MDCEV <- function(inputs) {
          for (j in 1:inputs$nAlt) {
            if (inputs$minX) {
              tmp <- inputs$continuousChoice[[j]] - (inputs$continuousChoice[[j]] >=
                                                       inputs$minConsumption[[j]]) * inputs$minConsumption[[j]]
              inputs$V[[j]] = inputs$V[[j]] + inputs$avail[[j]] *
                ((inputs$alpha[[j]] - 1) * log((tmp/inputs$gamma[[j]]) +
                                                 1) - log(inputs$cost[[j]]))
            }
            else {
              inputs$V[[j]] = inputs$V[[j]] + inputs$avail[[j]] *
                ((inputs$alpha[[j]] - 1) * log((inputs$continuousChoice[[j]]/inputs$gamma[[j]]) +
                                                 1) - log(inputs$cost[[j]]))
            }
          }
          term1 = (1 - inputs$totalChosen) * log(inputs$sigma)
          logfi = list()
          for (j in 1:inputs$nAlt) {
            if (inputs$minX) {
              tmp <- inputs$continuousChoice[[j]] - (inputs$continuousChoice[[j]] >=
                                                       inputs$minConsumption[[j]]) * inputs$minConsumption[[j]]
              logfi[[j]] = inputs$avail[[j]] * (log(1 -
                                                      inputs$alpha[[j]]) - log(tmp + inputs$gamma[[j]]))
            }
            else {
              logfi[[j]] = inputs$avail[[j]] * (log(1 -
                                                      inputs$alpha[[j]]) - log(inputs$continuousChoice[[j]] +
                                                                                 inputs$gamma[[j]]))
            }
          }
          term2 = 0
          for (j in 1:inputs$nAlt) term2 = term2 + inputs$avail[[j]] *
            (logfi[[j]] * inputs$discrete_choice[[j]])
          term3 = 0
          for (j in 1:inputs$nAlt) term3 = term3 + inputs$avail[[j]] *
            (inputs$cost[[j]]/exp(logfi[[j]]) * inputs$discrete_choice[[j]])
          term3 = log(term3)
          term4_1 = 0
          term4_2 = 0
          for (j in 1:inputs$nAlt) {
            term4_1 = term4_1 + inputs$avail[[j]] * (inputs$V[[j]]/inputs$sigma *
                                                       inputs$discrete_choice[[j]])
            term4_2 = term4_2 + inputs$avail[[j]] * exp(inputs$V[[j]]/inputs$sigma)
          }
          term4_2 = inputs$totalChosen * log(term4_2)
          term4 = term4_1 - term4_2
          rm(term4_1, term4_2)
          term5 = lfactorial(inputs$totalChosen - 1)
          P = exp(term1 + term2 + term3 + term4 + term5)
          rm(term1, term2, term3, term4, term5)
          if (any(inputs$chosenUnavail))
            P <- apollo_setRows(P, inputs$chosenUnavail,
                                0)
          return(P)
        }
      apollo_beta <- tryCatch(get("apollo_beta", envir = parent.frame(),
                                  inherits = TRUE), error = function(e) return(NULL))
      test <- !is.null(apollo_beta) && (functionality %in%
                                          c("preprocess", "gradient"))
      test <- test && all(sapply(mdcev_settings$V, is.function))
      test <- test && is.function(mdcev_settings$alpha)
      test <- test && is.function(mdcev_settings$gamma)
      test <- test && is.function(mdcev_settings$sigma)
      test <- test && apollo_inputs$apollo_control$analyticGrad
      mdcev_settings$gradient <- FALSE
      if (test) {
        mdcev_settings$dV <- apollo_dVdB(apollo_beta, apollo_inputs,
                                         mdcev_settings$V)
        mdcev_settings$dAlpha <- apollo_dVdB(apollo_beta,
                                             apollo_inputs, mdcev_settings$alpha)
        mdcev_settings$dGamma <- apollo_dVdB(apollo_beta,
                                             apollo_inputs, mdcev_settings$gamma)
        mdcev_settings$dSigma <- apollo_dVdB(apollo_beta,
                                             apollo_inputs, list(dSigma = mdcev_settings$sigma))[[1]]
      }
      rm(test)
      if (functionality == "preprocess") {
        mdcev_settings$V <- NULL
        mdcev_settings$alpha <- NULL
        mdcev_settings$gamma <- NULL
        mdcev_settings$sigma <- NULL
        return(mdcev_settings)
      }
    }
    testV <- any(sapply(mdcev_settings$V, is.function))
    testA <- any(sapply(mdcev_settings$alpha, is.function))
    testG <- any(sapply(mdcev_settings$gamma, is.function))
    testS <- is.function(mdcev_settings$sigma)
    if (testV)
      mdcev_settings$V = lapply(mdcev_settings$V, function(f) if (is.function(f))
        f()
        else f)
    if (testA)
      mdcev_settings$alpha = lapply(mdcev_settings$alpha, function(f) if (is.function(f))
        f()
        else f)
    if (testG)
      mdcev_settings$gamma = lapply(mdcev_settings$gamma, function(f) if (is.function(f))
        f()
        else f)
    if (testS)
      mdcev_settings$sigma = mdcev_settings$sigma()
    rm(testV, testA, testG, testS)
    mdcev_settings$V <- lapply(mdcev_settings$V, function(v) if (is.matrix(v) &&
                                                                 ncol(v) == 1)
      as.vector(v)
      else v)
    mdcev_settings$alpha <- lapply(mdcev_settings$alpha, function(a) if (is.matrix(a) &&
                                                                         ncol(a) == 1)
      as.vector(a)
      else a)
    mdcev_settings$gamma <- lapply(mdcev_settings$gamma, function(g) if (is.matrix(g) &&
                                                                         ncol(g) == 1)
      as.vector(g)
      else g)
    if (is.matrix(mdcev_settings$sigma) && ncol(mdcev_settings$sigma) ==
        1)
      mdcev_settings$sigma <- as.vector(mdcev_settings$sigma)
    if (mdcev_settings$hasOutside) {
      if (is.null(mdcev_settings$gamma[[mdcev_settings$outside]]))
        mdcev_settings$gamma[[mdcev_settings$outside]] <- 1
      tmp <- which(names(mdcev_settings$V) == mdcev_settings$outside)
      if (length(tmp) > 0)
        names(mdcev_settings$V)[tmp] <- "outside"
      tmp <- which(names(mdcev_settings$alpha) == mdcev_settings$outside)
      if (length(tmp) > 0)
        names(mdcev_settings$alpha)[tmp] <- "outside"
      tmp <- which(names(mdcev_settings$gamma) == mdcev_settings$outside)
      if (length(tmp) > 0)
        names(mdcev_settings$gamma)[tmp] <- "outside"
      rm(tmp)
    }
    if (any(mdcev_settings$alternatives != names(mdcev_settings$V)))
      mdcev_settings$V <- mdcev_settings$V[mdcev_settings$alternatives]
    if (any(mdcev_settings$alternatives != names(mdcev_settings$alpha)))
      mdcev_settings$alpha <- mdcev_settings$alpha[mdcev_settings$alternatives]
    if (any(mdcev_settings$alternatives != names(mdcev_settings$gamma)))
      mdcev_settings$gamma <- mdcev_settings$gamma[mdcev_settings$alternatives]
    if (!all(mdcev_settings$rows)) {
      mdcev_settings$V <- lapply(mdcev_settings$V, apollo_keepRows,
                                 r = mdcev_settings$rows)
      mdcev_settings$alpha <- lapply(mdcev_settings$alpha,
                                     apollo_keepRows, r = mdcev_settings$rows)
      mdcev_settings$gamma <- lapply(mdcev_settings$gamma,
                                     apollo_keepRows, r = mdcev_settings$rows)
      mdcev_settings$sigma <- apollo_keepRows(mdcev_settings$sigma,
                                              r = mdcev_settings$rows)
    }
  
    if (functionality %in% c("estimate", "conditionals",
                             "output", "components")) {
      L <- mdcev_settings$probs_MDCEV(mdcev_settings)
      if (any(!mdcev_settings$rows))
        L <- apollo_insertRows(L, mdcev_settings$rows, 1)
      return(L)
    }
    if (functionality %in% c("prediction", "raw")) {
      s <- mdcev_settings
      rm(mdcev_settings)
      if (!is.vector(s$sigma))
        stop("Forecasting not available for MDCEV with random sigma")
      if (!is.null(apollo_inputs$apollo_control$seed))
        seed <- apollo_inputs$apollo_control$seed + 5
      else seed <- 13 + 5
      set.seed(seed)
      tmp1 <- -log(-log(apollo_mlhs(s$nRep, s$nAlt, s$nObs)))
      epsL <- vector(mode = "list", length = s$nRep)
      tmp2 <- (0:(s$nObs - 1)) * s$nRep
      for (r in 1:s$nRep) epsL[[r]] <- tmp1[r + tmp2, ]
      rm(tmp1, tmp2)
      if (s$rawPrediction)
        XX <- array(NA, dim = c(s$nObs, s$nAlt, s$nRep),
                    dimnames = list(NULL, s$alternatives, NULL))
      else {
        Xm <- matrix(0, nrow = s$nObs, ncol = s$nAlt)
        Xv <- matrix(0, nrow = s$nObs, ncol = s$nAlt)
        Mm <- matrix(0, nrow = s$nObs, ncol = s$nAlt)
        Mv <- matrix(0, nrow = s$nObs, ncol = s$nAlt)
        Em <- matrix(0, nrow = s$nObs, ncol = s$nAlt)
        Ev <- matrix(0, nrow = s$nObs, ncol = s$nAlt)
      }
      X <- matrix(0, nrow = s$nObs, ncol = s$nAlt)
      extractDraw <- function(b, iInter, iIntra) {
        if (!is.list(b)) {
          if (is.vector(b) && length(b) == 1)
            return(rep(b, s$nObs))
          if (is.vector(b))
            return(b)
          if (is.matrix(b))
            return(b[, iInter])
          if (is.array(b) && length(dim(b)) == 3)
            return(b[, iInter, iIntra])
        }
        else {
          ans <- lapply(b, extractDraw, iInter = iInter,
                        iIntra = iIntra)
          ans <- do.call(cbind, ans)
          return(ans)
        }
      }
      if (anyNA(apollo_inputs$apollo_draws)) {
        nInter <- 0
        nIntra <- 0
      }
      else {
        nInter <- apollo_inputs$apollo_draws$interNDraws
        nIntra <- apollo_inputs$apollo_draws$intraNDraws
      }
      if (nInter == 0)
        nInter <- 1
      if (nIntra == 0)
        nIntra <- 1
      step1 <- ceiling(s$nRep/10)
      step2 <- ceiling(s$nRep/2)
      avail <- sapply(s$avail, function(a) if (length(a) ==
                                               1)
        rep(a, s$nObs)
        else a)
      cost <- sapply(s$cost, function(x) if (length(x) == 1)
        rep(x, s$nObs)
        else x)
      if (length(s$budget) == 1)
        budget <- rep(s$budget, s$nObs)
      else budget <- s$budget
      if (!is.null(apollo_inputs$silent))
        silent <- apollo_inputs$silent
      else silent <- FALSE
      if (!silent)
        cat("0%")
      if (s$hasOutside) {
        for (r in 1:s$nRep) {
          X <- 0 * X
          for (iInter in 1:nInter) {
            Xintra <- 0 * X
            for (iIntra in 1:nIntra) {
              phiP <- extractDraw(s$sigma, iInter, iIntra) *
                epsL[[r]]
              phiP <- exp(extractDraw(s$V, iInter, iIntra) +
                            phiP) * avail/cost
              gamma <- extractDraw(s$gamma, iInter, iIntra)
              alpha <- extractDraw(s$alpha, iInter, iIntra)
              for (i in 1:s$nObs) {
                p <- cost[i, 2:s$nAlt]
                b <- budget[i]
                g <- gamma[i, 2:s$nAlt]
                a0 <- alpha[i, 1]
                ak <- alpha[i, 2:s$nAlt]
                ph0 <- phiP[i, 1]
                phk <- phiP[i, 2:s$nAlt]
                orderofV = rank(-phk)
                M = 1
                stopping = FALSE
                while (!stopping) {
                  use = orderofV < M
                  lambda_1 = b + sum(p * g * use)
                  lambda_21 = ph0^(1/(1 - a0))
                  lambda_22 = sum(p * g * use * phk^(1/(1 -
                                                          a0)))
                  lambda_2 = lambda_21 + lambda_22
                  lambda = (lambda_1/lambda_2)^(a0 - 1)
                  if (M > sum(phk > lambda) || M > sum(phk >
                                                       0)) {
                    x0_1 = lambda_21 * lambda_1
                    Xintra[i, 1] = Xintra[i, 1] + x0_1/lambda_2
                    xk_1 = phk^(1/(1 - ak)) * lambda_1
                    Xintra[i, 2:s$nAlt] = Xintra[i, 2:s$nAlt] +
                      use * (xk_1/lambda_2 - 1) * g
                    stopping = TRUE
                  }
                  else M <- M + 1
                }
              }
            }
            X <- X + Xintra/nIntra
          }
          X <- X/nInter
          if (s$rawPrediction)
            XX[, , r] <- X
          else {
            Xm <- Xm + X/s$nRep
            Mm <- Mm + (X > 0)/s$nRep
            Em <- Em + X * cost/s$nRep
            Xv <- Xv + apply(X, MARGIN = 2, function(v) (v -
                                                           mean(v))^2)/s$nRep
            Mv <- Mv + apply(X > 0, MARGIN = 2, function(m) (m -
                                                               mean(m))^2)/s$nRep
            Ev <- Ev + apply(X * cost, MARGIN = 2, function(e) (e -
                                                                  mean(e))^2)/s$nRep
          }
          if (!silent) {
            if (r%%step2 == 0)
              cat(round(100 * r/s$nRep, 0), "%",
                  sep = "")
            else {
              if (r%%step1 == 0)
                cat(".")
            }
          }
        }
      }
      else {
        for (r in 1:s$nRep) {
          X <- 0 * X
          for (iInter in 1:nInter) {
            Xintra <- 0 * X
            for (iIntra in 1:nIntra) {
              phiP <- extractDraw(s$sigma, iInter, iIntra) *
                epsL[[r]]
              phiP <- exp(extractDraw(s$V, iInter, iIntra) +
                            phiP) * avail/cost
              gamma <- extractDraw(s$gamma, iInter, iIntra)
              alpha <- extractDraw(s$alpha, iInter, iIntra)
              for (i in 1:s$nObs) {
                p <- cost[i, ]
                b <- budget[i]
                g <- gamma[i, ]
                ak <- alpha[i, ]
                phk <- phiP[i, ]
                orderofV = rank(-phk)
                M = 1
                stopping = FALSE
                while (!stopping) {
                  use = orderofV < M
                  lambda_1 = b + sum(p * g * use)
                  lambda_2 = sum(p * g * use * phk^(1/(1 -
                                                         ak)))
                  lambda = (lambda_1/lambda_2)^(ak - 1)
                  if (M > sum(phk > lambda) || M > sum(phk >
                                                       0)) {
                    xk_1 = phk^(1/(1 - ak)) * lambda_1
                    Xintra[i, ] = Xintra[i, ] + use * (xk_1/lambda_2 -
                                                         1) * g
                    stopping = TRUE
                  }
                  else M <- M + 1
                }
              }
            }
            X <- X + Xintra/nIntra
          }
          X <- X/nInter
          if (s$rawPrediction)
            XX[, , r] <- X
          else {
            Xm <- Xm + X/s$nRep
            Mm <- Mm + (X > 0)/s$nRep
            Em <- Em + X * cost/s$nRep
            Xv <- Xv + apply(X, MARGIN = 2, function(v) (v -
                                                           mean(v))^2)/s$nRep
            Mv <- Mv + apply(X > 0, MARGIN = 2, function(m) (m -
                                                               mean(m))^2)/s$nRep
            Ev <- Ev + apply(X * cost, MARGIN = 2, function(e) (e -
                                                                  mean(e))^2)/s$nRep
          }
          if (!silent) {
            if (r%%step2 == 0)
              cat(round(100 * r/s$nRep, 0), "%",
                  sep = "")
            else {
              if (r%%step1 == 0)
                cat(".")
            }
          }
        }
      }
      if (!silent)
        cat("\n")
      if (s$rawPrediction)
        out <- XX
      else {
        out <- cbind(Xm, sqrt(Xv), Mm, sqrt(Mv), Em, sqrt(Ev))
        out <- apollo_insertRows(out, s$rows, NA)
        colN <- c("cont_mean", "cont_sd", "disc_mean",
                  "disc_sd", "expe_mean", "expe_sd")
        colnames(out) <- paste(names(s$continuousChoice),
                               rep(colN, each = s$nAlt), sep = "_")
      }
      return(out)
    }
  
  }
  
  
  #######################Segment 1  
  apollo_beta = apollo_beta1
  
  apollo_probabilities1=function(apollo_beta = apollo_beta1, apollo_inputs, functionality="estimate"){
    
    ### Attach inputs and detach after function exit
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    
    ### Create list of probabilities P
    P = list()
    
    ### Define individual alternatives
    alternatives = c("v1", 
                     "v2", 
                     "v3", 
                     "v4", 
                     "v5", 
                     "v6", 
                     "v7", 
                     "v8" 
    )
    
    ### Define availabilities
    avail = list(v1  = availnew1,   
                 v2     = availnew2,
                 v3   = availnew3,
                 v4 = availnew4,
                 v5 = availnew5,
                 v6   = availnew6,
                 v7 = availnew7,
                 v8 = availnew8 
    )
    
    ### Define continuous consumption for individual alternatives
    continuousChoice = list(v1  =clicknew1 * availnew1,
                            v2     =clicknew2 * availnew2,
                            v3   =clicknew3 * availnew3,
                            v4 =clicknew4 * availnew4,
                            v5 =clicknew5 * availnew5,
                            v6   =clicknew6 * availnew6,
                            v7 =clicknew7 * availnew7,
                            v8 =clicknew8 * availnew8#,
    )
    
    
    ### Define alpha parameters
    alpha = list(v1  = 1e-3 , 
                 v2     = 1e-3 , 
                 v3   = 1e-3 , 
                 v4 = 1e-3 , 
                 v5 = 1e-3 ,
                 v6   = 1e-3 ,
                 v7  = 1e-3 , 
                 v8 = 1e-3 #, 
    )
    
    
    ### Define costs for individual alternatives
    cost = list(v1      = 1, 
                v2         = 1,
                v3       = 1,
                v4     = 1,
                v5 = 1,
                v6       = 1, 
                v7      = 1,
                v8     = 1 #,
    )
    
    ### Define budget
    budget = clicktotal
    
    emdc_settings <- list(continuousChoice = continuousChoice, 
                          avail            = avail,
                          budget           = budget,
                          sigma            = 0.99, 
                          cost             = cost)
    
    
    ### ### Compute class-specific utilities
    V = list()
    
    V[["v1"    ]] =   delta_v1*(rootcat1==1) +  delta_v2*(rootcat1==2) +  delta_v3*(rootcat1==3) +  
      delta_v4*(rootcat1==4) +  delta_v5*(rootcat1==5) + 
      bweekend*weekend + bnpv*npvnew1 + bnclick*nclicknew1 + bnbuy*nbuynew1 + 
      nrccpv*rccpvnew1 + nrccclick*rccclicknew1 + nrccbuy*rccbuynew1 + badstockpsai*ad11 + bresidual*res1 + bvariety*Lag_Distinct
    
    V[["v2"    ]] = delta_v1*(rootcat2==1) +  delta_v2*(rootcat2==2) +  delta_v3*(rootcat2==3) +  
      delta_v4*(rootcat2==4) +  delta_v5*(rootcat2==5) + 
      bweekend*weekend + bnpv*npvnew2 + bnclick*nclicknew2 + bnbuy*nbuynew2 + 
      nrccpv*rccpvnew2 + nrccclick*rccclicknew2 + nrccbuy*rccbuynew2 + badstockpsai*ad21 + bresidual*res2 + bvariety*Lag_Distinct
    
    V[["v3"  ]] = delta_v1*(rootcat3==1) +  delta_v2*(rootcat3==2) +  delta_v3*(rootcat3==3) +  
      delta_v4*(rootcat3==4) +  delta_v5*(rootcat3==5) + 
      bweekend*weekend + bnpv*npvnew3 + bnclick*nclicknew3 + bnbuy*nbuynew3 + 
      nrccpv*rccpvnew3 + nrccclick*rccclicknew3 + nrccbuy*rccbuynew3 + badstockpsai*ad31 + bresidual*res3 + bvariety*Lag_Distinct
    
    V[["v4"]] = delta_v1*(rootcat4==1) +  delta_v2*(rootcat4==2) +  delta_v3*(rootcat4==3) +  
      delta_v4*(rootcat4==4) +  delta_v5*(rootcat4==5) + 
      bweekend*weekend + bnpv*npvnew4 + bnclick*nclicknew4 + bnbuy*nbuynew4 + 
      nrccpv*rccpvnew4 + nrccclick*rccclicknew4 + nrccbuy*rccbuynew4 + badstockpsai*ad41 + bresidual*res4 + bvariety*Lag_Distinct
    
    V[["v5"]] = delta_v1*(rootcat5==1) +  delta_v2*(rootcat5==2) +  delta_v3*(rootcat5==3) +  
      delta_v4*(rootcat5==4) +  delta_v5*(rootcat5==5) + 
      bweekend*weekend + bnpv*npvnew5 + bnclick*nclicknew5 + bnbuy*nbuynew5 + 
      nrccpv*rccpvnew5 + nrccclick*rccclicknew5 + nrccbuy*rccbuynew5 + badstockpsai*ad51 + bresidual*res5 + bvariety*Lag_Distinct
    
    V[["v6"  ]] = delta_v1*(rootcat6==1) +  delta_v2*(rootcat6==2) +  delta_v3*(rootcat6==3) +  
      delta_v4*(rootcat6==4) +  delta_v5*(rootcat6==5) + 
      bweekend*weekend + bnpv*npvnew6 + bnclick*nclicknew6 + bnbuy*nbuynew6 + 
      nrccpv*rccpvnew6 + nrccclick*rccclicknew6 + nrccbuy*rccbuynew6 + badstockpsai*ad61 + bresidual*res6 + bvariety*Lag_Distinct
    
    V[["v7"]] = delta_v1*(rootcat7==1) +  delta_v2*(rootcat7==2) +  delta_v3*(rootcat7==3) +  
      delta_v4*(rootcat7==4) +  delta_v5*(rootcat7==5) + 
      bweekend*weekend + bnpv*npvnew7 + bnclick*nclicknew7 + bnbuy*nbuynew7 + 
      nrccpv*rccpvnew7 + nrccclick*rccclicknew7 + nrccbuy*rccbuynew7 + badstockpsai*ad71 + bresidual*res7 + bvariety*Lag_Distinct
    
    V[["v8"]] = delta_v1*(rootcat8==1) +  delta_v2*(rootcat8==2) +  delta_v3*(rootcat8==3) +  
      delta_v4*(rootcat8==4) +  delta_v5*(rootcat8==5) + 
      bweekend*weekend + bnpv*npvnew8 + bnclick*nclicknew8 + bnbuy*nbuynew8 + 
      nrccpv*rccpvnew8 + nrccclick*rccclicknew8 + nrccbuy*rccbuynew8 + badstockpsai*ad81 + bresidual*res8 + bvariety*Lag_Distinct
    
    
    
    ### Define gamma parameters
    gamma = list(v1      = exp( gamma_v1*(rootcat1==1) +  gamma_v2*(rootcat1==2) +  gamma_v3*(rootcat1==3) +  
                                  gamma_v4*(rootcat1==4) +  gamma_v5*(rootcat1==5) ),
                 
                 v2      = exp( gamma_v1*(rootcat2==1) +  gamma_v2*(rootcat2==2) +  gamma_v3*(rootcat2==3) +  
                                  gamma_v4*(rootcat2==4) +  gamma_v5*(rootcat2==5) ), 
                 
                 
                 v3      = exp( gamma_v1*(rootcat3==1) +  gamma_v2*(rootcat3==2) +  gamma_v3*(rootcat3==3) +  
                                  gamma_v4*(rootcat3==4) +  gamma_v5*(rootcat3==5) ),
                 
                 v4      = exp( gamma_v1*(rootcat4==1) +  gamma_v2*(rootcat4==2) +  gamma_v3*(rootcat4==3) +  
                                  gamma_v4*(rootcat4==4) +  gamma_v5*(rootcat4==5) ), 
                 
                 v5      = exp( gamma_v1*(rootcat5==1) +  gamma_v2*(rootcat5==2) +  gamma_v3*(rootcat5==3) +  
                                  gamma_v4*(rootcat5==4) +  gamma_v5*(rootcat5==5) ),
                 
                 v6      = exp( gamma_v1*(rootcat6==1) +  gamma_v2*(rootcat6==2) +  gamma_v3*(rootcat6==3) +  
                                  gamma_v4*(rootcat6==4) +  gamma_v5*(rootcat6==5) ),
                 
                 v7      = exp( gamma_v1*(rootcat7==1) +  gamma_v2*(rootcat7==2) +  gamma_v3*(rootcat7==3) +  
                                  gamma_v4*(rootcat7==4) +  gamma_v5*(rootcat7==5) ), 
                 
                 v8      = exp( gamma_v1*(rootcat8==1) +  gamma_v2*(rootcat8==2) +  gamma_v3*(rootcat8==3) +  
                                  gamma_v4*(rootcat8==4) +  gamma_v5*(rootcat8==5) )
                 
    )
    
    delta = list( list(0,0,0,0,0,0,0 ,0 ) , 
                  list(   d12 * ( (rootcat2 == 1) * (rootcat1 == 2) + (rootcat2 == 2) * (rootcat1 == 1) )    +
                            d13 * ( (rootcat2 == 1) * (rootcat1 == 3) + (rootcat2 == 3) * (rootcat1 == 1) )  +
                            d14 * ( (rootcat2 == 1) * (rootcat1 == 4) + (rootcat2 == 4) * (rootcat1 == 1) )  +
                            d15 * ( (rootcat2 == 1) * (rootcat1 == 5) + (rootcat2 == 5) * (rootcat1 == 1) )  +
                            d23 * ( (rootcat2 == 2) * (rootcat1 == 3) + (rootcat2 == 3) * (rootcat1 == 2) )  +
                            d24 * ( (rootcat2 == 2) * (rootcat1 == 4) + (rootcat2 == 4) * (rootcat1 == 2) )  +
                            d25 * ( (rootcat2 == 2) * (rootcat1 == 5) + (rootcat2 == 5) * (rootcat1 == 2) )  +
                            d34 * ( (rootcat2 == 3) * (rootcat1 == 4) + (rootcat2 == 4) * (rootcat1 == 3) )  +
                            d35 * ( (rootcat2 == 3) * (rootcat1 == 5) + (rootcat2 == 5) * (rootcat1 == 3) )  +
                            d45 * ( (rootcat2 == 4) * (rootcat1 == 5) + (rootcat2 == 5) * (rootcat1 == 4) )  ,
                          0,0,0,0,0,0,0 ) ,
                  
                  list(   d12 * ( (rootcat3 == 1) * (rootcat1 == 2) + (rootcat3 == 2) * (rootcat1 == 1) )    +
                            d13 * ( (rootcat3 == 1) * (rootcat1 == 3) + (rootcat3 == 3) * (rootcat1 == 1) )  +
                            d14 * ( (rootcat3 == 1) * (rootcat1 == 4) + (rootcat3 == 4) * (rootcat1 == 1) )  +
                            d15 * ( (rootcat3 == 1) * (rootcat1 == 5) + (rootcat3 == 5) * (rootcat1 == 1) )  +
                            d23 * ( (rootcat3 == 2) * (rootcat1 == 3) + (rootcat3 == 3) * (rootcat1 == 2) )  +
                            d24 * ( (rootcat3 == 2) * (rootcat1 == 4) + (rootcat3 == 4) * (rootcat1 == 2) )  +
                            d25 * ( (rootcat3 == 2) * (rootcat1 == 5) + (rootcat3 == 5) * (rootcat1 == 2) )  +
                            d34 * ( (rootcat3 == 3) * (rootcat1 == 4) + (rootcat3 == 4) * (rootcat1 == 3) )  +
                            d35 * ( (rootcat3 == 3) * (rootcat1 == 5) + (rootcat3 == 5) * (rootcat1 == 3) )  +
                            d45 * ( (rootcat3 == 4) * (rootcat1 == 5) + (rootcat3 == 5) * (rootcat1 == 4) )  ,
                          
                          d12 * ( (rootcat3 == 1) * (rootcat2 == 2) + (rootcat3 == 2) * (rootcat2 == 1) )    +
                            d13 * ( (rootcat3 == 1) * (rootcat2 == 3) + (rootcat3 == 3) * (rootcat2 == 1) )  +
                            d14 * ( (rootcat3 == 1) * (rootcat2 == 4) + (rootcat3 == 4) * (rootcat2 == 1) )  +
                            d15 * ( (rootcat3 == 1) * (rootcat2 == 5) + (rootcat3 == 5) * (rootcat2 == 1) )  +
                            d23 * ( (rootcat3 == 2) * (rootcat2 == 3) + (rootcat3 == 3) * (rootcat2 == 2) )  +
                            d24 * ( (rootcat3 == 2) * (rootcat2 == 4) + (rootcat3 == 4) * (rootcat2 == 2) )  +
                            d25 * ( (rootcat3 == 2) * (rootcat2 == 5) + (rootcat3 == 5) * (rootcat2 == 2) )  +
                            d34 * ( (rootcat3 == 3) * (rootcat2 == 4) + (rootcat3 == 4) * (rootcat2 == 3) )  +
                            d35 * ( (rootcat3 == 3) * (rootcat2 == 5) + (rootcat3 == 5) * (rootcat2 == 3) )  +
                            d45 * ( (rootcat3 == 4) * (rootcat2 == 5) + (rootcat3 == 5) * (rootcat2 == 4) )  ,
                          0,0,0,0,0,0 ) , 
                  
                  list(     d12 * ( (rootcat4 == 1) * (rootcat1 == 2) + (rootcat4 == 2) * (rootcat1 == 1) )    +
                              d13 * ( (rootcat4 == 1) * (rootcat1 == 3) + (rootcat4 == 3) * (rootcat1 == 1) )  +
                              d14 * ( (rootcat4 == 1) * (rootcat1 == 4) + (rootcat4 == 4) * (rootcat1 == 1) )  +
                              d15 * ( (rootcat4 == 1) * (rootcat1 == 5) + (rootcat4 == 5) * (rootcat1 == 1) )  +
                              d23 * ( (rootcat4 == 2) * (rootcat1 == 3) + (rootcat4 == 3) * (rootcat1 == 2) )  +
                              d24 * ( (rootcat4 == 2) * (rootcat1 == 4) + (rootcat4 == 4) * (rootcat1 == 2) )  +
                              d25 * ( (rootcat4 == 2) * (rootcat1 == 5) + (rootcat4 == 5) * (rootcat1 == 2) )  +
                              d34 * ( (rootcat4 == 3) * (rootcat1 == 4) + (rootcat4 == 4) * (rootcat1 == 3) )  +
                              d35 * ( (rootcat4 == 3) * (rootcat1 == 5) + (rootcat4 == 5) * (rootcat1 == 3) )  +
                              d45 * ( (rootcat4 == 4) * (rootcat1 == 5) + (rootcat4 == 5) * (rootcat1 == 4) )  ,
                            
                            d12 * ( (rootcat4 == 1) * (rootcat2 == 2) + (rootcat4 == 2) * (rootcat2 == 1) )    +
                              d13 * ( (rootcat4 == 1) * (rootcat2 == 3) + (rootcat4 == 3) * (rootcat2 == 1) )  +
                              d14 * ( (rootcat4 == 1) * (rootcat2 == 4) + (rootcat4 == 4) * (rootcat2 == 1) )  +
                              d15 * ( (rootcat4 == 1) * (rootcat2 == 5) + (rootcat4 == 5) * (rootcat2 == 1) )  +
                              d23 * ( (rootcat4 == 2) * (rootcat2 == 3) + (rootcat4 == 3) * (rootcat2 == 2) )  +
                              d24 * ( (rootcat4 == 2) * (rootcat2 == 4) + (rootcat4 == 4) * (rootcat2 == 2) )  +
                              d25 * ( (rootcat4 == 2) * (rootcat2 == 5) + (rootcat4 == 5) * (rootcat2 == 2) )  +
                              d34 * ( (rootcat4 == 3) * (rootcat2 == 4) + (rootcat4 == 4) * (rootcat2 == 3) )  +
                              d35 * ( (rootcat4 == 3) * (rootcat2 == 5) + (rootcat4 == 5) * (rootcat2 == 3) )  +
                              d45 * ( (rootcat4 == 4) * (rootcat2 == 5) + (rootcat4 == 5) * (rootcat2 == 4) )  ,
                            
                            d12 * ( (rootcat4 == 1) * (rootcat3 == 2) + (rootcat4 == 2) * (rootcat3 == 1) )    +
                              d13 * ( (rootcat4 == 1) * (rootcat3 == 3) + (rootcat4 == 3) * (rootcat3 == 1) )  +
                              d14 * ( (rootcat4 == 1) * (rootcat3 == 4) + (rootcat4 == 4) * (rootcat3 == 1) )  +
                              d15 * ( (rootcat4 == 1) * (rootcat3 == 5) + (rootcat4 == 5) * (rootcat3 == 1) )  +
                              d23 * ( (rootcat4 == 2) * (rootcat3 == 3) + (rootcat4 == 3) * (rootcat3 == 2) )  +
                              d24 * ( (rootcat4 == 2) * (rootcat3 == 4) + (rootcat4 == 4) * (rootcat3 == 2) )  +
                              d25 * ( (rootcat4 == 2) * (rootcat3 == 5) + (rootcat4 == 5) * (rootcat3 == 2) )  +
                              d34 * ( (rootcat4 == 3) * (rootcat3 == 4) + (rootcat4 == 4) * (rootcat3 == 3) )  +
                              d35 * ( (rootcat4 == 3) * (rootcat3 == 5) + (rootcat4 == 5) * (rootcat3 == 3) )  +
                              d45 * ( (rootcat4 == 4) * (rootcat3 == 5) + (rootcat4 == 5) * (rootcat3 == 4) )  ,
                            0,0,0,0,0 ) ,
                  
                  list(       d12 * ( (rootcat5 == 1) * (rootcat1 == 2) + (rootcat5 == 2) * (rootcat1 == 1) )    +
                                d13 * ( (rootcat5 == 1) * (rootcat1 == 3) + (rootcat5 == 3) * (rootcat1 == 1) )  +
                                d14 * ( (rootcat5 == 1) * (rootcat1 == 4) + (rootcat5 == 4) * (rootcat1 == 1) )  +
                                d15 * ( (rootcat5 == 1) * (rootcat1 == 5) + (rootcat5 == 5) * (rootcat1 == 1) )  +
                                d23 * ( (rootcat5 == 2) * (rootcat1 == 3) + (rootcat5 == 3) * (rootcat1 == 2) )  +
                                d24 * ( (rootcat5 == 2) * (rootcat1 == 4) + (rootcat5 == 4) * (rootcat1 == 2) )  +
                                d25 * ( (rootcat5 == 2) * (rootcat1 == 5) + (rootcat5 == 5) * (rootcat1 == 2) )  +
                                d34 * ( (rootcat5 == 3) * (rootcat1 == 4) + (rootcat5 == 4) * (rootcat1 == 3) )  +
                                d35 * ( (rootcat5 == 3) * (rootcat1 == 5) + (rootcat5 == 5) * (rootcat1 == 3) )  +
                                d45 * ( (rootcat5 == 4) * (rootcat1 == 5) + (rootcat5 == 5) * (rootcat1 == 4) )  ,
                              
                              d12 * ( (rootcat5 == 1) * (rootcat2 == 2) + (rootcat5 == 2) * (rootcat2 == 1) )    +
                                d13 * ( (rootcat5 == 1) * (rootcat2 == 3) + (rootcat5 == 3) * (rootcat2 == 1) )  +
                                d14 * ( (rootcat5 == 1) * (rootcat2 == 4) + (rootcat5 == 4) * (rootcat2 == 1) )  +
                                d15 * ( (rootcat5 == 1) * (rootcat2 == 5) + (rootcat5 == 5) * (rootcat2 == 1) )  +
                                d23 * ( (rootcat5 == 2) * (rootcat2 == 3) + (rootcat5 == 3) * (rootcat2 == 2) )  +
                                d24 * ( (rootcat5 == 2) * (rootcat2 == 4) + (rootcat5 == 4) * (rootcat2 == 2) )  +
                                d25 * ( (rootcat5 == 2) * (rootcat2 == 5) + (rootcat5 == 5) * (rootcat2 == 2) )  +
                                d34 * ( (rootcat5 == 3) * (rootcat2 == 4) + (rootcat5 == 4) * (rootcat2 == 3) )  +
                                d35 * ( (rootcat5 == 3) * (rootcat2 == 5) + (rootcat5 == 5) * (rootcat2 == 3) )  +
                                d45 * ( (rootcat5 == 4) * (rootcat2 == 5) + (rootcat5 == 5) * (rootcat2 == 4) )  ,
                              
                              d12 * ( (rootcat5 == 1) * (rootcat3 == 2) + (rootcat5 == 2) * (rootcat3 == 1) )    +
                                d13 * ( (rootcat5 == 1) * (rootcat3 == 3) + (rootcat5 == 3) * (rootcat3 == 1) )  +
                                d14 * ( (rootcat5 == 1) * (rootcat3 == 4) + (rootcat5 == 4) * (rootcat3 == 1) )  +
                                d15 * ( (rootcat5 == 1) * (rootcat3 == 5) + (rootcat5 == 5) * (rootcat3 == 1) )  +
                                d23 * ( (rootcat5 == 2) * (rootcat3 == 3) + (rootcat5 == 3) * (rootcat3 == 2) )  +
                                d24 * ( (rootcat5 == 2) * (rootcat3 == 4) + (rootcat5 == 4) * (rootcat3 == 2) )  +
                                d25 * ( (rootcat5 == 2) * (rootcat3 == 5) + (rootcat5 == 5) * (rootcat3 == 2) )  +
                                d34 * ( (rootcat5 == 3) * (rootcat3 == 4) + (rootcat5 == 4) * (rootcat3 == 3) )  +
                                d35 * ( (rootcat5 == 3) * (rootcat3 == 5) + (rootcat5 == 5) * (rootcat3 == 3) )  +
                                d45 * ( (rootcat5 == 4) * (rootcat3 == 5) + (rootcat5 == 5) * (rootcat3 == 4) )  ,
                              
                              d12 * ( (rootcat5 == 1) * (rootcat4 == 2) + (rootcat5 == 2) * (rootcat4 == 1) )    +
                                d13 * ( (rootcat5 == 1) * (rootcat4 == 3) + (rootcat5 == 3) * (rootcat4 == 1) )  +
                                d14 * ( (rootcat5 == 1) * (rootcat4 == 4) + (rootcat5 == 4) * (rootcat4 == 1) )  +
                                d15 * ( (rootcat5 == 1) * (rootcat4 == 5) + (rootcat5 == 5) * (rootcat4 == 1) )  +
                                d23 * ( (rootcat5 == 2) * (rootcat4 == 3) + (rootcat5 == 3) * (rootcat4 == 2) )  +
                                d24 * ( (rootcat5 == 2) * (rootcat4 == 4) + (rootcat5 == 4) * (rootcat4 == 2) )  +
                                d25 * ( (rootcat5 == 2) * (rootcat4 == 5) + (rootcat5 == 5) * (rootcat4 == 2) )  +
                                d34 * ( (rootcat5 == 3) * (rootcat4 == 4) + (rootcat5 == 4) * (rootcat4 == 3) )  +
                                d35 * ( (rootcat5 == 3) * (rootcat4 == 5) + (rootcat5 == 5) * (rootcat4 == 3) )  +
                                d45 * ( (rootcat5 == 4) * (rootcat4 == 5) + (rootcat5 == 5) * (rootcat4 == 4) )  ,
                              
                              0,0,0,0) ,
                  
                  list(       d12 * ( (rootcat6 == 1) * (rootcat1 == 2) + (rootcat6 == 2) * (rootcat1 == 1) )    +
                                d13 * ( (rootcat6 == 1) * (rootcat1 == 3) + (rootcat6 == 3) * (rootcat1 == 1) )  +
                                d14 * ( (rootcat6 == 1) * (rootcat1 == 4) + (rootcat6 == 4) * (rootcat1 == 1) )  +
                                d15 * ( (rootcat6 == 1) * (rootcat1 == 5) + (rootcat6 == 5) * (rootcat1 == 1) )  +
                                d23 * ( (rootcat6 == 2) * (rootcat1 == 3) + (rootcat6 == 3) * (rootcat1 == 2) )  +
                                d24 * ( (rootcat6 == 2) * (rootcat1 == 4) + (rootcat6 == 4) * (rootcat1 == 2) )  +
                                d25 * ( (rootcat6 == 2) * (rootcat1 == 5) + (rootcat6 == 5) * (rootcat1 == 2) )  +
                                d34 * ( (rootcat6 == 3) * (rootcat1 == 4) + (rootcat6 == 4) * (rootcat1 == 3) )  +
                                d35 * ( (rootcat6 == 3) * (rootcat1 == 5) + (rootcat6 == 5) * (rootcat1 == 3) )  +
                                d45 * ( (rootcat6 == 4) * (rootcat1 == 5) + (rootcat6 == 5) * (rootcat1 == 4) )  ,
                              
                              d12 * ( (rootcat6 == 1) * (rootcat2 == 2) + (rootcat6 == 2) * (rootcat2 == 1) )    +
                                d13 * ( (rootcat6 == 1) * (rootcat2 == 3) + (rootcat6 == 3) * (rootcat2 == 1) )  +
                                d14 * ( (rootcat6 == 1) * (rootcat2 == 4) + (rootcat6 == 4) * (rootcat2 == 1) )  +
                                d15 * ( (rootcat6 == 1) * (rootcat2 == 5) + (rootcat6 == 5) * (rootcat2 == 1) )  +
                                d23 * ( (rootcat6 == 2) * (rootcat2 == 3) + (rootcat6 == 3) * (rootcat2 == 2) )  +
                                d24 * ( (rootcat6 == 2) * (rootcat2 == 4) + (rootcat6 == 4) * (rootcat2 == 2) )  +
                                d25 * ( (rootcat6 == 2) * (rootcat2 == 5) + (rootcat6 == 5) * (rootcat2 == 2) )  +
                                d34 * ( (rootcat6 == 3) * (rootcat2 == 4) + (rootcat6 == 4) * (rootcat2 == 3) )  +
                                d35 * ( (rootcat6 == 3) * (rootcat2 == 5) + (rootcat6 == 5) * (rootcat2 == 3) )  +
                                d45 * ( (rootcat6 == 4) * (rootcat2 == 5) + (rootcat6 == 5) * (rootcat2 == 4) )  ,
                              
                              d12 * ( (rootcat6 == 1) * (rootcat3 == 2) + (rootcat6 == 2) * (rootcat3 == 1) )    +
                                d13 * ( (rootcat6 == 1) * (rootcat3 == 3) + (rootcat6 == 3) * (rootcat3 == 1) )  +
                                d14 * ( (rootcat6 == 1) * (rootcat3 == 4) + (rootcat6 == 4) * (rootcat3 == 1) )  +
                                d15 * ( (rootcat6 == 1) * (rootcat3 == 5) + (rootcat6 == 5) * (rootcat3 == 1) )  +
                                d23 * ( (rootcat6 == 2) * (rootcat3 == 3) + (rootcat6 == 3) * (rootcat3 == 2) )  +
                                d24 * ( (rootcat6 == 2) * (rootcat3 == 4) + (rootcat6 == 4) * (rootcat3 == 2) )  +
                                d25 * ( (rootcat6 == 2) * (rootcat3 == 5) + (rootcat6 == 5) * (rootcat3 == 2) )  +
                                d34 * ( (rootcat6 == 3) * (rootcat3 == 4) + (rootcat6 == 4) * (rootcat3 == 3) )  +
                                d35 * ( (rootcat6 == 3) * (rootcat3 == 5) + (rootcat6 == 5) * (rootcat3 == 3) )  +
                                d45 * ( (rootcat6 == 4) * (rootcat3 == 5) + (rootcat6 == 5) * (rootcat3 == 4) )  ,
                              
                              d12 * ( (rootcat6 == 1) * (rootcat4 == 2) + (rootcat6 == 2) * (rootcat4 == 1) )    +
                                d13 * ( (rootcat6 == 1) * (rootcat4 == 3) + (rootcat6 == 3) * (rootcat4 == 1) )  +
                                d14 * ( (rootcat6 == 1) * (rootcat4 == 4) + (rootcat6 == 4) * (rootcat4 == 1) )  +
                                d15 * ( (rootcat6 == 1) * (rootcat4 == 5) + (rootcat6 == 5) * (rootcat4 == 1) )  +
                                d23 * ( (rootcat6 == 2) * (rootcat4 == 3) + (rootcat6 == 3) * (rootcat4 == 2) )  +
                                d24 * ( (rootcat6 == 2) * (rootcat4 == 4) + (rootcat6 == 4) * (rootcat4 == 2) )  +
                                d25 * ( (rootcat6 == 2) * (rootcat4 == 5) + (rootcat6 == 5) * (rootcat4 == 2) )  +
                                d34 * ( (rootcat6 == 3) * (rootcat4 == 4) + (rootcat6 == 4) * (rootcat4 == 3) )  +
                                d35 * ( (rootcat6 == 3) * (rootcat4 == 5) + (rootcat6 == 5) * (rootcat4 == 3) )  +
                                d45 * ( (rootcat6 == 4) * (rootcat4 == 5) + (rootcat6 == 5) * (rootcat4 == 4) )  ,
                              
                              d12 * ( (rootcat6 == 1) * (rootcat5 == 2) + (rootcat6 == 2) * (rootcat5 == 1) )    +
                                d13 * ( (rootcat6 == 1) * (rootcat5 == 3) + (rootcat6 == 3) * (rootcat5 == 1) )  +
                                d14 * ( (rootcat6 == 1) * (rootcat5 == 4) + (rootcat6 == 4) * (rootcat5 == 1) )  +
                                d15 * ( (rootcat6 == 1) * (rootcat5 == 5) + (rootcat6 == 5) * (rootcat5 == 1) )  +
                                d23 * ( (rootcat6 == 2) * (rootcat5 == 3) + (rootcat6 == 3) * (rootcat5 == 2) )  +
                                d24 * ( (rootcat6 == 2) * (rootcat5 == 4) + (rootcat6 == 4) * (rootcat5 == 2) )  +
                                d25 * ( (rootcat6 == 2) * (rootcat5 == 5) + (rootcat6 == 5) * (rootcat5 == 2) )  +
                                d34 * ( (rootcat6 == 3) * (rootcat5 == 4) + (rootcat6 == 4) * (rootcat5 == 3) )  +
                                d35 * ( (rootcat6 == 3) * (rootcat5 == 5) + (rootcat6 == 5) * (rootcat5 == 3) )  +
                                d45 * ( (rootcat6 == 4) * (rootcat5 == 5) + (rootcat6 == 5) * (rootcat5 == 4) )  ,
                              
                              0,0,0) , 
                  
                  list(       d12 * ( (rootcat7 == 1) * (rootcat1 == 2) + (rootcat7 == 2) * (rootcat1 == 1) )    +
                                d13 * ( (rootcat7 == 1) * (rootcat1 == 3) + (rootcat7 == 3) * (rootcat1 == 1) )  +
                                d14 * ( (rootcat7 == 1) * (rootcat1 == 4) + (rootcat7 == 4) * (rootcat1 == 1) )  +
                                d15 * ( (rootcat7 == 1) * (rootcat1 == 5) + (rootcat7 == 5) * (rootcat1 == 1) )  +
                                d23 * ( (rootcat7 == 2) * (rootcat1 == 3) + (rootcat7 == 3) * (rootcat1 == 2) )  +
                                d24 * ( (rootcat7 == 2) * (rootcat1 == 4) + (rootcat7 == 4) * (rootcat1 == 2) )  +
                                d25 * ( (rootcat7 == 2) * (rootcat1 == 5) + (rootcat7 == 5) * (rootcat1 == 2) )  +
                                d34 * ( (rootcat7 == 3) * (rootcat1 == 4) + (rootcat7 == 4) * (rootcat1 == 3) )  +
                                d35 * ( (rootcat7 == 3) * (rootcat1 == 5) + (rootcat7 == 5) * (rootcat1 == 3) )  +
                                d45 * ( (rootcat7 == 4) * (rootcat1 == 5) + (rootcat7 == 5) * (rootcat1 == 4) )  ,
                              
                              d12 * ( (rootcat7 == 1) * (rootcat2 == 2) + (rootcat7 == 2) * (rootcat2 == 1) )    +
                                d13 * ( (rootcat7 == 1) * (rootcat2 == 3) + (rootcat7 == 3) * (rootcat2 == 1) )  +
                                d14 * ( (rootcat7 == 1) * (rootcat2 == 4) + (rootcat7 == 4) * (rootcat2 == 1) )  +
                                d15 * ( (rootcat7 == 1) * (rootcat2 == 5) + (rootcat7 == 5) * (rootcat2 == 1) )  +
                                d23 * ( (rootcat7 == 2) * (rootcat2 == 3) + (rootcat7 == 3) * (rootcat2 == 2) )  +
                                d24 * ( (rootcat7 == 2) * (rootcat2 == 4) + (rootcat7 == 4) * (rootcat2 == 2) )  +
                                d25 * ( (rootcat7 == 2) * (rootcat2 == 5) + (rootcat7 == 5) * (rootcat2 == 2) )  +
                                d34 * ( (rootcat7 == 3) * (rootcat2 == 4) + (rootcat7 == 4) * (rootcat2 == 3) )  +
                                d35 * ( (rootcat7 == 3) * (rootcat2 == 5) + (rootcat7 == 5) * (rootcat2 == 3) )  +
                                d45 * ( (rootcat7 == 4) * (rootcat2 == 5) + (rootcat7 == 5) * (rootcat2 == 4) )  ,
                              
                              d12 * ( (rootcat7 == 1) * (rootcat3 == 2) + (rootcat7 == 2) * (rootcat3 == 1) )    +
                                d13 * ( (rootcat7 == 1) * (rootcat3 == 3) + (rootcat7 == 3) * (rootcat3 == 1) )  +
                                d14 * ( (rootcat7 == 1) * (rootcat3 == 4) + (rootcat7 == 4) * (rootcat3 == 1) )  +
                                d15 * ( (rootcat7 == 1) * (rootcat3 == 5) + (rootcat7 == 5) * (rootcat3 == 1) )  +
                                d23 * ( (rootcat7 == 2) * (rootcat3 == 3) + (rootcat7 == 3) * (rootcat3 == 2) )  +
                                d24 * ( (rootcat7 == 2) * (rootcat3 == 4) + (rootcat7 == 4) * (rootcat3 == 2) )  +
                                d25 * ( (rootcat7 == 2) * (rootcat3 == 5) + (rootcat7 == 5) * (rootcat3 == 2) )  +
                                d34 * ( (rootcat7 == 3) * (rootcat3 == 4) + (rootcat7 == 4) * (rootcat3 == 3) )  +
                                d35 * ( (rootcat7 == 3) * (rootcat3 == 5) + (rootcat7 == 5) * (rootcat3 == 3) )  +
                                d45 * ( (rootcat7 == 4) * (rootcat3 == 5) + (rootcat7 == 5) * (rootcat3 == 4) )  ,
                              
                              d12 * ( (rootcat7 == 1) * (rootcat4 == 2) + (rootcat7 == 2) * (rootcat4 == 1) )    +
                                d13 * ( (rootcat7 == 1) * (rootcat4 == 3) + (rootcat7 == 3) * (rootcat4 == 1) )  +
                                d14 * ( (rootcat7 == 1) * (rootcat4 == 4) + (rootcat7 == 4) * (rootcat4 == 1) )  +
                                d15 * ( (rootcat7 == 1) * (rootcat4 == 5) + (rootcat7 == 5) * (rootcat4 == 1) )  +
                                d23 * ( (rootcat7 == 2) * (rootcat4 == 3) + (rootcat7 == 3) * (rootcat4 == 2) )  +
                                d24 * ( (rootcat7 == 2) * (rootcat4 == 4) + (rootcat7 == 4) * (rootcat4 == 2) )  +
                                d25 * ( (rootcat7 == 2) * (rootcat4 == 5) + (rootcat7 == 5) * (rootcat4 == 2) )  +
                                d34 * ( (rootcat7 == 3) * (rootcat4 == 4) + (rootcat7 == 4) * (rootcat4 == 3) )  +
                                d35 * ( (rootcat7 == 3) * (rootcat4 == 5) + (rootcat7 == 5) * (rootcat4 == 3) )  +
                                d45 * ( (rootcat7 == 4) * (rootcat4 == 5) + (rootcat7 == 5) * (rootcat4 == 4) )  ,
                              
                              d12 * ( (rootcat7 == 1) * (rootcat5 == 2) + (rootcat7 == 2) * (rootcat5 == 1) )    +
                                d13 * ( (rootcat7 == 1) * (rootcat5 == 3) + (rootcat7 == 3) * (rootcat5 == 1) )  +
                                d14 * ( (rootcat7 == 1) * (rootcat5 == 4) + (rootcat7 == 4) * (rootcat5 == 1) )  +
                                d15 * ( (rootcat7 == 1) * (rootcat5 == 5) + (rootcat7 == 5) * (rootcat5 == 1) )  +
                                d23 * ( (rootcat7 == 2) * (rootcat5 == 3) + (rootcat7 == 3) * (rootcat5 == 2) )  +
                                d24 * ( (rootcat7 == 2) * (rootcat5 == 4) + (rootcat7 == 4) * (rootcat5 == 2) )  +
                                d25 * ( (rootcat7 == 2) * (rootcat5 == 5) + (rootcat7 == 5) * (rootcat5 == 2) )  +
                                d34 * ( (rootcat7 == 3) * (rootcat5 == 4) + (rootcat7 == 4) * (rootcat5 == 3) )  +
                                d35 * ( (rootcat7 == 3) * (rootcat5 == 5) + (rootcat7 == 5) * (rootcat5 == 3) )  +
                                d45 * ( (rootcat7 == 4) * (rootcat5 == 5) + (rootcat7 == 5) * (rootcat5 == 4) )  ,
                              
                              d12 * ( (rootcat7 == 1) * (rootcat6 == 2) + (rootcat7 == 2) * (rootcat6 == 1) )    +
                                d13 * ( (rootcat7 == 1) * (rootcat6 == 3) + (rootcat7 == 3) * (rootcat6 == 1) )  +
                                d14 * ( (rootcat7 == 1) * (rootcat6 == 4) + (rootcat7 == 4) * (rootcat6 == 1) )  +
                                d15 * ( (rootcat7 == 1) * (rootcat6 == 5) + (rootcat7 == 5) * (rootcat6 == 1) )  +
                                d23 * ( (rootcat7 == 2) * (rootcat6 == 3) + (rootcat7 == 3) * (rootcat6 == 2) )  +
                                d24 * ( (rootcat7 == 2) * (rootcat6 == 4) + (rootcat7 == 4) * (rootcat6 == 2) )  +
                                d25 * ( (rootcat7 == 2) * (rootcat6 == 5) + (rootcat7 == 5) * (rootcat6 == 2) )  +
                                d34 * ( (rootcat7 == 3) * (rootcat6 == 4) + (rootcat7 == 4) * (rootcat6 == 3) )  +
                                d35 * ( (rootcat7 == 3) * (rootcat6 == 5) + (rootcat7 == 5) * (rootcat6 == 3) )  +
                                d45 * ( (rootcat7 == 4) * (rootcat6 == 5) + (rootcat7 == 5) * (rootcat6 == 4) )  ,
                              
                              0,0 ) , 
                  
                  list(       d12 * ( (rootcat8 == 1) * (rootcat1 == 2) + (rootcat8 == 2) * (rootcat1 == 1) )    +
                                d13 * ( (rootcat8 == 1) * (rootcat1 == 3) + (rootcat8 == 3) * (rootcat1 == 1) )  +
                                d14 * ( (rootcat8 == 1) * (rootcat1 == 4) + (rootcat8 == 4) * (rootcat1 == 1) )  +
                                d15 * ( (rootcat8 == 1) * (rootcat1 == 5) + (rootcat8 == 5) * (rootcat1 == 1) )  +
                                d23 * ( (rootcat8 == 2) * (rootcat1 == 3) + (rootcat8 == 3) * (rootcat1 == 2) )  +
                                d24 * ( (rootcat8 == 2) * (rootcat1 == 4) + (rootcat8 == 4) * (rootcat1 == 2) )  +
                                d25 * ( (rootcat8 == 2) * (rootcat1 == 5) + (rootcat8 == 5) * (rootcat1 == 2) )  +
                                d34 * ( (rootcat8 == 3) * (rootcat1 == 4) + (rootcat8 == 4) * (rootcat1 == 3) )  +
                                d35 * ( (rootcat8 == 3) * (rootcat1 == 5) + (rootcat8 == 5) * (rootcat1 == 3) )  +
                                d45 * ( (rootcat8 == 4) * (rootcat1 == 5) + (rootcat8 == 5) * (rootcat1 == 4) )  ,
                              
                              d12 * ( (rootcat8 == 1) * (rootcat2 == 2) + (rootcat8 == 2) * (rootcat2 == 1) )    +
                                d13 * ( (rootcat8 == 1) * (rootcat2 == 3) + (rootcat8 == 3) * (rootcat2 == 1) )  +
                                d14 * ( (rootcat8 == 1) * (rootcat2 == 4) + (rootcat8 == 4) * (rootcat2 == 1) )  +
                                d15 * ( (rootcat8 == 1) * (rootcat2 == 5) + (rootcat8 == 5) * (rootcat2 == 1) )  +
                                d23 * ( (rootcat8 == 2) * (rootcat2 == 3) + (rootcat8 == 3) * (rootcat2 == 2) )  +
                                d24 * ( (rootcat8 == 2) * (rootcat2 == 4) + (rootcat8 == 4) * (rootcat2 == 2) )  +
                                d25 * ( (rootcat8 == 2) * (rootcat2 == 5) + (rootcat8 == 5) * (rootcat2 == 2) )  +
                                d34 * ( (rootcat8 == 3) * (rootcat2 == 4) + (rootcat8 == 4) * (rootcat2 == 3) )  +
                                d35 * ( (rootcat8 == 3) * (rootcat2 == 5) + (rootcat8 == 5) * (rootcat2 == 3) )  +
                                d45 * ( (rootcat8 == 4) * (rootcat2 == 5) + (rootcat8 == 5) * (rootcat2 == 4) )  ,
                              
                              d12 * ( (rootcat8 == 1) * (rootcat3 == 2) + (rootcat8 == 2) * (rootcat3 == 1) )    +
                                d13 * ( (rootcat8 == 1) * (rootcat3 == 3) + (rootcat8 == 3) * (rootcat3 == 1) )  +
                                d14 * ( (rootcat8 == 1) * (rootcat3 == 4) + (rootcat8 == 4) * (rootcat3 == 1) )  +
                                d15 * ( (rootcat8 == 1) * (rootcat3 == 5) + (rootcat8 == 5) * (rootcat3 == 1) )  +
                                d23 * ( (rootcat8 == 2) * (rootcat3 == 3) + (rootcat8 == 3) * (rootcat3 == 2) )  +
                                d24 * ( (rootcat8 == 2) * (rootcat3 == 4) + (rootcat8 == 4) * (rootcat3 == 2) )  +
                                d25 * ( (rootcat8 == 2) * (rootcat3 == 5) + (rootcat8 == 5) * (rootcat3 == 2) )  +
                                d34 * ( (rootcat8 == 3) * (rootcat3 == 4) + (rootcat8 == 4) * (rootcat3 == 3) )  +
                                d35 * ( (rootcat8 == 3) * (rootcat3 == 5) + (rootcat8 == 5) * (rootcat3 == 3) )  +
                                d45 * ( (rootcat8 == 4) * (rootcat3 == 5) + (rootcat8 == 5) * (rootcat3 == 4) )  ,
                              
                              d12 * ( (rootcat8 == 1) * (rootcat4 == 2) + (rootcat8 == 2) * (rootcat4 == 1) )    +
                                d13 * ( (rootcat8 == 1) * (rootcat4 == 3) + (rootcat8 == 3) * (rootcat4 == 1) )  +
                                d14 * ( (rootcat8 == 1) * (rootcat4 == 4) + (rootcat8 == 4) * (rootcat4 == 1) )  +
                                d15 * ( (rootcat8 == 1) * (rootcat4 == 5) + (rootcat8 == 5) * (rootcat4 == 1) )  +
                                d23 * ( (rootcat8 == 2) * (rootcat4 == 3) + (rootcat8 == 3) * (rootcat4 == 2) )  +
                                d24 * ( (rootcat8 == 2) * (rootcat4 == 4) + (rootcat8 == 4) * (rootcat4 == 2) )  +
                                d25 * ( (rootcat8 == 2) * (rootcat4 == 5) + (rootcat8 == 5) * (rootcat4 == 2) )  +
                                d34 * ( (rootcat8 == 3) * (rootcat4 == 4) + (rootcat8 == 4) * (rootcat4 == 3) )  +
                                d35 * ( (rootcat8 == 3) * (rootcat4 == 5) + (rootcat8 == 5) * (rootcat4 == 3) )  +
                                d45 * ( (rootcat8 == 4) * (rootcat4 == 5) + (rootcat8 == 5) * (rootcat4 == 4) )  ,
                              
                              d12 * ( (rootcat8 == 1) * (rootcat5 == 2) + (rootcat8 == 2) * (rootcat5 == 1) )    +
                                d13 * ( (rootcat8 == 1) * (rootcat5 == 3) + (rootcat8 == 3) * (rootcat5 == 1) )  +
                                d14 * ( (rootcat8 == 1) * (rootcat5 == 4) + (rootcat8 == 4) * (rootcat5 == 1) )  +
                                d15 * ( (rootcat8 == 1) * (rootcat5 == 5) + (rootcat8 == 5) * (rootcat5 == 1) )  +
                                d23 * ( (rootcat8 == 2) * (rootcat5 == 3) + (rootcat8 == 3) * (rootcat5 == 2) )  +
                                d24 * ( (rootcat8 == 2) * (rootcat5 == 4) + (rootcat8 == 4) * (rootcat5 == 2) )  +
                                d25 * ( (rootcat8 == 2) * (rootcat5 == 5) + (rootcat8 == 5) * (rootcat5 == 2) )  +
                                d34 * ( (rootcat8 == 3) * (rootcat5 == 4) + (rootcat8 == 4) * (rootcat5 == 3) )  +
                                d35 * ( (rootcat8 == 3) * (rootcat5 == 5) + (rootcat8 == 5) * (rootcat5 == 3) )  +
                                d45 * ( (rootcat8 == 4) * (rootcat5 == 5) + (rootcat8 == 5) * (rootcat5 == 4) )  ,
                              
                              d12 * ( (rootcat8 == 1) * (rootcat6 == 2) + (rootcat8 == 2) * (rootcat6 == 1) )    +
                                d13 * ( (rootcat8 == 1) * (rootcat6 == 3) + (rootcat8 == 3) * (rootcat6 == 1) )  +
                                d14 * ( (rootcat8 == 1) * (rootcat6 == 4) + (rootcat8 == 4) * (rootcat6 == 1) )  +
                                d15 * ( (rootcat8 == 1) * (rootcat6 == 5) + (rootcat8 == 5) * (rootcat6 == 1) )  +
                                d23 * ( (rootcat8 == 2) * (rootcat6 == 3) + (rootcat8 == 3) * (rootcat6 == 2) )  +
                                d24 * ( (rootcat8 == 2) * (rootcat6 == 4) + (rootcat8 == 4) * (rootcat6 == 2) )  +
                                d25 * ( (rootcat8 == 2) * (rootcat6 == 5) + (rootcat8 == 5) * (rootcat6 == 2) )  +
                                d34 * ( (rootcat8 == 3) * (rootcat6 == 4) + (rootcat8 == 4) * (rootcat6 == 3) )  +
                                d35 * ( (rootcat8 == 3) * (rootcat6 == 5) + (rootcat8 == 5) * (rootcat6 == 3) )  +
                                d45 * ( (rootcat8 == 4) * (rootcat6 == 5) + (rootcat8 == 5) * (rootcat6 == 4) )  ,
                              
                              d12 * ( (rootcat8 == 1) * (rootcat7 == 2) + (rootcat8 == 2) * (rootcat7 == 1) )    +
                                d13 * ( (rootcat8 == 1) * (rootcat7 == 3) + (rootcat8 == 3) * (rootcat7 == 1) )  +
                                d14 * ( (rootcat8 == 1) * (rootcat7 == 4) + (rootcat8 == 4) * (rootcat7 == 1) )  +
                                d15 * ( (rootcat8 == 1) * (rootcat7 == 5) + (rootcat8 == 5) * (rootcat7 == 1) )  +
                                d23 * ( (rootcat8 == 2) * (rootcat7 == 3) + (rootcat8 == 3) * (rootcat7 == 2) )  +
                                d24 * ( (rootcat8 == 2) * (rootcat7 == 4) + (rootcat8 == 4) * (rootcat7 == 2) )  +
                                d25 * ( (rootcat8 == 2) * (rootcat7 == 5) + (rootcat8 == 5) * (rootcat7 == 2) )  +
                                d34 * ( (rootcat8 == 3) * (rootcat7 == 4) + (rootcat8 == 4) * (rootcat7 == 3) )  +
                                d35 * ( (rootcat8 == 3) * (rootcat7 == 5) + (rootcat8 == 5) * (rootcat7 == 3) )  +
                                d45 * ( (rootcat8 == 4) * (rootcat7 == 5) + (rootcat8 == 5) * (rootcat7 == 4) )  ,
                              0 ) 
    )
    
    
    
    emdc_settings$utilityOutside = delta_outside
    emdc_settings$utilities = V
    emdc_settings$gamma = gamma
    emdc_settings$delta = delta

    # 
    P[["model"]] = apollo_emdc(emdc_settings, functionality)
    
    ### Take product across observation for same individual
    P = apollo_panelProd(P, apollo_inputs ,functionality)
    
    P = apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
  }
  
  apollo_prediction = function (model, apollo_probabilities = apollo_probabilities1, apollo_inputs, prediction_settings = list(), modelComponent = NA) {
    prediction_settings = list()
    modelComponent = NA
    
    if (is.null(prediction_settings$modelComponent)) {
      if (exists("modelComponent")) 
        prediction_settings$modelComponent = modelComponent
      else prediction_settings$modelComponent = NA
    }
    if (is.null(prediction_settings$runs))   prediction_settings$runs = 1
    if (is.null(prediction_settings$silent))   prediction_settings$silent = FALSE
    silent = prediction_settings$silent
    if (!is.null(apollo_inputs$silent) && apollo_inputs$silent) silent = TRUE
    
    if (is.null(prediction_settings$nRep)) prediction_settings$nRep <- 100L
    
    if (is.null(prediction_settings$summary)) prediction_settings$summary <- TRUE
    
    modelComponent = prediction_settings$modelComponent
    runs = prediction_settings$runs
    apollo_inputs$nRep <- prediction_settings$nRep
    
    apollo_compareInputs(apollo_inputs)
    apollo_randCoeff = apollo_inputs[["apollo_randCoeff"]]
    apollo_lcPars = apollo_inputs[["apollo_lcPars"]]
    apollo_checkArguments(apollo_probabilities, apollo_randCoeff, 
                          apollo_lcPars)
    if (!silent)     apollo_print("Running predictions from model using parameter estimates...")
    
    
    predictions = apollo_probabilities(apollo_beta, apollo_inputs, 
                                       functionality = "prediction")
    
    predictionamount = predictions$model[,1:9]  
    
    return(predictionamount)
  }
  
  apollo_inputs = apollo_validateInputs(apollo_beta = apollo_beta1)
  
  predictionprob1 = apollo_prediction(model, apollo_probabilities = apollo_probabilities1, apollo_inputs, prediction_settings=list(runs=1), modelComponent = NA)
  
  
  
  
  #######################Segment 2
  
  apollo_beta = apollo_beta2
  
  apollo_probabilities2=function(apollo_beta = apollo_beta2, apollo_inputs, functionality="estimate"){
    
    ### Attach inputs and detach after function exit
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    
    ### Create list of probabilities P
    P = list()
    
    ### Define individual alternatives
    alternatives = c("v1", 
                     "v2", 
                     "v3", 
                     "v4", 
                     "v5", 
                     "v6", 
                     "v7", 
                     "v8" #, 
                     # "outside"
    )
    
    ### Define availabilities
    avail = list(v1  = availnew1,    ######Need to check whether it works here
                 v2     = availnew2,
                 v3   = availnew3,
                 v4 = availnew4,
                 v5 = availnew5,
                 v6   = availnew6,
                 v7 = availnew7,
                 v8 = availnew8 #,
                 #outside  = availnew9
    )
    
    ### Define continuous consumption for individual alternatives
    continuousChoice = list(v1  =clicknew1 * availnew1,
                            v2     =clicknew2 * availnew2,
                            v3   =clicknew3 * availnew3,
                            v4 =clicknew4 * availnew4,
                            v5 =clicknew5 * availnew5,
                            v6   =clicknew6 * availnew6,
                            v7 =clicknew7 * availnew7,
                            v8 =clicknew8 * availnew8#,
                            #outside  =outsideclick
    )
    
    
    
    ### Define alpha parameters
    alpha = list(v1  = 1e-3 , 
                 v2     = 1e-3 , 
                 v3   = 1e-3 , 
                 v4 = 1e-3 , 
                 v5 = 1e-3 ,
                 v6   = 1e-3 ,
                 v7  = 1e-3 , 
                 v8 = 1e-3 #, 
                 #outside  = 1e-3 
    )
    
    
    ### Define costs for individual alternatives
    cost = list(v1      = 1, 
                v2         = 1,
                v3       = 1,
                v4     = 1,
                v5 = 1,
                v6       = 1, 
                v7      = 1,
                v8     = 1 #,
                #outside       = 1
    )
    
    ### Define budget
    budget = clicktotal
    
    emdc_settings <- list(continuousChoice = continuousChoice, 
                          avail            = avail,
                          budget           = budget,
                          sigma            = 0.99, 
                          cost             = cost)
    
    ### ### Compute class-specific utilities
    V = list()
    
    V[["v1"    ]] =   delta_v1*(rootcat1==1) +  delta_v2*(rootcat1==2) +  delta_v3*(rootcat1==3) +  
      delta_v4*(rootcat1==4) +  delta_v5*(rootcat1==5) + 
      bweekend*weekend + bnpv*npvnew1 + bnclick*nclicknew1 + bnbuy*nbuynew1 + 
      nrccpv*rccpvnew1 + nrccclick*rccclicknew1 + nrccbuy*rccbuynew1 + badstockpsai*ad11 + bresidual*res1 + bvariety*Lag_Distinct
    
    V[["v2"    ]] = delta_v1*(rootcat2==1) +  delta_v2*(rootcat2==2) +  delta_v3*(rootcat2==3) +  
      delta_v4*(rootcat2==4) +  delta_v5*(rootcat2==5) + 
      bweekend*weekend + bnpv*npvnew2 + bnclick*nclicknew2 + bnbuy*nbuynew2 + 
      nrccpv*rccpvnew2 + nrccclick*rccclicknew2 + nrccbuy*rccbuynew2 + badstockpsai*ad21 + bresidual*res2 + bvariety*Lag_Distinct
    
    V[["v3"  ]] = delta_v1*(rootcat3==1) +  delta_v2*(rootcat3==2) +  delta_v3*(rootcat3==3) +  
      delta_v4*(rootcat3==4) +  delta_v5*(rootcat3==5) + 
      bweekend*weekend + bnpv*npvnew3 + bnclick*nclicknew3 + bnbuy*nbuynew3 + 
      nrccpv*rccpvnew3 + nrccclick*rccclicknew3 + nrccbuy*rccbuynew3 + badstockpsai*ad31 + bresidual*res3 + bvariety*Lag_Distinct
    
    V[["v4"]] = delta_v1*(rootcat4==1) +  delta_v2*(rootcat4==2) +  delta_v3*(rootcat4==3) +  
      delta_v4*(rootcat4==4) +  delta_v5*(rootcat4==5) + 
      bweekend*weekend + bnpv*npvnew4 + bnclick*nclicknew4 + bnbuy*nbuynew4 + 
      nrccpv*rccpvnew4 + nrccclick*rccclicknew4 + nrccbuy*rccbuynew4 + badstockpsai*ad41 + bresidual*res4 + bvariety*Lag_Distinct
    
    V[["v5"]] = delta_v1*(rootcat5==1) +  delta_v2*(rootcat5==2) +  delta_v3*(rootcat5==3) +  
      delta_v4*(rootcat5==4) +  delta_v5*(rootcat5==5) + 
      bweekend*weekend + bnpv*npvnew5 + bnclick*nclicknew5 + bnbuy*nbuynew5 + 
      nrccpv*rccpvnew5 + nrccclick*rccclicknew5 + nrccbuy*rccbuynew5 + badstockpsai*ad51 + bresidual*res5 + bvariety*Lag_Distinct
    
    V[["v6"  ]] = delta_v1*(rootcat6==1) +  delta_v2*(rootcat6==2) +  delta_v3*(rootcat6==3) +  
      delta_v4*(rootcat6==4) +  delta_v5*(rootcat6==5) + 
      bweekend*weekend + bnpv*npvnew6 + bnclick*nclicknew6 + bnbuy*nbuynew6 + 
      nrccpv*rccpvnew6 + nrccclick*rccclicknew6 + nrccbuy*rccbuynew6 + badstockpsai*ad61 + bresidual*res6 + bvariety*Lag_Distinct
    
    V[["v7"]] = delta_v1*(rootcat7==1) +  delta_v2*(rootcat7==2) +  delta_v3*(rootcat7==3) +  
      delta_v4*(rootcat7==4) +  delta_v5*(rootcat7==5) + 
      bweekend*weekend + bnpv*npvnew7 + bnclick*nclicknew7 + bnbuy*nbuynew7 + 
      nrccpv*rccpvnew7 + nrccclick*rccclicknew7 + nrccbuy*rccbuynew7 + badstockpsai*ad71 + bresidual*res7 + bvariety*Lag_Distinct
    
    V[["v8"]] = delta_v1*(rootcat8==1) +  delta_v2*(rootcat8==2) +  delta_v3*(rootcat8==3) +  
      delta_v4*(rootcat8==4) +  delta_v5*(rootcat8==5) + 
      bweekend*weekend + bnpv*npvnew8 + bnclick*nclicknew8 + bnbuy*nbuynew8 + 
      nrccpv*rccpvnew8 + nrccclick*rccclicknew8 + nrccbuy*rccbuynew8 + badstockpsai*ad81 + bresidual*res8 + bvariety*Lag_Distinct
    
    ### Define gamma parameters
    gamma = list(v1      = exp( gamma_v1*(rootcat1==1) +  gamma_v2*(rootcat1==2) +  gamma_v3*(rootcat1==3) +  
                                  gamma_v4*(rootcat1==4) +  gamma_v5*(rootcat1==5) ), 
                 
                 v2      = exp( gamma_v1*(rootcat2==1) +  gamma_v2*(rootcat2==2) +  gamma_v3*(rootcat2==3) +  
                                  gamma_v4*(rootcat2==4) +  gamma_v5*(rootcat2==5) ), 
                 
                 
                 v3      = exp( gamma_v1*(rootcat3==1) +  gamma_v2*(rootcat3==2) +  gamma_v3*(rootcat3==3) +  
                                  gamma_v4*(rootcat3==4) +  gamma_v5*(rootcat3==5) ), 
                 
                 v4      = exp( gamma_v1*(rootcat4==1) +  gamma_v2*(rootcat4==2) +  gamma_v3*(rootcat4==3) +  
                                  gamma_v4*(rootcat4==4) +  gamma_v5*(rootcat4==5) ), 
                 
                 v5      = exp( gamma_v1*(rootcat5==1) +  gamma_v2*(rootcat5==2) +  gamma_v3*(rootcat5==3) +  
                                  gamma_v4*(rootcat5==4) +  gamma_v5*(rootcat5==5) ),
                 
                 v6      = exp( gamma_v1*(rootcat6==1) +  gamma_v2*(rootcat6==2) +  gamma_v3*(rootcat6==3) +  
                                  gamma_v4*(rootcat6==4) +  gamma_v5*(rootcat6==5) ), 
                 
                 v7      = exp( gamma_v1*(rootcat7==1) +  gamma_v2*(rootcat7==2) +  gamma_v3*(rootcat7==3) +  
                                  gamma_v4*(rootcat7==4) +  gamma_v5*(rootcat7==5) ), 
                 
                 v8      = exp( gamma_v1*(rootcat8==1) +  gamma_v2*(rootcat8==2) +  gamma_v3*(rootcat8==3) +  
                                  gamma_v4*(rootcat8==4) +  gamma_v5*(rootcat8==5) )
    )
    
    delta = list( list(0,0,0,0,0,0,0 ,0 ) , 
                  list(   d12 * ( (rootcat2 == 1) * (rootcat1 == 2) + (rootcat2 == 2) * (rootcat1 == 1) )    +
                            d13 * ( (rootcat2 == 1) * (rootcat1 == 3) + (rootcat2 == 3) * (rootcat1 == 1) )  +
                            d14 * ( (rootcat2 == 1) * (rootcat1 == 4) + (rootcat2 == 4) * (rootcat1 == 1) )  +
                            d15 * ( (rootcat2 == 1) * (rootcat1 == 5) + (rootcat2 == 5) * (rootcat1 == 1) )  +
                            d23 * ( (rootcat2 == 2) * (rootcat1 == 3) + (rootcat2 == 3) * (rootcat1 == 2) )  +
                            d24 * ( (rootcat2 == 2) * (rootcat1 == 4) + (rootcat2 == 4) * (rootcat1 == 2) )  +
                            d25 * ( (rootcat2 == 2) * (rootcat1 == 5) + (rootcat2 == 5) * (rootcat1 == 2) )  +
                            d34 * ( (rootcat2 == 3) * (rootcat1 == 4) + (rootcat2 == 4) * (rootcat1 == 3) )  +
                            d35 * ( (rootcat2 == 3) * (rootcat1 == 5) + (rootcat2 == 5) * (rootcat1 == 3) )  +
                            d45 * ( (rootcat2 == 4) * (rootcat1 == 5) + (rootcat2 == 5) * (rootcat1 == 4) )  ,
                          0,0,0,0,0,0,0 ) ,
                  
                  list(   d12 * ( (rootcat3 == 1) * (rootcat1 == 2) + (rootcat3 == 2) * (rootcat1 == 1) )    +
                            d13 * ( (rootcat3 == 1) * (rootcat1 == 3) + (rootcat3 == 3) * (rootcat1 == 1) )  +
                            d14 * ( (rootcat3 == 1) * (rootcat1 == 4) + (rootcat3 == 4) * (rootcat1 == 1) )  +
                            d15 * ( (rootcat3 == 1) * (rootcat1 == 5) + (rootcat3 == 5) * (rootcat1 == 1) )  +
                            d23 * ( (rootcat3 == 2) * (rootcat1 == 3) + (rootcat3 == 3) * (rootcat1 == 2) )  +
                            d24 * ( (rootcat3 == 2) * (rootcat1 == 4) + (rootcat3 == 4) * (rootcat1 == 2) )  +
                            d25 * ( (rootcat3 == 2) * (rootcat1 == 5) + (rootcat3 == 5) * (rootcat1 == 2) )  +
                            d34 * ( (rootcat3 == 3) * (rootcat1 == 4) + (rootcat3 == 4) * (rootcat1 == 3) )  +
                            d35 * ( (rootcat3 == 3) * (rootcat1 == 5) + (rootcat3 == 5) * (rootcat1 == 3) )  +
                            d45 * ( (rootcat3 == 4) * (rootcat1 == 5) + (rootcat3 == 5) * (rootcat1 == 4) )  ,
                          
                          d12 * ( (rootcat3 == 1) * (rootcat2 == 2) + (rootcat3 == 2) * (rootcat2 == 1) )    +
                            d13 * ( (rootcat3 == 1) * (rootcat2 == 3) + (rootcat3 == 3) * (rootcat2 == 1) )  +
                            d14 * ( (rootcat3 == 1) * (rootcat2 == 4) + (rootcat3 == 4) * (rootcat2 == 1) )  +
                            d15 * ( (rootcat3 == 1) * (rootcat2 == 5) + (rootcat3 == 5) * (rootcat2 == 1) )  +
                            d23 * ( (rootcat3 == 2) * (rootcat2 == 3) + (rootcat3 == 3) * (rootcat2 == 2) )  +
                            d24 * ( (rootcat3 == 2) * (rootcat2 == 4) + (rootcat3 == 4) * (rootcat2 == 2) )  +
                            d25 * ( (rootcat3 == 2) * (rootcat2 == 5) + (rootcat3 == 5) * (rootcat2 == 2) )  +
                            d34 * ( (rootcat3 == 3) * (rootcat2 == 4) + (rootcat3 == 4) * (rootcat2 == 3) )  +
                            d35 * ( (rootcat3 == 3) * (rootcat2 == 5) + (rootcat3 == 5) * (rootcat2 == 3) )  +
                            d45 * ( (rootcat3 == 4) * (rootcat2 == 5) + (rootcat3 == 5) * (rootcat2 == 4) )  ,
                          0,0,0,0,0,0 ) , 
                  
                  list(     d12 * ( (rootcat4 == 1) * (rootcat1 == 2) + (rootcat4 == 2) * (rootcat1 == 1) )    +
                              d13 * ( (rootcat4 == 1) * (rootcat1 == 3) + (rootcat4 == 3) * (rootcat1 == 1) )  +
                              d14 * ( (rootcat4 == 1) * (rootcat1 == 4) + (rootcat4 == 4) * (rootcat1 == 1) )  +
                              d15 * ( (rootcat4 == 1) * (rootcat1 == 5) + (rootcat4 == 5) * (rootcat1 == 1) )  +
                              d23 * ( (rootcat4 == 2) * (rootcat1 == 3) + (rootcat4 == 3) * (rootcat1 == 2) )  +
                              d24 * ( (rootcat4 == 2) * (rootcat1 == 4) + (rootcat4 == 4) * (rootcat1 == 2) )  +
                              d25 * ( (rootcat4 == 2) * (rootcat1 == 5) + (rootcat4 == 5) * (rootcat1 == 2) )  +
                              d34 * ( (rootcat4 == 3) * (rootcat1 == 4) + (rootcat4 == 4) * (rootcat1 == 3) )  +
                              d35 * ( (rootcat4 == 3) * (rootcat1 == 5) + (rootcat4 == 5) * (rootcat1 == 3) )  +
                              d45 * ( (rootcat4 == 4) * (rootcat1 == 5) + (rootcat4 == 5) * (rootcat1 == 4) )  ,
                            
                            d12 * ( (rootcat4 == 1) * (rootcat2 == 2) + (rootcat4 == 2) * (rootcat2 == 1) )    +
                              d13 * ( (rootcat4 == 1) * (rootcat2 == 3) + (rootcat4 == 3) * (rootcat2 == 1) )  +
                              d14 * ( (rootcat4 == 1) * (rootcat2 == 4) + (rootcat4 == 4) * (rootcat2 == 1) )  +
                              d15 * ( (rootcat4 == 1) * (rootcat2 == 5) + (rootcat4 == 5) * (rootcat2 == 1) )  +
                              d23 * ( (rootcat4 == 2) * (rootcat2 == 3) + (rootcat4 == 3) * (rootcat2 == 2) )  +
                              d24 * ( (rootcat4 == 2) * (rootcat2 == 4) + (rootcat4 == 4) * (rootcat2 == 2) )  +
                              d25 * ( (rootcat4 == 2) * (rootcat2 == 5) + (rootcat4 == 5) * (rootcat2 == 2) )  +
                              d34 * ( (rootcat4 == 3) * (rootcat2 == 4) + (rootcat4 == 4) * (rootcat2 == 3) )  +
                              d35 * ( (rootcat4 == 3) * (rootcat2 == 5) + (rootcat4 == 5) * (rootcat2 == 3) )  +
                              d45 * ( (rootcat4 == 4) * (rootcat2 == 5) + (rootcat4 == 5) * (rootcat2 == 4) )  ,
                            
                            d12 * ( (rootcat4 == 1) * (rootcat3 == 2) + (rootcat4 == 2) * (rootcat3 == 1) )    +
                              d13 * ( (rootcat4 == 1) * (rootcat3 == 3) + (rootcat4 == 3) * (rootcat3 == 1) )  +
                              d14 * ( (rootcat4 == 1) * (rootcat3 == 4) + (rootcat4 == 4) * (rootcat3 == 1) )  +
                              d15 * ( (rootcat4 == 1) * (rootcat3 == 5) + (rootcat4 == 5) * (rootcat3 == 1) )  +
                              d23 * ( (rootcat4 == 2) * (rootcat3 == 3) + (rootcat4 == 3) * (rootcat3 == 2) )  +
                              d24 * ( (rootcat4 == 2) * (rootcat3 == 4) + (rootcat4 == 4) * (rootcat3 == 2) )  +
                              d25 * ( (rootcat4 == 2) * (rootcat3 == 5) + (rootcat4 == 5) * (rootcat3 == 2) )  +
                              d34 * ( (rootcat4 == 3) * (rootcat3 == 4) + (rootcat4 == 4) * (rootcat3 == 3) )  +
                              d35 * ( (rootcat4 == 3) * (rootcat3 == 5) + (rootcat4 == 5) * (rootcat3 == 3) )  +
                              d45 * ( (rootcat4 == 4) * (rootcat3 == 5) + (rootcat4 == 5) * (rootcat3 == 4) )  ,
                            0,0,0,0,0 ) ,
                  
                  list(       d12 * ( (rootcat5 == 1) * (rootcat1 == 2) + (rootcat5 == 2) * (rootcat1 == 1) )    +
                                d13 * ( (rootcat5 == 1) * (rootcat1 == 3) + (rootcat5 == 3) * (rootcat1 == 1) )  +
                                d14 * ( (rootcat5 == 1) * (rootcat1 == 4) + (rootcat5 == 4) * (rootcat1 == 1) )  +
                                d15 * ( (rootcat5 == 1) * (rootcat1 == 5) + (rootcat5 == 5) * (rootcat1 == 1) )  +
                                d23 * ( (rootcat5 == 2) * (rootcat1 == 3) + (rootcat5 == 3) * (rootcat1 == 2) )  +
                                d24 * ( (rootcat5 == 2) * (rootcat1 == 4) + (rootcat5 == 4) * (rootcat1 == 2) )  +
                                d25 * ( (rootcat5 == 2) * (rootcat1 == 5) + (rootcat5 == 5) * (rootcat1 == 2) )  +
                                d34 * ( (rootcat5 == 3) * (rootcat1 == 4) + (rootcat5 == 4) * (rootcat1 == 3) )  +
                                d35 * ( (rootcat5 == 3) * (rootcat1 == 5) + (rootcat5 == 5) * (rootcat1 == 3) )  +
                                d45 * ( (rootcat5 == 4) * (rootcat1 == 5) + (rootcat5 == 5) * (rootcat1 == 4) )  ,
                              
                              d12 * ( (rootcat5 == 1) * (rootcat2 == 2) + (rootcat5 == 2) * (rootcat2 == 1) )    +
                                d13 * ( (rootcat5 == 1) * (rootcat2 == 3) + (rootcat5 == 3) * (rootcat2 == 1) )  +
                                d14 * ( (rootcat5 == 1) * (rootcat2 == 4) + (rootcat5 == 4) * (rootcat2 == 1) )  +
                                d15 * ( (rootcat5 == 1) * (rootcat2 == 5) + (rootcat5 == 5) * (rootcat2 == 1) )  +
                                d23 * ( (rootcat5 == 2) * (rootcat2 == 3) + (rootcat5 == 3) * (rootcat2 == 2) )  +
                                d24 * ( (rootcat5 == 2) * (rootcat2 == 4) + (rootcat5 == 4) * (rootcat2 == 2) )  +
                                d25 * ( (rootcat5 == 2) * (rootcat2 == 5) + (rootcat5 == 5) * (rootcat2 == 2) )  +
                                d34 * ( (rootcat5 == 3) * (rootcat2 == 4) + (rootcat5 == 4) * (rootcat2 == 3) )  +
                                d35 * ( (rootcat5 == 3) * (rootcat2 == 5) + (rootcat5 == 5) * (rootcat2 == 3) )  +
                                d45 * ( (rootcat5 == 4) * (rootcat2 == 5) + (rootcat5 == 5) * (rootcat2 == 4) )  ,
                              
                              d12 * ( (rootcat5 == 1) * (rootcat3 == 2) + (rootcat5 == 2) * (rootcat3 == 1) )    +
                                d13 * ( (rootcat5 == 1) * (rootcat3 == 3) + (rootcat5 == 3) * (rootcat3 == 1) )  +
                                d14 * ( (rootcat5 == 1) * (rootcat3 == 4) + (rootcat5 == 4) * (rootcat3 == 1) )  +
                                d15 * ( (rootcat5 == 1) * (rootcat3 == 5) + (rootcat5 == 5) * (rootcat3 == 1) )  +
                                d23 * ( (rootcat5 == 2) * (rootcat3 == 3) + (rootcat5 == 3) * (rootcat3 == 2) )  +
                                d24 * ( (rootcat5 == 2) * (rootcat3 == 4) + (rootcat5 == 4) * (rootcat3 == 2) )  +
                                d25 * ( (rootcat5 == 2) * (rootcat3 == 5) + (rootcat5 == 5) * (rootcat3 == 2) )  +
                                d34 * ( (rootcat5 == 3) * (rootcat3 == 4) + (rootcat5 == 4) * (rootcat3 == 3) )  +
                                d35 * ( (rootcat5 == 3) * (rootcat3 == 5) + (rootcat5 == 5) * (rootcat3 == 3) )  +
                                d45 * ( (rootcat5 == 4) * (rootcat3 == 5) + (rootcat5 == 5) * (rootcat3 == 4) )  ,
                              
                              d12 * ( (rootcat5 == 1) * (rootcat4 == 2) + (rootcat5 == 2) * (rootcat4 == 1) )    +
                                d13 * ( (rootcat5 == 1) * (rootcat4 == 3) + (rootcat5 == 3) * (rootcat4 == 1) )  +
                                d14 * ( (rootcat5 == 1) * (rootcat4 == 4) + (rootcat5 == 4) * (rootcat4 == 1) )  +
                                d15 * ( (rootcat5 == 1) * (rootcat4 == 5) + (rootcat5 == 5) * (rootcat4 == 1) )  +
                                d23 * ( (rootcat5 == 2) * (rootcat4 == 3) + (rootcat5 == 3) * (rootcat4 == 2) )  +
                                d24 * ( (rootcat5 == 2) * (rootcat4 == 4) + (rootcat5 == 4) * (rootcat4 == 2) )  +
                                d25 * ( (rootcat5 == 2) * (rootcat4 == 5) + (rootcat5 == 5) * (rootcat4 == 2) )  +
                                d34 * ( (rootcat5 == 3) * (rootcat4 == 4) + (rootcat5 == 4) * (rootcat4 == 3) )  +
                                d35 * ( (rootcat5 == 3) * (rootcat4 == 5) + (rootcat5 == 5) * (rootcat4 == 3) )  +
                                d45 * ( (rootcat5 == 4) * (rootcat4 == 5) + (rootcat5 == 5) * (rootcat4 == 4) )  ,
                              
                              0,0,0,0) ,
                  
                  list(       d12 * ( (rootcat6 == 1) * (rootcat1 == 2) + (rootcat6 == 2) * (rootcat1 == 1) )    +
                                d13 * ( (rootcat6 == 1) * (rootcat1 == 3) + (rootcat6 == 3) * (rootcat1 == 1) )  +
                                d14 * ( (rootcat6 == 1) * (rootcat1 == 4) + (rootcat6 == 4) * (rootcat1 == 1) )  +
                                d15 * ( (rootcat6 == 1) * (rootcat1 == 5) + (rootcat6 == 5) * (rootcat1 == 1) )  +
                                d23 * ( (rootcat6 == 2) * (rootcat1 == 3) + (rootcat6 == 3) * (rootcat1 == 2) )  +
                                d24 * ( (rootcat6 == 2) * (rootcat1 == 4) + (rootcat6 == 4) * (rootcat1 == 2) )  +
                                d25 * ( (rootcat6 == 2) * (rootcat1 == 5) + (rootcat6 == 5) * (rootcat1 == 2) )  +
                                d34 * ( (rootcat6 == 3) * (rootcat1 == 4) + (rootcat6 == 4) * (rootcat1 == 3) )  +
                                d35 * ( (rootcat6 == 3) * (rootcat1 == 5) + (rootcat6 == 5) * (rootcat1 == 3) )  +
                                d45 * ( (rootcat6 == 4) * (rootcat1 == 5) + (rootcat6 == 5) * (rootcat1 == 4) )  ,
                              
                              d12 * ( (rootcat6 == 1) * (rootcat2 == 2) + (rootcat6 == 2) * (rootcat2 == 1) )    +
                                d13 * ( (rootcat6 == 1) * (rootcat2 == 3) + (rootcat6 == 3) * (rootcat2 == 1) )  +
                                d14 * ( (rootcat6 == 1) * (rootcat2 == 4) + (rootcat6 == 4) * (rootcat2 == 1) )  +
                                d15 * ( (rootcat6 == 1) * (rootcat2 == 5) + (rootcat6 == 5) * (rootcat2 == 1) )  +
                                d23 * ( (rootcat6 == 2) * (rootcat2 == 3) + (rootcat6 == 3) * (rootcat2 == 2) )  +
                                d24 * ( (rootcat6 == 2) * (rootcat2 == 4) + (rootcat6 == 4) * (rootcat2 == 2) )  +
                                d25 * ( (rootcat6 == 2) * (rootcat2 == 5) + (rootcat6 == 5) * (rootcat2 == 2) )  +
                                d34 * ( (rootcat6 == 3) * (rootcat2 == 4) + (rootcat6 == 4) * (rootcat2 == 3) )  +
                                d35 * ( (rootcat6 == 3) * (rootcat2 == 5) + (rootcat6 == 5) * (rootcat2 == 3) )  +
                                d45 * ( (rootcat6 == 4) * (rootcat2 == 5) + (rootcat6 == 5) * (rootcat2 == 4) )  ,
                              
                              d12 * ( (rootcat6 == 1) * (rootcat3 == 2) + (rootcat6 == 2) * (rootcat3 == 1) )    +
                                d13 * ( (rootcat6 == 1) * (rootcat3 == 3) + (rootcat6 == 3) * (rootcat3 == 1) )  +
                                d14 * ( (rootcat6 == 1) * (rootcat3 == 4) + (rootcat6 == 4) * (rootcat3 == 1) )  +
                                d15 * ( (rootcat6 == 1) * (rootcat3 == 5) + (rootcat6 == 5) * (rootcat3 == 1) )  +
                                d23 * ( (rootcat6 == 2) * (rootcat3 == 3) + (rootcat6 == 3) * (rootcat3 == 2) )  +
                                d24 * ( (rootcat6 == 2) * (rootcat3 == 4) + (rootcat6 == 4) * (rootcat3 == 2) )  +
                                d25 * ( (rootcat6 == 2) * (rootcat3 == 5) + (rootcat6 == 5) * (rootcat3 == 2) )  +
                                d34 * ( (rootcat6 == 3) * (rootcat3 == 4) + (rootcat6 == 4) * (rootcat3 == 3) )  +
                                d35 * ( (rootcat6 == 3) * (rootcat3 == 5) + (rootcat6 == 5) * (rootcat3 == 3) )  +
                                d45 * ( (rootcat6 == 4) * (rootcat3 == 5) + (rootcat6 == 5) * (rootcat3 == 4) )  ,
                              
                              d12 * ( (rootcat6 == 1) * (rootcat4 == 2) + (rootcat6 == 2) * (rootcat4 == 1) )    +
                                d13 * ( (rootcat6 == 1) * (rootcat4 == 3) + (rootcat6 == 3) * (rootcat4 == 1) )  +
                                d14 * ( (rootcat6 == 1) * (rootcat4 == 4) + (rootcat6 == 4) * (rootcat4 == 1) )  +
                                d15 * ( (rootcat6 == 1) * (rootcat4 == 5) + (rootcat6 == 5) * (rootcat4 == 1) )  +
                                d23 * ( (rootcat6 == 2) * (rootcat4 == 3) + (rootcat6 == 3) * (rootcat4 == 2) )  +
                                d24 * ( (rootcat6 == 2) * (rootcat4 == 4) + (rootcat6 == 4) * (rootcat4 == 2) )  +
                                d25 * ( (rootcat6 == 2) * (rootcat4 == 5) + (rootcat6 == 5) * (rootcat4 == 2) )  +
                                d34 * ( (rootcat6 == 3) * (rootcat4 == 4) + (rootcat6 == 4) * (rootcat4 == 3) )  +
                                d35 * ( (rootcat6 == 3) * (rootcat4 == 5) + (rootcat6 == 5) * (rootcat4 == 3) )  +
                                d45 * ( (rootcat6 == 4) * (rootcat4 == 5) + (rootcat6 == 5) * (rootcat4 == 4) )  ,
                              
                              d12 * ( (rootcat6 == 1) * (rootcat5 == 2) + (rootcat6 == 2) * (rootcat5 == 1) )    +
                                d13 * ( (rootcat6 == 1) * (rootcat5 == 3) + (rootcat6 == 3) * (rootcat5 == 1) )  +
                                d14 * ( (rootcat6 == 1) * (rootcat5 == 4) + (rootcat6 == 4) * (rootcat5 == 1) )  +
                                d15 * ( (rootcat6 == 1) * (rootcat5 == 5) + (rootcat6 == 5) * (rootcat5 == 1) )  +
                                d23 * ( (rootcat6 == 2) * (rootcat5 == 3) + (rootcat6 == 3) * (rootcat5 == 2) )  +
                                d24 * ( (rootcat6 == 2) * (rootcat5 == 4) + (rootcat6 == 4) * (rootcat5 == 2) )  +
                                d25 * ( (rootcat6 == 2) * (rootcat5 == 5) + (rootcat6 == 5) * (rootcat5 == 2) )  +
                                d34 * ( (rootcat6 == 3) * (rootcat5 == 4) + (rootcat6 == 4) * (rootcat5 == 3) )  +
                                d35 * ( (rootcat6 == 3) * (rootcat5 == 5) + (rootcat6 == 5) * (rootcat5 == 3) )  +
                                d45 * ( (rootcat6 == 4) * (rootcat5 == 5) + (rootcat6 == 5) * (rootcat5 == 4) )  ,
                              
                              0,0,0) , 
                  
                  list(       d12 * ( (rootcat7 == 1) * (rootcat1 == 2) + (rootcat7 == 2) * (rootcat1 == 1) )    +
                                d13 * ( (rootcat7 == 1) * (rootcat1 == 3) + (rootcat7 == 3) * (rootcat1 == 1) )  +
                                d14 * ( (rootcat7 == 1) * (rootcat1 == 4) + (rootcat7 == 4) * (rootcat1 == 1) )  +
                                d15 * ( (rootcat7 == 1) * (rootcat1 == 5) + (rootcat7 == 5) * (rootcat1 == 1) )  +
                                d23 * ( (rootcat7 == 2) * (rootcat1 == 3) + (rootcat7 == 3) * (rootcat1 == 2) )  +
                                d24 * ( (rootcat7 == 2) * (rootcat1 == 4) + (rootcat7 == 4) * (rootcat1 == 2) )  +
                                d25 * ( (rootcat7 == 2) * (rootcat1 == 5) + (rootcat7 == 5) * (rootcat1 == 2) )  +
                                d34 * ( (rootcat7 == 3) * (rootcat1 == 4) + (rootcat7 == 4) * (rootcat1 == 3) )  +
                                d35 * ( (rootcat7 == 3) * (rootcat1 == 5) + (rootcat7 == 5) * (rootcat1 == 3) )  +
                                d45 * ( (rootcat7 == 4) * (rootcat1 == 5) + (rootcat7 == 5) * (rootcat1 == 4) )  ,
                              
                              d12 * ( (rootcat7 == 1) * (rootcat2 == 2) + (rootcat7 == 2) * (rootcat2 == 1) )    +
                                d13 * ( (rootcat7 == 1) * (rootcat2 == 3) + (rootcat7 == 3) * (rootcat2 == 1) )  +
                                d14 * ( (rootcat7 == 1) * (rootcat2 == 4) + (rootcat7 == 4) * (rootcat2 == 1) )  +
                                d15 * ( (rootcat7 == 1) * (rootcat2 == 5) + (rootcat7 == 5) * (rootcat2 == 1) )  +
                                d23 * ( (rootcat7 == 2) * (rootcat2 == 3) + (rootcat7 == 3) * (rootcat2 == 2) )  +
                                d24 * ( (rootcat7 == 2) * (rootcat2 == 4) + (rootcat7 == 4) * (rootcat2 == 2) )  +
                                d25 * ( (rootcat7 == 2) * (rootcat2 == 5) + (rootcat7 == 5) * (rootcat2 == 2) )  +
                                d34 * ( (rootcat7 == 3) * (rootcat2 == 4) + (rootcat7 == 4) * (rootcat2 == 3) )  +
                                d35 * ( (rootcat7 == 3) * (rootcat2 == 5) + (rootcat7 == 5) * (rootcat2 == 3) )  +
                                d45 * ( (rootcat7 == 4) * (rootcat2 == 5) + (rootcat7 == 5) * (rootcat2 == 4) )  ,
                              
                              d12 * ( (rootcat7 == 1) * (rootcat3 == 2) + (rootcat7 == 2) * (rootcat3 == 1) )    +
                                d13 * ( (rootcat7 == 1) * (rootcat3 == 3) + (rootcat7 == 3) * (rootcat3 == 1) )  +
                                d14 * ( (rootcat7 == 1) * (rootcat3 == 4) + (rootcat7 == 4) * (rootcat3 == 1) )  +
                                d15 * ( (rootcat7 == 1) * (rootcat3 == 5) + (rootcat7 == 5) * (rootcat3 == 1) )  +
                                d23 * ( (rootcat7 == 2) * (rootcat3 == 3) + (rootcat7 == 3) * (rootcat3 == 2) )  +
                                d24 * ( (rootcat7 == 2) * (rootcat3 == 4) + (rootcat7 == 4) * (rootcat3 == 2) )  +
                                d25 * ( (rootcat7 == 2) * (rootcat3 == 5) + (rootcat7 == 5) * (rootcat3 == 2) )  +
                                d34 * ( (rootcat7 == 3) * (rootcat3 == 4) + (rootcat7 == 4) * (rootcat3 == 3) )  +
                                d35 * ( (rootcat7 == 3) * (rootcat3 == 5) + (rootcat7 == 5) * (rootcat3 == 3) )  +
                                d45 * ( (rootcat7 == 4) * (rootcat3 == 5) + (rootcat7 == 5) * (rootcat3 == 4) )  ,
                              
                              d12 * ( (rootcat7 == 1) * (rootcat4 == 2) + (rootcat7 == 2) * (rootcat4 == 1) )    +
                                d13 * ( (rootcat7 == 1) * (rootcat4 == 3) + (rootcat7 == 3) * (rootcat4 == 1) )  +
                                d14 * ( (rootcat7 == 1) * (rootcat4 == 4) + (rootcat7 == 4) * (rootcat4 == 1) )  +
                                d15 * ( (rootcat7 == 1) * (rootcat4 == 5) + (rootcat7 == 5) * (rootcat4 == 1) )  +
                                d23 * ( (rootcat7 == 2) * (rootcat4 == 3) + (rootcat7 == 3) * (rootcat4 == 2) )  +
                                d24 * ( (rootcat7 == 2) * (rootcat4 == 4) + (rootcat7 == 4) * (rootcat4 == 2) )  +
                                d25 * ( (rootcat7 == 2) * (rootcat4 == 5) + (rootcat7 == 5) * (rootcat4 == 2) )  +
                                d34 * ( (rootcat7 == 3) * (rootcat4 == 4) + (rootcat7 == 4) * (rootcat4 == 3) )  +
                                d35 * ( (rootcat7 == 3) * (rootcat4 == 5) + (rootcat7 == 5) * (rootcat4 == 3) )  +
                                d45 * ( (rootcat7 == 4) * (rootcat4 == 5) + (rootcat7 == 5) * (rootcat4 == 4) )  ,
                              
                              d12 * ( (rootcat7 == 1) * (rootcat5 == 2) + (rootcat7 == 2) * (rootcat5 == 1) )    +
                                d13 * ( (rootcat7 == 1) * (rootcat5 == 3) + (rootcat7 == 3) * (rootcat5 == 1) )  +
                                d14 * ( (rootcat7 == 1) * (rootcat5 == 4) + (rootcat7 == 4) * (rootcat5 == 1) )  +
                                d15 * ( (rootcat7 == 1) * (rootcat5 == 5) + (rootcat7 == 5) * (rootcat5 == 1) )  +
                                d23 * ( (rootcat7 == 2) * (rootcat5 == 3) + (rootcat7 == 3) * (rootcat5 == 2) )  +
                                d24 * ( (rootcat7 == 2) * (rootcat5 == 4) + (rootcat7 == 4) * (rootcat5 == 2) )  +
                                d25 * ( (rootcat7 == 2) * (rootcat5 == 5) + (rootcat7 == 5) * (rootcat5 == 2) )  +
                                d34 * ( (rootcat7 == 3) * (rootcat5 == 4) + (rootcat7 == 4) * (rootcat5 == 3) )  +
                                d35 * ( (rootcat7 == 3) * (rootcat5 == 5) + (rootcat7 == 5) * (rootcat5 == 3) )  +
                                d45 * ( (rootcat7 == 4) * (rootcat5 == 5) + (rootcat7 == 5) * (rootcat5 == 4) )  ,
                              
                              d12 * ( (rootcat7 == 1) * (rootcat6 == 2) + (rootcat7 == 2) * (rootcat6 == 1) )    +
                                d13 * ( (rootcat7 == 1) * (rootcat6 == 3) + (rootcat7 == 3) * (rootcat6 == 1) )  +
                                d14 * ( (rootcat7 == 1) * (rootcat6 == 4) + (rootcat7 == 4) * (rootcat6 == 1) )  +
                                d15 * ( (rootcat7 == 1) * (rootcat6 == 5) + (rootcat7 == 5) * (rootcat6 == 1) )  +
                                d23 * ( (rootcat7 == 2) * (rootcat6 == 3) + (rootcat7 == 3) * (rootcat6 == 2) )  +
                                d24 * ( (rootcat7 == 2) * (rootcat6 == 4) + (rootcat7 == 4) * (rootcat6 == 2) )  +
                                d25 * ( (rootcat7 == 2) * (rootcat6 == 5) + (rootcat7 == 5) * (rootcat6 == 2) )  +
                                d34 * ( (rootcat7 == 3) * (rootcat6 == 4) + (rootcat7 == 4) * (rootcat6 == 3) )  +
                                d35 * ( (rootcat7 == 3) * (rootcat6 == 5) + (rootcat7 == 5) * (rootcat6 == 3) )  +
                                d45 * ( (rootcat7 == 4) * (rootcat6 == 5) + (rootcat7 == 5) * (rootcat6 == 4) )  ,
                              
                              0,0 ) , 
                  
                  list(       d12 * ( (rootcat8 == 1) * (rootcat1 == 2) + (rootcat8 == 2) * (rootcat1 == 1) )    +
                                d13 * ( (rootcat8 == 1) * (rootcat1 == 3) + (rootcat8 == 3) * (rootcat1 == 1) )  +
                                d14 * ( (rootcat8 == 1) * (rootcat1 == 4) + (rootcat8 == 4) * (rootcat1 == 1) )  +
                                d15 * ( (rootcat8 == 1) * (rootcat1 == 5) + (rootcat8 == 5) * (rootcat1 == 1) )  +
                                d23 * ( (rootcat8 == 2) * (rootcat1 == 3) + (rootcat8 == 3) * (rootcat1 == 2) )  +
                                d24 * ( (rootcat8 == 2) * (rootcat1 == 4) + (rootcat8 == 4) * (rootcat1 == 2) )  +
                                d25 * ( (rootcat8 == 2) * (rootcat1 == 5) + (rootcat8 == 5) * (rootcat1 == 2) )  +
                                d34 * ( (rootcat8 == 3) * (rootcat1 == 4) + (rootcat8 == 4) * (rootcat1 == 3) )  +
                                d35 * ( (rootcat8 == 3) * (rootcat1 == 5) + (rootcat8 == 5) * (rootcat1 == 3) )  +
                                d45 * ( (rootcat8 == 4) * (rootcat1 == 5) + (rootcat8 == 5) * (rootcat1 == 4) )  ,
                              
                              d12 * ( (rootcat8 == 1) * (rootcat2 == 2) + (rootcat8 == 2) * (rootcat2 == 1) )    +
                                d13 * ( (rootcat8 == 1) * (rootcat2 == 3) + (rootcat8 == 3) * (rootcat2 == 1) )  +
                                d14 * ( (rootcat8 == 1) * (rootcat2 == 4) + (rootcat8 == 4) * (rootcat2 == 1) )  +
                                d15 * ( (rootcat8 == 1) * (rootcat2 == 5) + (rootcat8 == 5) * (rootcat2 == 1) )  +
                                d23 * ( (rootcat8 == 2) * (rootcat2 == 3) + (rootcat8 == 3) * (rootcat2 == 2) )  +
                                d24 * ( (rootcat8 == 2) * (rootcat2 == 4) + (rootcat8 == 4) * (rootcat2 == 2) )  +
                                d25 * ( (rootcat8 == 2) * (rootcat2 == 5) + (rootcat8 == 5) * (rootcat2 == 2) )  +
                                d34 * ( (rootcat8 == 3) * (rootcat2 == 4) + (rootcat8 == 4) * (rootcat2 == 3) )  +
                                d35 * ( (rootcat8 == 3) * (rootcat2 == 5) + (rootcat8 == 5) * (rootcat2 == 3) )  +
                                d45 * ( (rootcat8 == 4) * (rootcat2 == 5) + (rootcat8 == 5) * (rootcat2 == 4) )  ,
                              
                              d12 * ( (rootcat8 == 1) * (rootcat3 == 2) + (rootcat8 == 2) * (rootcat3 == 1) )    +
                                d13 * ( (rootcat8 == 1) * (rootcat3 == 3) + (rootcat8 == 3) * (rootcat3 == 1) )  +
                                d14 * ( (rootcat8 == 1) * (rootcat3 == 4) + (rootcat8 == 4) * (rootcat3 == 1) )  +
                                d15 * ( (rootcat8 == 1) * (rootcat3 == 5) + (rootcat8 == 5) * (rootcat3 == 1) )  +
                                d23 * ( (rootcat8 == 2) * (rootcat3 == 3) + (rootcat8 == 3) * (rootcat3 == 2) )  +
                                d24 * ( (rootcat8 == 2) * (rootcat3 == 4) + (rootcat8 == 4) * (rootcat3 == 2) )  +
                                d25 * ( (rootcat8 == 2) * (rootcat3 == 5) + (rootcat8 == 5) * (rootcat3 == 2) )  +
                                d34 * ( (rootcat8 == 3) * (rootcat3 == 4) + (rootcat8 == 4) * (rootcat3 == 3) )  +
                                d35 * ( (rootcat8 == 3) * (rootcat3 == 5) + (rootcat8 == 5) * (rootcat3 == 3) )  +
                                d45 * ( (rootcat8 == 4) * (rootcat3 == 5) + (rootcat8 == 5) * (rootcat3 == 4) )  ,
                              
                              d12 * ( (rootcat8 == 1) * (rootcat4 == 2) + (rootcat8 == 2) * (rootcat4 == 1) )    +
                                d13 * ( (rootcat8 == 1) * (rootcat4 == 3) + (rootcat8 == 3) * (rootcat4 == 1) )  +
                                d14 * ( (rootcat8 == 1) * (rootcat4 == 4) + (rootcat8 == 4) * (rootcat4 == 1) )  +
                                d15 * ( (rootcat8 == 1) * (rootcat4 == 5) + (rootcat8 == 5) * (rootcat4 == 1) )  +
                                d23 * ( (rootcat8 == 2) * (rootcat4 == 3) + (rootcat8 == 3) * (rootcat4 == 2) )  +
                                d24 * ( (rootcat8 == 2) * (rootcat4 == 4) + (rootcat8 == 4) * (rootcat4 == 2) )  +
                                d25 * ( (rootcat8 == 2) * (rootcat4 == 5) + (rootcat8 == 5) * (rootcat4 == 2) )  +
                                d34 * ( (rootcat8 == 3) * (rootcat4 == 4) + (rootcat8 == 4) * (rootcat4 == 3) )  +
                                d35 * ( (rootcat8 == 3) * (rootcat4 == 5) + (rootcat8 == 5) * (rootcat4 == 3) )  +
                                d45 * ( (rootcat8 == 4) * (rootcat4 == 5) + (rootcat8 == 5) * (rootcat4 == 4) )  ,
                              
                              d12 * ( (rootcat8 == 1) * (rootcat5 == 2) + (rootcat8 == 2) * (rootcat5 == 1) )    +
                                d13 * ( (rootcat8 == 1) * (rootcat5 == 3) + (rootcat8 == 3) * (rootcat5 == 1) )  +
                                d14 * ( (rootcat8 == 1) * (rootcat5 == 4) + (rootcat8 == 4) * (rootcat5 == 1) )  +
                                d15 * ( (rootcat8 == 1) * (rootcat5 == 5) + (rootcat8 == 5) * (rootcat5 == 1) )  +
                                d23 * ( (rootcat8 == 2) * (rootcat5 == 3) + (rootcat8 == 3) * (rootcat5 == 2) )  +
                                d24 * ( (rootcat8 == 2) * (rootcat5 == 4) + (rootcat8 == 4) * (rootcat5 == 2) )  +
                                d25 * ( (rootcat8 == 2) * (rootcat5 == 5) + (rootcat8 == 5) * (rootcat5 == 2) )  +
                                d34 * ( (rootcat8 == 3) * (rootcat5 == 4) + (rootcat8 == 4) * (rootcat5 == 3) )  +
                                d35 * ( (rootcat8 == 3) * (rootcat5 == 5) + (rootcat8 == 5) * (rootcat5 == 3) )  +
                                d45 * ( (rootcat8 == 4) * (rootcat5 == 5) + (rootcat8 == 5) * (rootcat5 == 4) )  ,
                              
                              d12 * ( (rootcat8 == 1) * (rootcat6 == 2) + (rootcat8 == 2) * (rootcat6 == 1) )    +
                                d13 * ( (rootcat8 == 1) * (rootcat6 == 3) + (rootcat8 == 3) * (rootcat6 == 1) )  +
                                d14 * ( (rootcat8 == 1) * (rootcat6 == 4) + (rootcat8 == 4) * (rootcat6 == 1) )  +
                                d15 * ( (rootcat8 == 1) * (rootcat6 == 5) + (rootcat8 == 5) * (rootcat6 == 1) )  +
                                d23 * ( (rootcat8 == 2) * (rootcat6 == 3) + (rootcat8 == 3) * (rootcat6 == 2) )  +
                                d24 * ( (rootcat8 == 2) * (rootcat6 == 4) + (rootcat8 == 4) * (rootcat6 == 2) )  +
                                d25 * ( (rootcat8 == 2) * (rootcat6 == 5) + (rootcat8 == 5) * (rootcat6 == 2) )  +
                                d34 * ( (rootcat8 == 3) * (rootcat6 == 4) + (rootcat8 == 4) * (rootcat6 == 3) )  +
                                d35 * ( (rootcat8 == 3) * (rootcat6 == 5) + (rootcat8 == 5) * (rootcat6 == 3) )  +
                                d45 * ( (rootcat8 == 4) * (rootcat6 == 5) + (rootcat8 == 5) * (rootcat6 == 4) )  ,
                              
                              d12 * ( (rootcat8 == 1) * (rootcat7 == 2) + (rootcat8 == 2) * (rootcat7 == 1) )    +
                                d13 * ( (rootcat8 == 1) * (rootcat7 == 3) + (rootcat8 == 3) * (rootcat7 == 1) )  +
                                d14 * ( (rootcat8 == 1) * (rootcat7 == 4) + (rootcat8 == 4) * (rootcat7 == 1) )  +
                                d15 * ( (rootcat8 == 1) * (rootcat7 == 5) + (rootcat8 == 5) * (rootcat7 == 1) )  +
                                d23 * ( (rootcat8 == 2) * (rootcat7 == 3) + (rootcat8 == 3) * (rootcat7 == 2) )  +
                                d24 * ( (rootcat8 == 2) * (rootcat7 == 4) + (rootcat8 == 4) * (rootcat7 == 2) )  +
                                d25 * ( (rootcat8 == 2) * (rootcat7 == 5) + (rootcat8 == 5) * (rootcat7 == 2) )  +
                                d34 * ( (rootcat8 == 3) * (rootcat7 == 4) + (rootcat8 == 4) * (rootcat7 == 3) )  +
                                d35 * ( (rootcat8 == 3) * (rootcat7 == 5) + (rootcat8 == 5) * (rootcat7 == 3) )  +
                                d45 * ( (rootcat8 == 4) * (rootcat7 == 5) + (rootcat8 == 5) * (rootcat7 == 4) )  ,
                              0 ) 
    )
    
    
    
    emdc_settings$utilityOutside = delta_outside
    emdc_settings$utilities = V
    emdc_settings$gamma = gamma
    emdc_settings$delta = delta
    
    # 

    P[["model"]] = apollo_emdc(emdc_settings, functionality)
    
    ### Take product across observation for same individual
    P = apollo_panelProd(P, apollo_inputs ,functionality)
    
    P = apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
  }
  
  apollo_prediction = function (model, apollo_probabilities = apollo_probabilities2, apollo_inputs, prediction_settings = list(), modelComponent = NA) {
    prediction_settings = list()
    modelComponent = NA
    
    if (is.null(prediction_settings$modelComponent)) {
      if (exists("modelComponent")) 
        prediction_settings$modelComponent = modelComponent
      else prediction_settings$modelComponent = NA
    }
    if (is.null(prediction_settings$runs))   prediction_settings$runs = 1
    if (is.null(prediction_settings$silent))   prediction_settings$silent = FALSE
    silent = prediction_settings$silent
    if (!is.null(apollo_inputs$silent) && apollo_inputs$silent) silent = TRUE
    
    if (is.null(prediction_settings$nRep)) prediction_settings$nRep <- 100L
    
    if (is.null(prediction_settings$summary)) prediction_settings$summary <- TRUE
    
    modelComponent = prediction_settings$modelComponent
    runs = prediction_settings$runs
    apollo_inputs$nRep <- prediction_settings$nRep
    
    apollo_compareInputs(apollo_inputs)
    apollo_randCoeff = apollo_inputs[["apollo_randCoeff"]]
    apollo_lcPars = apollo_inputs[["apollo_lcPars"]]
    apollo_checkArguments(apollo_probabilities, apollo_randCoeff, 
                          apollo_lcPars)
    if (!silent)     apollo_print("Running predictions from model using parameter estimates...")
    
    
    ###################HERE IT REQUIRES APOLLO_PROBABILITIES
    predictions = apollo_probabilities(apollo_beta, apollo_inputs, 
                                       functionality = "prediction")
    
    predictionamount = predictions$model[,1:9] 
    
    return(predictionamount)
  }
  
  apollo_inputs = apollo_validateInputs(apollo_beta = apollo_beta2)
  
  predictionprob2 = apollo_prediction(model, apollo_probabilities = apollo_probabilities2, apollo_inputs, prediction_settings=list(runs=1), modelComponent = NA)
  
  
  finalpredictionprob = predictionprob1 * membership + predictionprob2 * (1 - membership )

  
  pricema = database[,103:110]
  
  adclickma = exp(-0.130*pricema)/(1+exp(-0.130*pricema))

  
  finalclickprob = t( apply (cbind(database[,15:22],finalpredictionprob[,2:9], adclickma), 1, add_restrict_imputevalue) )
  
  return(finalclickprob)
  
}

  
  
  
# ########################################################################################
# # we use demand model code to predict expected number of clicks given datam and predicted n_ad_vec under S_vec
# ########################################################################################
# 

database1 = read.csv("cleaned data.csv",header=TRUE)


i=1

while(i<=nrow(database1)){
  pv_id_temp = database1$pv_id2[i]
  # find all index where datam$pv_id == pv_id_temp
  index_temp = which(datam$pv_id2 == pv_id_temp)
  if(length(index_temp)==0){
    i=i+1
  }else{
    # create a 8 by 1 vector recording the subcat_id3 with positive n_ad for this pv
    index_temp2 = which(n_ad_vec[index_temp]>0)
    n_ad_temp = n_ad_vec[index_temp][index_temp2]
    
    # update subcat_id
    database1[i, 16:23] = rep(datam$subcat_id3[index_temp[index_temp2]], times = n_ad_temp)
    
    # update npv, click, buy
    index_temp3 = rep(index_temp[index_temp2], times = n_ad_temp)
    database1[i, 32:39] = log(datam$npv14[index_temp3]+1)
    database1[i, 41:48] = log(datam$click14[index_temp3]+1)
    database1[i, 50:57] = datam$buy14[index_temp3]
    
    # update rcc_pv, rcc_click, rcc_buy
    database1[i, 59:66] = log(datam$rcc_pv[index_temp3]+1)
    database1[i, 68:75] = log(datam$rcc_click[index_temp3]+1)
    database1[i, 77:84] = log(datam$rcc_buy[index_temp3]+1)
    
    # update availnew1 to availnew9
    temp = database1[i, 17:23] - database1[i, 16:22]
    temp = c(as.numeric(database1[i, 17]), as.vector(temp))
    temp2 = rep(0, 8)
    temp2[temp>0] =1
    database1[i, 86:93] = temp2
    
    if(i%%500==0){
      cat(i, fill=TRUE)
    }
    i=i+1
  }
}

database = database1
database0 = database1



apollo_beta1 = c(
  gamma_v1      = 1.301, 
  gamma_v2      = 0.600, 
  gamma_v3      = .501, 
  gamma_v4      = .477, 
  gamma_v5      = .664, 
  delta_v1      = -3.967,
  delta_v2         = -3.849,
  delta_v3       = -3.897,
  delta_v4     = -3.357,
  delta_v5 = -3.348,
  delta_outside      = 0,
  
  
  
  
  bweekend = 0.088, bnpv = 0.026,  bnclick = 0.109,  bnbuy = 	-0.080  , nrccpv=0.010, nrccclick=0.348, nrccbuy=0.025, 
  badstockpsai = 0.637 ,  bresidual = 0.139 , bvariety = 0.043,
  # Compl/subst
  d12= 0.005, d13= 0.014, 
  d14= 0.022, d15= 0.003, 
  d23= 0.015, d24= 0.017, 
  d25= 0.014, d34= 0.023, 
  d35= 0.042, d45= -0.002
  
  
  
)

apollo_beta2 = c(
  
  gamma_v1      = 0.903	, 
  gamma_v2      = 0.557, 
  gamma_v3      = 0.489, 
  gamma_v4      = 0.397, 
  gamma_v5      = 0.197, 
  delta_v1      = -2.536,
  delta_v2         = -2.465,
  delta_v3       = -2.394,
  delta_v4     = -2.529,
  delta_v5 = -2.412,
  delta_outside      = 0,
  
  bweekend = -0.029	, bnpv = 0.024,  bnclick = 0.062,  bnbuy = -0.727 , nrccpv=0.037, nrccclick=0.196, nrccbuy=	-0.476, 
  badstockpsai = 0.305	 ,  bresidual = 0.137 , bvariety = 0.026	,
  d12= 0.009, d13= 0.001, 
  d14= 0.006, d15= -0.006, 
  d23= 0.005, d24= 0.008, 
  d25= -0.007, d34= 0.023, 
  d35= 0.013, d45= -0.008
)

membership0 = read.csv("membership.csv",header=TRUE)
membership = membership0[,2]


apollo_fixed = c( "delta_outside", "gamma_outside", "sigma"    ) #, "delta_v1_b"



set.seed(11)

#####################Prediction results
t1<-Sys.time()
t1

pred_click = choiceprobpred(database1)

summary(pred_click)

t2<-Sys.time()

timeprint = t2-t1
timeprint

# ########################################################################################








########################################################################################
########################################################################################
# function of calculating expected profit, payment, cpc, ctr,  
# for the auction for one individual impression
# Here the CPC is based on the second highest bid within each candidate subcat, which is different from previous code
EPROFIT=function(n_ad_disp, n_pot_ad, ad_index, ctrvec,
                 b, QS, v, tao, u){
  # Arguments
  # n_ad_disp: a N by 1 vector which records the assigned ad quota for each potential category for the focal consumer,
  #            N=n_pot_cat
  # n_pot_ad: a vector records the number of competing ads in each potential category.it has the same length as S
  # ad_index: a M by 1 vector records the index of all potential ads in datab. it has the length of sum(n_pot_ad)=M
  # pv_id: a N by 1 vector records the index of pv in datam
  # ctrvec: a 8 by 1 vector recording the CTR for the subcat of the 8 ads chosen for display (a function of n_ad_disp)
  
  # b: a M by 1 vector of bid from potential ads
  # QS: a M by 1 vector of quality score of potential ads
  # v: a M by 1 vector of value-per-click of potential ads
  # tao: variance of error term in ranking model
  # u: a M by 1 vector, where each element is drawn from rnorm(0,1)
 
  # Output
  # profit: total profit of advertisers
  # payment: total payment of advertisers
  # cpc: average cpc in this auction
  # ctr: average CTR in this auction
  
  N=length(n_ad_disp)
  M=length(ad_index)
  
  # define a M by 1 vector mark, for those selected ads, mark=1
  mark=rep(0, M)
  
  # define a M by 1 vector mark2, for all potential ads from subcats with positive n_ad_disp
  mark2=rep(0, M)
  
  # define a 8 by 1 vector, recording the index of ad
  index_chosen=rep(0, 8)
  
  # define cpc as a 8 by 1 vector
  cpc = rep(0, 8)
  
  i=1
  j=1
  l=1
  
  while(i<=N){
    if(n_ad_disp[i]>0){
      ntemp=min(n_ad_disp[i], n_pot_ad[i]) # deal with the case where n_pot_ad is too small
      npottemp=n_pot_ad[i]
      indextemp=j:(j+npottemp-1)
      mark2[indextemp]=1
      
      wb_temp=b[indextemp]*QS[indextemp]*exp(u[indextemp]*tao^0.5)
      # attache a 0 at the end of wbtemp for the calculation of cpc in case all 8 ads were chosen from one subcat
      wb_temp=c(wb_temp, 0)
      QS_temp=c(QS[indextemp], 1)
      error_temp=c(exp(u[indextemp]*tao^0.5), 1)
      
      # sort wb_temp
      out=sort(wb_temp, decreasing=TRUE, index.return=TRUE)
      mark[indextemp[out$ix[1:ntemp]]]=1
      index_chosen[l:(l+ntemp-1)]=indextemp[out$ix[1:ntemp]]
      
      # calculate cpc
      cpc[l:(l+ntemp-1)] = wb_temp[out$ix[1:ntemp+1]]/(QS_temp[out$ix[1:ntemp]]*error_temp[out$ix[1:ntemp]])
      
      l=l+ntemp
    }
    j=j+n_pot_ad[i]
    i=i+1
  }
  
  index_chosen2=index_chosen[which(index_chosen>0)]
  ctrvec2=ctrvec[which(index_chosen>0)]
  
  # # Next we calculate CPC for 8 chosen ads across all potential categories
  # wb_all=b*QS*exp(u*tao^0.5)
  # # we only focus on weighted bids from those ads from subcats with positive n_ad_disp
  # # We set the weighted bid of ads from non-chosen subcats to be relatively small
  # wb_all[which(mark2==0)]=min(wb_all)-1
  # 
  # out=sort(wb_all, decreasing=TRUE, index.return=TRUE)
  # 
  # # define the index of second-highest
  # index_2nd=mapply(min, match(index_chosen, out$ix)+1, length(wb_all))
  # index_2nd=out$ix[index_2nd]
  # 
  # # define CPC for chosen ads
  # cpc=wb_all[index_2nd]/(QS[index_chosen]*exp(u[index_chosen]*tao^0.5))
  
  
  # define profit: from clicks
  profit=ctrvec2*(v[index_chosen2]-cpc[which(index_chosen>0)])
  
  # define payment
  pay=ctrvec2*cpc[which(index_chosen>0)]
  
  # define avgcpc
  if(sum(ctrvec2)>0){
    avgcpc=sum(ctrvec2*cpc[which(index_chosen>0)])/sum(ctrvec2)
  }else{
    avgcpc=mean(cpc[which(index_chosen>0)])
  }
  # avgcpc=mean(b[index_chosen])
  
  # define avgvpc
  if(sum(ctrvec2)>0){
    avgvpc=sum(ctrvec2*v[index_chosen])/sum(ctrvec2)
  }else{
    avgvpc=mean(v[index_chosen])
  }
  
  return(list(profit=sum(profit), pay=sum(pay), ctr=mean(ctrvec2), cpc=avgcpc, vpc=avgvpc)) 
}



########################################################################################


datam=read.csv("Delivery_Match_1day.csv", header=T)
datab=read.csv("VPC_Reg_ratio=0.csv", header=T)

T=max(datab$thedate) # we only use the first 12 days in estimation
datam=datam[datam$thedate<=T, ]


load("Supply 1 data.Rdata")

load("bideqm_1day_20231205.RData")

########################################################################################
# We first modify the date_subcat_id and N_ad in datam
# First, we check how many date_subcat_id in datam has no correspondence in datab
N_datecat_all=length(unique(datam$date_subcat_id))

# For each date_subcat_id in datam, 
# if it doesnt appear in datab, we check if the corresponding subcat has other date_subcat appearing in datab
# If yes, we replace original date_subcat_id with other date_subcat_id associated with this subcat
# Otherwise, we leave the oringal date_subcat_id

# we define a matrix N_datecat_all by 3, the first column records if this date_subcat_id needs a change
# the second column records the new date_subcat_id it changes to
# the third column records the new N_ad
datecatmat=matrix(0, N_datecat_all, 3)

i=1
set.seed(123)
while(i<=N_datecat_all){
  index=which(datab$date_subcat_id==i)
  if(length(index)==0){
    subcat_id=datam$subcat_id3[match(i, datam$date_subcat_id)]
    # find all date_subcat_id of this subcat_id in datab
    temp=unique(datab$date_subcat_id[which(datab$subcat_id3==subcat_id)])
    # check if temp is empty
    if(length(temp)>0){
      # then we randomly select one date_cat_id in temp for the old date_subcat_id==i
      id=sample(1:(length(temp)), 1)
      datecatmat[i, 1]=1
      datecatmat[i, 2]=temp[id]
      datecatmat[i, 3]=length(which(datab$date_subcat_id==datecatmat[i, 2]))
    }
  }
  i=i+1
  if(i%%100==0){
    cat(i, fill=TRUE)
  }
}
rm(.Random.seed)

# Next we modify date_subcat_id and N_ad in datam
i=1
while(i<=N_datecat_all){
  if(datecatmat[i, 1]==1){
    index1=which(datam$date_subcat_id==i)
    index2=match(datecatmat[i, 2], datam$date_subcat_id)
    
    # modify
    datam$date_subcat_id[index1]=datecatmat[i, 2]
    datam$N_ad[index1]=datecatmat[i, 3]
  }
  i=i+1
  if(i%%100==0){
    cat(i, fill=TRUE)
  }
}

summary(datam$date_subcat_id)
summary(datam$N_ad)

# Finally we deal with those date_subcat_id==0
index=which(datam$date_subcat_id==0)

i=1
set.seed(123)
while(i<=length(index)){
  j=index[i]
  
  subcat_id=datam$subcat_id3[j]
  # find all date_subcat_id of this subcat_id in datab
  temp=unique(datab$date_subcat_id[which(datab$subcat_id3==subcat_id)])
  # check if temp is empty
  if(length(temp)>0){
    # then we randomly select one date_cat_id in temp for the old date_subcat_id==i
    id=sample(1:(length(temp)), 1)
    datecatidnew=temp[id]
    
    # modify
    index2=match(datecatidnew, datam$date_subcat_id)
    datam$date_subcat_id[j]=datecatidnew
    datam$N_ad[j]=length(which(datab$date_subcat_id==datecatidnew))
  }
  i=i+1
  if(i%%1000==0){
    cat(i, fill=TRUE)
  }
}
rm(.Random.seed)
length(which(datam$date_subcat_id==0))/nrow(datam)

summary(datam$date_subcat_id)
summary(datam$N_ad)

########################################################################################
# take the mean of beta_u_mat
beta_u_vec=apply(beta_u_mat, 2, mean)
beta_s_vec=apply(beta_s_mat, 2, mean)


# we add one row to datab, this row refers to a hypothetical adid which belongs to all datesubcat that are missing in datab
datab=rbind(datab, 0)
datab$adid2[nrow(datab)]=max(datab$adid2[-nrow(datab)])+1
datab$rprice[nrow(datab)]=median(datab$rprice[-nrow(datab)])
datab$qs_post[nrow(datab)]=median(datab$qs_post[-nrow(datab)])



########################################################################################
# Set parameters
T=max(datab$thedate)
N_user=length(unique(datam$user_id3))
N_ad=8 # no of ads for display
N_cat_all=length(unique(datam$subcat_id3))
N_userroot=length(unique(datam$usroot_id)) # here the userroot refers to <user, rootcat>
N_ad_all=length(unique(datab$adid2))
N_datecat_all=max(unique(datab$date_subcat_id))

# no of pvs
N_pv=length(unique(datam$pv_id2))

# some useful index vectors
pv_id=datam$pv_id2
cat_id_s=datam$subcat_id3
userroot_id_s=datam$usroot_id
date_cat_id_s=datam$date_subcat_id
ad_id_b=datab$adid2
date_cat_id_b=datab$date_subcat_id

QS_post=datab$qs_post

vpc=datab$vpc
bid=interout[, 2]
vpc=c(vpc, median(vpc))
bid=c(bid, median(bid))
summary(vpc)
summary(datab$vpc)
summary(bid)
summary(datab$bid)

# define the vector of interest score
S_vec=datam$S_vec
# define tao, the variance of error term in ranking model
tao=1.208



########################################################################################
# create a list, which records the following information for each pv_id
# n_pot_subcat: a number of competing subcats
# pv_id: the corresponding rows of pvs in datam
# user_id: a number associated with the current pv
# date_subcat: a n_pot_subcat by 1 vector recording those date_subcat_id for this pv
# n_ad_vec: a n_pot_subcat by 1 vector recording the n_ad for each candidate subcat
# ad_index: a M by 1 vecotr recording the rows in datab for all potential ads asscoiated with all potential subcat
# if the datesubcat_id does not have correspondence in datab, we refer it to the last row in datab

pv_list=list()

i=1
count=0 # count records those date_subcat that has to refer to hypothetical advertisers in the last row of datab
while(i<=N_pv){
  pv_id_temp=which(datam$pv_id2==i)
  j=pv_id_temp[1]
  n_pot_subcat=datam$N_pot_subcat[j]
  user_id=datam$user_id3[j]
  
  n_ad_vec=ceiling(datam$N_ad[pv_id_temp]) # some N_ad might not be integers
  date_subcat=date_cat_id_s[pv_id_temp]
  ad_index=rep(0, sum(n_ad_vec))
  j=1
  k=1
  while(j<=length(date_subcat)){
    indextemp=which((date_cat_id_b %in% date_subcat[j])==TRUE)
    if(length(indextemp)==0 | date_subcat[j]==0){
      count=count+1
      ad_index[k:(k+n_ad_vec[j]-1)]=nrow(datab)
    }else{
      if(length(indextemp)!=n_ad_vec[j]){
        cat("pv_id: ", i, fill=TRUE)
        cat("n_pot_ad is incorrect for the jth subcat ", j, fill=TRUE)
        cat("\n")
      }else{
        ad_index[k:(k+n_ad_vec[j]-1)]=indextemp
      }
    }
    k=k+n_ad_vec[j]
    j=j+1
  }
  
  
  pv_list[[i]]=list(n_pot_subcat=n_pot_subcat, pv_id=pv_id_temp, user_id=user_id, date_subcat=date_subcat,
                    n_ad_vec=n_ad_vec, ad_index=ad_index)
  
  if(i%%100==0){
    cat(i, fill=TRUE)
  }
  i=i+1
}
count



t2<-Sys.time()
t2

timeprint = t2-t1
timeprint



########################################################################################
# transfer S_vec into a length(S_vec) by M matrix, where each element is the corresponding n_ad
# generate a vec Sim_vec, N_pv by 1, each element is the corresponding n_related2 for an impression
n_ad_vec=rep(0, nrow(datam))
# Sim_vec=rep(0, N_pv)

i=1
j=1
set.seed(123)

while(i<=nrow(datam)){
  n_pot=datam$N_pot_subcat[i]
  S_temp=S_vec[i:(i+n_pot-1)]
  prob=exp(S_temp)/sum(exp(S_temp))
  n_ad_temp=rmultinom(1, 8, prob)
  n_ad_vec[i:(i+n_pot-1)]=n_ad_temp
  
  
  i=i+n_pot
  j=j+1
  if(j%%100==0){
    cat(j, fill=TRUE)
  }
}
rm(.Random.seed)




t4<-Sys.time()
t4

########################################################################################
# Start formal simulations
ndraw=100

# define four matrix
profitmat=matrix(0, N_pv, ndraw)
paymat=matrix(0, N_pv, ndraw)
clickmat=matrix(0, N_pv, ndraw)
cpcmat=matrix(0, N_pv, ndraw)
vpcmat=matrix(0, N_pv, ndraw)
n_subcat_disp_mat=matrix(0, N_pv, ndraw)
avg_n_ad_disp_mat=matrix(0, N_pv, ndraw)

# define n_ad_vec_mat
n_ad_vec_mat = matrix(0, nrow(datam), ndraw)

m=1

set.seed(123)

while(m<=ndraw){
  # for each m, simulate n_ad_vec and Sim_vec first
  k=1 # pointer for n_ad_vec
  n_ad_vec=rep(0, nrow(datam))
  
  # S_vec_draw=S_vec+rnorm(length(S_vec), 0, sd=sd_s)
  S_vec_draw=S_vec
  while(k<=nrow(datam)){
    n_pot=datam$N_pot_subcat[k]
    S_temp=S_vec_draw[k:(k+n_pot-1)]
    prob=exp(S_temp)/sum(exp(S_temp))
    n_ad_temp=rmultinom(1, 8, prob)
    n_ad_vec[k:(k+n_pot-1)]=n_ad_temp
    k=k+n_pot
  }
  n_ad_vec_mat[, m] = n_ad_vec
  m=m+1
}
rm(.Random.seed)

 summary(interout[, 2]/datab$bid[1:nrow(interout)])
 index=1:nrow(interout)
 length(index)/nrow(interout)
 median(interout[index, 2]/datab$bid[index])
 bid=interout[, 2]
 bid[-index]=datab$bid[-c(index, nrow(datab))]*(median(interout[index, 2]/datab$bid[index]))
 summary(bid/datab$bid[1:nrow(interout)])
 bid=c(bid, median(bid))


begin=proc.time()[3]

m=1

while(m<=ndraw){
  n_ad_vec =n_ad_vec_mat[, m]
  
 
  # Given n_ad_vec and Sim_vec, predict outcome for each pv
  i=1

  while(i<=N_pv){ #i<=N_pv
    n_ad_disp=n_ad_vec[pv_list[[i]]$pv_id]
    n_subcat_disp=length(n_ad_disp[n_ad_disp>0])
    avg_n_ad_disp=mean(n_ad_disp[n_ad_disp>0])

    # define n_pot_ad
    n_pot_ad=pv_list[[i]]$n_ad_vec

    # define user_id
    user_id=pv_list[[i]]$user_id

    # define ad_index
    ad_index=pv_list[[i]]$ad_index

    ######################################################
    # we next calculate ctr using the code written by Yao
    # match pv_id with datam
    database1 = database0[database0$pv_id2==i, ]
    
    # change clicknew1-clicknew7 to 0
    database1[, 6:13] = 0
    database1[, 103:110] = 0
    
    # update other variables in database1, index_temp is index where datam$pv_id == i
    index_temp=pv_list[[i]]$pv_id
    # create a 8 by 1 vector recording the subcat_id3 with positive n_ad for this pv
    index_temp2 = which(n_ad_disp>0)
    n_ad_temp = n_ad_vec[index_temp][index_temp2]
    
    # update subcat_id
    database1[, 15:22] = rep(datam$subcat_id3[index_temp[index_temp2]], times = n_ad_temp)
    
    # update npv, click, buy
    index_temp3 = rep(index_temp[index_temp2], times = n_ad_temp)
    database1[, 24:31] = log(datam$npv14[index_temp3]+1)
    database1[, 33:40] = log(datam$nclick14[index_temp3]+1)
    database1[, 42:49] = datam$nbuy14[index_temp3]
    
    # update rcc_pv, rcc_click, rcc_buy
    database1[, 51:58] = log(datam$rcc_pv[index_temp3]+1)
    database1[, 60:67] = log(datam$rcc_click[index_temp3]+1)
    database1[, 69:76] = log(datam$rcc_buy[index_temp3]+1)
    
    # update availnew1 to availnew9
    temp = database1[, 16:22] - database1[, 15:21]
    temp = c(as.numeric(database1[, 16]), as.vector(temp))
    temp2 = rep(0, 8)
    temp2[temp!=0] = 1 #original was temp2[temp>0] =1
    database1[, 78:85] = temp2
    
    
    # we need to replicate database1 to have two same rows to run the code
    database1 = rbind(database1, database1)
    database = database1
    
    membership = membership0[membership0[,1] == database1$user_id3 , 2]
    
    pred_click = choiceprobpred(database1)
    # pred_click
    ctrvec = pred_click[1, ]
    

    ######################################################
    # define bid, vpc, qscore and u
    v=vpc[ad_index]
    b=bid[ad_index]
    QS=QS_post[ad_index]

    # draw random numbers
    u=rnorm(sum(n_pot_ad))
    # rm(.Random.seed)

    # get profit, pay, click, and cpc
    out=EPROFIT(n_ad_disp, n_pot_ad, ad_index, ctrvec,
                b, QS, v, tao, u)

    # record results
    profitmat[i, m]=out$profit
    paymat[i, m]=out$pay
    clickmat[i, m]=out$ctr
    cpcmat[i, m]=out$cpc
    vpcmat[i, m]=out$vpc
    n_subcat_disp_mat[i, m]=n_subcat_disp
    avg_n_ad_disp_mat[i, m]=avg_n_ad_disp


    i=i+1
    
    if (i%%1000 == 0) {
      cat(i, fill=TRUE)
      print(Sys.time())
    }
    
    
  }

  if(m %% 1==0){
    cat("iteration of i, m", c(i,m), "\n")
    output_1day=list(profit=profitmat, pay=paymat, cpc=cpcmat, click=clickmat, vpc=vpcmat, n_subcat_disp=n_subcat_disp_mat, avg_n_ad_disp=avg_n_ad_disp_mat)
    save(output_1day,
         file = "results_from_CF_aucsim_1day.RData")

  }
  
  m=m+1
}


end=proc.time()[3] 
cat("Time used:", end-begin, fill=TRUE)


t5<-Sys.time()
t5


timeprint = t5-t4
timeprint

