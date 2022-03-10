#' Compute difference of monadic coefficients across blocks
#' 
#' The function computes, for a given monadic predictor name, the pair-wise
#' differences in coefficients associated with that predictor across blocks.
#'
#' @param model An object inhereting from class \code{mmsbm}.
#' @param .var Character string with exact name of predictor. Must match 
#'             one of the Monadic predictor names reported by \code{summary()}
#'             function. 
#' @param family For bipartite networks: equal to \code{1} for family 1,
#'               or equal to \code{2} for family 2. Deafults to \code{1}.  
#' @param hmm.state For dynamic models: integer value indicating the HMM state
#'                  used for selecting the coefficient matrix. 
#'
#' @return Table of coefficients differences and their p-values. 
#' @export
#'
#' @examples
coef_delta <- function(model, .var, family = 1, hmm.state = 1){
  vc <- vcov(model)
  if(family == 1){
    vc <- vc$MonadCoef1
    coefs <- model$MonadCoef1
  } else {
    vc <- vc$MonadCoef2
    coefs <- model$MonadCoef2
  }
  vc_idx <- grep(.var, colnames(vc)) 
  vc_sub <- vc[vc_idx, vc_idx]
  all_vc <- diag(vc_sub)
  c_idx <- grep(.var, rownames(coefs))
  n_diff <- sum(lower.tri(vc_sub))
  diffs <- array(NA, n_diff)
  se_diffs <- array(NA, n_diff)
  compare_idx <- which(lower.tri(vc_sub), arr.ind = TRUE)
  for(i in 1:nrow(compare_idx)){
    diffs[i] <- coefs[c_idx, compare_idx[i, 1], hmm.state] - coefs[c_idx, compare_idx[i, 2], hmm.state]
    se_diffs[i] <- sqrt(all_vc[compare_idx[i, 1]] + all_vc[compare_idx[i, 2]] - 2*vc_sub[compare_idx[i,, drop=FALSE]])
  }
  res <- cbind(round(diffs, 4), round(se_diffs,2), 
               round(diffs/se_diffs,3), 
               round(pnorm(abs(diffs/se_diffs), lower.tail = FALSE)*2,3))
  colnames(res) <- c("Est. Difference", "Std. Error", "z value","Pr(>|z|)")
  rownames(res) <- paste("Group",compare_idx[,1], "vs. Group",compare_idx[,2])
  return(res)
}
