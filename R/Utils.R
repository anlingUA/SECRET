#' cell type proportion results evaluation
#'
#' @param trueP true proportion
#' @param estP estimated proportion
#'
#' @return RMSD, mAD, R
#' @export
#'
eval_prop <- function(trueP, estP){
  if (all(sapply(lapply(estP,FUN=rownames), FUN = identical, rownames(trueP)))){
    evals<-lapply(estP, function(a1){
      RMSD <- round(sqrt(mean(as.matrix((trueP - a1)^2), na.rm = T)), digits = 3)
      mAD <- round(mean(as.matrix(abs(trueP - a1)), na.rm = T), digits = 3)
      pearson <- round(cor(c(as.matrix(trueP)), c(as.matrix(a1))), digits = 3)

      return(list(estP = a1, RMSD = RMSD, mAD = mAD, pearson = abs(pearson)))
    })
    est_overall <- NULL
    for (l in 1:length(estP)){
      est_overall <- rbind(est_overall, c(evals[[l]]$RMSD, evals[[l]]$mAD, evals[[l]]$pearson))
    }
    colnames(est_overall) <- c("RMSD","mAD","R")
    rownames(est_overall) <- names(estP)
  } else {print("make sure rownames are in the same order!")}
  est_overall
}


