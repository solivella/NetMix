cLL <- function(x, select_States = FALSE, select_Blocks = FALSE){
    LB <- x$LowerBound
    if(select_States){
      LB <- LB - lfactorial(x$n_states)
      
    }
    if(select_Blocks){
      LB <- LB - lfactorial(x$n_blocks)
    }
    return(LB)
}
