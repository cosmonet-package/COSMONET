#' @title Split dataset in training and testing set
#'
#' @description This function splits data in training and testing set.
#' @param x input matrix \eqn{n1xp}.
#' @param q percentage of the sample size.
#'
#' @return Training and testing set.
#' @export
SplitData <- function(x, q){
  
  # q percentage of the sample size
  sample_size <- floor(q * nrow(x))
  # set the seed to make your partition reproducible
  set.seed(2020)
  ind_training <- sample(seq_len(nrow(x)), size = sample_size)
  training <- x[ind_training, ]
  testing <- x[-ind_training, ]
  
  return(list(training=training, testing=testing))
}