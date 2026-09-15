get_median <- function(data, indices) {
  return(median(data[indices],na.mr=TRUE))
}