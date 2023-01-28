
#' Logarithmic transformation
#'
#' @param y The original value
#'
#' @return Log transformed value
#' @export
#'
#' @examples
#' ylog = slog(5)
slog <- function(y) {
  sy <- (y > 0) * y + (y <= 0) * 1
  return(log(sy))
}

#' Power transformation
#'
#' @param y Input value
#' @param pw Parameter of power transformation
#'
#' @return Power transformed value
#' @export
#'
#' @examples
#' ypow = spower(2,0.5)
spower <- function(y, pw) {
  sy <- (y > 0) * y + (y <= 0) * 1
  return(sy^pw)
}



#'  Yeo-Johnson transformation
#'
#' @description Yeo-Johnson transformation is used to transform a continuous (numeric) variable
#'              so that the resulting variable looks more normally distributed.
#'              See the details in the paper by
#'              Yeo, I. K., and Johnson, R. A. (2000). A new family of power transformations to improve normality or symmetry, Biometrika,87, 954--959.
#'
#' @param y numeric variable
#' @param theta parameter of the transformation
#'
#' @return Transformed value
#' @export
#'
#' @examples
#' ytrans <-YJtrans(5,0.2)
#' ytrans <-YJtrans(-5,0.5)
YJtrans <- function(y, theta) {
  sg <- y >= 0
  if (theta == 0) {
    temp <- slog(y + 1) * sg + (1 - sg) * (0.5 - 0.5 * (y - 1)^2)
  }
  if (theta == 2) {
    temp <- sg * (-0.5 + 0.5 * (y + 1)^2) - slog(-y + 1) * (1 - sg)
  }
  if ((theta != 0) & (theta != 2)) {
    temp <- sg * (spower(y + 1, theta) - 1)/theta + (1 - sg) * (1 - spower(-y + 1, 2 - theta))/(2 - theta)
  }
  return(temp)
}


#' Inverse of Yeo-Johnson transformation
#'
#' @param y Numeric value
#' @param theta parameter for the inverse transformation
#'
#' @return Inverse transformed value
#' @export
#' @description It computes the inverse of Yeo and Johnson (2000) transformation.
#'              See the paper by Yeo and Johnson (2000) for details.
#' @examples
#' Itrans = IYJtrans(-2,0.5)
IYJtrans <- function(y, theta) {
  sg <- y >= 0
  if (theta == 0) {
    temp <- (exp(y) - 1) * sg + (1 - sg) * (1 - spower(-2 * y + 1, 0.5))
  }
  if (theta == 2) {
    temp <- sg * (-1 + spower(2 * y + 1, 0.5)) + (1 - exp(-y)) * (1 - sg)
  }
  if ((theta != 0) & (theta != 2)) {
    temp <- sg * (spower(abs(theta) * y + 1, 1/theta) - 1) + (1 - sg) * (1 - spower(1 - (2 - theta) * y, 1/(2 - theta)))
  }
  return(temp)
}


#' Derivative of Yeo-Johnson transformation
#'
#' @param y Input numeric value
#' @param theta Parameter of YJ transformation
#'
#' @return derivative of YJ transform at specific value
#' @export
#'
#' @examples
#' dtrans <- DYJtrans(2,0.5)
DYJtrans <- function(y, theta) {
  sg <- y >= 0
  temp <- spower(y + 1, theta - 1) * sg + spower(-y + 1, 1 - theta) * (1 - sg)
  return(temp)
}
