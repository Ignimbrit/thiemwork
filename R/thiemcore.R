#' Use Sichardt's equation to calculate the range of a well funnel cone
#'
#' @description \code{sichardt} needs the drawdown at the well and the
#' hydraulic conductiviry of the aquifer in which it is placed and will give
#' you a fair estimate of the dimension of the ensuing well funnel cone range.
#'
#' @param dh the drawdown of the well in m
#' @param kf the hydraulic conductivity of the aquifer in m/s
#'
#' @return the range in m
#'
#' @export
#'
#' @examples
#' # first example
#' sichardt(dh = 2, kf = 10^-4)
#'
#' # second example
#' aquifer_cond <- 2.5*10^-3
#' drawdown_steps <- seq(0.2, 5, 0.2)
#'
#' ranges <- vapply(
#' X = drawdown_steps,
#' FUN = sichardt,
#' FUN.VALUE = vector("double", length = 1),
#' kf = aquifer_cond)#' plot(ranges~drawdown_steps)
#'

sichardt <- function(dh, kf){
  3000*dh*sqrt(kf)
}

#' Use Kussakin's equation to calculate the range of a well funnel cone
#'
#' @description \code{kussakin} needs the drawdown at the well, the
#' hydraulic potential of the undisturbed gw-surface as meters above aquitard
#' and the hydraulic conductiviry of the aquifer in which it is placed-
#' It will in return give you a fair estimate of the dimension of the ensuing
#' well funnel cone range.
#'
#' @param dh the drawdown of the well in m
#' @param kf the hydraulic conductivity of the aquifer in m/s
#' @param hm the height of the water column above the quitard in m
#'
#' @return the range in m
#'
#' @export
#'
#' @examples
#' # first example
#' kussakin(dh = 2, kf = 10^-4, hm = 10)
#'
#' # second example
#' aquifer_cond <- 2.5*10^-3 # m/s
#' aquifer_height <- 8 # m
#' drawdown_steps <- seq(0.2, 5, 0.2) # m
#'
#' ranges <- vapply(
#' X = drawdown_steps,
#' FUN = kussakin,
#' FUN.VALUE = vector("double", length = 1),
#' kf = aquifer_cond, hm = aquifer_height)#' plot(ranges~drawdown_steps)
#'

kussakin <- function(dh, kf, hm){
  575*dh*sqrt((kf*hm))
}

#' Use the Dupuit-Thiem equation to calculate the pumping rate of a well
#'
#' @description \code{thiem_Q} is an implementation of the widely used
#' well equation by Dupuit-Thiem FOR UNCONFINED AQUIFERS.
#'
#' @param h1 the height of the water column above the aquitard at monitoring well 1 in m
#' @param h2 the height of the water column above the aquitard at monitoring well 2 in m
#' @param r1 the distance from the pumping well to the monitoring well No 1 in m
#' @param r2 the distance from the pumping well to the monitoring well No 2 in m
#' @param kf the hydraulic conductivity of the aquifer in m/s
#'
#' @return the pumping rate in m³/s
#'
#' @export
#'
#' @examples
#' # first example
#' thiem_Q(h1 = 10, h2 = 8, r1 = 50, r2 = 0.1, kf = 10^-3)
#'
#' # second example
#' drawdown_steps <- seq(0.2, 5, 0.2) # m
#' aquifer_cond <- 2.5*10^-3 # m/s
#' monitoring_well_potential = 10 # m water column above aquitard
#' monitoring_well_distance = 30 # m
#' pumping_well_radius = 0.08 # m == DN160
#'
#' rates <- vapply(
#' X = drawdown_steps,
#' FUN = function(h2, h1, r1, r2, kf){ # rearranging input order for vapply
#'   thiem_Q(h1 = h1, h2 = h2, r1 = r1, r2 = r2, kf = kf)
#' },
#' FUN.VALUE = vector("double", length = 1),
#' h1 = monitoring_well_potential,
#' r1 = monitoring_well_distance,
#' r2 = pumping_well_radius,
#' kf = aquifer_cond)
#'
#' plot(drawdown_steps~rates)
#'

thiem_Q <- function(h1, h2, r1, r2, kf){
  pi*kf*((h2^2-h1^2)/(log((r2/r1))))
}

#' Use the Dupuit-Thiem equation to calculate the water column level in a well
#'
#' @description \code{thiem_h} uses a heuristic to determine the height of
#' the water level above aquitard left in a well that is pumping at rate Q. To
#' get the often more useful drawdown of a well see \code{\link{thiem_h}}.
#'
#' @param Q the pumping rate of the well in m³/s
#' @param h0 the natural height of the water column above the aquitard (in absence of a well) in m
#' @param r_well the radius of the well
#' @param kf the hydraulic conductivity of the aquifer in m/s
#' @param dh_search_min the minimum possible drawdown of the well to be considered in m. Defaults to 0.001m == 1mm
#' @param dh_search_max the minimum possible drawdown of the well to be considered in m. Defaults to \code{h0}
#'
#' @return the water column level in the well in m above aquitard
#'
#' @export
#'
#' @examples
#' # first example
#' thiem_h(Q = 0.001, h0 = 10, r_well = 0.08, kf = 10^-4)
#'

thiem_h <- function(Q, h0, r_well, kf, dh_search_min = 0.001, dh_search_max = h0){
  
  speculative_Q <- function(dh, h0, r_well, kf){
    pi*kf*((h0^2-(h0-dh)^2)/(log(((sichardt(dh = dh, kf = kf))/r_well))))
  }
  
  q_speculation_quality <- function(dh, Q, h0, r_well, kf){
    Q_spec <- speculative_Q(dh = dh, h0 = h0, r_well = r_well, kf = kf)
    sqrt(((Q-Q_spec)^2))
  }
  
  ompitimization_res <- stats::optimize(
    q_speculation_quality,
    interval = c(dh_search_min, dh_search_max),
    Q = Q, h0 = h0, r_well = r_well, kf = kf
  )
  
  dh <- ompitimization_res$minimum
  
  h0-dh
}

#' Use the Dupuit-Thiem equation to calculate the drawdown in a well
#'
#' @description \code{thiem_dh} is a thin wrapper around
#' \code{\link{thiem_h}} that returns the drawdown of the well instead of
#' its water coloumn height, which is usually more useful.
#'
#' @param Q the pumping rate of the well in m³/s
#' @param h0 the natural height of the water column above the aquitard (in absence of a well) in m
#' @param r_well the radius of the well
#' @param kf the hydraulic conductivity of the aquifer in m/s
#' @param dh_search_min the minimum possible drawdown of the well to be considered in m. Defaults to 0.001m == 1mm
#' @param dh_search_max the minimum possible drawdown of the well to be considered in m. Defaults to \code{h0}
#'
#' @return the drawdown of the well in m
#'
#' @export
#'
#' @examples
#' # first example
#' thiem_dh(Q = 0.001, h0 = 10, r_well = 0.08, kf = 10^-4)
#'

thiem_dh <- function(Q, h0, r_well, kf, dh_search_min = 0.001, dh_search_max = h0){
  h <- thiem_h(Q = Q, h0 = h0, r_well = r_well, kf = kf, dh_search_min = dh_search_min, dh_search_max = dh_search_max)
  h0-h
}

#' Use the Dupuit-Thiem equation to calculate the cone shape of a well
#'
#' @description \code{thiem_coneshape} calculates the hydraulic potential
#' \code{h} at a given distance \code{x} from a pumping well
#'
#' @param x the distance in m from the well whre h is to be calculated
#' @param Q the pumping rate of the well in m³/s
#' @param h0 the natural height of the water column above the aquitard (in absence of a well) in m
#' @param r_well the radius of the well
#' @param kf the hydraulic conductivity of the aquifer in m/s
#' @param dh_search_min the minimum possible drawdown of the well to be considered in m. Defaults to 0.001m == 1mm
#' @param dh_search_max the minimum possible drawdown of the well to be considered in m. Defaults to \code{h0}
#'
#' @return hydraulic potential \code{h} as m above aquitard in distance \code{x} to the well
#'
#' @export
#'
#' @examples
#' # first example
#' x_range <- seq(0.5, 80, 0.5)
#' h_at_x <- vapply(
#' x_range,
#' thiem_coneshape,
#' vector("double", length = 1),
#' Q = 0.002, h0 = 10, r_well = 0.08, kf = 1*10^-4
#' )
#'
#' plot(h_at_x~x_range)
#'

thiem_coneshape <- function(x, Q, h0, r_well, kf, dh_search_min = 0.001, dh_search_max = h0){
  h <-  thiem_h(Q = Q, h0 = h0, r_well = r_well, kf = kf, dh_search_min = dh_search_min, dh_search_max = dh_search_max)
  dh <- h0-h
  R <- sichardt(dh = dh, kf = kf)
  
  if(x >= R){
    h0
  } else {
    sqrt(((Q*(2.3*(log10(x)-log10(r_well))))/(pi*kf))+h^2)
  }
}
