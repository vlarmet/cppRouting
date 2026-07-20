#' Default car profile for OSM PBF import
#'
#' Highway tag values kept for car routing, with a class id and a default speed (km/h) used unless
#' overridden by a \code{maxspeed} tag. Mirrors a typical osm2po car profile (\code{wtr.finalMask = car}):
#' motorway/trunk/primary/secondary (+ their \code{_link} variants), tertiary, residential, road and
#' unclassified. \code{service}, \code{living_street}, \code{track}, \code{path}, cycleways and footways are
#' intentionally excluded, as in that reference configuration.
#'
#' @return A data.frame with columns \code{highway}, \code{clazz} and \code{speed}, suitable as the
#' \code{profile} argument of \link{cpp_extract} for a \code{.pbf} source.
#' @export
osm_car_profile <- function(){
  data.frame(
    highway = c("motorway", "motorway_link", "trunk", "trunk_link",
                "primary", "primary_link", "secondary", "secondary_link",
                "tertiary", "tertiary_link", "residential", "road", "unclassified"),
    clazz   = c(11, 12, 13, 14, 15, 16, 21, 22, 31, 32, 41, 42, 43),
    speed   = c(120, 30, 90, 30, 70, 30, 60, 30, 40, 20, 40, 50, 50),
    stringsAsFactors = FALSE
  )
}
