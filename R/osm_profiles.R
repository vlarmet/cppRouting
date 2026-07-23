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

#' Default pedestrian profile for OSM PBF import
#'
#' Highway tag values kept for foot routing, with a class id and a flat assumed walking pace (km/h,
#' not a legal speed limit -- \code{use_maxspeed} has no effect for this mode, see
#' \link{graph_from_osmpbf}). Excludes \code{motorway}/\code{trunk} (+ \code{_link}), which are not
#' walkable; includes \code{footway}/\code{path}/\code{pedestrian}/\code{steps} in addition to ordinary
#' roads, since pedestrians can generally use any street. \code{clazz} ids reuse \link{osm_car_profile}'s
#' values for highway types shared with that profile, so the same real-world road has the same class id
#' regardless of mode.
#'
#' @return A data.frame with columns \code{highway}, \code{clazz} and \code{speed}, suitable as the
#' \code{profile} element of \link{cpp_extract}'s \code{osm} argument for a \code{.pbf} source.
#' @export
osm_foot_profile <- function(){
  data.frame(
    highway = c("primary", "primary_link", "secondary", "secondary_link",
                "tertiary", "tertiary_link", "residential", "road", "unclassified",
                "living_street", "service", "track",
                "footway", "path", "pedestrian", "steps"),
    clazz   = c(15, 16, 21, 22, 31, 32, 41, 42, 43, 44, 45, 46, 51, 52, 53, 54),
    speed   = rep(5, 16),
    stringsAsFactors = FALSE
  )
}

#' Default bicycle profile for OSM PBF import
#'
#' Highway tag values kept for bicycle routing, with a class id and a flat assumed cycling pace (km/h,
#' not a legal speed limit -- \code{use_maxspeed} has no effect for this mode, see
#' \link{graph_from_osmpbf}). Excludes \code{motorway}/\code{trunk} (+ \code{_link}) and pedestrian-only
#' ways (\code{footway}/\code{pedestrian}/\code{steps}); includes \code{cycleway}/\code{path}/\code{track}
#' in addition to ordinary roads. \code{clazz} ids reuse \link{osm_car_profile}'s values for highway types
#' shared with that profile.
#'
#' @return A data.frame with columns \code{highway}, \code{clazz} and \code{speed}, suitable as the
#' \code{profile} element of \link{cpp_extract}'s \code{osm} argument for a \code{.pbf} source.
#' @export
osm_bike_profile <- function(){
  data.frame(
    highway = c("primary", "primary_link", "secondary", "secondary_link",
                "tertiary", "tertiary_link", "residential", "road", "unclassified",
                "living_street", "service", "track", "path", "cycleway"),
    clazz   = c(15, 16, 21, 22, 31, 32, 41, 42, 43, 44, 45, 46, 52, 61),
    speed   = rep(15, 14),
    stringsAsFactors = FALSE
  )
}

#' Default heavy goods vehicle (truck) profile for OSM PBF import
#'
#' Identical highway/clazz/speed table to \link{osm_car_profile} (trucks use the same road classes and,
#' as a simplification, the same legal speed limits in this version) -- the actual difference from
#' \code{mode = "car"} is in access-tag resolution (\code{hgv=no}/\code{hgv=destination}/\code{hgv=yes}),
#' not the set of roads considered. Weight/height/width restrictions (\code{maxweight}, \code{maxheight},
#' ...) are not resolved against a vehicle's actual dimensions in this version.
#'
#' @return A data.frame with columns \code{highway}, \code{clazz} and \code{speed}, suitable as the
#' \code{profile} element of \link{cpp_extract}'s \code{osm} argument for a \code{.pbf} source.
#' @export
osm_hgv_profile <- function(){
  osm_car_profile()
}
