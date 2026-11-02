#' Longitudes for drawing only.
#'
#' When samples straddle the 180th meridian (e.g. Commander Islands at +166
#' and the Aleutians at -170), plotting raw longitudes stretches the map
#' across the whole globe. Shifting the positive ones by -360 puts them next
#' to their neighbours (+166 becomes -194). Applied only when the longitudes
#' span more than 180 degrees, so projects that do not cross the date line
#' are drawn exactly as before. Coordinates used for analysis are untouched.
plot_longitude <- function(lon) {
  finite <- lon[is.finite(lon)]
  if (length(finite) == 0 || diff(range(finite)) <= 180) return(lon)
  ifelse(is.finite(lon) & lon > 0, lon - 360, lon)
}

#' Adds a pie radius column, pie_r, for geom_scatterpie.
#'
#' Without an explicit radius scatterpie uses 1/50 of the samples' longitude
#' span, so projects from one region get specks and a single site gets no
#' pies, on a map that always shows all of California. This takes 1/50 of the
#' drawn map's width instead: the samples plus the outline.
with_pie_r <- function(d, map) {
  d$pie_r <- diff(range(c(map$long, d$Longitude), na.rm = TRUE)) / 50
  d
}
