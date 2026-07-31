#' Basemap layers sized to a set of sample coordinates.
#'
#' The workflow is not tied to any one region, so the map extent is derived from
#' the coordinates themselves rather than hardcoded. A world outline is always
#' drawn; US state boundaries are added when the extent overlaps the continental
#' US, since that is where finer internal borders help most and `maps` ships
#' them for free.
#'
#' @param coords data frame with numeric Longitude and Latitude columns
#' @param margin_frac fraction of the coordinate span to pad the extent by
#' @return list(layers, xlim, ylim, r) where `layers` goes straight into a
#'   ggplot and `r` is a sensible scatterpie radius in degrees for this extent.
make_basemap <- function(coords, margin_frac = 0.12) {
  lon <- coords$Longitude[is.finite(coords$Longitude)]
  lat <- coords$Latitude[is.finite(coords$Latitude)]

  if (length(lon) == 0 || length(lat) == 0) {
    stop("make_basemap: no finite coordinates to draw")
  }

  xr <- range(lon)
  yr <- range(lat)

  # A single locality (or a tight cluster) has ~zero span, so enforce a floor of
  # half a degree of padding to keep the map from collapsing to a point.
  pad_x <- max(diff(xr) * margin_frac, 0.5)
  pad_y <- max(diff(yr) * margin_frac, 0.5)

  xlim <- c(xr[1] - pad_x, xr[2] + pad_x)
  ylim <- c(yr[1] - pad_y, yr[2] + pad_y)

  world <- ggplot2::map_data("world")
  layers <- list(
    ggplot2::geom_polygon(
      data = world,
      ggplot2::aes(x = long, y = lat, group = group),
      fill = "grey97", color = "grey60", linewidth = 0.2
    )
  )

  overlaps_us <- xlim[1] < -66 && xlim[2] > -125 && ylim[1] < 50 && ylim[2] > 24
  if (overlaps_us) {
    layers <- c(layers, list(
      ggplot2::geom_polygon(
        data = ggplot2::map_data("state"),
        ggplot2::aes(x = long, y = lat, group = group),
        fill = NA, color = "grey60", linewidth = 0.2
      )
    ))
  }

  layers <- c(layers, list(
    ggplot2::coord_quickmap(xlim = xlim, ylim = ylim, expand = FALSE)
  ))

  list(
    layers = layers,
    xlim = xlim,
    ylim = ylim,
    # Pie radius in data units, scaled to the extent so the same code works for
    # a 2-degree island and a continent-wide sample.
    r = 0.025 * max(diff(xlim), diff(ylim))
  )
}

#' Discrete fill palette for ancestry components.
#'
#' Okabe-Ito, which is colourblind-safe, recycled if K exceeds its length.
ancestry_palette <- function(n) {
  base <- c(
    "#0072B2", "#E69F00", "#009E73", "#CC79A7",
    "#56B4E9", "#D55E00", "#F0E442", "#999999"
  )
  if (n <= length(base)) return(base[seq_len(n)])
  rep(base, length.out = n)
}
