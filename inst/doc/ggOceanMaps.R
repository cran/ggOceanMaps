## ---- include = FALSE---------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(collapse = TRUE, message = FALSE, warning = FALSE, comment = "#>")

if (!requireNamespace("ggOceanMapsData", quietly = TRUE)) {
utils::install.packages("ggOceanMapsData", repos = c("https://mikkovihtakari.github.io/drat", "https://cloud.r-project.org"))
}

## ----eval = FALSE-------------------------------------------------------------
#  install.packages("ggOceanMapsData", repos = c("https://mikkovihtakari.github.io/drat", "https://cloud.r-project.org"))

## ----setup--------------------------------------------------------------------
library(ggOceanMapsData)
library(ggOceanMaps)

## -----------------------------------------------------------------------------
library(ggOceanMaps)
basemap(limits = 60)

## -----------------------------------------------------------------------------
basemap(limits = c(-20, 20, 40, 59))

## ---- echo = FALSE, fig.width=10----------------------------------------------
start.lat <- data.frame(lon = seq(0, 360, 1), lat = 60)
end.lat <- data.frame(lon = seq(0, 360, 1), lat = 80)

start.lat.lon <- data.frame(lon = c(seq(120, 180, 1), seq(-180, -120, 1)), lat = 60)
end.lat.lon <- data.frame(lon = c(seq(120, 180, 1), seq(-180, -120, 1)), lat = 80)
lon.segments <- data.frame(lon = c(rep(-120, 2), rep(120, 2)), lat = rep(c(60, 80), 2))

real.lims <- auto_limits(expand.grid(lon = c(-120, 180, 120), lat = c(60, 60, 80)))
limits <- real.lims$projLimits

dt <- data.frame(lon = c(real.lims$projLimits[1:2]), lat = real.lims$projLimits[3:4])

labels.lat <- data.frame(lon = c(60, 60), lat = c(60, 80), label = c("limits[3]\n(Min. latitude)", "limits[4]\n(Max. latitude)"))
labels.lon <- data.frame(lon = c(120, -120), lat = c(rep(60, 2), rep(80, 2)), label = rep(c("limits[1]\n(Start longitude)", "limits[2]\n(End longitude)")))

p1 <- ggplot() +
  geom_spatial_path(data = start.lat, aes(x = lon, y = lat), crs = 4326, linetype = 2) +
  geom_spatial_path(data = end.lat, aes(x = lon, y = lat), crs = 4326, linetype = 2) +
  geom_spatial_path(data = start.lat.lon, aes(x = lon, y = lat), crs = 4326, color = "red") +
  geom_spatial_path(data = end.lat.lon, aes(x = lon, y = lat), crs = 4326, color = "red") +
  geom_spatial_path(data = lon.segments, aes(x = lon, y = lat, group = lon), crs = 4326, color = "red") +
  layer_spatial(data = real.lims$projBound, fill = NA, color = "red", size = 2) +
  scale_x_continuous("Longitude (decimal degrees)", breaks = seq(0, 360, 30)) +
  geom_label(data = data.frame(lon = 0, lat = 0, label = "The pole"), aes(x = lon, y = lat, label = label), size = FS(8)) +
  geom_spatial_label(data = labels.lat, aes(x = lon, y = lat, label = label), crs = 4326, size = FS(8)) +
  geom_spatial_label(data = labels.lon, aes(x = lon, y = lat, label = label), crs = 4326, color = "red", size = FS(8)) +
  ylab("Longitude (decimal degrees)") +
  coord_sf(crs = 3995, label_axes = "EEEE") +
  theme_map(grid.col = "grey70", grid.size = 0.1, base_size = 8)

labels.proj <- data.frame(lon = c(limits[1], limits[2], mean(limits[1:2]), mean(limits[1:2])), lat = c(mean(limits[3:4]), mean(limits[3:4]), limits[3], limits[4]), label = c("limits[1]\n(Min. longitude)", "limits[2]\n(Max. longitude)", "limits[3]\n(Min. latitude)", "limits[4]\n(Max. latitude)"))

p2 <- ggplot() +
  geom_vline(xintercept = real.lims$projLimits[1:2], linetype = 2) +
  geom_hline(yintercept = real.lims$projLimits[3:4], linetype = 2) +
  layer_spatial(data = real.lims$projBound, fill = NA, color = "red", size = 2) +
  geom_label(data = labels.proj, aes(x = lon, y = lat, label = label), color = "red", size = FS(8)) +
  geom_label(data = data.frame(lon = 0, lat = 0, label = "The pole"), aes(x = lon, y = lat, label = label), size = FS(8)) +
  geom_spatial_path(data = start.lat, aes(x = lon, y = lat), crs = 4326, color = NA) +
  labs(y = "Latitude (meters from the pole)", x = "Longitude (meters from the pole)") +
  coord_sf(crs = 3995, datum = 3995) +
  theme_map(grid.col = "grey70", grid.size = 0.1, base_size = 8) + 
  theme(plot.margin = margin(t = 20, unit = "pt"))

cowplot::plot_grid(p1, p2, labels = "AUTO", align = "hv")

## -----------------------------------------------------------------------------
dt <- data.frame(lon = c(160, 160, -160, -160), lat = c(60, 80, 80, 60))

basemap(limits = c(160, -160, 60, 80)) +
  geom_spatial_polygon(data = dt, aes(x = lon, y = lat), fill = NA, color = "red")

## -----------------------------------------------------------------------------
basemap(limits = 60, projection.grid = TRUE, grid.col = "red")

## -----------------------------------------------------------------------------
basemap(limits = c(-2e6, 1e6, 0, 3e6), shapefiles = "Arctic") 

## -----------------------------------------------------------------------------
dt <- expand.grid(lon = c(160, -160), lat = c(60, 80))

basemap(data = dt) +
  geom_spatial_point(data = dt, aes(x = lon, y = lat), color = "red")

## -----------------------------------------------------------------------------
basemap(limits = c(100, 160, -20, 30), bathymetry = TRUE)

## -----------------------------------------------------------------------------
basemap(limits = 60, glaciers = TRUE, bathymetry = TRUE)

## -----------------------------------------------------------------------------
dt <- data.frame(lon = c(seq(-180, 0, 30), seq(30, 180, 30)), lat = -70)
basemap(limits = -60, glaciers = TRUE) + geom_spatial_point(data = dt, aes(x = lon, y = lat), color = "red")

## -----------------------------------------------------------------------------
basemap(limits = -60, glaciers = TRUE) + 
  geom_point(data = transform_coord(dt), aes(x = lon, y = lat), color = "red")

## -----------------------------------------------------------------------------
dt <- data.frame(lon = c(-100, -80, -60), lat = c(10, 25, 40), var = c("a", "a", "b"))
basemap(data = dt) + geom_point(data = dt, aes(x = lon, y = lat), color = "red")

## -----------------------------------------------------------------------------
transform_coord(data.frame(lon = -80, lat = 25), bind = TRUE)

## -----------------------------------------------------------------------------
basemap(limits = c(-160, -80, 60, 85), rotate = TRUE)

## -----------------------------------------------------------------------------
qmap(dt, color = I("red")) # set color
qmap(dt, color = var) # map color

## -----------------------------------------------------------------------------
basemap(limits = c(0, 46, 70, 81), bathymetry = TRUE, bathy.style = "poly_greys")
basemap(limits = c(0, 46, 70, 81), bathymetry = TRUE, bathy.style = "contour_blues")
basemap(limits = c(0, 46, 70, 81), bathymetry = TRUE, bathy.style = "contour_grey")

## -----------------------------------------------------------------------------
basemap(limits = c(-140, -105, 20, 40), bathymetry = TRUE) + scale_fill_viridis_d("Water depth (m)")

## -----------------------------------------------------------------------------
basemap(limits = c(0, 60, 68, 80), bathymetry = TRUE, bathy.style = "contour_blues") + scale_color_hue()

## -----------------------------------------------------------------------------
basemap(limits = c(-20, 30, 55, 70), glaciers = TRUE, 
        bathymetry = TRUE, bathy.style = "poly_greys",
        land.col = "#eeeac4", gla.col = "cadetblue", 
        land.border.col = NA, gla.border.col = NA,
        grid.size = 0.05)

## -----------------------------------------------------------------------------
basemap(limits = c(124, 148, 31, 50), grid.col = NA) + labs(x = NULL, y = "Only latitude for you, ...")

## -----------------------------------------------------------------------------
basemap(limits = c(-75, -45, 62, 78), rotate = TRUE) + 
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "tr", which_north = "true") 

## -----------------------------------------------------------------------------
p <- basemap(-60)
attributes(p)

## -----------------------------------------------------------------------------
data(fishingAreasNor, package = "ggOceanMapsData")

basemap(limits = raster::extent(fishingAreasNor)[1:4]) + 
  annotation_spatial(fishingAreasNor, fill = NA) + 
  coord_sf(expand = FALSE)

## -----------------------------------------------------------------------------
labels <- sp::SpatialPointsDataFrame(rgeos::gCentroid(fishingAreasNor, byid=TRUE), 
                                     data = fishingAreasNor@data)
labels <- df_spatial(labels)

# To correct for a bug in ggspatial 1.1.1
if("x_without_geom[df_geom$feature_id, ]" %in% names(labels)) {
  names(labels)[names(labels) == "x_without_geom[df_geom$feature_id, ]"] <- "FID"
}
# Bug correction end.

fishingAreasNor@data$area <- raster::area(fishingAreasNor)/1e9 # calculate area in 1000 km2

p <- basemap(limits = raster::extent(fishingAreasNor)[1:4]) +
  annotation_spatial(fishingAreasNor, aes(fill = area)) +
  geom_spatial_text(data = labels, aes(x = x, y = y, label = FID), size = FS(8), fontface = 2) +
  scale_fill_distiller(name = "Area\n(1000 km2)",
                       palette = "Spectral", na.value = "white",
                       limits = c(0, 500), oob = scales::squish)

reorder_layers(p)

## -----------------------------------------------------------------------------
p <- reorder_layers(p)
tmp <- sapply(p$layers, function(k) !is.null(k$mapping$label)) # the layer with label mapping
p$layers <- c(p$layers[-which(tmp)], p$layers[which(tmp)])
p

## ----eval = FALSE-------------------------------------------------------------
#  etopoPath <- "" # Replace by the path to the folder where the ETOPO1 grd file is located.
#  lims <- c(-8, 65, 68, 82)
#  projection <- "+init=epsg:32636"
#  basemap(limits = lims)

## ----eval = FALSE-------------------------------------------------------------
#  rb <- raster_bathymetry(bathy = paste(etopoPath, "ETOPO1_Ice_g_gmt4.grd", sep = "/"),
#                          depths = c(50, 100, 200, 300, 500, 1000, 1500, 2000, 4000, 6000, 10000),
#                          proj.out = projection,
#                          boundary = lims
#  )

## ----eval = FALSE-------------------------------------------------------------
#  class(rb)
#  names(rb)
#  raster::plot(rb$raster)

## ----eval = FALSE-------------------------------------------------------------
#  bs_bathy <- vector_bathymetry(rb)
#  sp::plot(bs_bathy)

## ----eval = FALSE-------------------------------------------------------------
#  NEDPath <- "" # Natural Earth Data location
#  outPath <- "" # Data output location

## ----eval = FALSE-------------------------------------------------------------
#  world <- rgdal::readOGR(paste(NEDPath, "ne_10m_land/ne_10m_land.shp", sep = "/"))
#  islands <- rgdal::readOGR(paste(NEDPath, "ne_10m_minor_islands/ne_10m_minor_islands.shp", sep = "/"))
#  world <- rbind(world, islands)
#  
#  bs_land <- clip_shapefile(world, lims)
#  bs_land <- sp::spTransform(bs_land, CRSobj = sp::CRS(projection))
#  rgeos::gIsValid(bs_land) # Has to return TRUE, if not use rgeos::gBuffer
#  bs_land <- rgeos::gBuffer(bs_land, byid = TRUE, width = 0)
#  sp::plot(bs_land)

## ----eval = FALSE-------------------------------------------------------------
#  glaciers <- rgdal::readOGR(paste(NEDPath, "ne_10m_glaciated_areas/ne_10m_glaciated_areas.shp", sep = "/"))
#  rgeos::gIsValid(glaciers) # Needs buffering
#  glaciers <- rgeos::gBuffer(glaciers, byid = TRUE, width = 0)
#  
#  bs_glacier <- clip_shapefile(glaciers, lims)
#  bs_glacier <- sp::spTransform(bs_glacier, CRSobj = sp::CRS(projection))
#  rgeos::gIsValid(bs_glacier)
#  sp::plot(bs_glacier)

## ----eval = FALSE-------------------------------------------------------------
#  save(bs_bathy, bs_land, bs_glacier, file = paste(outPath, "bs_shapes.rda", sep = "/"), compress = "xz")

## ----echo = FALSE-------------------------------------------------------------
data(bs_shapes, package = "ggOceanMapsData")

## -----------------------------------------------------------------------------
basemap(shapefiles = list(land = bs_land, glacier = bs_glacier, bathy = bs_bathy), bathymetry = TRUE, glaciers = TRUE)

## -----------------------------------------------------------------------------
basemap(limits = c(10, 53, 70, 80), shapefiles = list(land = bs_land, glacier = bs_glacier, bathy = bs_bathy), bathymetry = TRUE, glaciers = TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  basemap(limits = c(160, -160, 0, 30)) # not evaluated

## -----------------------------------------------------------------------------
basemap(limits = c(160, -160, 30, 60), rotate = TRUE)

## -----------------------------------------------------------------------------
basemap(limits = c(-180, 180, -70, -60))

## ----eval = FALSE-------------------------------------------------------------
#  basemap(limits = raster::extent(fishingAreasNor)[1:4], bathymetry = TRUE) +
#    annotation_spatial(fishingAreasNor, aes(fill = area))
#  #> Error: Continuous value supplied to discrete scale

## ----eval = FALSE-------------------------------------------------------------
#  basemap(limits = raster::extent(fishingAreasNor)[1:4], bathymetry = TRUE,
#          bathy.style = "contour_blues", legends = FALSE) +
#    annotation_spatial(fishingAreasNor, aes(fill = area), alpha = 0.4) +
#    coord_sf(expand = FALSE)

## -----------------------------------------------------------------------------
citation("ggOceanMaps")

