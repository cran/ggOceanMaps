## ---- include = FALSE---------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(collapse = TRUE, message = FALSE, warning = FALSE, comment = "#>")


## ----setup--------------------------------------------------------------------
library(ggOceanMaps)

## -----------------------------------------------------------------------------
etopoPath <- ""
NEDPath <- ""
outPath <- ""

## ----eval = FALSE-------------------------------------------------------------
#  rb <- raster_bathymetry(bathy = paste(etopoPath, "ETOPO1_Ice_g_gmt4.grd", sep = "/"),
#                          depths = c(50, 300, 500, 1000, 1500, 2000, 4000, 6000, 10000),
#                          proj.out = "+init=epsg:4326",
#                          boundary = c(-180.0083, 180.0083, -90, 90),
#                          aggregation.factor = 6
#  )
#  
#  dd_bathy <- vector_bathymetry(rb, drop.crumbs = 50, remove.holes = 50)
#  
#  save(dd_bathy, file = paste(outPath, "dd_bathy.rda", sep = "/"), compress = "xz")

## ----eval = FALSE-------------------------------------------------------------
#  world <- rgdal::readOGR(paste(NEDPath, "ne_10m_land/ne_10m_land.shp", sep = "/"))
#  islands <- rgdal::readOGR(paste(NEDPath, "ne_10m_minor_islands/ne_10m_minor_islands.shp", sep = "/"))
#  world <- rbind(world, islands)
#  
#  dd_land <- clip_shapefile(world, c(-180, 180, -90, 90))
#  
#  save(dd_land, file = paste(outPath, "dd_land.rda", sep = "/"), compress = "xz")

## ----eval = FALSE-------------------------------------------------------------
#  glaciers <- rgdal::readOGR(paste(NEDPath, "ne_10m_glaciated_areas/ne_10m_glaciated_areas.shp", sep = "/"))
#  iceshelves <- rgdal::readOGR(paste(NEDPath, "ne_10m_antarctic_ice_shelves_polys/ne_10m_antarctic_ice_shelves_polys.shp", sep = "/"))
#  
#  glaciers <- rbind(glaciers, iceshelves)
#  glaciers <- rgeos::gBuffer(glaciers, byid = TRUE, width = 0)
#  
#  dd_glacier <- clip_shapefile(glaciers, c(-180, 180, -90, 90))
#  dd_glacier <- rgeos::gBuffer(dd_glacier, byid = FALSE, width = 0.1)
#  dd_glacier <- rgeos::gBuffer(dd_glacier, byid = FALSE, width = -0.1)
#  
#  save(dd_glacier, file = paste(outPath, "dd_glacier.rda", sep = "/"), compress = "xz")

## ----eval = FALSE-------------------------------------------------------------
#  rb <- raster_bathymetry(bathy = paste(etopoPath, "ETOPO1_Ice_g_gmt4.grd", sep = "/"),
#                          depths = c(50, 300, 500, 1000, 1500, 2000, 4000, 6000, 10000),
#                          proj.out = "+init=epsg:3995",
#                          boundary = c(-180.0083, 180.0083, 30, 90),
#                          aggregation.factor = 2
#  )
#  
#  arctic_bathy <- vector_bathymetry(rb)
#  
#  save(arctic_bathy, file = paste(outPath, "arctic_bathy.rda", sep = "/"), compress = "xz")

## ----eval = FALSE-------------------------------------------------------------
#  arctic_land <- clip_shapefile(world, c(-180, 180, 30, 90))
#  arctic_land <- sp::spTransform(arctic_land, sp::CRS(sp::proj4string(arctic_bathy)))
#  arctic_land <- rgeos::gBuffer(arctic_land, byid = TRUE, width = 0)
#  
#  save(arctic_land, file = paste(outPath, "arctic_land.rda", sep = "/"), compress = "xz")

## ----eval = FALSE-------------------------------------------------------------
#  arctic_glacier <- clip_shapefile(glaciers, c(-180, 180, 40, 90))
#  arctic_glacier <- sp::spTransform(arctic_glacier, sp::CRS(sp::proj4string(arctic_bathy)))
#  arctic_glacier <- rgeos::gBuffer(arctic_glacier, byid = FALSE, width = 1000)
#  arctic_glacier <- rgeos::gBuffer(arctic_glacier, byid = FALSE, width = -1000)
#  
#  save(arctic_glacier, file = paste(outPath, "arctic_glacier.rda", sep = "/"), compress = "xz")

## ----eval = FALSE-------------------------------------------------------------
#  rb <- raster_bathymetry(bathy = paste(etopoPath, "ETOPO1_Ice_g_gmt4.grd", sep = "/"),
#                          depths = c(50, 300, 500, 1000, 1500, 2000, 4000, 6000, 10000),
#                          proj.out = "+init=epsg:3031",
#                          boundary = c(-180.0083, 180.0083, -80, -30),
#                          aggregation.factor = 2
#  )
#  
#  antarctic_bathy <- vector_bathymetry(rb)
#  
#  save(antarctic_bathy, file = paste(outPath, "antarctic_bathy.rda", sep = "/"), compress = "xz")

## ----eval = FALSE-------------------------------------------------------------
#  antarctic_land <- clip_shapefile(world, c(-180, 180, -90, -30))
#  antarctic_land <- sp::spTransform(antarctic_land, sp::CRS(sp::proj4string(antarctic_bathy)))
#  antarctic_land <- rgeos::gBuffer(antarctic_land, byid = TRUE, width = 0)
#  
#  save(antarctic_land, file = paste(outPath, "antarctic_land.rda", sep = "/"), compress = "xz")

## ----eval = FALSE-------------------------------------------------------------
#  antarctic_glacier <- clip_shapefile(glaciers, c(-180, 180, -90, -30))
#  antarctic_glacier <- sp::spTransform(antarctic_glacier, sp::CRS(sp::proj4string(antarctic_bathy)))
#  antarctic_glacier <- rgeos::gBuffer(antarctic_glacier, byid = FALSE, width = 1000)
#  antarctic_glacier <- rgeos::gBuffer(antarctic_glacier, byid = FALSE, width = -1000)
#  
#  save(antarctic_glacier, file = paste(outPath, "antarctic_glacier.rda", sep = "/"), compress = "xz")

