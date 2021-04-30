######################################
#
# Reduce resolution of CMAP data for 
# 
#  Written by M. Fidino 12/19/2019 mdy
#
######################################

library(raster)
library(dplyr)
library(sf)
library(fasterize)

# read in the chicago data layer
my_raster_path <- 
	"D:/GIS/cmap/landcover_2010_chicagoregion.img"

# read it in
my_map <- raster::raster(
	my_raster_path
)

# crop it down to Chicago
city_outline <- sf::st_read(
	"D:/GIS/IL City Shape Files", 
	layer = "tl_2017_17_place"
) %>% 
	select("NAME")

# get Chicago
chicago <- city_outline[city_outline$NAME == "Chicago",]

# Transform to same projection as the raster
chicago <- sf::st_transform(
	chicago,
	crs = sf::st_crs(
		my_map
	)
)

# crop the Chicagoland region to just Chicago
chicago_map <- raster::crop(
	my_map,
	extent(
		chicago
	)
)

# The crop only get's us to the rectangular extent of the chicago
#   polygon. We want to NA the values outside of the Chicago polygon
#   Adding a buffer of 1000m to account for sites or conflicts
#   on the edge.
chicago_buffer <- sf::st_buffer(
	chicago,
	1000
)

chicago_map <- raster::mask(
	chicago_map,
	chicago_buffer
)

# Aggregate the land-use classes. Using a custom function
#   to make a table of each land-use class, calculate the
#   proportion, and then return the first value (which is canopy).
canopy <- raster::aggregate(chicago_map,
														500, 
														function(x,...){
										          tmp <- tabulate(x, 7)
	                            tmp <- tmp / sum(tmp)
	                            return(tmp[1])
	                          }
)

# do the same for grass cover, lulc category 2.
grass  <- raster::aggregate(chicago_map,
														500, 
														function(x,...){
															tmp <- tabulate(x, 7)
															tmp <- tmp / sum(tmp)
															return(tmp[2])
														}
)

# imperv is 3 classes together (5,6 and 7)
imperv <- raster::aggregate(chicago_map,
														500, 
														function(x,...){
															tmp <- tabulate(x, 7)
															tmp <- tmp / sum(tmp)
															return(sum(tmp[5:7]))
														}
)

# convert the three layers to a raster stack
chicago_stack <- do.call(
	raster::stack,
	list(
		canopy,
		grass,
		imperv
	)
)

# add column names
colnames(values(chicago_stack)) <- c(
	"canopy",
	"grass",
	"imperv"
)

# mask the edges again
chicago_stack <- raster::mask(
	chicago_stack,
	chicago
)

# save as RDS
saveRDS(
	chicago_stack,
	"./data/high_res_reduction.rds"
)

# save it as a multi-layer raster as well
raster::writeRaster(
	chicago_stack,
	filename="chicago_highres_500res.tif",
	options="INTERLEAVE=BAND",
	overwrite=TRUE
)

#################################

# do the same thing but with the housing density data
pop_data <- sf::st_read(
	"D:/GIS/housing_density", 
	layer = "il_blk10_Census_change_1990_2010_PLA2"
) %>% 
	dplyr::select("HU10")

# transform to raster crs

pop_data <- sf::st_transform(
	pop_data,
	sf::st_crs(
		my_map
	)
)

# fix any potential issues before cropping
pop_data <- sf::st_make_valid(
	pop_data
)

# crop it down to the chicago area
pop_data <- sf::st_crop(
	pop_data,
	raster::extent(
		chicago_stack
	)
) %>% 
	sf::st_buffer(.,0)

# make a temporary raster from the one we just generated
tmp <- dropLayer(
	chicago_stack,
	2:3
) %>% as(
	'SpatialPolygonsDataFrame'
	) %>% st_as_sf

# add a cell column to summarise by
tmp$cell <- stringr::str_pad(
	as.character(
		1:nrow(tmp)
	),
	width = 5,
	pad = "0"
)

# determine intersection
#   group by cell
#   sum the values
#   make it a POLYGON so we can rasterize it
hou_intersect <- sf::st_intersection(
	tmp,
	pop_data
) %>% 
	dplyr::group_by(
		cell
	) %>% 
	dplyr::summarise_at(
		.vars = "HU10",
		.funs = sum
	) %>% 
	sf::st_cast(
		"POLYGON"
	)

# make it a raster
hu10_raster <- fasterize::fasterize(
	hou_intersect,
	dropLayer(
		chicago_stack,
		2:3
	),
	"HU10",
	fun = "last"
)

# save the rds
saveRDS(
	hu10_raster,
	"./data/hu10_raster.rds"
)

# Bring in income 
income <- raster::raster(
	"./data/Final_Income_Raster.tif"
)

income <- raster::projectRaster(
	income, 
	chicago_stack
)

# mask values outside of Chicago
income <- raster::mask(
	income,
	chicago
)


# bring in vacancy
vacancy <- raster::raster(
	"./data/Vacancy_311.tif"
)

vacancy <- raster::projectRaster(
	vacancy, 
	chicago_stack
)

vacancy <- raster::mask(
	vacancy,
	chicago
)

# and distance to a natural water source
water <- raster::raster(
	"./data/Final_WaterDist_Raster.tif"
)

water <- projectRaster(
	water,
	chicago_stack
)


water <- raster::mask(
	water,
	chicago
)


# add it as another layer to the other rasters to make a final object
all_raw_layers <- raster::stack(
	chicago_stack,
	hu10_raster,
	income,
	vacancy,
	water
)

colnames(values(all_raw_layers))[4:7] <- c(
	"houses",
	"income",
	"vacancy",
	"dist2water"
)


saveRDS(
	all_raw_layers,
	"./data/all_raw_layers.rds"
)

raster::writeRaster(
	all_raw_layers,
	filename="chicago_variables_raster_500.tif",
	options="INTERLEAVE=BAND",
	overwrite=TRUE
)
