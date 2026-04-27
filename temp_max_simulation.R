library(BrazilMet)
library(dplyr)
library(sf)
library(ggplot2)
library(geobr)
library(terra)
library(gstat)
library(ggspatial)
library(scales)
library(gstat)

crs_geo  <- 4326
crs_proj <- 5880
res_m   <- 10000   
set.seed(123)

br_proj  <- read_country(year = 2020) |> st_transform(crs_proj)
ufs_proj <- read_state(year = 2020)   |> st_transform(crs_proj)
bb <- st_bbox(br_proj)

stations_operating <- see_stations_info() |>
  filter(
    situation_operation == "operating",
    !is.na(longitude_degrees),
    !is.na(latitude_degrees)
  ) |>
  st_as_sf(coords = c("longitude_degrees", "latitude_degrees"), crs = crs_geo)

coords_geo <- st_coordinates(stations_operating)

stations_temp <- stations_operating |>
  mutate(
    lon = coords_geo[, 1],
    lat = coords_geo[, 2],
    abs_lat = abs(lat),
    
    efeito_altitude = -0.006 * altitude_m,
    litoral = if_else(lon > -41.5, 1, 0),
    
    dominio = case_when(
      uf %in% c("AM", "RR", "AP", "PA", "AC", "RO") ~ "amazonico",
      uf %in% c("CE", "RN", "PB", "PE", "AL", "SE") ~ "semiarido_leste",
      uf %in% c("PI", "BA") & abs_lat < 16 ~ "semiarido_norte_bahia_piaui",
      uf %in% c("MT", "MS", "GO", "TO", "DF") ~ "cerrado",
      uf %in% c("PR", "SC", "RS") ~ "subtropical",
      uf %in% c("ES", "RJ", "SP") & litoral == 1 ~ "litoraneo_sudeste",
      uf == "BA" & litoral == 1 ~ "litoraneo_bahia",
      uf == "MA" & abs_lat < 8 ~ "transicao_amazonia",
      TRUE ~ "tropical_interior"
    ),
    
    tmax_base = case_when(
      dominio == "amazonico" ~ 33.5,
      dominio == "transicao_amazonia" ~ 33.0,
      dominio == "semiarido_leste" ~ 35.5,
      dominio == "semiarido_norte_bahia_piaui" ~ 36.0,
      dominio == "cerrado" ~ 33.0,
      dominio == "subtropical" ~ 28.5,
      dominio == "litoraneo_sudeste" ~ 35.0,
      dominio == "litoraneo_bahia" ~ 30.5,
      dominio == "tropical_interior" ~ 31.5,
      TRUE ~ 31.0
    ),
    
    ajuste_latitude = case_when(
      abs_lat <= 5  ~  1.2,
      abs_lat <= 10 ~  0.8,
      abs_lat <= 15 ~  0.4,
      abs_lat <= 20 ~  0.0,
      abs_lat <= 25 ~ -0.8,
      abs_lat <= 30 ~ -1.8,
      TRUE          ~ -2.8
    ),
    
    ajuste_litoral = if_else(litoral == 1, -1.2, 0),
    
    ajuste_regional = case_when(
      uf %in% c("PI", "CE", "RN") & altitude_m < 300 ~  1.0,
      uf %in% c("PE", "PB", "AL", "SE") & altitude_m > 700 ~ -1.0,
      uf == "MG" & altitude_m > 900 ~ -1.5,
      uf %in% c("SP", "RJ", "ES") & altitude_m > 800 ~ -1.3,
      uf %in% c("SC", "RS", "PR") & altitude_m > 900 ~ -1.8,
      uf %in% c("MT", "MS") & altitude_m < 250 ~  0.8,
      uf %in% c("AM", "PA", "AP", "RR", "RO", "AC") & altitude_m < 150 ~ 0.6,
      TRUE ~ 0
    ),
    
    temperatura_maxima = tmax_base +
      ajuste_latitude +
      ajuste_litoral +
      ajuste_regional +
      efeito_altitude +
      rnorm(n(), 0, 0.7),
    
    temperatura_maxima = round(pmin(pmax(temperatura_maxima, 20), 40), 1)
  ) |>
  select(
    -lon, -lat, -abs_lat, -efeito_altitude, -litoral,
    -dominio, -tmax_base, -ajuste_latitude,
    -ajuste_litoral, -ajuste_regional
  )

stations_proj <- st_transform(stations_temp, crs_proj)

grade <- rast(
  xmin = bb["xmin"], xmax = bb["xmax"],
  ymin = bb["ymin"], ymax = bb["ymax"],
  resolution = res_m,
  crs = st_crs(br_proj)$wkt
)

grade_df <- as.data.frame(grade, xy = TRUE, cells = TRUE, na.rm = FALSE)

grade_sf <- st_as_sf(
  grade_df,
  coords = c("x", "y"),
  crs = crs_proj
)

idw_sf <- gstat::idw(
  formula = temperatura_maxima ~ 1,
  locations = stations_proj,
  newdata = grade_sf,
  idp = 2.0,
  maxdist = 510000
)

idw_coords <- st_coordinates(idw_sf)

idw_df <- idw_sf |>
  st_drop_geometry() |>
  mutate(
    x = idw_coords[, 1],
    y = idw_coords[, 2]
  ) |>
  select(x, y, temperatura_maxima = var1.pred)

r_temp <- rast(idw_df, type = "xyz")
names(r_temp) <- "temperatura_maxima"
crs(r_temp) <- st_crs(br_proj)$wkt

br_vect <- vect(br_proj)

r_temp <- r_temp |>
  crop(br_vect) |>
  mask(br_vect)

df_temp <- as.data.frame(r_temp, xy = TRUE, na.rm = TRUE)

bbox_sf <- st_as_sfc(bb)

st_crs(bbox_sf) <- st_crs(br_proj)

mask_out <- st_difference(bbox_sf, st_union(br_proj))

ggplot() +
  geom_raster(
    data = df_temp,
    aes(x = x, y = y, fill = temperatura_maxima)
  ) +
  geom_sf(
    data = mask_out,
    fill = "white",
    color = NA
  ) +
  geom_sf(
    data = ufs_proj,
    fill = NA,
    color = "gray",
    linewidth = 0.3,
    alpha = 0.95
  ) +
  geom_sf(
    data = br_proj,
    fill = NA,
    color = "black",
    linewidth = 1
  ) +
  
  scale_fill_viridis_c(
    name = "Temperatura máxima (°C)",
    breaks = pretty(df_temp$temperatura_maxima, n = 6),
    labels = label_number(
      accuracy = 0.1,
      decimal.mark = ",",
      big.mark = "."
    )
  ) +
  coord_sf(
    xlim = c(bb["xmin"], bb["xmax"]),
    ylim = c(bb["ymin"], bb["ymax"]),
    expand = FALSE,
    clip = "on"
  ) +
  annotation_scale(
    location = "bl",
    width_hint = 0.25,
    text_cex = 0.8,
    line_width = 0.7
  ) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_fancy_orienteering,
    height = unit(1.0, "cm"),
    width  = unit(1.0, "cm")
  ) +
  labs(
    title = "Temperatura máxima plausível interpolada",
    x = NULL,
    y = NULL
  ) +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 13),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    plot.margin = margin(8, 12, 8, 8)
  )

