library(shiny)
library(DT)
library(BrazilMet)
library(dplyr)
library(sf)
library(geobr)
library(terra)
library(gstat)
library(leaflet)
library(scales)
library(htmltools)
library(ggplot2)
library(plotly)

crs_geo  <- 4326
crs_proj <- 5880
res_m    <- 10000

br_proj  <- read_country(year = 2020) |> st_transform(crs_proj)
ufs_proj <- read_state(year = 2020) |> st_transform(crs_proj)
bb <- st_bbox(br_proj)
br_vect <- vect(br_proj)

gerar_temperatura_plausivel <- function(df, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  coords <- st_coordinates(df)
  
  df |>
    mutate(
      lon = coords[, 1],
      lat = coords[, 2],
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
      temperatura_maxima = round(
        pmin(
          pmax(
            tmax_base + ajuste_latitude + ajuste_litoral +
              ajuste_regional + efeito_altitude + rnorm(n(), 0, 0.7),
            20
          ),
          40
        ),
        1
      )
    ) |>
    select(
      -lon, -lat, -abs_lat, -efeito_altitude, -litoral,
      -dominio, -tmax_base, -ajuste_latitude,
      -ajuste_litoral, -ajuste_regional
    )
}

interpolar_temperatura <- function(stations_sf) {
  stations_proj <- st_transform(stations_sf, crs_proj)
  
  grade <- rast(
    xmin = bb["xmin"], xmax = bb["xmax"],
    ymin = bb["ymin"], ymax = bb["ymax"],
    resolution = res_m,
    crs = st_crs(br_proj)$wkt
  )
  
  grade_sf <- st_as_sf(
    as.data.frame(grade, xy = TRUE, cells = TRUE, na.rm = FALSE),
    coords = c("x", "y"),
    crs = crs_proj
  )
  
  idw_sf <- gstat::idw(
    formula = temperatura_maxima ~ 1,
    locations = stations_proj,
    newdata = grade_sf,
    idp = 2,
    maxdist = 510000
  )
  
  xy <- st_coordinates(idw_sf)
  
  idw_df <- idw_sf |>
    st_drop_geometry() |>
    mutate(x = xy[, 1], y = xy[, 2]) |>
    select(x, y, temperatura_maxima = var1.pred)
  
  r_temp <- rast(idw_df, type = "xyz")
  names(r_temp) <- "temperatura_maxima"
  crs(r_temp) <- st_crs(br_proj)$wkt
  
  r_temp |>
    crop(br_vect) |>
    mask(br_vect)
}

stations_base <- see_stations_info() |>
  filter(
    situation_operation == "operating",
    !is.na(longitude_degrees),
    !is.na(latitude_degrees)
  ) |>
  st_as_sf(coords = c("longitude_degrees", "latitude_degrees"), crs = crs_geo) |>
  gerar_temperatura_plausivel(seed = 123)

ui <- fluidPage(
  titlePanel("Temperatura máxima plausível interpolada"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "estacao",
        "Estação",
        choices = stations_base$station_code,
        selected = stations_base$station_code[1]
      ),
      numericInput(
        "temp_manual",
        "Temperatura manual (°C)",
        value = 30,
        min = 0,
        max = 40,
        step = 0.1
      ),
      actionButton("aplicar_manual", "Aplicar na estação"),
      br(), br(),
      actionButton("random_estacao", "Gerar aleatória na estação"),
      actionButton("random_todas", "Gerar aleatória em todas"),
      actionButton("resetar", "Resetar"),
      br(), br(),
      checkboxInput("mostrar_pontos", "Mostrar pontos das estações", TRUE),
      sliderInput("tam_ponto", "Tamanho dos pontos", min = 1, max = 12, value = 3, step = 1),
      br(),
      tags$strong("Estação selecionada"),
      verbatimTextOutput("info_estacao")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Mapa", leafletOutput("mapa", height = "750px")),
        tabPanel(
          "Análises",
          br(),
          plotlyOutput("grafico_densidade", height = "420px"),
          br(),
          DTOutput("tabela_resumo"),
          br(),
          downloadButton("baixar_csv", "Baixar CSV da tabela"),
          br(), br(),
          DTOutput("tabela")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  stations_rv <- reactiveVal(stations_base)
  estacao_selecionada <- reactiveVal(stations_base$station_code[1])
  
  observeEvent(input$aplicar_manual, {
    obj <- stations_rv()
    idx <- which(obj$station_code == input$estacao)
    
    if (length(idx) == 1) {
      obj$temperatura_maxima[idx] <- round(input$temp_manual, 1)
      stations_rv(obj)
      estacao_selecionada(input$estacao)
    }
  })
  
  observeEvent(input$random_estacao, {
    obj <- stations_rv()
    idx <- which(obj$station_code == input$estacao)
    
    if (length(idx) == 1) {
      temp_nova <- gerar_temperatura_plausivel(obj[idx, ], seed = sample.int(1e6, 1))$temperatura_maxima
      obj$temperatura_maxima[idx] <- temp_nova
      stations_rv(obj)
      estacao_selecionada(input$estacao)
      updateNumericInput(session, "temp_manual", value = temp_nova)
    }
  })
  
  observeEvent(input$random_todas, {
    obj <- gerar_temperatura_plausivel(stations_rv(), seed = sample.int(1e6, 1))
    stations_rv(obj)
  })
  
  observeEvent(input$resetar, {
    stations_rv(stations_base)
    estacao_selecionada(stations_base$station_code[1])
    updateSelectInput(session, "estacao", selected = stations_base$station_code[1])
    updateNumericInput(session, "temp_manual", value = stations_base$temperatura_maxima[1])
  })
  
  observeEvent(input$estacao, {
    obj <- stations_rv()
    idx <- which(obj$station_code == input$estacao)
    
    if (length(idx) == 1) {
      estacao_selecionada(input$estacao)
      updateNumericInput(session, "temp_manual", value = obj$temperatura_maxima[idx])
    }
  })
  
  observeEvent(input$tabela_cell_edit, {
    info <- input$tabela_cell_edit
    obj <- stations_rv()
    df <- obj |> st_drop_geometry()
    
    col_real <- info$col + 1
    
    if (names(df)[col_real] == "temperatura_maxima") {
      novo_valor <- suppressWarnings(as.numeric(info$value))
      
      if (!is.na(novo_valor) && novo_valor >= 0 && novo_valor <= 40) {
        obj$temperatura_maxima[info$row] <- round(novo_valor, 1)
        stations_rv(obj)
        
        if (obj$station_code[info$row] == input$estacao) {
          updateNumericInput(session, "temp_manual", value = round(novo_valor, 1))
        }
      }
    }
  })
  
  raster_temp <- reactive({
    interpolar_temperatura(stations_rv())
  })
  
  raster_temp_ll <- reactive({
    project(raster_temp(), "EPSG:4326", method = "bilinear")
  })
  
  stations_ll <- reactive({
    st_transform(stations_rv(), 4326) |>
      mutate(
        popup = paste0(
          "<b>Código:</b> ", station_code, "<br/>",
          "<b>Município:</b> ", station_municipality, "<br/>",
          "<b>UF:</b> ", uf, "<br/>",
          "<b>Altitude:</b> ", round(altitude_m, 1), " m<br/>",
          "<b>Temperatura máxima:</b> ", format(temperatura_maxima, decimal.mark = ","), " °C"
        )
      )
  })
  
  pal_raster <- reactive({
    vals <- values(raster_temp_ll())
    colorNumeric(
      palette = "viridis",
      domain = vals[!is.na(vals)],
      na.color = "transparent"
    )
  })
  
  resumo_temperatura <- reactive({
    df <- stations_rv() |> st_drop_geometry()
    
    data.frame(
      Medida = c(
        "Temperatura máxima média",
        "Desvio padrão",
        "Mediana",
        "Maior observada",
        "Menor observada"
      ),
      Valor = c(
        mean(df$temperatura_maxima, na.rm = TRUE),
        sd(df$temperatura_maxima, na.rm = TRUE),
        median(df$temperatura_maxima, na.rm = TRUE),
        max(df$temperatura_maxima, na.rm = TRUE),
        min(df$temperatura_maxima, na.rm = TRUE)
      )
    ) |>
      mutate(Valor = round(Valor, 2))
  })
  
  tabela_exportacao <- reactive({
    stations_rv() |>
      st_drop_geometry() |>
      transmute(
        station_code,
        station_municipality,
        uf,
        altitude_m,
        temperatura_maxima
      )
  })
  
  output$baixar_csv <- downloadHandler(
    filename = function() {
      paste0("temperaturas_estacoes_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv2(
        tabela_exportacao(),
        file = file,
        row.names = FALSE
      )
    }
  )
  
  output$info_estacao <- renderText({
    obj <- stations_rv()
    idx <- which(obj$station_code == estacao_selecionada())
    req(length(idx) == 1)
    
    paste0(
      "Código: ", obj$station_code[idx], "\n",
      "Município: ", obj$station_municipality[idx], "\n",
      "UF: ", obj$uf[idx], "\n",
      "Altitude: ", round(obj$altitude_m[idx], 1), " m\n",
      "Temperatura máxima: ", format(obj$temperatura_maxima[idx], decimal.mark = ","), " °C"
    )
  })
  
  output$mapa <- renderLeaflet({
    leaflet(options = leafletOptions(zoomControl = TRUE)) |>
      addProviderTiles(providers$CartoDB.Positron) |>
      fitBounds(lng1 = -74, lat1 = -34, lng2 = -34, lat2 = 6)
  })
  
  observe({
    req(raster_temp_ll())
    
    mapa <- leafletProxy("mapa") |>
      clearImages() |>
      clearMarkers() |>
      clearControls() |>
      clearShapes() |>
      addRasterImage(
        raster_temp_ll(),
        colors = pal_raster(),
        opacity = 0.8,
        project = TRUE
      ) |>
      addPolylines(
        data = st_transform(ufs_proj, 4326),
        color = "gray50",
        weight = 0.6,
        opacity = 0.8
      ) |>
      addPolylines(
        data = st_transform(br_proj, 4326),
        color = "black",
        weight = 1.2,
        opacity = 1
      )
    
    if (isTRUE(input$mostrar_pontos)) {
      est <- stations_ll()
      cod_sel <- estacao_selecionada()
      
      mapa <- mapa |>
        addCircleMarkers(
          data = est |> filter(station_code != cod_sel),
          radius = input$tam_ponto,
          color = "red",
          stroke = TRUE,
          weight = 1,
          fillOpacity = 0.85,
          popup = ~popup,
          layerId = ~station_code
        )
      
      est_sel <- est |> filter(station_code == cod_sel)
      
      if (nrow(est_sel) == 1) {
        mapa <- mapa |>
          addCircleMarkers(
            data = est_sel,
            radius = input$tam_ponto + 2,
            color = "yellow",
            fillColor = "yellow",
            stroke = TRUE,
            weight = 2,
            fillOpacity = 1,
            popup = ~popup,
            layerId = ~station_code
          )
      }
    }
    
    vals <- values(raster_temp_ll())
    vals <- vals[!is.na(vals)]
    
    mapa |>
      addLegend(
        pal = pal_raster(),
        values = vals,
        title = "Temperatura máxima (°C)",
        position = "bottomright",
        labFormat = labelFormat(
          digits = 1,
          transform = function(x) round(x, 1)
        )
      )
  })
  
  observeEvent(input$mapa_marker_click, {
    click <- input$mapa_marker_click
    req(click$id)
    
    obj <- stations_rv()
    idx <- which(obj$station_code == click$id)
    
    if (length(idx) == 1) {
      estacao_selecionada(click$id)
      updateSelectInput(session, "estacao", selected = click$id)
      updateNumericInput(session, "temp_manual", value = obj$temperatura_maxima[idx])
    }
  })
  
  output$grafico_densidade <- renderPlotly({
    df <- stations_rv() |> st_drop_geometry()
    
    df <- df |>
      mutate(
        faixa = cut(
          temperatura_maxima,
          breaks = c(20, 24, 28, 32, 36, 40),
          include.lowest = TRUE,
          right = TRUE,
          labels = c("20–24", "24–28", "28–32", "32–36", "36–40")
        )
      )
    
    cores <- c(
      "20–24" = "#2c7bb6",
      "24–28" = "#74add1",
      "28–32" = "#abd9e9",
      "32–36" = "#fdae61",
      "36–40" = "#d7191c"
    )
    
    media <- mean(df$temperatura_maxima, na.rm = TRUE)
    mediana <- median(df$temperatura_maxima, na.rm = TRUE)
    
    dens <- density(df$temperatura_maxima, na.rm = TRUE)
    
    dens_df <- data.frame(
      x = dens$x,
      y = dens$y
    )
    
    ymax <- max(dens_df$y, na.rm = TRUE)
    
    p <- plot_ly()
    
    p <- p |>
      add_trace(
        data = dens_df,
        x = ~x,
        y = ~y,
        type = "scatter",
        mode = "lines",
        fill = "tozeroy",
        line = list(color = "black", width = 2),
        fillcolor = "rgba(248,118,109,0.55)",
        hovertemplate = paste(
          "Temperatura: %{x:.2f} °C",
          "<br>Densidade: %{y:.4f}",
          "<extra></extra>"
        ),
        showlegend = FALSE
      )
    
    for (fx in names(cores)) {
      sub <- df |> filter(faixa == fx)
      
      if (nrow(sub) > 0) {
        p <- p |>
          add_segments(
            data = sub,
            x = ~temperatura_maxima,
            xend = ~temperatura_maxima,
            y = -0.0035,
            yend = -0.0005,
            line = list(color = cores[fx], width = 2),
            opacity = 1,
            hovertemplate = paste0(
              "Temperatura: %{x:.1f} °C",
              "<br>Faixa térmica: ", fx,
              "<extra></extra>"
            ),
            showlegend = FALSE
          )
      }
    }
    
    p <- p |>
      add_segments(
        x = media, xend = media,
        y = 0, yend = ymax * 1.05,
        line = list(color = "black", width = 2, dash = "dash"),
        showlegend = FALSE,
        hoverinfo = "skip"
      ) |>
      add_segments(
        x = mediana, xend = mediana,
        y = 0, yend = ymax * 1.05,
        line = list(color = "black", width = 2, dash = "dot"),
        showlegend = FALSE,
        hoverinfo = "skip"
      ) |>
      add_annotations(
        x = media,
        y = ymax * 0.98,
        text = paste0("Média = ", round(media, 2)),
        showarrow = FALSE,
        xanchor = "left",
        yanchor = "bottom",
        font = list(size = 16, color = "black")
      ) |>
      add_annotations(
        x = mediana,
        y = ymax * 0.90,
        text = paste0("Mediana = ", round(mediana, 2)),
        showarrow = FALSE,
        xanchor = "left",
        yanchor = "bottom",
        font = list(size = 16, color = "black")
      ) |>
      layout(
        title = list(
          text = "<b>Distribuição das temperaturas máximas</b>",
          x = 0.02,
          xanchor = "left",
          font = list(size = 24)
        ),
        xaxis = list(
          title = "Temperatura máxima (°C)",
          zeroline = FALSE,
          range = c(
            min(df$temperatura_maxima, na.rm = TRUE) - 1,
            max(df$temperatura_maxima, na.rm = TRUE) + 1
          )
        ),
        yaxis = list(
          title = "Densidade",
          zeroline = FALSE,
          range = c(-0.005, ymax * 1.12)
        ),
        margin = list(t = 90, l = 70, r = 30, b = 70),
        plot_bgcolor = "white",
        paper_bgcolor = "white",
        showlegend = FALSE
      )
    
    p
  })
  
  output$tabela_resumo <- renderDT({
    resumo_temperatura() |>
      datatable(
        rownames = FALSE,
        options = list(
          dom = "t",
          ordering = FALSE,
          pageLength = 10
        )
      )
  })
  
  output$tabela <- renderDT({
    tabela_exportacao() |>
      datatable(
        editable = list(target = "cell", disable = list(columns = 0:3)),
        options = list(pageLength = 15, scrollX = TRUE)
      )
  })
}

shinyApp(ui, server)