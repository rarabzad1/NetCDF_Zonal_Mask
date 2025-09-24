# -------------------------------
# Dockerfile for NetCDF Zonal Mask Shiny App
# -------------------------------

# Use official R Shiny image as base
FROM rocker/shiny:latest

# System libraries required for spatial packages
RUN apt-get update && apt-get install -y \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    unzip \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install R packages required by the app
RUN R -e "install.packages(c(
    'shiny', 'leaflet', 'sf', 'DT', 'shinyWidgets', 'zip', 'shinyjs',
    'ncdf4', 'geosphere', 'dplyr', 'sp', 'lwgeom', 'rmapshaper', 'raster'
), repos = 'https://cloud.r-project.org')"

# Copy app files into Docker image
COPY app.R /srv/shiny-server/
# If you have other scripts (e.g., nc_spatial_mask.R), copy them too
# COPY nc_spatial_mask.R /srv/shiny-server/

# Make app.R the default app for shiny server
RUN chown -R shiny:shiny /srv/shiny-server/

# Expose default Shiny port
EXPOSE 3838

# Run Shiny Server
CMD ["/usr/bin/shiny-server"]
