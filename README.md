# NetCDF Zonal Mask Shiny App

[![R](https://img.shields.io/badge/R-4.2%2B-blue.svg)](https://www.r-project.org/) 
[![Shiny](https://img.shields.io/badge/Shiny-App-green.svg)](https://shiny.rstudio.com/)

## Overview

This Shiny application provides a **user-friendly interface to apply spatial masking to NetCDF files** using hydrological shapefiles (HRU or other polygon masks). It wraps the `nc_spatial_mask()` function to allow users to:

- Select spatio-temporal variables from a NetCDF file.
- Upload a shapefile (zipped) to define a spatial mask.
- Generate a masked NetCDF file.
- Optionally visualize time-averaged values of masked variables in a dedicated **Plot** tab.

This tool is ideal for hydrologists, climate modelers, and GIS analysts who need to **aggregate or extract NetCDF data over defined spatial zones**.

---

## Features

- Automatically detect **spatio-temporal variables** in uploaded NetCDF files.
- Upload shapefiles as `.zip` archives containing `.shp`, `.shx`, `.dbf`, `.prj` (and optional `.cpg`).
- Generate masked NetCDF files without specifying an output filename (handled internally).
- Download the masked NetCDF file directly from the app.
- Conditional **Plot** tab for time-averaged masked variable values, with dynamic subtabs per variable.
- Detailed logging of all operations for transparency and debugging.

---

## Installation

### Required Packages

The app checks for and installs missing R packages automatically:

- `shiny`, `leaflet`, `sf`, `DT`, `shinyWidgets`, `zip`, `shinyjs`, `ncdf4`, `geosphere`, `dplyr`, `sp`, `lwgeom`, `rmapshaper`, `raster`

### Clone the Repository

```bash
git clone https://github.com/rarabzad1/NetCDF_Zonal_Mask.git
cd NetCDF_Zonal_Mask
````

### Run Locally

```r
library(shiny)
runApp("app.R")
```

Or open `app.R` in RStudio and click **Run App**.

---

## Inputs

1. **NetCDF File (`.nc`)**

   * Main dataset containing spatio-temporal variables.
   * App detects variables with at least **one temporal and two spatial dimensions**.

2. **Shapefile ZIP (`.zip`)**

   * A zipped shapefile defining the spatial mask.
   * Must include `.shp`, `.shx`, `.dbf`, `.prj` files; optional `.cpg` supported.

3. **Variable Selection (Optional)**

   * Multiselect input of spatio-temporal variables to mask.
   * Defaults to all detected variables if none are selected.

4. **Plot Tab Visibility (Optional)**

   * Checkbox to hide or show the **Plot** tab.

---

## Outputs

1. **Masked NetCDF File**

   * Automatically named and saved internally by `nc_spatial_mask()`.
   * Downloadable from the app.

2. **Log Tab**

   * Displays detailed messages for every step (NetCDF inspection, shapefile extraction, masking, errors/warnings).

3. **Plot Tab (Optional)**

   * Shows **time-averaged values** of masked variables.
   * Organized into subtabs for each variable.
   * Provides visual inspection of spatial aggregation results.

---

## Usage

1. Upload a NetCDF file.
2. Upload a shapefile ZIP.
3. Optionally select the variables to mask.
4. Click **"Generate Masked netCDF"**.
5. Download the masked NetCDF once processing is complete.
6. Optionally view the **Plot** tab to inspect results.

---

## Notes

* `ncFileOut` is always `NULL`; the app manages output filenames internally.
* Temporary working directories are used for intermediate files; outputs are copied to a safe temp folder for download.
* Logging provides full traceability, especially helpful on hosted environments like **Posit Cloud**.
* File upload limit is set to **200 MB** (adjustable in `options(shiny.maxRequestSize)`).

---

## Repository

* GitHub: [NetCDF\_Zonal\_Mask](https://github.com/rarabzad1/NetCDF_Zonal_Mask)
* Hosted App (Posit Cloud): [Launch App](https://raventools-netcdf-zonal-mask.share.connect.posit.cloud)

---

## License

MIT License. You are free to use and modify for research or operational purposes.

---

## Contact

**Rezgar Arabzadeh**

* GitHub: [@rarabzad1](https://github.com/rarabzad1)
* Email: `rarabzad@uwaterloo.ca`
