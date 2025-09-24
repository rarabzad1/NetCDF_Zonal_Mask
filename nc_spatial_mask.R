#' @title crops a netcdf file by a mask
#'
#' @description
#' crops a sptio-temporal  netcdf file by a masking area
#' @param ncFile file path of the netcdf file to be cropped
#' @param maskFile file path of the shapefile of the masking area
#' @param ncFileOut (optional) file path of the output of the masked netcdf file. If missing, a file in the same location of the \code{ncFile} is created
#' @return a netcdf file
#' @export nc_spatial_mask
#' @importFrom raster crs shapefile
#' @importFrom ncdf4 nc_open nc_close ncvar_get ncdim_def ncvar_def ncvar_put nc_create
#' @importFrom sp spTransform
#' @importFrom sf st_buffer st_union st_as_sf st_transform st_contains
#' @examples
#' dir.create("c:/rdrs")
#' setwd("c:/rdrs")
#' # download data
#' download.file("https://github.com/rarabzad/RDRS/raw/main/data/data.zip","data.zip")
#' download.file("https://github.com/rarabzad/RDRS/raw/main/data/hru.zip","hru.zip")
#' unzip("data.zip")
#' unzip("hru.zip")
#' hru<-raster::shapefile("hru/finalcat_hru_info.shp")
#' hru<-hru[hru$SubId == 11004375,] # a subset of the basin
#' hru<-as_Spatial(st_sf(st_buffer(st_union(st_make_valid(st_as_sf(hru))),dist = 10),data.frame(id=1))) # simplifying the selected portion of the basin
#' writeOGR(obj = hru,dsn = getwd(),layer = "hru",driver = "ESRI Shapefile",overwrite_layer = T)
#' nc_spatial_mask(ncFile=list.files(pattern  = "*.nc")[1],
#'                                     maskFile = "hru.shp")    
#' @author Rezgar Arabzadeh, University of Waterloo, October 2023
nc_spatial_mask<-function(ncFile,
                            maskFile,
                            ncFileOut=NULL,
                            var=NULL)
{
  library(ncdf4)
  library(sf)
  library(sp)
  library(raster)
  library(dplyr)
  library(geosphere)
  library(lwgeom)
  library(rmapshaper)
  library(zip)
  find_spt_vars <- function(nc)
  {
    if (is.null(nc) || !"var" %in% names(nc)) stop("`nc` must be an already-opened ncdf4 object.")
    
    norm <- function(x) tolower(gsub("[^a-z0-9]", "", as.character(x)))
    time_words <- c("time","t","tim","date","day","hour","minute","second","record","step","tau")
    lat_words  <- c("lat","latitude","y","northing")
    lon_words  <- c("lon","longitude","x","easting")
    time_regex <- paste0("\\b(", paste(time_words, collapse="|"), ")\\b")
    lat_regex  <- paste0("\\b(", paste(lat_words, collapse="|"), ")\\b")
    lon_regex  <- paste0("\\b(", paste(lon_words, collapse="|"), ")\\b")
    
    units_indicate_time <- function(un) {
      if (is.null(un)) return(FALSE)
      un <- tolower(un)
      grepl("since", un) || grepl("\\b(sec|second|seconds|min|minute|minutes|hour|hr|day|doy|julian|year|month)\\b", un)
    }
    units_indicate_lat <- function(un) {
      if (is.null(un)) return(FALSE)
      un <- tolower(un)
      grepl("degree_north|degrees_north|deg_north|degree_n|degrees_n", un) ||
        (grepl("degree|deg", un) && grepl("north|n", un))
    }
    units_indicate_lon <- function(un) {
      if (is.null(un)) return(FALSE)
      un <- tolower(un)
      grepl("degree_east|degrees_east|deg_east|degree_e", un) ||
        (grepl("degree|deg", un) && grepl("east|e", un))
    }
    
    vars <- nc$var
    dims_info <- nc$dim
    found <- character(0)
    
    for (vname in names(vars)) {
      v <- vars[[vname]]
      dim_objs <- v$dim
      if (is.null(dim_objs) || length(dim_objs) == 0) next
      dim_names <- sapply(dim_objs, function(d) d$name)
      
      is_time <- is_lat <- is_lon <- logical(length(dim_names))
      
      # robustly fetch coordinates attribute (if present)
      coords_attr_raw <- try(ncatt_get(nc, vname, "coordinates"), silent = TRUE)
      coord_tokens <- NULL
      if (!inherits(coords_attr_raw, "try-error") && !is.null(coords_attr_raw)) {
        # ncatt_get returns a list; the attribute value is typically coords_attr_raw$value
        if (is.list(coords_attr_raw) && !is.null(coords_attr_raw$value) && is.character(coords_attr_raw$value)) {
          coord_tokens <- unlist(strsplit(coords_attr_raw$value, "\\s+"))
          coord_tokens <- tolower(coord_tokens[nzchar(coord_tokens)])
        } else if (is.character(coords_attr_raw) && length(coords_attr_raw) == 1) {
          coord_tokens <- unlist(strsplit(coords_attr_raw, "\\s+"))
          coord_tokens <- tolower(coord_tokens[nzchar(coord_tokens)])
        }
      }
      
      for (i in seq_along(dim_names)) {
        dn <- dim_names[i]
        if (is.na(dn) || nchar(dn) == 0) next
        dn_norm <- norm(dn)
        if (grepl(time_regex, dn_norm, perl = TRUE)) is_time[i] <- TRUE
        if (grepl(lat_regex, dn_norm, perl = TRUE))  is_lat[i]  <- TRUE
        if (grepl(lon_regex, dn_norm, perl = TRUE))  is_lon[i]  <- TRUE
        
        dimobj <- if (!is.null(dims_info[[dn]])) dims_info[[dn]] else NULL
        if (!is.null(dimobj) && !is.null(dimobj$units)) {
          if (units_indicate_time(dimobj$units)) is_time[i] <- TRUE
          if (units_indicate_lat(dimobj$units))  is_lat[i]  <- TRUE
          if (units_indicate_lon(dimobj$units))  is_lon[i]  <- TRUE
        }
        
        coordvar <- NULL
        if (!is.null(vars[[dn]])) coordvar <- vars[[dn]]
        if (!is.null(coordvar) && !is.null(coordvar$units)) {
          if (units_indicate_time(coordvar$units)) is_time[i] <- TRUE
          if (units_indicate_lat(coordvar$units))  is_lat[i]  <- TRUE
          if (units_indicate_lon(coordvar$units))  is_lon[i]  <- TRUE
        }
        
        # use coordinates attribute tokens as hints if available
        if (!is.null(coord_tokens) && dn %in% coord_tokens) {
          if (any(grepl("lat", coord_tokens))) is_lat[i] <- TRUE
          if (any(grepl("lon", coord_tokens))) is_lon[i] <- TRUE
        }
        
        # safe attempt to read coordinate variable values and infer ranges
        vals_try <- try(ncvar_get(nc, dn), silent = TRUE)
        if (!inherits(vals_try, "try-error") && is.numeric(vals_try)) {
          rng <- range(as.numeric(vals_try), na.rm = TRUE)
          if (is.finite(rng[1]) && is.finite(rng[2])) {
            if (rng[1] >= -90 && rng[2] <= 90) is_lat[i] <- TRUE
            if ((rng[1] >= -180 && rng[2] <= 180) || (rng[1] >= 0 && rng[2] <= 360)) is_lon[i] <- TRUE
          }
        }
      } # end dims loop
      
      n_time <- sum(is_time, na.rm = TRUE)
      n_lat  <- sum(is_lat, na.rm = TRUE)
      n_lon  <- sum(is_lon, na.rm = TRUE)
      
      if ((n_time >= 1) && ((n_lat + n_lon) >= 2)) found <- c(found, vname)
    }
    
    unique(found)
  }
  if(!file.exists(ncFile))   stop("provided netcdf file doesn't exist!")
  if(!file.exists(maskFile)) stop("provided mask file doesn't exist!")
  if(!is.null(ncFileOut)) if(sub(".*\\.", "", basename(ncFileOut)) != "nc") stop ("wrong 'ncFileOut' extension specified. only '*.nc' file are accepted!")
  boundary<-st_read(maskFile)
  if(is.na(crs(boundary))) stop("provided mask file has no projection system!")
  nc<-nc_open(ncFile)
  vars<-find_spt_vars(nc)
  if(!is.null(var))
  {
    if(all(!var %in% vars)) stop("wrong variables specified!")
    if(any(!var %in% vars))
    {
      cat("The following variables are not available in the NetCDF file:\n")
      cat(paste0(var[!var %in% vars],"\n"))
    }
    vars<-var[var %in% vars]
  }
  boundary<-st_transform(boundary,crs("+proj=aea +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")) 
  boundary_buffered<-st_union(st_buffer(st_as_sf(boundary),dist = 10000)) # more than half of the rdrs grid cell sizes
  lat<-ncvar_get(nc,"lat")
  lon<-ncvar_get(nc,"lon")
  latlon<-st_as_sf(data.frame(lon=c(lon),lat=c(lat)),
                   coords = c("lon", "lat"),
                   crs = st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  latlon<-st_transform(latlon,st_crs(boundary_buffered))
  mask<-matrix(st_contains(boundary_buffered, latlon,sparse = F),nrow(lat),ncol(lat))
  if(!any(mask)) stop("the provided mask and the netcdf files don't overlap!")
  id<-which(mask,arr.ind = T)
  if(sum(mask)>1)
  {
    frame<-rbind(apply(id,2,min),apply(id,2,max))
    frameR<-frame[1,1]:frame[2,1]
    frameC<-frame[1,2]:frame[2,2]
  }else{
    frame<-id
    frameR<-frame[1,1,drop=F]
    frameC<-frame[1,2,drop=F]
  }
  maskRC<-mask[frameR,frameC,drop=F]
  maskRC<-ifelse(maskRC,1,NA)
  latRC<-lat[frameR,frameC,drop=F]
  lonRC<-lon[frameR,frameC,drop=F]
  varsData<-list()
  for(j in 1:length(vars))
  {
    cat(paste0("(",j,"/",length(vars),")", " #### Masking: ",vars[j]," ####\n"))
    var_dims<-unlist(lapply(nc$var[[vars[j]]]$dim,function(x) x$name))
    start<-c(frame[1,1],frame[1,2])
    count<-c(length(frameR),length(frameC))
    id_time<-which(!(var_dims %in% c("rlon","rlat")))
    start[id_time]<-1
    count[id_time]<-nc$var[[vars[j]]]$dim[[id_time]]$len
    names(start)<-names(count)<-var_dims
    subset<-ncvar_get(nc    = nc,
                      varid = vars[j],
                      start = start,
                      count = count,collapse_degen = F)
    if(any(is.na(maskRC)))
    {
      if(length(start)>2)
      {
        for(i in 1:dim(subset)[3])  subset[,,i]<-subset[,,i]*maskRC
      }else{
        subset<-subset*maskRC
      }
    }
    if(any(vars[j]==c("RDRS_v2.1_P_GZ_SFC","CaSR_v3.1_P_GZ_SFC")))
    {
      if(length(dim(subset))==3) subset<-apply(subset,c(1,2),mean)
      geo2ele<-function(gph)  (gph*10*9.81)*6371000/(9.81*6371000-gph*10*9.81)
      subset<-geo2ele(subset)
    }
    varsData[[j]]<-subset
  }
  rlon_dim  <- nc$dim$rlon
  rlat_dim  <- nc$dim$rlat
  time_dim  <- nc$dim$time
  ncVars<-list()
  ncVars[[1]] <- nc$var$lon
  ncVars[[2]] <- nc$var$lat
  for(j in 1:length(vars))
  {
    ncVars[[2+j]] <- nc$var[[vars[j]]]
    if(any(vars[j]==c("RDRS_v2.1_P_GZ_SFC","CaSR_v3.1_P_GZ_SFC")))
    {
      spatial_dim_id<-which(unlist(lapply(nc$var[[vars[j]]]$dim,function(x) x$name)) != "time")
      ncVars[[2+j]]<-
      ncvar_def(name = ncVars[[2+j]]$name,
                units = ncVars[[2+j]]$units,
                dim = ncVars[[2+j]]$dim[spatial_dim_id],
                missval = ncVars[[2+j]]$missval,
                longname = ncVars[[2+j]]$longname,
                prec = ncVars[[2+j]]$prec,
                compression = ncVars[[2+j]]$compression,
                chunksizes = ncVars[[2+j]]$chunksizes[spatial_dim_id])
      ncVars[[2+j]]$name <- "Geopotential_Elevation"
      ncVars[[2+j]]$units <- "m"
      ncVars[[2+j]]$longname <- "Analysis: Geopotential Elevation (MASL) converted from Geopotential Height using nc_spatial_mask tool"
    }
  }
  for(j in 1:length(ncVars))
  {
    ncVars[[j]]$dim[[1]]$vals<-ncVars[[j]]$dim[[1]]$vals[frameR]
    ncVars[[j]]$dim[[2]]$vals<-ncVars[[j]]$dim[[2]]$vals[frameC]
    ncVars[[j]]$dim[[1]]$len<-length(frameR)
    ncVars[[j]]$dim[[2]]$len<-length(frameC)
    ncVars[[j]]$size[1:2]<-c(length(frameR),length(frameC))
    ncVars[[j]]$varsize[1:2]<-c(length(frameR),length(frameC))
    ncVars[[j]]$chunksizes<-unlist(lapply(ncVars[[j]]$dim,function(x) x$len))
  }
  if(is.null(ncFileOut))
  {
    ncFileOut<-paste0(sub("\\.\\w+$", "", basename(ncFile)),"_","masked",".nc")
    if(dirname(ncFile)==".")
    {
      ncFileOut<-file.path(getwd(),ncFileOut)
    }else{
      ncFileOut<-file.path(dirname(ncFile),ncFileOut)
    }
  }
  ncnew  <- nc_create(filename = ncFileOut,
                      vars = ncVars,
                      force_v4 = T)
  ncvar_put( ncnew, ncVars[[1]],  lonRC)
  ncvar_put( ncnew, ncVars[[2]],  latRC)
  for(j in 1:length(vars))
  {
    start<-rep(1,length(dim(varsData[[j]])))
    count<-dim(varsData[[j]])
    ncvar_put(ncnew,
              ncVars[[j+2]],
              varsData[[j]])
  }
  nc_close(ncnew)
  nc_close(nc)
}
