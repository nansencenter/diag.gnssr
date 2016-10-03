#=
 = Create a Sentinel-1 subdomain in netCDF following a conversion of the full scene from
 = SAFE format to netCDF.  A short text file (hdr) that partly duplicates information in
 = manifest.dim (SAFE XML file) is also created.  The sigma-naught of a calibrated netCDF
 = file is assumed to be available (using SNAP toolbox or Nansat) for the full acquisition
 = domain.  An intermediate netCDF subdomain is first obtained using the netCDF Operators
 = (nco.sourceforge.net); data is then copied to a final netCDF file.  Include least-squares
 = mappings between the native SAR and the arbitrary netCDF grids - RD September 2016.
 =#

using My, LightXML, NetCDF, Optim

if (argc = length(ARGS)) != 6
  print("\nUsage: jjj $(basename(@__FILE__)) S1A_EW_GRDM_1SDH_20151203T160430_20151203T160530_008880_00CB22_53F6.zip 60.050 25.900 50000 50000 53\n")
  print("       jjj $(basename(@__FILE__)) SAFE_format.zip                                                     reflat reflon merdis zondis far_range_pixels_to_omit\n\n")
  exit(1)
end
xmlfil = split(ARGS[1], ".")[1] * ".SAFE/manifest.dim"
xyzfil = split(ARGS[1], ".")[1] * ".xyzzy.nc"
calfil = split(ARGS[1], ".")[1] * "_Cal.nc"
netfil = split(ARGS[1], ".")[1] * ".nc"
stmfil = split(ARGS[1], ".")[1]
reflat = parse(Float64, ARGS[2])
reflon = parse(Float64, ARGS[3])
merdis = parse(Float64, ARGS[4])
zondis = parse(Float64, ARGS[5])
farpix = parse(  Int64, ARGS[6])

timbeg = join(split(split(ARGS[1], "_")[5], 'T'))                             # associate the scene with a single
timend = join(split(split(ARGS[1], "_")[6], 'T'))                             # midpoint time
timdif = datesous(timbeg, timend, "sc") / 2.0
timcen = dateadd( timbeg, timdif, "sc")
datcen = timcen[1: 4] * "-" * timcen[ 5: 6] * "-" * timcen[ 7: 8] * "-" * timcen[9:10]
hmscen = timcen[9:10] * ":" * timcen[11:12] * ":" * timcen[13:14]
stem   = timcen[1: 4] * "-" * timcen[ 5: 6] * "-" * timcen[ 7: 8] * "-" * timcen[9:14]

xmldat = split(string(get_elements_by_tagname(root(parse_file(xmlfil)), "Dataset_Sources")[1]), "\n")
function xmlparse(xmldat::Array{SubString{ASCIIString},1}, key::ASCIIString)
  line = filter(x -> contains(x, key), xmldat)
  split(split(line[1], "<")[2], ">")[2]
end

nlines        = parse(Int64,   xmlparse(xmldat,        "Raster height"))      # and telescope through the data to get
npixels       = parse(Int64,   xmlparse(xmldat,         "Raster width"))      # the location closest to that desired
pixel_spacing = parse(Float64, xmlparse(xmldat, "Range sample spacing"))      # (assuming range = azimuth resolution)
mindis = minlat = minlon = 9e9
minlin = div( nlines, 2)
minpix = div(npixels, 2)
for a = 25:-1:0
  for b = minlin - 2^a : 2^a : minlin + 2^a
    if b >= 0 && b < nlines
      for c = minpix - 2^a : 2^a : minpix + 2^a
        if c >= 0 && c < npixels
          gdlatitude  = ncread(calfil,  "latitude", start=[c+1,b+1], count=[1,1])[1]
          gdlongitude = ncread(calfil, "longitude", start=[c+1,b+1], count=[1,1])[1]
          tmpdis = (reflat - gdlatitude)^2 + (reflon - gdlongitude)^2
          if tmpdis < mindis
            mindis = tmpdis
            minlat = gdlatitude
            minlon = gdlongitude
            minlin = b
            minpix = c
#           @printf("a b c are (%8d %8d %8d) and lat lon mindis are %lf %lf %e ***\n", a, b, c, gdlatitude, gdlongitude, mindis)
#         else
#           @printf("a b c are (%8d %8d %8d) and lat lon tmpdis are %lf %lf %e\n", a, b, c, gdlatitude, gdlongitude, tmpdis)
          end
        end
      end
    end
  end
end
print("\nreflat reflon are $reflat $reflon\n")
print("minlat minlon are $minlat $minlon at line pixel $minlin $minpix\n")
print("allowable distance between them is about $(10 * pixel_spacing) m\n")
print("actual    distance between them is about $(111000 * mindis^0.5) m\n")
if 111000 * mindis^0.5 > 10 * pixel_spacing ; exit(0) ; end

reduce = 0                                                                    # calculate the desired number of
gridint = 1                                                                   # subscenes in azimuth and range
subscenesize = 1
resol = pixel_spacing
startline_valid  = minlin - convert(Int64, div(merdis, resol))
stopline_valid   = minlin + convert(Int64, div(merdis, resol))
startpixel_valid = minpix - convert(Int64, div(zondis, resol))
stopixel_valid   = minpix + convert(Int64, div(zondis, resol))
if startline_valid  < 1       ; startline_valid  = 1       ; end
if startpixel_valid < 1       ; startpixel_valid = 1       ; end
if stopline_valid   > nlines  ; stopline_valid   = nlines  ; end
if stopixel_valid   > npixels - farpix ; stopixel_valid   = npixels -farpix ; end
numlines_valid   = stopline_valid -  startline_valid + 1
numpixels_valid  = stopixel_valid - startpixel_valid + 1
aziblocks        = numlines_valid
rngblocks        = numpixels_valid
nlines2skipbeg   = startline_valid - 1
npixels2skipbeg  = startpixel_valid - 1
nlines2skipend   = nlines - stopline_valid
npixels2skipend  = npixels - stopixel_valid

@printf("reduce power      %d\n", reduce)                                     # create the intermediate netCDF
@printf("new grid int      %d\n", gridint)                                    # subdomain using netCDF Operators
@printf("subscenesize      %d\n", subscenesize)
@printf("aziblocks         %d\n", aziblocks)
@printf("rngblocks         %d\n", rngblocks)
@printf("nlines2skipbeg   %10d   npixels2skipbeg  %10d\n", nlines2skipbeg, npixels2skipbeg)
@printf("startline_valid  %10d   startpixel_valid %10d\n", startline_valid, startpixel_valid)
@printf("numlines_valid   %10d   numpixels_valid  %10d\n", numlines_valid, numpixels_valid)
@printf("stopline_valid   %10d   stopixel_valid   %10d\n", stopline_valid, stopixel_valid)
@printf("nlines2skipend   %10d   npixels2skipend  %10d\n", nlines2skipend, npixels2skipend)
@printf("nlines           %10d   npixels          %10d\n\n", nlines, npixels)
print("ncks -d y,$(startline_valid-1),$(stopline_valid-1) -d x,$(startpixel_valid-1),$(stopixel_valid-1) $calfil $xyzfil\n")
  run(`ncks -d y,$(startline_valid-1),$(stopline_valid-1) -d x,$(startpixel_valid-1),$(stopixel_valid-1) $calfil $xyzfil`)

if     resol >= 10000 ; resolution = @sprintf(    "%.0f", resol)              # then create an output file
elseif resol >= 1000  ; resolution = @sprintf(   "0%.0f", resol)
elseif resol >= 100   ; resolution = @sprintf(  "00%.0f", resol)
elseif resol >= 10    ; resolution = @sprintf( "000%.0f", resol)
elseif resol >= 1     ; resolution = @sprintf("0000%.0f", resol)
else                  ; resolution = "00000" ; end
polars = [xmlparse(xmldat, "mds1_tx_rx_polar"), xmlparse(xmldat, "mds2_tx_rx_polar"),
          xmlparse(xmldat, "mds3_tx_rx_polar"), xmlparse(xmldat, "mds4_tx_rx_polar")]
if "HH" in polars                                                       ; pol = "HH" ; end
if                                     "VV" in polars                   ; pol = "VV" ; end
if "HH" in polars && "HV" in polars                                     ; pol = "DH" ; end
if                                     "VV" in polars && "VH" in polars ; pol = "DV" ; end
if "HH" in polars && "HV" in polars && "VV" in polars && "VH" in polars ; pol = "QP" ; end
stem *= pol ; outfil = @sprintf("%s.%s.sar.nc", stem, resolution)
print("\n/home/ricani/bin/nc.template.sar $outfil $datcen $aziblocks $rngblocks\n")
    run(`/home/ricani/bin/nc.template.sar $outfil $datcen $aziblocks $rngblocks`)

print("reading sigo, angl, lats, and lons from $xyzfil\n")                    # read and orient the subdomain grids
#nc = NetCDF.open(xyzfil, mode=NC_NOWRITE, readdimvar=false)
sigo = ncread(xyzfil,      "Sigma0_HH", start=[1,1], count=[-1,-1])
angl = ncread(xyzfil, "incident_angle", start=[1,1], count=[-1,-1])
lats = ncread(xyzfil,       "latitude", start=[1,1], count=[-1,-1])
lons = ncread(xyzfil,      "longitude", start=[1,1], count=[-1,-1])
temp = Array(Float64, rngblocks, aziblocks)
ncclose(xyzfil)
rm(xyzfil)

if     lats[1,1] > lats[1,aziblocks] && lons[1,1] < lons[rngblocks,1]
  comment = "with no lat/lon switch"
elseif lats[1,1] < lats[1,aziblocks] && lons[1,1] < lons[rngblocks,1]
  for a = 1:aziblocks, b = 1:rngblocks ; temp[b,a] = sigo[            b,aziblocks+1-a] ; end ; for a = 1:aziblocks, b = 1:rngblocks ; sigo[b,a] = temp[b,a] ; end
  for a = 1:aziblocks, b = 1:rngblocks ; temp[b,a] = angl[            b,aziblocks+1-a] ; end ; for a = 1:aziblocks, b = 1:rngblocks ; angl[b,a] = temp[b,a] ; end
  for a = 1:aziblocks, b = 1:rngblocks ; temp[b,a] = lats[            b,aziblocks+1-a] ; end ; for a = 1:aziblocks, b = 1:rngblocks ; lats[b,a] = temp[b,a] ; end
  for a = 1:aziblocks, b = 1:rngblocks ; temp[b,a] = lons[            b,aziblocks+1-a] ; end ; for a = 1:aziblocks, b = 1:rngblocks ; lons[b,a] = temp[b,a] ; end
  comment = "with lat switch"
elseif lats[1,1] > lats[1,aziblocks] && lons[1,1] > lons[rngblocks,1]
  for a = 1:aziblocks, b = 1:rngblocks ; temp[b,a] = sigo[rngblocks+1-b,            a] ; end ; for a = 1:aziblocks, b = 1:rngblocks ; sigo[b,a] = temp[b,a] ; end
  for a = 1:aziblocks, b = 1:rngblocks ; temp[b,a] = angl[rngblocks+1-b,            a] ; end ; for a = 1:aziblocks, b = 1:rngblocks ; angl[b,a] = temp[b,a] ; end
  for a = 1:aziblocks, b = 1:rngblocks ; temp[b,a] = lats[rngblocks+1-b,            a] ; end ; for a = 1:aziblocks, b = 1:rngblocks ; lats[b,a] = temp[b,a] ; end
  for a = 1:aziblocks, b = 1:rngblocks ; temp[b,a] = lons[rngblocks+1-b,            a] ; end ; for a = 1:aziblocks, b = 1:rngblocks ; lons[b,a] = temp[b,a] ; end
  comment = "with lon switch"
else
  for a = 1:aziblocks, b = 1:rngblocks ; temp[b,a] = sigo[rngblocks+1-b,aziblocks+1-a] ; end ; for a = 1:aziblocks, b = 1:rngblocks ; sigo[b,a] = temp[b,a] ; end
  for a = 1:aziblocks, b = 1:rngblocks ; temp[b,a] = angl[rngblocks+1-b,aziblocks+1-a] ; end ; for a = 1:aziblocks, b = 1:rngblocks ; angl[b,a] = temp[b,a] ; end
  for a = 1:aziblocks, b = 1:rngblocks ; temp[b,a] = lats[rngblocks+1-b,aziblocks+1-a] ; end ; for a = 1:aziblocks, b = 1:rngblocks ; lats[b,a] = temp[b,a] ; end
  for a = 1:aziblocks, b = 1:rngblocks ; temp[b,a] = lons[rngblocks+1-b,aziblocks+1-a] ; end ; for a = 1:aziblocks, b = 1:rngblocks ; lons[b,a] = temp[b,a] ; end
  comment = "with lat+lon switch"
end
ullat = lats[1,1]                 ; ullon = lons[1,1]
urlat = lats[rngblocks,1]         ; urlon = lons[rngblocks,1]
lllat = lats[1,aziblocks]         ; lllon = lons[1,aziblocks]
lrlat = lats[rngblocks,aziblocks] ; lrlon = lons[rngblocks,aziblocks]

print("writing sigo, angl, lats, and lons   to $outfil $comment\n")           # and store the (reoriented) fields
ncwrite(sigo, outfil, "sigo", start=[1,1,1], count=[-1,-1,-1])
ncwrite(angl, outfil, "angl", start=[1,1,1], count=[-1,-1,-1])                # also identify a subset of locations
ncwrite(lats, outfil, "lats", start=[1,1,1], count=[-1,-1,-1])                # for mapping between the native SAR
ncwrite(lons, outfil, "lons", start=[1,1,1], count=[-1,-1,-1])                # and arbitrary netCDF grids
ncclose(outfil)

aziskip = convert(Int64, div(aziblocks, 25)) ; azisteps = 0 ; for b = 1:aziskip:aziblocks ; azisteps += 1 ; end
rngskip = convert(Int64, div(rngblocks, 25)) ; rngsteps = 0 ; for b = 1:rngskip:rngblocks ; rngsteps += 1 ; end
print("mapping between $azisteps X $rngsteps native SAR grid locations and corresponding netCDF \n")

varmas = Array(Float64, azisteps * rngsteps, 4)                               # and define a least squares metric
varcol = Array(Float64, azisteps * rngsteps)                                  # (this sqerror closure requires data
function sqerror(coef::Array{Float64,1})                                      #  arrays in global scope)
  err = 0.0
  for i in 1:azisteps * rngsteps
    res  = coef[1] * varmas[i,1] + coef[2] * varmas[i,2] + coef[3] * varmas[i,3] + coef[4] * varmas[i,4]
    err += (varcol[i] - res)^2
  end
  return err
end

a = 1 ; for b = 1:aziskip:aziblocks, c = 1:rngskip:rngblocks                  # first get a least-squares mapping from
  varmas[a,1] = 1.0                                                           # indecies of lat/lon to actual lat/lon
  varmas[a,2] = c
  varmas[a,3] = b
  varmas[a,4] = b * c
  a += 1
end
a = 1 ; for b = 1:aziskip:aziblocks, c = 1:rngskip:rngblocks ; varcol[a] =               lats[c,b] ; a += 1 ; end ; linreslat = optimize(sqerror, [0.0, 0.0, 0.0, 0.0], iterations = 10000)
a = 1 ; for b = 1:aziskip:aziblocks, c = 1:rngskip:rngblocks ; varcol[a] =               lons[c,b] ; a += 1 ; end ; linreslon = optimize(sqerror, [0.0, 0.0, 0.0, 0.0], iterations = 10000)

a = 1 ; for b = 1:aziskip:aziblocks, c = 1:rngskip:rngblocks                  # then get a reverse mapping from actual
  varmas[a,1] = 1.0                                                           # lat/lon to arbitrary netCDF lat/lon
  varmas[a,2] = lons[c,b]                                                     # (used to plot track data in the native
  varmas[a,3] = lats[c,b]                                                     #  SAR coordinate system)
  varmas[a,4] = lats[c,b] * lons[c,b]
  a += 1
end
a = 1 ; for b = 1:aziskip:aziblocks, c = 1:rngskip:rngblocks ; varcol[a] =  50.0 - (b - 1) * 0.001 ; a += 1 ; end ; tarreslat = optimize(sqerror, [0.0, 0.0, 0.0, 0.0], iterations = 100000)
a = 1 ; for b = 1:aziskip:aziblocks, c = 1:rngskip:rngblocks ; varcol[a] = 280.0 + (c - 1) * 0.001 ; a += 1 ; end ; tarreslon = optimize(sqerror, [0.0, 0.0, 0.0, 0.0], iterations = 100000)

degtorad = pi / 180.0                                                         # finally get a track angle and save such
trackval = cos((lats[1,aziblocks] + lats[1,1]) / 2.0 * degtorad)              # subdomain values to a hdr file along
trackyal =      lats[1,aziblocks] - lats[1,1]                                 # with duplicate parts of the manifest
trackxal =     (lons[1,aziblocks] - lons[1,1]) * trackval
trackdir = (450.0 - atan2(trackyal, trackxal) / degtorad) % 360.0

hdrfil = @sprintf("%s.%s.hdr", stem, resolution)
print("\n") ; fpa = My.ouvre(hdrfil, "w") ; print("\n")
form = @sprintf("%-23s%s\n",                "Product",                                                           xmlparse(xmldat,                 "Product name")) ; write(fpa, form)
form = @sprintf("%-23s%-18s%s\n",           "ProductType",          xmlparse(xmldat,       "Satellite mission"), xmlparse(xmldat,                 "Product type")) ; write(fpa, form)
form = @sprintf("%-23s%s\n",                "SystemID",                                                          xmlparse(xmldat, "Processing system identifier")) ; write(fpa, form)
form = @sprintf("%-23s%-18s%-18s%-18s%s\n", "Polarizations",        xmlparse(xmldat,        "mds1_tx_rx_polar"), xmlparse(xmldat,             "mds2_tx_rx_polar"),
                                                                    xmlparse(xmldat,        "mds3_tx_rx_polar"), xmlparse(xmldat,             "mds4_tx_rx_polar")) ; write(fpa, form)
form = @sprintf("%-23s%-18s%s\n",           "AscDescMode",          xmlparse(xmldat, "ASCENDING or DESCENDING"), xmlparse(xmldat,         "Right or left facing")) ; write(fpa, form)
form = @sprintf("%-23s%-18s%s\n",           "RngAziPixelSpacing",   xmlparse(xmldat,    "Range sample spacing"), xmlparse(xmldat,       "Azimuth sample spacing")) ; write(fpa, form)
form = @sprintf("%-23s%s\n",                "AntennaBeams",                                                      xmlparse(xmldat,             "Acquisition mode")) ; write(fpa, form)
form = @sprintf("%-23s%d\n",                "NumberLines",                                                                                                 nlines) ; write(fpa, form)
form = @sprintf("%-23s%d\n",                "NumberPixels",                                                                                               npixels) ; write(fpa, form)
form = @sprintf("%-23s%-18s%-18s%s\n",      "SceneCentreTime",                                                                             datcen, hmscen, datcen) ; write(fpa, form)
form = @sprintf("%-23s%.6f\n",              "UpperLeft_Latitude",                                                                                           ullat) ; write(fpa, form)
form = @sprintf("%-23s%.6f\n",              "UpperLeft_Longitude",                                                                                          ullon) ; write(fpa, form)
form = @sprintf("%-23s%.6f\n",              "UpperRight_Latitude",                                                                                          urlat) ; write(fpa, form)
form = @sprintf("%-23s%.6f\n",              "UpperRight_Longitude",                                                                                         urlon) ; write(fpa, form)
form = @sprintf("%-23s%.6f\n",              "LowerLeft_Latitude",                                                                                           lllat) ; write(fpa, form)
form = @sprintf("%-23s%.6f\n",              "LowerLeft_Longitude",                                                                                          lllon) ; write(fpa, form)
form = @sprintf("%-23s%.6f\n",              "LowerRight_Latitude",                                                                                          lrlat) ; write(fpa, form)
form = @sprintf("%-23s%.6f\n",              "LowerRight_Longitude",                                                                                         lrlon) ; write(fpa, form)
form = @sprintf("%-23s%.6f\n",              "CalculatedTrackAngle",                                                                                      trackdir) ; write(fpa, form)
form = @sprintf("%-23s%s\n",                "ScalingOffsetFromLUT",                                                                                    "0.000000") ; write(fpa, form)
form = @sprintf("%-23s%d\n",                "StartLineValid",                                                                                     startline_valid) ; write(fpa, form)
form = @sprintf("%-23s%d\n",                "StopLineValid",                                                                                       stopline_valid) ; write(fpa, form)
form = @sprintf("%-23s%d\n",                "StartPixelValid",                                                                                   startpixel_valid) ; write(fpa, form)
form = @sprintf("%-23s%d\n",                "StopPixelValid",                                                                                      stopixel_valid) ; write(fpa, form)
form = @sprintf("%-23s%d\n",                "NumAzimuthSubScenes",                                                                                      aziblocks) ; write(fpa, form)
form = @sprintf("%-23s%d\n",                "NumRangeSubScenes",                                                                                        rngblocks) ; write(fpa, form)
form = @sprintf("%-23s%d\n",                "SubSceneSize",                                                                                                     1) ; write(fpa, form)
form = @sprintf("%-23s %18.10f %18.10f %18.10f %18.10f\n",        "lat_plane_coefs", linreslat.minimum[1], linreslat.minimum[2], linreslat.minimum[3], linreslat.minimum[4]) ; write(fpa, form)
form = @sprintf("%-23s %18.10f %18.10f %18.10f %18.10f\n",        "lon_plane_coefs", linreslon.minimum[1], linreslon.minimum[2], linreslon.minimum[3], linreslon.minimum[4]) ; write(fpa, form)
form = @sprintf("%-23s %18.10f %18.10f %18.10f %18.10f\n", "target_lat_plane_coefs", tarreslat.minimum[1], tarreslat.minimum[2], tarreslat.minimum[3], tarreslat.minimum[4]) ; write(fpa, form)
form = @sprintf("%-23s %18.10f %18.10f %18.10f %18.10f\n", "target_lon_plane_coefs", tarreslon.minimum[1], tarreslon.minimum[2], tarreslon.minimum[3], tarreslon.minimum[4]) ; write(fpa, form)
close(fpa)
exit(0)


#=
saffil = stmfil[1] * ".SAFE/manifest.safe" ; safdat = string(get_elements_by_tagname(root(parse_file(saffil)), "metadataSection")[1])
dimfil = stmfil[1] * ".SAFE/manifest.dim"  ; dimdat = string(get_elements_by_tagname(root(parse_file(dimfil)), "metadataSection")[1])
coords = split(safdat, "gml:coordinates")[2][2:end-2]
timsta = split(safdat,  "safe:startTime")[2][2:end-2]
timend = split(safdat,   "safe:stopTime")[2][2:end-2]
ascdes = split(safdat,         "s1:pass")[2][2:end-2]
coords = split(xmldat, "gml:coordinates")[2][2:end-2]
coords = split(xmldat,  "safe:startTime")[2][2:end-2]
=#
