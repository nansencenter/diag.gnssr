#=
 = Identify and count the location of all observations in the input file,
 = where locations are defined at the resolution of a grid (for subsequent
 = collocation) and only for elevations below sea level (i.e., excluding
 = inland waters).  Ignored for the moment is where ECMWF Interim sea ice
 = is nearby, so that the corresponding count can be given as a negative
 = number - RD November 2016.
 =#

using My, NetCDF
const CUTOFF           = -44.0                          # SST minimum value

if (argc = length(ARGS)) != 2
  print("\nUsage: jjj $(basename(@__FILE__)) gnss.tds1.txt.all /home/ricani/data/ecmwf/ECMWF_Interim_invariants_0.125.nc\n\n")
  exit(1)
end

tats = collect(90.000:-0.125:-90.000)                                         # first read ECMWF geopotential (height
tons = collect( 0.000: 0.125:359.875)                                         # relative to sea level) and land-sea mask
topo =   ncread(ARGS[2],   "z", start=[1,1,1], count=[-1,-1,-1])
land =   ncread(ARGS[2], "lsm", start=[1,1,1], count=[-1,-1,-1])
scla = ncgetatt(ARGS[2],   "z", "scale_factor") ; adda = ncgetatt(ARGS[2],   "z", "add_offset")
sclb = ncgetatt(ARGS[2], "lsm", "scale_factor") ; addb = ncgetatt(ARGS[2], "lsm", "add_offset")
typeof(scla) != Void && typeof(adda) != Void && (topo = topo * scla + adda)
typeof(sclb) != Void && typeof(addb) != Void && (land = land * sclb + addb)

lats = collect(-90.000:0.125: 90.000)                                         # then define the collocation grid and
lons = collect(  0.000:0.125:359.875)                                         # initialize the count and SST mask
subs = Set(Array(Tuple{Float64, Float64}, 0))
numb = zeros(length(lons), length(lats))
mask =  ones(length(lons), length(lats))

fpa = My.ouvre(ARGS[1],"r")                                                   # identify and count the collocations
for line in readlines(fpa)                                                    # (as long as we are away from land)
  vals = split(line)
  lat = float(vals[2])
  lon = float(vals[3])
# sst = float(vals[14])
  dellat, indlat = findmin(abs(tats - lat))
  dellon, indlon = findmin(abs(tons - lon))
  if topo[indlon,indlat,1] < 50 && land[indlon,indlat,1] < 0.5
    dellat, indlat = findmin(abs(lats - lat))
    dellon, indlon = findmin(abs(lons - lon))
    push!(subs, (lats[indlat], lons[indlon]))
                     numb[indlon,indlat] += 1.0
#   if sst < CUTOFF  mask[indlon,indlat] = -1.0  end
  end
end

fpa = My.ouvre("$(ARGS[1]).locate", "w")                                      # and save them
for loc in subs
  (lat, lon) = loc
  indlat = findfirst(lats, lat)
  indlon = findfirst(lons, lon)
  line = @sprintf("%9.3f %9.3f %8.0f\n", lat, lon, numb[indlon,indlat] * mask[indlon,indlat])
  write(fpa, line)
end
close(fpa)
exit(0)
