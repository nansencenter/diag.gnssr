#=
 = Loop through ECMWF Interim analyses and extract wind speed at the
 = position of a set of GNSS-R observations - RD April, November 2016.
 =#

using My, NetCDF
const TIMS             = 1336                           # number of 6-h intervals from 2015060100 to 2016042918
const RESOL            = 0.125                          # resolution of the ECMWF Interim analysis
const UCUR             = 1                              # indecies of the wind components
const VCUR             = 2
const PARAMS           = 2
const MISS             = -9999.0                        # generic missing value

if (argc = length(ARGS)) != 2
  print("\nUsage: jjj $(basename(@__FILE__)) gnss.tds1.txt.all.locate_1.0_calib.sort /net/sverdrup-1/vol/sat_auxdata/model/ecmwf/0.125-deg\n")
  print(  "   or: jjj $(basename(@__FILE__)) gnss.tds1.txt.all.locate_1.0_calib.sort   /mnt/10.11.12.232/sat_auxdata/model/ecmwf/0.125-deg\n\n")
  exit(1)
end

gnss = Set(Array(Tuple{Float64, Float64}, 0))                                 # read the GNSS-R locations
fpa = My.ouvre(ARGS[1], "r") ; lines = readlines(fpa) ; close(fpa)
for line in lines
  vals = split(line)
  push!(gnss, (float(vals[2]), float(vals[1])))
end

file = Array(Any, 11, 3)
file[ 1,1] = ARGS[2] * "/interim_2015-06_oper.nc" ; file[ 1,2] =    1: 120 ; file[ 1,3] = 120
file[ 2,1] = ARGS[2] * "/interim_2015-07_oper.nc" ; file[ 2,2] =  121: 244 ; file[ 2,3] = 124
file[ 3,1] = ARGS[2] * "/interim_2015-08_oper.nc" ; file[ 3,2] =  245: 368 ; file[ 3,3] = 124
file[ 4,1] = ARGS[2] * "/interim_2015-09_oper.nc" ; file[ 4,2] =  369: 488 ; file[ 4,3] = 120
file[ 5,1] = ARGS[2] * "/interim_2015-10_oper.nc" ; file[ 5,2] =  489: 612 ; file[ 5,3] = 124
file[ 6,1] = ARGS[2] * "/interim_2015-11_oper.nc" ; file[ 6,2] =  613: 732 ; file[ 6,3] = 120
file[ 7,1] = ARGS[2] * "/interim_2015-12_oper.nc" ; file[ 7,2] =  733: 856 ; file[ 7,3] = 124
file[ 8,1] = ARGS[2] * "/interim_2015-01_oper.nc" ; file[ 8,2] =  857: 980 ; file[ 8,3] = 124
file[ 9,1] = ARGS[2] * "/interim_2015-02_oper.nc" ; file[ 9,2] =  981:1096 ; file[ 9,3] = 116
file[10,1] = ARGS[2] * "/interim_2015-03_oper.nc" ; file[10,2] = 1097:1220 ; file[10,3] = 124
file[11,1] = ARGS[2] * "/interim_2015-04_oper.nc" ; file[11,2] = 1221:1336 ; file[11,3] = 116

wcom = Array(  Int16, TIMS, 2)                                                # and for each location read
wspd = Array(Float64, TIMS, 3)                                                # the wind component timeseries
for lonlat in gnss
  wspd[:,:] = MISS

  (lon, lat) = lonlat
  latind = Int((90.0 - lat) / RESOL + 1)
  lonind = Int(        lon  / RESOL + 1)
  for a = 1:11
    x = ncread(file[a,1], "u10", start=[lonind,latind,1], count=[1,1,file[a,3]]) ; wcom[file[a,2],1] = x
    x = ncread(file[a,1], "v10", start=[lonind,latind,1], count=[1,1,file[a,3]]) ; wcom[file[a,2],2] = x
    addu = ncgetatt(file[a,1], "u10", "add_offset") ; sclu = ncgetatt(file[a,1], "u10", "scale_factor")
    addv = ncgetatt(file[a,1], "v10", "add_offset") ; sclv = ncgetatt(file[a,1], "v10", "scale_factor")
    for b = file[a,2]
      if wcom[b,1] == -32767 || wcom[b,2] == -32767
        wspd[b,1] = wspd[b,2] = wspd[b,3] = MISS
      else
        utmp = wcom[b,1] * sclu + addu
        vtmp = wcom[b,2] * sclv + addv
        wspd[b,1] = (utmp^2 + vtmp^2)^0.5
        wspd[b,2] =  utmp
        wspd[b,3] =           vtmp
      end
    end
  end

  tmp = @sprintf("%9.3f.%9.3f", lat, lon) ; tmq = replace(tmp, " ", ".")      # then save wind speed
  tmr = "ecmwf/ecmwf.$tmq" ; fpa = My.ouvre(tmr, "w")
  date = "2015060100"
  for a = 1:TIMS
    form = @sprintf("%10s0000 %11.4f %11.4f %11.4f %11.4f %11.4f -9999.0000  -9999.0000  -9999.0000  -9999.0000\n",
           date, lat, lon, wspd[a,1], wspd[a,2], wspd[a,3])
    write(fpa, form)
    date = My.dateadd(date, 6, "hr")
  end
  close(fpa)
end
exit(0)

#=
lats =  -79.875:0.25:79.875  ; latn = length(lats)
lons = -179.875:0.25:179.875 ; lonn = length(lons)
  lons =   0.125:0.25:359.875 ; lonn = length(lons)                           # locs[   1,  1] = (-180.00,-90.00)
end                                                                           # locs[   1,720] = (-180.00, 89.75)
locs = [   (x, y)        for x=lons, y=lats]                                  # locs[1440,  1] = ( 179.75,-90.00)
mask = [in((x, y), buoy) for x=lons, y=lats]                                  # locs[1440,720] = ( 179.75, 89.75)
fend = locs[mask .== true] ; (locn,) = size(fend)                             # and define all filename endings
=#
