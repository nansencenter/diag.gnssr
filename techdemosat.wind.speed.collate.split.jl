#=
 = Split the 6-h 0.125-deg observations into calibration and validation groups, where the
 = calibration subset is taken to be at locations within each LORES-degree box that contains
 = the largest number of available observations - RD March, November 2016.
 =#

using My
const LORES            = 1.0                            # resolution of the calib/valid subset
const HIRES            = 0.125                          # resolution of the gridded observations
const CUTOFF           = 0                              # minimum number of obs at calib location
const LAT              = 1
const LON              = 2
const NUM              = 3
const PARAMS           = 3

if (argc = length(ARGS)) != 1
  print("\nUsage: jjj $(basename(@__FILE__)) gnss.tds1.txt.all.locate\n\n")
  exit(1)
end

ratio = LORES / HIRES                                                         # define the output grids at low res
lats = convert(Int64, 180 / LORES) + 1
lons = convert(Int64, 360 / LORES)
data = zeros(lons, lats, PARAMS)

fpa = My.ouvre(ARGS[1], "r") ; lines = readlines(fpa) ; close(fpa)            # loop through the available obs and
tmp = @sprintf("%s_%.1f_calib", ARGS[1], LORES) ; fpb = My.ouvre(tmp, "w")    # retain the largest number in memory
tmp = @sprintf("%s_%.1f_valid", ARGS[1], LORES) ; fpc = My.ouvre(tmp, "w")    # and write the rest directly to file

for line in lines
  tmp = split(line)
  lat = float(tmp[1])
  lon = float(tmp[2])
  num = float(tmp[3])
  hilatind = (lat + 90.0) / HIRES
  hilonind =  lon         / HIRES
  lolatind = trunc(Int, hilatind / ratio) + 1
  lolonind = trunc(Int, hilonind / ratio) + 1

  if num > data[lolonind,lolatind,NUM] && num > CUTOFF
    if data[lolonind,lolatind,NUM] > 0
      line = @sprintf("%9.3f %9.3f %8.0f\n", data[lolonind,lolatind,LAT], data[lolonind,lolatind,LON], data[lolonind,lolatind,NUM])
      write(fpc, line)
    end
    data[lolonind,lolatind,LAT] = lat
    data[lolonind,lolatind,LON] = lon
    data[lolonind,lolatind,NUM] = num
  else
    line = @sprintf("%9.3f %9.3f %8.0f\n", lat, lon, num)
    write(fpc, line)
  end
end

for a = 1:lats, b = 1:lons
  if data[b,a,NUM] > 0
    line = @sprintf("%9.3f %9.3f %8.0f\n", data[b,a,LAT], data[b,a,LON], data[b,a,NUM])
    write(fpb, line)
  end
end

close(fpb)
close(fpc)
exit(0)

#=
minlons = collect(0.0:RESOL:360.0-RESOL)                                      # define a search grid for speed
midlons = minlons + RESOL / 2.0
maxlons = minlons + RESOL
minlats = collect(-90.0:RESOL:90.0-RESOL)
midlats = minlats + RESOL / 2.0
maxlats = minlats + RESOL
numlons = length(minlons)
numlats = length(minlats)

valn = length(lines)
vals = Array(Float64, valn, PARAMS)
vals = SharedArray(Float64, (valn, PARAMS))

locs = Array(Tuple{Float64, Float64}, 0)                                      # find the calib locations
for a = 1:numlons                                                             # (largest number of daily obs
  for b = 1:numlats                                                           #  in each gridbox if available)
#   @printf("%8.2f %8.2f\n", midlons[a], midlats[b])
    maxlat = maxlon = maxnum = -1.0
    for c = 1:valn
      if vals[c,NUM] > maxnum && minlons[a] <= vals[c,LON] < maxlons[a] &&
                                 minlats[b] <= vals[c,LAT] < maxlats[b]
        maxlat = vals[c,LAT]
        maxlon = vals[c,LON]
        maxnum = vals[c,NUM]
      end
    end
    if maxnum > CUTOFF
      push!(locs, (maxlat, maxlon))
    end
  end
end

function locate(minlat::Float64, maxlat::Float64, minlon::Float64, maxlon::Float64)
  maxlat = maxlon = maxnum = -1.0
  for c = 1:valn
    if vals[c,NUM] > maxnum && minlon <= vals[c,LON] < maxlon && minlat <= vals[c,LAT] < maxlat
      maxlat = vals[c,LAT]
      maxlon = vals[c,LON]
      maxnum = vals[c,NUM]
    end
  end
  (maxlat, maxlon, maxnum)
end

locs = Array(Tuple{Float64, Float64}, 0)                                      # find the calib locations
for a = 1:1#numlons                                                             # (largest number of daily obs
  for b = 1:numlats                                                           #  in each gridbox if available)
#   @printf("%8.2f %8.2f\n", midlons[a], midlats[b])
    (maxlat, maxlon, maxnum) = locate(minlats[b], maxlats[b], minlons[a], maxlons[a])
    if maxnum > CUTOFF
      push!(locs, (maxlat, maxlon))
    end
  end
end
=#
