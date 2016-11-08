#=
 = Collate TDS-1 observations into averages at six-hourly 0.125-degree resolution
 = (corresponding to that of the ECMWF Interim grid) including all input variables
 = (speed, incang, elo, azo, ela, aza, and gain) over one day - RD October 2016.
 =#

using My
const RDAT             = "20140701000000"               # reference date (for date calculations)

if (argc = length(ARGS)) != 2
  print("\nUsage: jjj $(basename(@__FILE__)) gnss_obs_file days_since_$(RDAT[1:8])\n")
  print(" e.g.: jjj $(basename(@__FILE__)) gnss.tds1.txt 200\n\n")
  exit(1)
end                                                                           # calculate the 6-h index of the beginning
dnum = parse(Int64, ARGS[2]) * 4 + 1                                          # of the day of interest (relative to RDAT)

inds = [4 5 6 7 8 9 10]                                                       # define the indecies of the six-hourly averages,
lats = collect(-90.000:0.125:90.000)                                          # the output grid, the days in the month, and
lons = collect(  0.000:0.125:359.875)                                         # then initialize the (large) averaging arrays

dats = Array(UTF8String, 0)
datn = Array(   Float64, 0)
      date  = RDAT ;            push!(dats, date) ; push!(datn, My.datesous(RDAT, date, "hr")) ; date = My.dateadd(date, 6, "hr")
while date != "20170101000000"  push!(dats, date) ; push!(datn, My.datesous(RDAT, date, "hr")) ; date = My.dateadd(date, 6, "hr")  end
sums = zeros(length(lons), length(lats), 4, length(inds))
numb = zeros(length(lons), length(lats), 4)

tot = false                                                                   # find the nearest gridbox in time and
fpa = My.ouvre(ARGS[1], "r")                                                  # space and sum any values found
for line in readlines(fpa)
  vals = split(line)
  spd = float(vals[4])
  if spd > -333.0 && spd < 333.0
    dat = My.datesous(RDAT, vals[1], "hr")
    deldat, inddat = findmin(abs(datn - dat))
    inddat += 1 - dnum

    if in(inddat, [1, 2, 3, 4])
# @show line
      lat = float(vals[2])
      lon = float(vals[3]) ; lon < 0 && (lon += 360) ; lon >= 360 && (lon -= 360)
      dellat, indlat = findmin(abs(lats - lat))
      dellon, indlon = findmin(abs(lons - lon))
      for (a, vind) in enumerate(inds)
        sums[indlon,indlat,inddat,a] += float(vals[vind])
      end
      numb[indlon,indlat,inddat] += 1.0
      tot = true
    end
  end
end
close(fpa)

if tot                                                                        # and save any valid averages
  fpb = My.ouvre(ARGS[1] * "." * dats[dnum][1:8], "w")
  for a = 1:4
    for (indlat, lat) in enumerate(lats)
      for (indlon, lon) in enumerate(lons)
        if numb[indlon,indlat,a] > 0
          sums[indlon,indlat,a,:] /= numb[indlon,indlat,a]
          form = @sprintf("%s %11.4f %11.4f %11.4f %11.4f %11.4f %11.4f %11.4f %11.4f %11.4f\n",
            dats[dnum+a], lat, lon,  sums[indlon,indlat,a,1], sums[indlon,indlat,a,2], sums[indlon,indlat,a,3],
            sums[indlon,indlat,a,4], sums[indlon,indlat,a,5], sums[indlon,indlat,a,6], sums[indlon,indlat,a,7])
          write(fpb, form)
        end
      end
    end
  end
  close(fpb)
end

exit(0)

#=
form = @sprintf("%s %11.4f %11.4f %11.4f %11.4f %11.4f %11.4f %11.4f %11.4f %11.4f\n", dat, lat, lon, spd, inc, elo, azo, ela, aza, gai)
       20150603152622    -51.6240    109.3445      4.8000    153.6559     66.2709    -56.8312     66.2296    -56.4421      0.2260
=#
