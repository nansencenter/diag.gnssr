#=
 = Loop through all locations of interest and assemble valid collocations
 = from the in situ and analysis subdirs - RD May, November 2016.
 =#

using My
const WSPD             = 1                              # indecies of all data variables
const UWND             = 2
const VWND             = 3
const PARS             = 3

const OBS              = 1                              # indecies of data source
const BEF              = 2                              # (OBS = GNSS-R and BEF,NOW,AFT = ERA Interim)
const NOW              = 3
const AFT              = 4
const SRCS             = 4

const TIMS             = 1336                           # number in timeseries
const MISS             = -9999.0                        # generic missing value

if (argc = length(ARGS)) != 1
  print("\nUsage: jjj $(basename(@__FILE__)) gnss.tds1.txt.all.locate_1.0_calib_obs\n\n")
  exit(1)
end

function read_nth_line(fn::AbstractString, ln::Int64)
  stream = open(fn, "r") 
  for i = 1:ln-1  readline(stream)  end
  line =          readline(stream)
  close(stream)
  return line
end

fpa = My.ouvre(ARGS[1],           "r")
fpb = My.ouvre(ARGS[1] * ".comb", "w")

for line in eachline(fpa)                                                     # loop through the insitu locations
  data = fill(MISS, PARS, SRCS)                                               # and append analysis to insitu data
  vals = split(line)
  dat            =       vals[1]  ; datind = Int(4 * My.datesous("20150531180000", dat, "dy"))
  lat            = float(vals[2])
  lon            = float(vals[3])
  data[WSPD,OBS] = float(vals[4])
  out = @sprintf("%s %9.3f %9.3f", dat, lat, lon)
  tmp = @sprintf("%9.3f.%9.3f", lat, lon) ; tail = replace(tmp, " ", ".")

  tmp = split(read_nth_line("ecmwf/ecmwf.$tail.bef", datind))
  data[WSPD,BEF] = float(tmp[4]) ; data[UWND,BEF] = float(tmp[5]) ; data[VWND,BEF] = float(tmp[6])
  newdat = tmp[1] ; if dat != newdat  println("ERROR : $dat != $newdat") ; exit(-1)  end
  tmp = split(read_nth_line("ecmwf/ecmwf.$tail",     datind))
  data[WSPD,NOW] = float(tmp[4]) ; data[UWND,NOW] = float(tmp[5]) ; data[VWND,NOW] = float(tmp[6])
  newdat = tmp[1] ; if dat != newdat  println("ERROR : $dat != $newdat") ; exit(-1)  end
  tmp = split(read_nth_line("ecmwf/ecmwf.$tail.aft", datind))
  data[WSPD,AFT] = float(tmp[4]) ; data[UWND,AFT] = float(tmp[5]) ; data[VWND,AFT] = float(tmp[6])
  newdat = tmp[1] ; if dat != newdat  println("ERROR : $dat != $newdat") ; exit(-1)  end

  if data[WSPD,OBS] > -333.0 && data[WSPD,BEF] > -333.0 && data[WSPD,NOW] > -333.0 && data[WSPD,AFT] > -333.0 &&
     data[WSPD,OBS] <  333.0 && data[WSPD,BEF] <  333.0 && data[WSPD,NOW] <  333.0 && data[WSPD,AFT] <  333.0
    data[UWND,OBS] = data[WSPD,OBS] * data[UWND,NOW] / data[WSPD,NOW]
    data[VWND,OBS] = data[WSPD,OBS] * data[VWND,NOW] / data[WSPD,NOW]
    tmp = @sprintf(" %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f",
          data[WSPD,OBS], data[UWND,OBS], data[VWND,OBS], data[WSPD,BEF], data[UWND,BEF], data[VWND,BEF],
          data[WSPD,NOW], data[UWND,NOW], data[VWND,NOW], data[WSPD,AFT], data[UWND,AFT], data[VWND,AFT])
    out *= tmp ; out *= line[52:end] ; write(fpb, out)
  end
end

close(fpa)
close(fpb)
exit(0)

#=
all/gnss.tds1.txt.all.locate_1.0_calib_obs
20150604000000    -54.8750     40.5000      8.9100    174.1286     84.6298    -91.5800     84.5951    -91.3217     12.4745
form = @sprintf("%s %11.4f %11.4f %11.4f %11.4f %11.4f %11.4f %11.4f %11.4f %11.4f\n", dat, lat, lon, spd, inc, elo, azo, ela, aza, gai)
ecmwf/ecmwf...-54.500...337.750.bef
20150604000000    -54.5000    337.7500      5.5934  -9999.0000  -9999.0000  -9999.0000  -9999.0000  -9999.0000  -9999.0000
gnss.tds1.txt.all.locate_1.0_calib_obs.comb
dat, lat, lon, obsw, obsu, obsv, befw, befu, befv, noww, nowu, nowv, aftw, aftu, aftv, inc, elo, azo, ela, aza, gai
=#
