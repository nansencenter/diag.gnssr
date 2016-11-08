#=
 = Loop through all locations of interest and assemble valid collocations
 = from the in situ and analysis subdirs - RD May, November 2016.
 =#

using My
const WSPD             = 4                              # indecies of all data variables
const MISS             = -9999.0                        # generic missing value

if (argc = length(ARGS)) != 1
  print("\nUsage: jjj $(basename(@__FILE__)) gnss.tds1.txt.all.locate_1.0_calib_obs\n\n")
  exit(1)
end

wind = WSPD                                                                   # identify the output variable
dirs = ["ecmwf"]
dirn = length(dirs)

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
  vals = split(line)
  dat =       vals[   1]  ; datind = Int(4 * My.datesous("20150531180000", dat, "dy"))
  lat = float(vals[   2])
  lon = float(vals[   3])
  cur = float(vals[wind])
  out = @sprintf("%s %9.3f %9.3f %9.3f", dat, lat, lon, cur)
  tmp = @sprintf("%9.3f.%9.3f", lat, lon) ; tail = replace(tmp, " ", ".")

  bef = fill(MISS, dirn)                                                      # add analysis bef/now/aft to insitu data
  now = fill(MISS, dirn)
  aft = fill(MISS, dirn)
  flag = true
  for (a, dira) in enumerate(dirs)
    tmp = split(read_nth_line("$dira/$dira.$tail.bef", datind)) ; bef[a] = float(tmp[wind])
    newdat = tmp[1] ; if dat != newdat  println("ERROR : $dat != $newdat") ; exit(-1)  end
    tmp = split(read_nth_line("$dira/$dira.$tail",     datind)) ; now[a] = float(tmp[wind])
    newdat = tmp[1] ; if dat != newdat  println("ERROR : $dat != $newdat") ; exit(-1)  end
    tmp = split(read_nth_line("$dira/$dira.$tail.aft", datind)) ; aft[a] = float(tmp[wind])
    newdat = tmp[1] ; if dat != newdat  println("ERROR : $dat != $newdat") ; exit(-1)  end
    if bef[a] < -333.0 || bef[a] > 333.0 || now[a] < -333.0 || now[a] > 333.0 ||
       aft[a] < -333.0 || aft[a] > 333.0  flag = false  end
  end

  if flag                                                                     # and store the line if all values exist
    for (a, dira) in enumerate(dirs)
      tmp = @sprintf(" %9.3f %9.3f %9.3f", bef[a], now[a], aft[a]) ; out *= tmp
    end
    out *= line[52:end] ; write(fpb, out)
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
dat, lat, lon, spd, bef, now, aft, inc, elo, azo, ela, aza, gai
=#
