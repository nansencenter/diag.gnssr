#=
 = Loop through the six-hourly ECMWF timeseries and create the corresponding forward and
 = backward extrapolated timeseries for the two (u,v) wind components - RD April, November 2016.
 =#

using My, Interpolations
const UWND             = 1                              # indecies of all data variables
const VWND             = 2
const WSPD             = 3
const PARS             = 3

const BEF              = 1                              # indecies of the source of extrapolations
const NOW              = 2
const AFT              = 3
const SRCS             = 3

const EXTRA            = 9                              # number of points used for extrapolation
const TIMS             = 1336                           # number in timeseries
const MISS             = -9999.0                        # generic missing value

if (argc = length(ARGS)) != 2
  print("\nUsage: jjj $(basename(@__FILE__)) ecmwf z.listah\n\n")
  exit(1)
end

inner = div(EXTRA - 1, 2)
outer = div(EXTRA + 1, 2)
dats = Array(String,  TIMS)
data = Array(Float64, TIMS, SRCS, PARS)

fpa = My.ouvre("$(ARGS[1])/$(ARGS[2])", "r")                                  # loop through the list of locations
files = readlines(fpa) ; close(fpa)                                           # and process each timeseries
for fila in files
  fila = strip(fila)
  if isfile("$(ARGS[1])/$fila.bef") && isfile("$(ARGS[1])/$fila.aft")  continue  end

  fpa = My.ouvre("$(ARGS[1])/$fila", "r", false)
  lines = readlines(fpa) ; close(fpa)
  for (a, line) in enumerate(lines)
    vals = split(line)
    data[a,NOW,UWND] = float(vals[5])
    data[a,NOW,VWND] = float(vals[6])
    dats[a]          =       vals[1]
  end

  for a = 1:PARS-1                                                            # set to missing the first few BEF
    for b = 1:EXTRA+1                                                         # extrapolations (not defined below)
      data[b,BEF,a] = MISS
    end
  end

  for a = 1:PARS-1                                                            # simultaneously extrap from BEF and AFT
    for b = 1+outer:TIMS-outer
      tmp = vec(data[b-inner:b+inner,NOW,a])
      if all(-333 .< tmp .< 333)
        tmpmax = maximum(tmp)
        tmpmin = minimum(tmp)
        itp = interpolate(tmp, BSpline(Quadratic(Line())), OnCell())
        tmpbef = itp[10] ; tmpbef > tmpmax && (tmpbef = tmpmax) ; tmpbef < tmpmin && (tmpbef = tmpmin)
        tmpaft = itp[ 0] ; tmpaft > tmpmax && (tmpaft = tmpmax) ; tmpaft < tmpmin && (tmpaft = tmpmin)
        data[b+outer,BEF,a] = tmpbef
        data[b-outer,AFT,a] = tmpaft
      else
        data[b+outer,BEF,a] = data[b-outer,AFT,a] = MISS
      end
    end
  end

  for a = 1:PARS-1                                                            # set to missing the last few AFT
    for b = 0:EXTRA                                                           # extrapolations (not defined above)
      data[TIMS-b,AFT,a] = MISS
    end
  end

  for a = 1:TIMS
    if data[a,BEF,UWND] > MISS && data[a,BEF,VWND] > MISS
      data[a,BEF,WSPD] = (data[a,BEF,UWND]^2 + data[a,BEF,VWND]^2)^0.5
    else
      data[a,BEF,WSPD] =  data[a,BEF,UWND]   = data[a,BEF,VWND] = MISS
    end
    if data[a,AFT,UWND] > MISS && data[a,AFT,VWND] > MISS
      data[a,AFT,WSPD] = (data[a,AFT,UWND]^2 + data[a,AFT,VWND]^2)^0.5
    else
      data[a,AFT,WSPD] =  data[a,AFT,UWND]   = data[a,AFT,VWND] = MISS
    end
  end

  filb = "$fila.bef"                                                          # then save all extrapolations
  filc = "$fila.aft"
  fpb = My.ouvre("$(ARGS[1])/$filb", "w", false)
  fpc = My.ouvre("$(ARGS[1])/$filc", "w", false)
  (lll, lat, lon) = split(replace(fila, r"[\.]{2,}", " "))
  for a = 1:TIMS
    formb = @sprintf("%s %11.4f %11.4f %11.4f %11.4f %11.4f  -9999.0000  -9999.0000  -9999.0000  -9999.0000\n",
      dats[a], float(lat), float(lon), data[a,BEF,WSPD], data[a,BEF,UWND], data[a,BEF,VWND])
    formc = @sprintf("%s %11.4f %11.4f %11.4f %11.4f %11.4f  -9999.0000  -9999.0000  -9999.0000  -9999.0000\n",
      dats[a], float(lat), float(lon), data[a,AFT,WSPD], data[a,AFT,UWND], data[a,AFT,VWND])
    write(fpb, formb)
    write(fpc, formc)
  end
  close(fpb)
  close(fpc)
end
exit(0)

#=
20150604060000    -13.0000    123.8750     10.1032  -9999.0000  -9999.0000  -9999.0000  -9999.0000  -9999.0000  -9999.0000
=#
