#=
 = Plot an initial inspection timeseries of parameters on multiple panels
 = - RD August, September 2016.
 =#

using My, NetCDF, Winston

if (argc = length(ARGS)) != 1
  print("\nUsage: jjj $(basename(@__FILE__)) JA2_GPN_2PdP273_092_20151203_155445_20151203_165058.nc\n\n")
  exit(1)
end
altfil = ARGS[1]

altfil == "JA2_GPN_2PdP273_092_20151203_155445_20151203_165058.nc" && (timsta = "2015-12-03-160035" ; timfin = "2015-12-03-160053")
altfil == "JA2_GPN_2PdP273_111_20151204_094249_20151204_103902.nc" && (timsta = "2015-12-04-103300" ; timfin = "2015-12-04-103313")
timbeg = My.datesous("2000-01-01-000000", timsta, "sc")
timend = My.datesous("2000-01-01-000000", timfin, "sc")
timall = ncread(altfil, "time")
timmsk = timbeg .< timall .< timend
timlen = length(timall[timmsk])
print("found $timlen points along the track\n")

timstr = Array(UTF8String, timlen)                                            # parse the dates as string labels
for a = 1:timlen
  timstr[a] = My.dateadd("2000-01-01-000000", round(timall[timmsk][a] * 1000) / 1000, "sc")
end

pars = ["lat", "lon", "rad_distance_to_land", "swh_ku", "swh_c", "wind_speed_alt", "range_ku", "range_c", "alt"]
parn = length(pars)
data = Array(Float64, timlen, parn)

for (a, para) in enumerate(pars)                                              # read all data and rescale if needed
  vara = ncread(  altfil, para)
  scla = ncgetatt(altfil, para, "scale_factor")
  typeof(scla) != Void && (vara = vara * scla)
  data[:,a] = vara[timmsk]
end

data[:,7] = data[:,9] - data[:,7]
data[:,8] = data[:,9] - data[:,8]
fpa = My.ouvre(altfil * ".txt", "w")                                          # echo all data to an ASCII file for
for a = 1:timlen, b = 1:parn                                                  # later geospatial plots
  form  = @sprintf("%15.8lf", data[a,b])
  form *= b == parn ? "\n" : " "
  write(fpa, form)
# print(form)
end
close(fpa)

ppp = Winston.Table(parn-1,1) ; setattr(ppp, "cellpadding", -0.5)             # and make timeseries plot
for (a, para) in enumerate(pars)
  if a < parn
#   print("adding $para\n")
    tpos = (a, 1)
    tmp = Winston.FramedPlot()
    ppp[tpos...] = Winston.add(tmp)
    tmp = Winston.Curve(collect(1:timlen), data[:,a])
          Winston.add(ppp[tpos...], tmp)
    tmp = Winston.PlotLabel(.53, .92, para, "size", 3.4)
          Winston.add(ppp[tpos...], tmp)
  end
end

xyzzy = ARGS[1]*".png"
print("writing $xyzzy\n")
Winston.savefig(ppp, xyzzy, "width", 1000, "height", 1700)
exit(0)
