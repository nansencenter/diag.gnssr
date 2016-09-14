#=
 = Plot all parameters on multiple panels using a command like
 =
 = jjj plot.gnssr.track.altimetry.jl JA2_GPN_2PdP273_092_20151203_155445_20151203_165058.nc
 = jjj plot.gnssr.track.altimetry.jl JA2_GPN_2PdP273_111_20151204_094249_20151204_103902.nc
 =
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

for (a, para) in enumerate(pars)
  vara = ncread(  altfil, para)
  scla = ncgetatt(altfil, para, "scale_factor")
  typeof(scla) != Void && (vara = vara * scla)
  data[:,a] = vara[timmsk]
end

data[:,7] = data[:,9] - data[:,7]
data[:,8] = data[:,9] - data[:,8]

fpa = My.ouvre(altfil * ".txt", "w")
for a = 1:timlen, b = 1:parn
  form  = @sprintf("%15.8lf", data[a,b])
  form *= b == parn ? "\n" : " "
  write(fpa, form)
  print(form)
end
close(fpa)

ppp = Winston.Table(parn-1,1) ; setattr(ppp, "cellpadding", -0.5)
for (a, para) in enumerate(pars)
  if a < parn
    print("adding $para\n")
    tpos = (a, 1)
    tmp = Winston.FramedPlot() #title=para, xlabel="Time (s)")
#         title="$varname Spectra (dB)", xlog = true,
#         xlabel="Timescale (days)", xrange = (1/1000,1/2), yrange = (ymin,ymax))
#         xlog = true, xrange = (1,timlen)) #, yrange = (ymin,ymax))
# setattr(tmp.x1, "ticks",          xposa) ; setattr(tmp.x2, "ticks",          xposb) ; setattr(tmp.y1, "ticks",          yposa)
# setattr(tmp.x1, "tickdir",            1) ; setattr(tmp.x2, "tickdir",           -1) ; setattr(tmp.y1, "tickdir",            1)
# setattr(tmp.x1, "ticklabels_offset",  0) ; setattr(tmp.x2, "ticklabels_offset", -7) ; setattr(tmp.y1, "ticklabels_offset",  0)
# setattr(tmp.x1, "ticklabels",     xlaba) ; setattr(tmp.x2, "ticklabels",     xlabb) ; setattr(tmp.y1, "ticklabels",     ylaba)
# setattr(tmp.x1, "draw_subticks",  false) ; setattr(tmp.x2, "draw_subticks",   true) ; setattr(tmp.y1, "draw_subticks",   true)
# tpos[1] <= 2 && setattr(tmp.x1, :ticklabels_style, Dict{Symbol, Any}(:color => "transparent"))
# tpos[1] >= 2 && setattr(tmp.x2, :ticklabels_style, Dict{Symbol, Any}(:color => "transparent"))
# tpos[2] == 1 && setattr(tmp.y1, :ticklabels_style, Dict{Symbol, Any}(:color => "black"))
# tpos[2] == 2 && setattr(tmp.y1, :ticklabels_style, Dict{Symbol, Any}(:color => "transparent"))
    ppp[tpos...] = Winston.add(tmp)
    tmp = Winston.Curve(collect(1:timlen), data[:,a])
#       style(tmp, kind = kynd[b])
#       setattr(tmp, label = "($(specstr[a,z])) $(dirs[a])")
          Winston.add(ppp[tpos...], tmp)
    tmp = Winston.PlotLabel(.53, .92, para, "size", 3.4)
          Winston.add(ppp[tpos...], tmp)

# ump = Array(Any, 8)
# cols = [  "red",  "blue", "green", "orange",    "red",   "blue",  "green", "orange"]
# kynd = ["solid", "solid", "solid",  "solid", "dashed", "dashed", "dashed", "dashed"]

# b = 1
# for a in order
#   ump[a] = Winston.Curve(specval[1:end,z], spectra[a,1:end,z], "color", parse(Winston.Colorant, cols[b]))
#            style(ump[a], kind = kynd[b])
#            setattr(ump[a], label = "($(specstr[a,z])) $(dirs[a])")
#            Winston.add(ppp[tpos...], ump[a])
#   b += 1
# end
# tmp = Winston.Legend(.23, .92, Any[ump[order[1]], ump[order[2]], ump[order[3]], ump[order[4]]])
#       Winston.add(ppp[tpos...], tmp)
# tmp = Winston.Legend(.58, .92, Any[ump[order[5]], ump[order[6]], ump[order[7]], ump[order[8]]])
#       Winston.add(ppp[tpos...], tmp)
# tmp = Winston.DataLabel(0.0012, ymin + 0.12 * (ymax - ymin), varname, "texthalign", "left", "size", 1.4)
#       Winston.add(ppp[tpos...], tmp)
  end
end

xyzzy = ARGS[1]*".png"
print("writing $xyzzy\n")
Winston.savefig(ppp, xyzzy, "width", 1000, "height", 1700)
exit(0)


  for a in order
    contains(ARGS[z], "shfx") && (specstr[a,z] = @sprintf("%.0lf", specvar[a,z]))
    contains(ARGS[z], "lhfx") && (specstr[a,z] = @sprintf("%.0lf", specvar[a,z]))
    contains(ARGS[z], "wspd") && (specstr[a,z] = @sprintf("%.1lf", specvar[a,z]))
    contains(ARGS[z], "airt") && (specstr[a,z] = @sprintf("%.1lf", specvar[a,z]))
    contains(ARGS[z], "sstt") && (specstr[a,z] = @sprintf("%.1lf", specvar[a,z]))
    contains(ARGS[z], "shum") && (specstr[a,z] = @sprintf("%.1lf", specvar[a,z]))
    print("$(pars[z]) variance $(specstr[a,z]) in $(dirs[a])\n")
  end

  xlaba = ["1000", "365", "182", "100", "60", "30", "14", "7", "3.5", "2"]
  xlabb = ["0.002", "0.3"]
  xposa = float(xlaba).^-1
  xposb = float(xlabb)

  contains(ARGS[z], "shfx") && (varname = "a) Sensible Heat Flux" ; ymin = -15 ; ymax = 40 ; tpos = (1,1))
  contains(ARGS[z], "lhfx") && (varname = "b) Latent Heat Flux"   ; ymin = -15 ; ymax = 40 ; tpos = (1,2))
  contains(ARGS[z], "wspd") && (varname = "c) Wind Speed"         ; ymin = -30 ; ymax = 15 ; tpos = (2,1))
  contains(ARGS[z], "shum") && (varname = "d) Specific Humidity"  ; ymin = -30 ; ymax = 15 ; tpos = (2,2))
  contains(ARGS[z], "sstt") && (varname = "e) Sea Surface Temp"   ; ymin = -30 ; ymax = 15 ; tpos = (3,1))
  contains(ARGS[z], "airt") && (varname = "f) Air Temperature"    ; ymin = -30 ; ymax = 15 ; tpos = (3,2))
  tpos[1] == 1 && (ylaba = ["-10", "0", "10", "20", "30", "40"])
  tpos[1] != 1 && (ylaba = ["-30", "-20", "-10", "0", "10"])
  yposa = float(ylaba)

  tmp = Winston.FramedPlot(
#         title="$varname Spectra (dB)", xlog = true,
#         xlabel="Timescale (days)", xrange = (1/1000,1/2), yrange = (ymin,ymax))
          xlog = true, xrange = (1/1000,1/2), yrange = (ymin,ymax))
  setattr(tmp.x1, "ticks",          xposa) ; setattr(tmp.x2, "ticks",          xposb) ; setattr(tmp.y1, "ticks",          yposa)
  setattr(tmp.x1, "tickdir",            1) ; setattr(tmp.x2, "tickdir",           -1) ; setattr(tmp.y1, "tickdir",            1)
  setattr(tmp.x1, "ticklabels_offset",  0) ; setattr(tmp.x2, "ticklabels_offset", -7) ; setattr(tmp.y1, "ticklabels_offset",  0)
  setattr(tmp.x1, "ticklabels",     xlaba) ; setattr(tmp.x2, "ticklabels",     xlabb) ; setattr(tmp.y1, "ticklabels",     ylaba)
  setattr(tmp.x1, "draw_subticks",  false) ; setattr(tmp.x2, "draw_subticks",   true) ; setattr(tmp.y1, "draw_subticks",   true)
  tpos[1] <= 2 && setattr(tmp.x1, :ticklabels_style, Dict{Symbol, Any}(:color => "transparent"))
  tpos[1] >= 2 && setattr(tmp.x2, :ticklabels_style, Dict{Symbol, Any}(:color => "transparent"))
  tpos[2] == 1 && setattr(tmp.y1, :ticklabels_style, Dict{Symbol, Any}(:color => "black"))
  tpos[2] == 2 && setattr(tmp.y1, :ticklabels_style, Dict{Symbol, Any}(:color => "transparent"))
  ppp[tpos...] = Winston.add(tmp)

  ump = Array(Any, 8)
  cols = [  "red",  "blue", "green", "orange",    "red",   "blue",  "green", "orange"]
  kynd = ["solid", "solid", "solid",  "solid", "dashed", "dashed", "dashed", "dashed"]

  b = 1
  for a in order
    ump[a] = Winston.Curve(specval[1:end,z], spectra[a,1:end,z], "color", parse(Winston.Colorant, cols[b]))
             style(ump[a], kind = kynd[b])
             setattr(ump[a], label = "($(specstr[a,z])) $(dirs[a])")
             Winston.add(ppp[tpos...], ump[a])
    b += 1
  end
  tmp = Winston.Legend(.23, .92, Any[ump[order[1]], ump[order[2]], ump[order[3]], ump[order[4]]])
        Winston.add(ppp[tpos...], tmp)
  tmp = Winston.Legend(.58, .92, Any[ump[order[5]], ump[order[6]], ump[order[7]], ump[order[8]]])
        Winston.add(ppp[tpos...], tmp)
  tmp = Winston.DataLabel(0.0012, ymin + 0.12 * (ymax - ymin), varname, "texthalign", "left", "size", 1.4)
        Winston.add(ppp[tpos...], tmp)
