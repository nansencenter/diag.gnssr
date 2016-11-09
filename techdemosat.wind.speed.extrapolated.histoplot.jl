#=
 = Loop through the available collocations and plot the corresponding forward and
 = backward extrapolations relative to the actual (uninterpolated) values.  Note that
 = BEF refers to an extrapolation using analysis data from before the target value and
 = AFT refers to an extrapolation using data from afterward - RD June 2016.
 =#

using My, Winston
const ODAT             = 1                              # identify indecies of the input data:
const OLAT             = 2                              # date/lat/lon on the collocation grid
const OLON             = 3
const OCUR             = 4                              # then five buoy parameters

const MINAVG           = 10                             # minimum number of samples for an average
const MISS             = -9999.0                        # generic missing value

if (argc = length(ARGS)) != 2
  print("\nUsage: jjj $(basename(@__FILE__)) gnss.tds1.txt.all.locate_1.0_extra_obs.comb uwnd\n\n")
  exit(1)
end

step = 0.1 ; bound = collect(-40.0:step:40.0)

#=
grid = zeros(length(bound), length(bound), 2)                                 # read the before and after grids
orid = zeros(length(bound), length(bound), 2)
fname = ARGS[1] * "." * ARGS[2] * ".extra.dat"
fpa = My.ouvre(fname, "r")
for (a, vala) in enumerate(bound)
  for (b, valb) in enumerate(bound)
    line = readline(fpa)
    (grid[b,a,1], grid[b,a,2], orid[b,a,1], orid[b,a,2]) = float(split(line))
  end
end
close(fpa)

gric = zeros(length(bound), length(bound), 2)                                 # and the calibrated before and after grids
oric = zeros(length(bound), length(bound), 2)
fname = ARGS[1] * "." * ARGS[2] * ".extra.dau"
fpa = My.ouvre(fname, "r")
for (a, vala) in enumerate(bound)
  for (b, valb) in enumerate(bound)
    line = readline(fpa)
    (gric[b,a,1], gric[b,a,2], oric[b,a,1], oric[b,a,2]) = float(split(line))
  end
end
close(fpa)

sumb = zeros(length(bound)) ; meanb = zeros(length(bound))                    # as well as the corresponding count and sum
suma = zeros(length(bound)) ; meana = zeros(length(bound))                    # of extrapolation values in each TOTN interval
oumb = zeros(length(bound)) ; oeanb = zeros(length(bound))
ouma = zeros(length(bound)) ; oeana = zeros(length(bound))
fname = ARGS[1] * "." * ARGS[2] * ".extra.sum"
fpa = My.ouvre(fname, "r")
for (a, vala) in enumerate(bound)
  line = readline(fpa)
  (sumb[a], meanb[a], suma[a], meana[a], oumb[a], oeanb[a], ouma[a], oeana[a]) = float(split(line))
end
close(fpa)

sucb = zeros(length(bound)) ; meacb = zeros(length(bound))                    # and the corresponding calibrated values
suca = zeros(length(bound)) ; meaca = zeros(length(bound))
oucb = zeros(length(bound)) ; oeacb = zeros(length(bound))
ouca = zeros(length(bound)) ; oeaca = zeros(length(bound))
fname = ARGS[1] * "." * ARGS[2] * ".extra.sun"
fpa = My.ouvre(fname, "r")
for (a, vala) in enumerate(bound)
  line = readline(fpa)
  (sucb[a], meacb[a], suca[a], meaca[a], oucb[a], oeacb[a], ouca[a], oeaca[a]) = float(split(line))
end
close(fpa)

tumb = zeros(length(bound)) ; teanb = zeros(length(bound))                    # and allow a smoothing of these bin averages
tuma = zeros(length(bound)) ; teana = zeros(length(bound))
for (a, vala) in enumerate(bound)
  if a == 1 || a == length(bound)
    tumb[a] = sumb[a] ; teanb[a] = meanb[a]
    tuma[a] = suma[a] ; teana[a] = meana[a]
  else
    tumb[a] = sumb[a-1] + sumb[a] + sumb[a+1] ; teanb[a] = meanb[a-1] + meanb[a] + meanb[a+1]
    tuma[a] = suma[a-1] + suma[a] + suma[a+1] ; teana[a] = meana[a-1] + meana[a] + meana[a+1]
  end
end
maskb = tumb .>= MINAVG ; avgb = teanb[maskb] ./ tumb[maskb] ; bbndb = bound[maskb]
maska = tuma .>= MINAVG ; avga = teana[maska] ./ tuma[maska] ; bbnda = bound[maska]

tumb = zeros(length(bound)) ; teanb = zeros(length(bound))
tuma = zeros(length(bound)) ; teana = zeros(length(bound))
for (a, vala) in enumerate(bound)
  if a == 1 || a == length(bound)
    tumb[a] = oumb[a] ; teanb[a] = oeanb[a]
    tuma[a] = ouma[a] ; teana[a] = oeana[a]
  else
    tumb[a] = oumb[a-1] + oumb[a] + oumb[a+1] ; teanb[a] = oeanb[a-1] + oeanb[a] + oeanb[a+1]
    tuma[a] = ouma[a-1] + ouma[a] + ouma[a+1] ; teana[a] = oeana[a-1] + oeana[a] + oeana[a+1]
  end
end
maskb = tumb .>= MINAVG ; ovgb = teanb[maskb] ./ tumb[maskb] ; obndb = bound[maskb]
maska = tuma .>= MINAVG ; ovga = teana[maska] ./ tuma[maska] ; obnda = bound[maska]

tumb = zeros(length(bound)) ; teanb = zeros(length(bound))                    # and the corresponding calibrated bin averages
tuma = zeros(length(bound)) ; teana = zeros(length(bound))
for (a, vala) in enumerate(bound)
  if a == 1 || a == length(bound)
    tumb[a] = sucb[a] ; teanb[a] = meacb[a]
    tuma[a] = suca[a] ; teana[a] = meaca[a]
  else
    tumb[a] = sucb[a-1] + sucb[a] + sucb[a+1] ; teanb[a] = meacb[a-1] + meacb[a] + meacb[a+1]
    tuma[a] = suca[a-1] + suca[a] + suca[a+1] ; teana[a] = meaca[a-1] + meaca[a] + meaca[a+1]
  end
end
maskb = tumb .>= MINAVG ; avcb  = teanb[maskb] ./ tumb[maskb] ; bbncb = bound[maskb]
maska = tuma .>= MINAVG ; avca  = teana[maska] ./ tuma[maska] ; bbnca = bound[maska]

tumb = zeros(length(bound)) ; teanb = zeros(length(bound))
tuma = zeros(length(bound)) ; teana = zeros(length(bound))
for (a, vala) in enumerate(bound)
  if a == 1 || a == length(bound)
    tumb[a] = oucb[a] ; teanb[a] = oeacb[a]
    tuma[a] = ouca[a] ; teana[a] = oeaca[a]
  else
    tumb[a] = oucb[a-1] + oucb[a] + oucb[a+1] ; teanb[a] = oeacb[a-1] + oeacb[a] + oeacb[a+1]
    tuma[a] = ouca[a-1] + ouca[a] + ouca[a+1] ; teana[a] = oeaca[a-1] + oeaca[a] + oeaca[a+1]
  end
end
maskb = tumb .>= MINAVG ; ovcb  = teanb[maskb] ./ tumb[maskb] ; obncb = bound[maskb]
maska = tuma .>= MINAVG ; ovca  = teana[maska] ./ tuma[maska] ; obnca = bound[maska]
=#
fname = ARGS[1] * ".extra.reg"                                                # finally read the regression coefficient pairs
fpa = My.ouvre(fname, "r")
line = readline(fpa)
(intb, slob, inta, sloa, ontb, olob, onta, oloa) = float(split(line))
close(fpa)

fname = ARGS[1] * ".extra.reh"                                                # and the corresponding calibrated pairs
fpa = My.ouvre(fname, "r")
line = readline(fpa)
(incb, slcb, inca, slca, oncb, olcb, onca, olca) = float(split(line))
close(fpa)

varname = "Wind Speed (ms^{-1})"                                              # define the plot title and limits
analysis = "ERA Interim"
plotitle = analysis * " " * varname

function point(bound::Array{Float64,1}, grid::Array{Float64,3}, plotind::Int64)
  xpts = Array(Float64, 0)
  ypts = Array(Float64, 0)
  zpts = Array(Float64, 0)
  for (a, vala) in enumerate(bound)
    for (b, valb) in enumerate(bound)
      if grid[b,a,plotind] > 0
        push!(xpts, vala)
        push!(ypts, valb)
        push!(zpts, float(grid[b,a,plotind]))
      end
    end
  end
  return(xpts, ypts, zpts)
end

xmin = -5 ; xmax = 45 ; ymin = -5 ; ymax = 45
#=
contains(ARGS[1], "shfx") && (xmin = -300 ; xmax =  900 ; ymin =  -600 ; ymax = 1200)
contains(ARGS[1], "lhfx") && (xmin = -100 ; xmax = 1200 ; ymin = -1000 ; ymax = 2000)
contains(ARGS[1], "wspd") && (xmin =   -5 ; xmax =   45 ; ymin =   -35 ; ymax =   70)
contains(ARGS[1], "airt") && (xmin =  -20 ; xmax =   40 ; ymin =   -25 ; ymax =   45)
contains(ARGS[1], "sstt") && (xmin =  -10 ; xmax =   40 ; ymin =   -10 ; ymax =   40)
contains(ARGS[1], "shum") && (xmin =    0 ; xmax =   30 ; ymin =   -15 ; ymax =   40)
=#

ppp = Winston.Table(2,4) #; setattr(ppp, "cellpadding", -0.5)                  # and then create the before and after
for z = 1:2                                                                   # scatterplots, respectively (make xpts
# (xpts, ypts, zpts) = point(bound, grid, z)                                  # and ypts refer to grid midpoints)
#  xpts += 0.5 * step
#  ypts += 0.5 * step

  z == 1 && (varname = "Extrapolation from before vs\n" * plotitle)
  z == 2 && (varname = "Extrapolation from after vs\n" * plotitle)
  tmp = Winston.FramedPlot(xrange = (xmin,xmax), yrange = (ymin,ymax)) #, title = varname)
# setattr(tmp.x1, "tickdir",           -1) ; setattr(tmp.x2, "tickdir",           -1) ; setattr(tmp.y1, "tickdir",            1) ; setattr(tmp.y2, "tickdir",           -1)
# setattr(tmp.x1, "ticklabels_offset",  0) ; setattr(tmp.x2, "ticklabels_offset", -7) ; setattr(tmp.y1, "ticklabels_offset",  0) ; setattr(tmp.y2, "ticklabels_offset", -7)
# setattr(tmp.x1, "draw_subticks",  false) ; setattr(tmp.x2, "draw_subticks",   true) ; setattr(tmp.y1, "draw_subticks",  false) ; setattr(tmp.y1, "draw_subticks",   true)
#           setattr(tmp.x1, :ticklabels_style, Dict{Symbol, Any}(:color => "antiquewhite"))
# z == 1 && setattr(tmp.y1, :ticklabels_style, Dict{Symbol, Any}(:color => "black"))
# z == 2 && setattr(tmp.y1, :ticklabels_style, Dict{Symbol, Any}(:color => "antiquewhite"))
#           setattr(tmp.x1, "draw_subticks", false)
# z == 1 && setattr(tmp.y1, "draw_subticks",  true)
# z == 2 && setattr(tmp.y1, "draw_subticks", false)
  setattr(tmp.x1, "draw_subticks", false) ; setattr(tmp.x2, "draw_subticks", false)
  setattr(tmp.y1, "draw_subticks", false) ; setattr(tmp.y2, "draw_subticks", false)
  setattr(tmp.x1, :ticklabels_style, Dict{Symbol, Any}(:color => "white"))
  setattr(tmp.y1, :ticklabels_style, Dict{Symbol, Any}(:color => "white"))
  ppp[1,z] = Winston.add(tmp)
# if z == 2
# tmp = Winston.PlotLabel(0.50, 0.90, plotitle, "texthalign", "center", "size", 2.0)
#       Winston.add(ppp[1,z], tmp)
# end

  z == 1 && (inter = intb ; slope = slob)
  z == 2 && (inter = inta ; slope = sloa)
  tmp = Winston.Slope(slope, (0, inter), kind = "solid", "linewidth", 1, "color", parse(Winston.Colorant, "green"))
        Winston.add(ppp[1,z], tmp)
  tmp = Winston.Slope(     1, (0,    0), kind = "solid")
        Winston.add(ppp[1,z], tmp)
#=
  cols = ["red",  "blue", "green", "orange"]
  lims = [    1,      10,     100,     1000]
  for (a, color) in enumerate(cols)
    mask = zpts .>= lims[a]
    tmp = Winston.Points(xpts[mask], ypts[mask], kind = "filled circle", "color", parse(Winston.Colorant, cols[a]), symbolsize = 0.1)
          Winston.add(ppp[1,z], tmp)
    if z == 1
      tmp = Winston.Curve(bbndb, avgb, kind = "solid")
            Winston.add(ppp[1,z], tmp)
#     tmp = Winston.PlotLabel(0.08, 1.00 - a * 0.07, "<span foreground=\"$(cols[length(cols) - a + 1])\">\\geq $(lims[length(cols) - a + 1])</span>", "texthalign", "left", "size", 3.0)
#           Winston.add(ppp[1,z], tmp)
    end
    if z == 2
      tmp = Winston.Curve(bbnda, avga, kind = "solid")
            Winston.add(ppp[1,z], tmp)
    end
  end
=#
end

for z = 1:2                                                                   # scatterplots, respectively (make xpts
# (xpts, ypts, zpts) = point(bound, gric, z)                                  # and ypts refer to grid midpoints)
#  xpts += 0.5 * step
#  ypts += 0.5 * step

  z == 1 && (varname = "Calibrated extrapolation from\nbefore vs " * analysis)
  z == 2 && (varname = "Calibrated extrapolation from\nafter vs " * analysis)
  tmp = Winston.FramedPlot(xrange = (xmin,xmax), yrange = (ymin,ymax)) #, title = varname)
# setattr(tmp.x1, "tickdir",           -1) ; setattr(tmp.x2, "tickdir",           -1) ; setattr(tmp.y1, "tickdir",            1) ; setattr(tmp.y2, "tickdir",           -1)
# setattr(tmp.x1, "ticklabels_offset",  0) ; setattr(tmp.x2, "ticklabels_offset", -7) ; setattr(tmp.y1, "ticklabels_offset",  0) ; setattr(tmp.y2, "ticklabels_offset", -7)
# setattr(tmp.x1, "draw_subticks",  false) ; setattr(tmp.x2, "draw_subticks",   true) ; setattr(tmp.y1, "draw_subticks",  false) ; setattr(tmp.y1, "draw_subticks",   true)
  setattr(tmp.x1, "draw_subticks", false) ; setattr(tmp.x2, "draw_subticks", false)
  setattr(tmp.y1, "draw_subticks", false) ; setattr(tmp.y2, "draw_subticks", false)
  setattr(tmp.x1, :ticklabels_style, Dict{Symbol, Any}(:color => "white"))
  setattr(tmp.y1, :ticklabels_style, Dict{Symbol, Any}(:color => "white"))
#           setattr(tmp.x1, :ticklabels_style, Dict{Symbol, Any}(:color => "black"))
# z == 1 && setattr(tmp.y1, :ticklabels_style, Dict{Symbol, Any}(:color => "black"))
# z == 2 && setattr(tmp.y1, :ticklabels_style, Dict{Symbol, Any}(:color => "white"))
#           setattr(tmp.x1, "draw_subticks", false)
#           setattr(tmp.y1, "draw_subticks",  true)
  ppp[2,z] = Winston.add(tmp)
# tmp = Winston.PlotLabel(0.50, 0.90, plotitle, "texthalign", "center", "size", 2.0)
#       Winston.add(ppp[2,z], tmp)

  z == 1 && (inter = incb ; slope = slcb)
  z == 2 && (inter = inca ; slope = slca)
  tmp = Winston.Slope(slope, (0, inter), kind = "solid", "linewidth", 1, "color", parse(Winston.Colorant, "green"))
        Winston.add(ppp[2,z], tmp)
  tmp = Winston.Slope(     1, (0,    0), kind = "solid")
        Winston.add(ppp[2,z], tmp)
#=
  cols = ["red",  "blue", "green", "orange"]
  lims = [    1,      10,     100,     1000]
  for (a, color) in enumerate(cols)
    mask = zpts .>= lims[a]
    tmp = Winston.Points(xpts[mask], ypts[mask], kind = "filled circle", "color", parse(Winston.Colorant, cols[a]), symbolsize = 0.1)
          Winston.add(ppp[2,z], tmp)
    if z == 1
      tmp = Winston.Curve(bbncb, avcb, kind = "solid")
            Winston.add(ppp[2,z], tmp)
    end
    if z == 2
      tmp = Winston.Curve(bbnca, avca, kind = "solid")
            Winston.add(ppp[2,z], tmp)
    end
  end
=#
end

for z = 1:2                                                                   # scatterplots, respectively (make xpts
# (xpts, ypts, zpts) = point(bound, orid, z)                                  # and ypts refer to grid midpoints)
#  xpts += 0.5 * step
#  ypts += 0.5 * step

  z == 1 && (varname = "Extrapolation from\nbefore vs obs")
  z == 2 && (varname = "Extrapolation from\nafter vs obs")
  tmp = Winston.FramedPlot(xrange = (xmin,xmax), yrange = (ymin,ymax)) #, title = varname)
  setattr(tmp.x1, "draw_subticks", false) ; setattr(tmp.x2, "draw_subticks", false)
  setattr(tmp.y1, "draw_subticks", false) ; setattr(tmp.y2, "draw_subticks", false)
  setattr(tmp.x1, :ticklabels_style, Dict{Symbol, Any}(:color => "white"))
  setattr(tmp.y1, :ticklabels_style, Dict{Symbol, Any}(:color => "white"))
  ppp[1,z+2] = Winston.add(tmp)
# tmp = Winston.PlotLabel(0.50, 0.90, plotitle, "texthalign", "center", "size", 2.0)
#       Winston.add(ppp[1,z+2], tmp)

  z == 1 && (inter = ontb ; slope = olob)
  z == 2 && (inter = onta ; slope = oloa)
  tmp = Winston.Slope(slope, (0, inter), kind = "solid", "linewidth", 1, "color", parse(Winston.Colorant, "green"))
        Winston.add(ppp[1,z+2], tmp)
  tmp = Winston.Slope(     1, (0,    0), kind = "solid")
        Winston.add(ppp[1,z+2], tmp)
#=
  cols = ["red",  "blue", "green", "orange"]
  lims = [    1,      10,     100,     1000]
  for (a, color) in enumerate(cols)
    mask = zpts .>= lims[a]
    tmp = Winston.Points(xpts[mask], ypts[mask], kind = "filled circle", "color", parse(Winston.Colorant, cols[a]), symbolsize = 0.1)
          Winston.add(ppp[1,z+2], tmp)
    if z == 1
      tmp = Winston.Curve(obndb, ovgb, kind = "solid") 
            Winston.add(ppp[1,z+2], tmp)
    end
    if z == 2
      tmp = Winston.Curve(obnda, ovga, kind = "solid") 
            Winston.add(ppp[1,z+2], tmp)
    end
  end
=#
end

for z = 1:2                                                                   # scatterplots, respectively (make xpts
# (xpts, ypts, zpts) = point(bound, oric, z)                                  # and ypts refer to grid midpoints)
#  xpts += 0.5 * step
#  ypts += 0.5 * step

  z == 1 && (varname = "Calibrated extrapolation from\nbefore vs obs")
  z == 2 && (varname = "Calibrated extrapolation from\nafter vs obs")
  tmp = Winston.FramedPlot(xrange = (xmin,xmax), yrange = (ymin,ymax)) #, title = varname)
  setattr(tmp.x1, "draw_subticks", false) ; setattr(tmp.x2, "draw_subticks", false)
  setattr(tmp.y1, "draw_subticks", false) ; setattr(tmp.y2, "draw_subticks", false)
  setattr(tmp.x1, :ticklabels_style, Dict{Symbol, Any}(:color => "white"))
  setattr(tmp.y1, :ticklabels_style, Dict{Symbol, Any}(:color => "white"))
  ppp[2,z+2] = Winston.add(tmp)
# tmp = Winston.PlotLabel(0.50, 0.90, plotitle, "texthalign", "center", "size", 2.0)
#       Winston.add(ppp[2,z+2], tmp)

  z == 1 && (inter = oncb ; slope = olcb)
  z == 2 && (inter = onca ; slope = olca)
  tmp = Winston.Slope(slope, (0, inter), kind = "solid", "linewidth", 1, "color", parse(Winston.Colorant, "green"))
        Winston.add(ppp[2,z+2], tmp)
  tmp = Winston.Slope(     1, (0,    0), kind = "solid")
        Winston.add(ppp[2,z+2], tmp)
#=
  cols = ["red",  "blue", "green", "orange"]
  lims = [    1,      10,     100,     1000]
  for (a, color) in enumerate(cols)
    mask = zpts .>= lims[a]
    tmp = Winston.Points(xpts[mask], ypts[mask], kind = "filled circle", "color", parse(Winston.Colorant, cols[a]), symbolsize = 0.1)
          Winston.add(ppp[2,z+2], tmp)
    if z == 1
      tmp = Winston.Curve(obncb, ovcb, kind = "solid")
            Winston.add(ppp[2,z+2], tmp)
    end
    if z == 2
      tmp = Winston.Curve(obnca, ovca, kind = "solid")
            Winston.add(ppp[2,z+2], tmp)
    end
  end
=#
end

xyzzy = ARGS[1] * ".extra.png"
print("writing $xyzzy\n")
Winston.savefig(ppp, xyzzy, "width", 1700, "height", 1000)
exit(0)


#=
    ump = Array(Any, length(cols))
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

          title="$varname Spectra (dB)", xlog = true,
          xlabel="Timescale (days)", xrange = (1/1000,1/2), yrange = (ymin,ymax))
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
=#
