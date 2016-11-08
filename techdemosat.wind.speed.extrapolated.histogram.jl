#=
 = Loop through the available collocations and grid the corresponding forward and
 = backward extrapolations relative to the actual (uninterpolated) values.  Note that
 = BEF refers to an extrapolation using analysis data from before the target value and
 = AFT refers to an extrapolation using data from afterward.  Relative to the values at
 = the extrapolation target time (TOTN), both local (binwise) and global regressions of
 = the two separate extrapolations (from before and after) are also saved.  The global
 = regression is then applied to the data and the same files are stored - RD June,
 = November 2016.
 =#

using My
const ODAT             = 1                              # identify indecies of the input data:
const OLAT             = 2                              # date/lat/lon on the collocation grid
const OLON             = 3
const OCUR             = 4                              # plus the GNSS-R retrieval and three
const TOTB             = 5                              # ECMWF Interim indecies (b/n/a)
const TOTN             = 6
const TOTA             = 7
const RANGE            = 2.5:0.5:17.5

const CUTOFF           = 400                            # number of collocations in a sample
const MISS             = -9999.0                        # generic missing value

if (argc = length(ARGS)) != 2
  print("\nUsage: jjj $(basename(@__FILE__)) gnss.tds1.txt.all.locate_1.0_extra_obs.comb\n\n")
  exit(1)
end

step = 0.1 ; bound = collect(-40.0:step:80.0)
gridb = zeros(length(bound), length(bound)) ; meanb = zeros(length(bound))
grida = zeros(length(bound), length(bound)) ; meana = zeros(length(bound))
oridb = zeros(length(bound), length(bound)) ; oeanb = zeros(length(bound))
orida = zeros(length(bound), length(bound)) ; oeana = zeros(length(bound))
regb = Array(Float64, 0)
regn = Array(Float64, 0)
rega = Array(Float64, 0)
oegb = Array(Float64, 0)
oego = Array(Float64, 0)
oega = Array(Float64, 0)

fpa = My.ouvre(ARGS[1], "r") ; tinea = readlines(fpa) ; close(fpa)            # grid the collocations and save values
tinuma = length(tinea)                                                        # for a regression versus TOTN and OCUR
for a = 1:tinuma
  vals = float(split(tinea[a]))
  if vals[TOTB] > -333 && vals[TOTB] < 333 &&
     vals[TOTN] > -333 && vals[TOTN] < 333 &&
     vals[TOTA] > -333 && vals[TOTA] < 333 &&
     vals[OCUR] > -333 && vals[OCUR] < 333
    delbef, indbef = findmin(abs(bound - vals[TOTB])) ; bound[indbef] > vals[TOTB] && indbef > 1 && (indbef -= 1)
    delnow, indnow = findmin(abs(bound - vals[TOTN])) ; bound[indnow] > vals[TOTN] && indnow > 1 && (indnow -= 1)
    delaft, indaft = findmin(abs(bound - vals[TOTA])) ; bound[indaft] > vals[TOTA] && indaft > 1 && (indaft -= 1)
    delobs, indobs = findmin(abs(bound - vals[OCUR])) ; bound[indobs] > vals[OCUR] && indobs > 1 && (indobs -= 1)
    gridb[indbef,indnow] += 1 ; meanb[indnow] += vals[TOTB]
    grida[indaft,indnow] += 1 ; meana[indnow] += vals[TOTA]
    oridb[indbef,indobs] += 1 ; oeanb[indobs] += vals[TOTB]
    orida[indaft,indobs] += 1 ; oeana[indobs] += vals[TOTA]
    push!(regb, vals[TOTB])
    push!(regn, vals[TOTN])
    push!(rega, vals[TOTA])
    push!(oegb, vals[TOTB])
    push!(oego, vals[OCUR])
    push!(oega, vals[TOTA])
  end
end

fname = ARGS[1] * ".extra.dat"
fpa = My.ouvre(fname, "w")                                                    # and save the grids
for (a, vala) in enumerate(bound)
  for (b, valb) in enumerate(bound)
    @printf(fpa, "%15.8f %15.8f %15.8f %15.8f\n", gridb[b,a], grida[b,a], oridb[b,a], orida[b,a])
  end
end
close(fpa)

sumb = sum(gridb, 1)                                                          # as well as the corresponding count and sum
suma = sum(grida, 1)                                                          # of extrapolation values in each TOTN interval
oumb = sum(oridb, 1)
ouma = sum(orida, 1)
fname = ARGS[1] * ".extra.sum"
fpa = My.ouvre(fname, "w")
for (a, vala) in enumerate(bound)
  @printf(fpa, "%22.0f %33.11f %22.0f %33.11f %22.0f %33.11f %22.0f %33.11f\n",
    sumb[a], meanb[a], suma[a], meana[a], oumb[a], oeanb[a], ouma[a], oeana[a])
end
close(fpa)

fname = ARGS[1] * ".extra.reg"                                                # and finally save the regression coefficients
fpa = My.ouvre(fname, "w")
(intb, slob) = linreg(regn, regb)
(inta, sloa) = linreg(regn, rega)
(ontb, olob) = linreg(oego, oegb)
(onta, oloa) = linreg(oego, oega)
@printf(fpa, "%33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f\n",
  intb, slob, inta, sloa, ontb, olob, onta, oloa)
close(fpa)

gridb = zeros(length(bound), length(bound)) ; meanb = zeros(length(bound))    # reinitialize the variables for calibration
grida = zeros(length(bound), length(bound)) ; meana = zeros(length(bound))
oridb = zeros(length(bound), length(bound)) ; oeanb = zeros(length(bound))
orida = zeros(length(bound), length(bound)) ; oeana = zeros(length(bound))
regb = Array(Float64, 0)
regn = Array(Float64, 0)
rega = Array(Float64, 0)
oegb = Array(Float64, 0)
oego = Array(Float64, 0)
oega = Array(Float64, 0)

for a = 1:tinuma                                                              # now calibrate the same collocations using
  vals = float(split(tinea[a]))                                               # the global regressions above and regrid
  if vals[TOTB] > -333 && vals[TOTB] < 333 &&
     vals[TOTN] > -333 && vals[TOTN] < 333 &&
     vals[TOTA] > -333 && vals[TOTA] < 333 &&
     vals[OCUR] > -333 && vals[OCUR] < 333
    tmpb = (vals[TOTB] - intb) / slob
    tmpn =  vals[TOTN]
    tmpa = (vals[TOTA] - inta) / sloa
    tmpo =  vals[OCUR]
    delbef, indbef = findmin(abs(bound - tmpb)) ; bound[indbef] > tmpb && indbef > 1 && (indbef -= 1)
    delnow, indnow = findmin(abs(bound - tmpn)) ; bound[indnow] > tmpn && indnow > 1 && (indnow -= 1)
    delaft, indaft = findmin(abs(bound - tmpa)) ; bound[indaft] > tmpa && indaft > 1 && (indaft -= 1)
    delobs, indobs = findmin(abs(bound - tmpo)) ; bound[indobs] > tmpo && indobs > 1 && (indobs -= 1)
    gridb[indbef,indnow] += 1 ; meanb[indnow] += tmpb
    grida[indaft,indnow] += 1 ; meana[indnow] += tmpa
    oridb[indbef,indobs] += 1 ; oeanb[indobs] += tmpb
    orida[indaft,indobs] += 1 ; oeana[indobs] += tmpa
    push!(regb, tmpb)
    push!(regn, tmpn)
    push!(rega, tmpa)
    push!(oegb, tmpb)
    push!(oego, tmpo)
    push!(oega, tmpa)
  end
end

fname = ARGS[1] * ".extra.dau"
fpa = My.ouvre(fname, "w")                                                    # and save the grids
for (a, vala) in enumerate(bound)
  for (b, valb) in enumerate(bound)
    @printf(fpa, "%15.8f %15.8f %15.8f %15.8f\n", gridb[b,a], grida[b,a], oridb[b,a], orida[b,a])
  end
end
close(fpa)

sumb = sum(gridb, 1)                                                          # as well as the corresponding count and sum
suma = sum(grida, 1)                                                          # of extrapolation values in each TOTN interval
oumb = sum(oridb, 1)
ouma = sum(orida, 1)
fname = ARGS[1] * ".extra.sun"
fpa = My.ouvre(fname, "w")
for (a, vala) in enumerate(bound)
  @printf(fpa, "%22.0f %33.11f %22.0f %33.11f %22.0f %33.11f %22.0f %33.11f\n",
    sumb[a], meanb[a], suma[a], meana[a], oumb[a], oeanb[a], ouma[a], oeana[a])
end
close(fpa)

fname = ARGS[1] * ".extra.reh"                                                # and finally save the regression coefficients
fpa = My.ouvre(fname, "w")
(intb, slob) = linreg(regn, regb)
(inta, sloa) = linreg(regn, rega)
(ontb, olob) = linreg(oego, oegb)
(onta, oloa) = linreg(oego, oega)
@printf(fpa, "%33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f %33.11f\n",
  intb, slob, inta, sloa, ontb, olob, onta, oloa)
close(fpa)
exit(0)


#=
count(shfs, shfn, data[:,SHFX,BEF], data[:,SHFX,NOW], data[:,SHFX,AFT])
function count(bound::Array{Float64,1}, grid::Array{Float64,3}, bef::Array{Float64,1}, now::Array{Float64,1}, aft::Array{Float64,1})
  for a = 1:dirn
    delbef, indbef = findmin(abs(bound - bef[a])) ; bound[indbef] > bef[a] && indbef > 1 && (indbef -= 1)
    delnow, indnow = findmin(abs(bound - now[a])) ; bound[indnow] > now[a] && indnow > 1 && (indnow -= 1)
    delaft, indaft = findmin(abs(bound - aft[a])) ; bound[indaft] > aft[a] && indaft > 1 && (indaft -= 1)
    grid[indbef,indnow,a] += 1 ; grid[indaft,indnow,a] += 1
  end
end
=#
