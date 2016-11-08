#=
 = Perform a pair of analysis calibration and performance estimates, followed by
 = recalibration using the OLS slope and intercept from the opposing collocations,
 = followed by another calibration and performance estimate.  The important issue
 = of data outlier detection and removal is treated using the DetMCD algorithm in
 = R.  For morphing experiments, the linearity of the bijective mapping into and
 = out of Gaussian-anamorphosis space is also checked (given that the in situ
 = slope and intercept remain the reference) - RD June, July, October, November
 = 2016.
 =#

using My, Morph, RCall ; R"library(DetMCD)"
const ODAT             = 1                              # identify indecies of the input data:
const OLAT             = 2                              # date/lat/lon on the collocation grid
const OLON             = 3
const OCUR             = 4                              # plus the GNSS-R retrieval and three
const TOTB             = 5                              # ECMWF Interim indecies (b/n/a)
const TOTN             = 6
const TOTA             = 7

const MISS             = -9999.0                        # generic missing value
const MORPH            = false                          # perform Gaussian anamorphosis
const EXTRA            = true                           # recalibrate the extrapolated data using extra collocations
const MCDTRIM          = 0.5                            # Minimum Covariance Determinant trimming (nonoutlier percent)

if (argc = length(ARGS)) != 2
  print("\nUsage: jjj $(basename(@__FILE__)) gnss.tds1.txt.all.locate_1.0_calib_obs.comb 1.0\n\n")
  exit(1)
end
rescale = float(ARGS[2])

#=
 = Function returning triple collocation cal/val measures for a group of analyses, following McColl
 = et al. (2014).  Inputs are an array of collocated values and stats are returned for a collocation
 = set, where it is assumed that extrapolation from before and after is done using the same analysis,
 = so no consideration of relative effective resolution is necessary (cf. Vogelzang et al. 2011)
 =#

function triple(curr::Array{Float64,3})
#=mask = masquepourcent(curr[1, :,2], MCDTRIM) &                              # get the parametric center of mass
         masquepourcent(curr[1, :,1], MCDTRIM) &                              # after trimming extreme values first
         masquepourcent(curr[2, :,1], MCDTRIM)
  sampsitu =      curr[1,mask,2]
  samprefa =      curr[1,mask,1]
  samprefb =      curr[2,mask,1]
  mass     = mean(curr[2,mask,2])
# @show length(mask) length(mask[mask])

  mass     = mean(curr[2,:,2])
  sampsitu =      curr[1,:,2] #- mass
  samprefa =      curr[1,:,1] #- mass
  samprefb =      curr[2,:,1] #- mass

  avg1 = mean(sampsitu)                                                       # and use a robust calculation of covariance
  avg2 = mean(samprefa)                                                       # (two-pass here, but more algorithms are at
  avg3 = mean(samprefb)                                                       # en.wikipedia.org/wiki/Algorithms_for_calculating_variance)
  cv11 = mean((sampsitu - avg1) .* (sampsitu - avg1))
  cv12 = mean((sampsitu - avg1) .* (samprefa - avg2))
  cv13 = mean((sampsitu - avg1) .* (samprefb - avg3))
  cv22 = mean((samprefa - avg2) .* (samprefa - avg2))
  cv23 = mean((samprefa - avg2) .* (samprefb - avg3))
  cv33 = mean((samprefb - avg3) .* (samprefb - avg3))
=#
  mass = mean(curr[2,:,2])
# curr[1,:,2] -= mass
# curr[1,:,1] -= mass
# curr[2,:,1] -= mass
  temp = [curr[1,:,2]' curr[1,:,1]' curr[2,:,1]']
  remp = rcopy(R"DetMCD($temp, alpha = $MCDTRIM)")
  mask = falses(length(temp[:,1])) ; for a in remp[:Hsubsets]  mask[a] = true  end
# mass = mean(curr[2,mask,2])

  avg1 = remp[:center][1]
  avg2 = remp[:center][2]
  avg3 = remp[:center][3]
  cv11 = remp[:cov][1,1]
  cv12 = remp[:cov][1,2]
  cv13 = remp[:cov][1,3]
  cv22 = remp[:cov][2,2]
  cv23 = remp[:cov][2,3]
  cv33 = remp[:cov][3,3]

  bet1 = 0.5 * (cv12 + cv13) / cv23
  alp1 = avg1 - bet1 * 0.5 * (avg2 + avg3) #+ mass * (1.0 - bet1)
# bet1 = cv13 / cv23
# alp1 = avg1 - bet1 * avg3
  bet2 = bet3 = 1.0 ; alp2 = alp3 = 0.0

# bet2 = cv23 / cv13
# bet3 = cv23 / cv12
# alp2 = avg2 - bet2 * avg1 + mass * (1.0 - bet2)
# alp3 = avg3 - bet3 * avg1 + mass * (1.0 - bet3)

  tmpval = cv11 - cv12 * cv13 / cv23 ; sig1 = tmpval > 0 ? sqrt(tmpval) : 0.0
  tmpval = cv22 - cv12 * cv23 / cv13 ; sig2 = tmpval > 0 ? sqrt(tmpval) : 0.0
  tmpval = cv33 - cv13 * cv23 / cv12 ; sig3 = tmpval > 0 ? sqrt(tmpval) : 0.0
  tmpval = cv12 * cv13 / cv11 / cv23 ; cor1 = tmpval > 0 ? sqrt(tmpval) : 0.0
  tmpval = cv12 * cv23 / cv22 / cv13 ; cor2 = tmpval > 0 ? sqrt(tmpval) : 0.0
  tmpval = cv13 * cv23 / cv33 / cv12 ; cor3 = tmpval > 0 ? sqrt(tmpval) : 0.0

  return(mass, alp1, bet1, sig1, cor1, alp2, bet2, sig2, cor2, alp3, bet3, sig3, cor3)    # then return all statistics
end

#=
 = main program
 =#

const MEMO             = 1                              # center-of-mass parameter
const MEMB             = 2                              # error model x = ALPH + BETA * truth + error
const MEMA             = 3                              # error model x = ALPH + BETA * truth + error
const MEMS             = 3                              # number of triple collocation members

const MASS             = 1                              # center-of-mass parameter (again...)
const ALPH             = 2                              # error model x = ALPH + BETA * truth + error
const BETA             = 3                              # error model x = ALPH + BETA * truth + error
const SIGM             = 4                              # triple coll RMSE
const CORR             = 5                              # triple coll correlation coefficient
const PARS             = 5                              # number of triple collocation parameters
const CUTOFF           = 170.0                            # gain cutoff

ARGS222 = replace(ARGS[1], "calib", "valid")                                  # read both sets of collocations
fpa = My.ouvre(ARGS[1], "r") ; minea = readlines(fpa) ; close(fpa)            # dat, lat, lon, spd, bef, now,
fpb = My.ouvre(ARGS222, "r") ; mineb = readlines(fpb) ; close(fpb)            # aft, inc, elo, azo, ela, aza, gai
minuma = length(minea) ; print("$(ARGS[1]) has $minuma lines\n")
minumb = length(mineb) ;   print("$ARGS222 has $minumb lines\n")

maska  = Array(Bool, minuma)                                                  # filter data by passing high gain
maskb  = Array(Bool, minumb)
#for a = 1:minuma  vala = float(split(minea[a])) ; maska[a] = true ; vala[8] > CUTOFF && (maska[a] = false)  end
#for a = 1:minumb  valb = float(split(mineb[a])) ; maskb[a] = true ; valb[8] > CUTOFF && (maskb[a] = false)  end
#for a = 1:minuma  vala = float(split(minea[a])) ; maska[a] = false ; abs(vala[4] - vala[6]) < 1 && (maska[a] = true)  end
#for a = 1:minumb  valb = float(split(mineb[a])) ; maskb[a] = false ; abs(valb[4] - valb[6]) < 1 && (maskb[a] = true)  end
maska = fill(true, minuma) ; maskb  = fill(true, minumb)
linea = minea[maska] ; linuma = length(linea) ; cura = zeros(2, linuma, 2) ; print("$(ARGS[1]) has $linuma good lines\n")
lineb = mineb[maskb] ; linumb = length(lineb) ; curb = zeros(2, linumb, 2) ;   print("$ARGS222 has $linumb good lines\n")

if rescale > 0                                                                # rescale the first set of input data values
  print("rescaling $(ARGS[1])\n")
  for a = 1:linuma
    tmpbef = float(linea[a][46:56]) + rescale * randn()
    tmpnow = float(linea[a][56:66])
    tmpaft = float(linea[a][66:76]) + rescale * randn()
    out = @sprintf("%s %9.3f %9.3f %9.3f %s", linea[a][1:44], tmpbef, tmpnow, tmpaft, linea[a][76:end])
    linea[a] = out
  end
  print("rescaling $ARGS222\n")                                               # rescale the second set of input data values
  for a = 1:linumb
    tmpbef = float(lineb[a][46:56]) + rescale * randn()
    tmpnow = float(lineb[a][56:66])                     
    tmpaft = float(lineb[a][66:76]) + rescale * randn()
    out = @sprintf("%s %9.3f %9.3f %9.3f %s", lineb[a][1:44], tmpbef, tmpnow, tmpaft, lineb[a][76:end])
    lineb[a] = out
  end
end

if EXTRA
  fname = replace(ARGS[1], "calib", "extra") * ".extra.reg"                   # also read the regression coefficient pairs
  fpa = My.ouvre(fname, "r") ; line = readline(fpa) ; close(fpa)              # for calibrating the extrapolations relative
  (intb, slob, inta, sloa) = float(split(line))                               # to the extra collocation target (TOTN)
# intz = 0.5 * (intb + inta) ; intb = inta = intz
# sloz = 0.5 * (slob + sloa) ; slob = sloa = sloz
end

refa = Array(Float64, linuma)                                                 # and calculate a pair of reference variables
refb = Array(Float64, linumb)                                                 # (either from observations or from analyses)
for a = 1:linuma
  vala = float(split(linea[a]))
  EXTRA && (vala[TOTB] = (vala[TOTB] - intb) / slob ;
            vala[TOTA] = (vala[TOTA] - inta) / sloa)
# refa[a] = vala[OCUR]
# refa[a] = vala[TOTN]
# refa[a] = 0.9 * vala[OCUR] + 0.1 * vala[TOTN]
  refa[a] = 0.5 * (vala[OCUR] + vala[TOTN])
end
for a = 1:linumb
  vala = float(split(lineb[a]))
  EXTRA && (vala[TOTB] = (vala[TOTB] - intb) / slob ;
            vala[TOTA] = (vala[TOTA] - inta) / sloa)
# refb[a] = vala[OCUR]
# refb[a] = vala[TOTN]
# refb[a] = 0.9 * vala[OCUR] + 0.1 * vala[TOTN]
  refb[a] = 0.5 * (vala[OCUR] + vala[TOTN])
end

statis = [MISS for a = 1:4, b = 1:MEMS, c = 1:PARS]                           # allocate a set of global cal/val arrays
allmas = [MISS for a = 1:4]
allalp = [MISS for a = 1:4]
allbet = [MISS for a = 1:4]
allsig = [MISS for a = 1:4]
allcor = [MISS for a = 1:4]

if MORPH                                                                      # define for both collocation sets a mapping
  tmpa = Array(Float64, linuma)                                               # (Gaussian anamorphosis) based on the in situ
  tmpb = Array(Float64, linumb)                                               # observations (which are taken to be unbiased)
  for a = 1:linuma  vals = float(split(linea[a])) ; tmpa[a] = vals[OCUR]  end # and apply each map to the opposite set
  for a = 1:linumb  vals = float(split(lineb[a])) ; tmpb[a] = vals[OCUR]  end
  shfa = minimum(tmpa) < 0 ? -2.0 * minimum(tmpa) : 0.0
  shfb = minimum(tmpb) < 0 ? -2.0 * minimum(tmpb) : 0.0
  rawa = sort(tmpa) + shfa
  rawb = sort(tmpb) + shfb

  R"library(RGeostats)"
  R"dba  = db.create(aa = c($rawa))"
  R"fia  = anam.fit(dba, name='aa', type='emp')"
  R"ana <- anam.z2y(fia, dba, 'aa')"
  cdfa = unique(rcopy(R"db.extract(ana, names = c('Gaussian.aa'))"))
  R"rm(dba, fia, ana)"
  rawa = unique(rawa - shfa)
  if length(rawa) != length(cdfa)
    print("\nERROR: length(rawa) $(length(rawa)) != $(length(cdfa)) length(cdfa)\n\n")
    exit(1)
  end

  R"dbb  = db.create(bb = c($rawb))"
  R"fib  = anam.fit(dbb, name='bb', type='emp')"
  R"anb <- anam.z2y(fib, dbb, 'bb')"
  cdfb = unique(rcopy(R"db.extract(anb, names = c('Gaussian.bb'))"))
  R"rm(dbb, fib, anb)"
  rawb = unique(rawb - shfb)
  if length(rawb) != length(cdfb)
    print("\nERROR: length(rawb) $(length(rawb)) != $(length(cdfb)) length(cdfb)\n\n")
    exit(1)
  end                                                                         # finally extend and equate the tails

  rawa[  1] > rawb[  1] && (rawa[  1] = rawb[  1]) ; rawb[  1] > rawa[  1] && (rawb[  1] = rawa[  1])
  rawa[end] < rawb[end] && (rawa[end] = rawb[end]) ; rawb[end] < rawa[end] && (rawb[end] = rawa[end])
  cdfa[  1] > cdfb[  1] && (cdfa[  1] = cdfb[  1]) ; cdfb[  1] > cdfa[  1] && (cdfb[  1] = cdfa[  1])
  cdfa[end] < cdfb[end] && (cdfa[end] = cdfb[end]) ; cdfb[end] < cdfa[end] && (cdfb[end] = cdfa[end])
  cdfa[  1] >        -9 && (cdfa[  1] =        -9) ; cdfb[  1] >        -9 && (cdfb[  1] =        -9)
  cdfa[end] <         9 && (cdfa[end] =         9) ; cdfb[end] <         9 && (cdfb[end] =         9)

# tmp      = Winston.FramedPlot(title="Empirical anamorphosis", xlabel="Gaussian Current", ylabel="Actual Current")
# ppp      = Winston.add(tmp)
# tmp      = Winston.Curve(cdfa, rawa, "color", parse(Winston.Colorant,   "black"))
#            Winston.add(ppp, tmp)
# tmp      = Winston.Curve(cdfb, rawb, "color", parse(Winston.Colorant,   "red"))
#            Winston.add(ppp, tmp)
# xyzzy = "gaussian_anamorphosis_test.png"
# print("writing $xyzzy\n")
# Winston.savefig(ppp, xyzzy, "width", 1700, "height", 1000)
end

for a = 1:linuma                                                              # report cal/val parameters for the first set
  vals = float(split(linea[a]))                                               # (in original units regardless of morphing)
  EXTRA && (vals[TOTB] = (vals[TOTB] - intb) / slob ;
            vals[TOTA] = (vals[TOTA] - inta) / sloa)
  cura[1,a,:] = [vals[TOTB] vals[OCUR]]
  cura[2,a,:] = [vals[TOTA] refa[a]   ]
end
a = 1 ; (mass, alp1, bet1, sig1, cor1, alp2, bet2, sig2, cor2, alp3, bet3, sig3, cor3) = triple(cura)
statis[a,MEMO,MASS] =        statis[a,MEMB,MASS] =        statis[a,MEMA,MASS] =        allmas[a] = mass
statis[a,MEMO,ALPH] = alp1 ; statis[a,MEMB,ALPH] = alp2 ; statis[a,MEMA,ALPH] = alp3 ; allalp[a] = alp1 #0.5 * (alp2 + alp3)
statis[a,MEMO,BETA] = bet1 ; statis[a,MEMB,BETA] = bet2 ; statis[a,MEMA,BETA] = bet3 ; allbet[a] = bet1 #0.5 * (bet2 + bet3)
statis[a,MEMO,SIGM] = sig1 ; statis[a,MEMB,SIGM] = sig2 ; statis[a,MEMA,SIGM] = sig3 ; allsig[a] = sig1 #0.5 * (sig2 + sig3)
statis[a,MEMO,CORR] = cor1 ; statis[a,MEMB,CORR] = cor2 ; statis[a,MEMA,CORR] = cor3 ; allcor[a] = cor1 #0.5 * (cor2 + cor3)
#a = 1 ; (allmas[a], allalp[a], allbet[a], allsig[a], allcor[a]) = triple(cura)

@printf("\nnumb = %15d for %s\n", linuma, ARGS[1])
@printf("cala = %15.8f mean(vals[TOTB]) = %15.8f\n",      allalp[a],  mean(cura[1,:,1]))
@printf("calb = %15.8f mean(vals[TOTA]) = %15.8f\n",      allbet[a],  mean(cura[2,:,1]))
@printf("mean = %15.8f mean(vals[OCUR]) = %15.8f\n", mean(allmas[a]), mean(cura[1,:,2]))
@printf("%33s %8s %8s %8s %8s\n", " ", "allalp", "allbet", "allsig", "allcor")
@printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", allalp[a], allbet[a], allsig[a], allcor[a])
@printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", statis[a,MEMO,ALPH], statis[a,MEMO,BETA], statis[a,MEMO,SIGM], statis[a,MEMO,CORR])
@printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", statis[a,MEMB,ALPH], statis[a,MEMB,BETA], statis[a,MEMB,SIGM], statis[a,MEMB,CORR])
@printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", statis[a,MEMA,ALPH], statis[a,MEMA,BETA], statis[a,MEMA,SIGM], statis[a,MEMA,CORR])

for a = 1:linumb                                                              # report cal/val parameters for the second set
  vals = float(split(lineb[a]))                                               # (in original units regardless of morphing)
  EXTRA && (vals[TOTB] = (vals[TOTB] - intb) / slob ;
            vals[TOTA] = (vals[TOTA] - inta) / sloa)
  curb[1,a,:] = [vals[TOTB] vals[OCUR]]
  curb[2,a,:] = [vals[TOTA] refb[a]   ]
end
a = 2 ; (mass, alp1, bet1, sig1, cor1, alp2, bet2, sig2, cor2, alp3, bet3, sig3, cor3) = triple(curb)
statis[a,MEMO,MASS] =        statis[a,MEMB,MASS] =        statis[a,MEMA,MASS] =        allmas[a] = mass
statis[a,MEMO,ALPH] = alp1 ; statis[a,MEMB,ALPH] = alp2 ; statis[a,MEMA,ALPH] = alp3 ; allalp[a] = alp1 #0.5 * (alp2 + alp3)
statis[a,MEMO,BETA] = bet1 ; statis[a,MEMB,BETA] = bet2 ; statis[a,MEMA,BETA] = bet3 ; allbet[a] = bet1 #0.5 * (bet2 + bet3)
statis[a,MEMO,SIGM] = sig1 ; statis[a,MEMB,SIGM] = sig2 ; statis[a,MEMA,SIGM] = sig3 ; allsig[a] = sig1 #0.5 * (sig2 + sig3)
statis[a,MEMO,CORR] = cor1 ; statis[a,MEMB,CORR] = cor2 ; statis[a,MEMA,CORR] = cor3 ; allcor[a] = cor1 #0.5 * (cor2 + cor3)
#a = 2 ; (allmas[a], allalp[a], allbet[a], allsig[a], allcor[a]) = triple(curb)

@printf("\nnumb = %15d for %s\n", linumb, ARGS222)
@printf("cala = %15.8f mean(vals[TOTB]) = %15.8f\n",      allalp[a],  mean(curb[1,:,1]))
@printf("calb = %15.8f mean(vals[TOTA]) = %15.8f\n",      allbet[a],  mean(curb[2,:,1]))
@printf("mean = %15.8f mean(vals[OCUR]) = %15.8f\n", mean(allmas[a]), mean(curb[1,:,2]))
@printf("%33s %8s %8s %8s %8s\n", " ", "allalp", "allbet", "allsig", "allcor")
@printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", allalp[a], allbet[a], allsig[a], allcor[a])
@printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", statis[a,MEMO,ALPH], statis[a,MEMO,BETA], statis[a,MEMO,SIGM], statis[a,MEMO,CORR])
@printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", statis[a,MEMB,ALPH], statis[a,MEMB,BETA], statis[a,MEMB,SIGM], statis[a,MEMB,CORR])
@printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", statis[a,MEMA,ALPH], statis[a,MEMA,BETA], statis[a,MEMA,SIGM], statis[a,MEMA,CORR])

if MORPH                                                                      # save linear mapping metrics
  rawavga11 = mean(cura[1,:,1]) ; rawavga21 = mean(cura[2,:,1]) ; rawavga22 = mean(cura[2,:,2])
  rawvara11 =  var(cura[1,:,1]) ; rawvara21 =  var(cura[2,:,1]) ; rawvara22 =  var(cura[2,:,2])
  rawavgb11 = mean(curb[1,:,1]) ; rawavgb21 = mean(curb[2,:,1]) ; rawavgb22 = mean(curb[2,:,2])
  rawvarb11 =  var(curb[1,:,1]) ; rawvarb21 =  var(curb[2,:,1]) ; rawvarb22 =  var(curb[2,:,2])
end

MORPH && (fpb = My.ouvre(ARGS[1] * ".cali.pair.morph", "w"))
MORPH || (fpb = My.ouvre(ARGS[1] * ".cali.pair",       "w"))
form = @sprintf("  mean param   MASS is %6.2f\n", mean(allmas[1]))
write(fpb, form)
form = @sprintf("  mean param   MASS is %6.2f\n", mean(allmas[2]))
write(fpb, form)
form = @sprintf("%77s %8s %8s %8s %8s\n", " ", "allalp", "allbet", "allsig", "allcor")
write(fpb, form)
form = @sprintf("%77s %8.3f %8.3f %8.3f %8.3f\n", ARGS[1], allalp[1], allbet[1], allsig[1], allcor[1])
write(fpb, form)
form = @sprintf("%77s %8.3f %8.3f %8.3f %8.3f\n", ARGS222, allalp[2], allbet[2], allsig[2], allcor[2])
write(fpb, form)
form = @sprintf("%77s %8.3f %8.3f %8.3f %8.3f\n", "obs", statis[1,MEMO,ALPH], statis[1,MEMO,BETA], statis[1,MEMO,SIGM], statis[1,MEMO,CORR])
write(fpb, form)
form = @sprintf("%77s %8.3f %8.3f %8.3f %8.3f\n", "obs", statis[2,MEMO,ALPH], statis[2,MEMO,BETA], statis[2,MEMO,SIGM], statis[2,MEMO,CORR])
write(fpb, form)
form = @sprintf("%77s %8.3f %8.3f %8.3f %8.3f\n", "bef", statis[1,MEMB,ALPH], statis[1,MEMB,BETA], statis[1,MEMB,SIGM], statis[1,MEMB,CORR])
write(fpb, form)
form = @sprintf("%77s %8.3f %8.3f %8.3f %8.3f\n", "bef", statis[2,MEMB,ALPH], statis[2,MEMB,BETA], statis[2,MEMB,SIGM], statis[2,MEMB,CORR])
write(fpb, form)
form = @sprintf("%77s %8.3f %8.3f %8.3f %8.3f\n", "aft", statis[1,MEMA,ALPH], statis[1,MEMA,BETA], statis[1,MEMA,SIGM], statis[1,MEMA,CORR])
write(fpb, form)
form = @sprintf("%77s %8.3f %8.3f %8.3f %8.3f\n", "aft", statis[2,MEMA,ALPH], statis[2,MEMA,BETA], statis[2,MEMA,SIGM], statis[2,MEMA,CORR])
write(fpb, form)
close(fpb)

fpb = My.ouvre(ARGS[1] * ".recalibrate", "w")                                 # and save these bias corrections
@printf(fpb, "%33.11f %33.11f %33.11f %33.11f\n",
allalp[1], allbet[1], allalp[2], allbet[2])
close(fpb)

if MORPH
  for a = 1:linuma                                                            # compute cal/val parameters for the first set
    vals = float(split(linea[a]))                                             # (in Gaussian units if morphing)
    EXTRA && (vals[TOTB] = (vals[TOTB] - intb) / slob ;
              vals[TOTA] = (vals[TOTA] - inta) / sloa)
    MORPH && (vals[TOTB] = rawtonorm(cdfb, rawb, vals[TOTB]) ;
              vals[TOTA] = rawtonorm(cdfb, rawb, vals[TOTA]) ;
              vals[OCUR] = rawtonorm(cdfb, rawb, vals[OCUR]) ;
              refa[a]    = rawtonorm(cdfb, rawb, refa[a]))
    cura[1,a,:] = [vals[TOTB] vals[OCUR]]
    cura[2,a,:] = [vals[TOTA] refa[a]   ]
  end
  a = 1 ; (mass, alp1, bet1, sig1, cor1, alp2, bet2, sig2, cor2, alp3, bet3, sig3, cor3) = triple(cura)
  statis[a,MEMO,MASS] =        statis[a,MEMB,MASS] =        statis[a,MEMA,MASS] =        allmas[a] = mass
  statis[a,MEMO,ALPH] = alp1 ; statis[a,MEMB,ALPH] = alp2 ; statis[a,MEMA,ALPH] = alp3 ; allalp[a] = alp1 #0.5 * (alp2 + alp3)
  statis[a,MEMO,BETA] = bet1 ; statis[a,MEMB,BETA] = bet2 ; statis[a,MEMA,BETA] = bet3 ; allbet[a] = bet1 #0.5 * (bet2 + bet3)
  statis[a,MEMO,SIGM] = sig1 ; statis[a,MEMB,SIGM] = sig2 ; statis[a,MEMA,SIGM] = sig3 ; allsig[a] = sig1 #0.5 * (sig2 + sig3)
  statis[a,MEMO,CORR] = cor1 ; statis[a,MEMB,CORR] = cor2 ; statis[a,MEMA,CORR] = cor3 ; allcor[a] = cor1 #0.5 * (cor2 + cor3)
# a = 1 ; (allmas[a], allalp[a], allbet[a], allsig[a], allcor[a]) = triple(cura)

  for a = 1:linumb                                                            # compute cal/val parameters for the second set
    vals = float(split(lineb[a]))                                             # (in Gaussian units if morphing)
    EXTRA && (vals[TOTB] = (vals[TOTB] - intb) / slob ;
              vals[TOTA] = (vals[TOTA] - inta) / sloa)
    MORPH && (vals[TOTB] = rawtonorm(cdfa, rawa, vals[TOTB]) ;
              vals[TOTA] = rawtonorm(cdfa, rawa, vals[TOTA]) ;
              vals[OCUR] = rawtonorm(cdfa, rawa, vals[OCUR]) ;
              refb[a]    = rawtonorm(cdfb, rawb, refb[a]))
    curb[1,a,:] = [vals[TOTB] vals[OCUR]]
    curb[2,a,:] = [vals[TOTA] refb[a]   ]
  end
  a = 2 ; (mass, alp1, bet1, sig1, cor1, alp2, bet2, sig2, cor2, alp3, bet3, sig3, cor3) = triple(curb)
  statis[a,MEMO,MASS] =        statis[a,MEMB,MASS] =        statis[a,MEMA,MASS] =        allmas[a] = mass
  statis[a,MEMO,ALPH] = alp1 ; statis[a,MEMB,ALPH] = alp2 ; statis[a,MEMA,ALPH] = alp3 ; allalp[a] = alp1 #0.5 * (alp2 + alp3)
  statis[a,MEMO,BETA] = bet1 ; statis[a,MEMB,BETA] = bet2 ; statis[a,MEMA,BETA] = bet3 ; allbet[a] = bet1 #0.5 * (bet2 + bet3)
  statis[a,MEMO,SIGM] = sig1 ; statis[a,MEMB,SIGM] = sig2 ; statis[a,MEMA,SIGM] = sig3 ; allsig[a] = cor1 #0.5 * (sig2 + sig3)
  statis[a,MEMO,CORR] = cor1 ; statis[a,MEMB,CORR] = cor2 ; statis[a,MEMA,CORR] = cor3 ; allcor[a] = sig1 #0.5 * (cor2 + cor3)
# a = 2 ; (allmas[a], allalp[a], allbet[a], allsig[a], allcor[a]) = triple(curb)
                                                                              # save linear mapping metrics
  gauavga11 = mean(cura[1,:,1]) ; gauavga21 = mean(cura[2,:,1]) ; gauavga22 = mean(cura[2,:,2])
  gauvara11 =  var(cura[1,:,1]) ; gauvara21 =  var(cura[2,:,1]) ; gauvara22 =  var(cura[2,:,2])
  gauavgb11 = mean(curb[1,:,1]) ; gauavgb21 = mean(curb[2,:,1]) ; gauavgb22 = mean(curb[2,:,2])
  gauvarb11 =  var(curb[1,:,1]) ; gauvarb21 =  var(curb[2,:,1]) ; gauvarb22 =  var(curb[2,:,2])
end

for a = 1:linuma                                                              # apply recalibration either in Gaussian or original units
  vals = float(split(linea[a]))                                               # using calibration parameters (and anamorphosis) from the
  EXTRA && (vals[TOTB] = (vals[TOTB] - intb) / slob ;                         # other set and report new cal/val parameters in original
            vals[TOTA] = (vals[TOTA] - inta) / sloa)                          # units
  MORPH && (vals[TOTB] = rawtonorm(cdfb, rawb, vals[TOTB]) ;
            vals[TOTA] = rawtonorm(cdfb, rawb, vals[TOTA]))
# vals[TOTB] = (vals[TOTB] - allalp[2]) / allbet[2]
# vals[TOTA] = (vals[TOTA] - allalp[2]) / allbet[2]
  vals[OCUR] = (vals[OCUR] - allalp[2]) / allbet[2]#; if vals[OCUR] < 0  vals[OCUR] = 0.0  end
# vals[OCUR] = (vals[OCUR] -       2.0) / 0.5
##vals[TOTB] = (vals[TOTB] - statis[2,MEMB,ALPH]) / statis[2,MEMB,ALPH]
##vals[TOTA] = (vals[TOTA] - statis[2,MEMA,ALPH]) / statis[2,MEMA,BETA]
  MORPH && (vals[TOTB] = normtoraw(cdfb, rawb, vals[TOTB]) ;
            vals[TOTA] = normtoraw(cdfb, rawb, vals[TOTA]))
  cura[1,a,:] = [vals[TOTB] vals[OCUR]]
  cura[2,a,:] = [vals[TOTA] refa[a]   ]
end
a = 3 ; (mass, alp1, bet1, sig1, cor1, alp2, bet2, sig2, cor2, alp3, bet3, sig3, cor3) = triple(cura)
statis[a,MEMO,MASS] =        statis[a,MEMB,MASS] =        statis[a,MEMA,MASS] =        allmas[a] = mass
statis[a,MEMO,ALPH] = alp1 ; statis[a,MEMB,ALPH] = alp2 ; statis[a,MEMA,ALPH] = alp3 ; allalp[a] = alp1 #0.5 * (alp2 + alp3)
statis[a,MEMO,BETA] = bet1 ; statis[a,MEMB,BETA] = bet2 ; statis[a,MEMA,BETA] = bet3 ; allbet[a] = bet1 #0.5 * (bet2 + bet3)
statis[a,MEMO,SIGM] = sig1 ; statis[a,MEMB,SIGM] = sig2 ; statis[a,MEMA,SIGM] = sig3 ; allsig[a] = sig1 #0.5 * (sig2 + sig3)
statis[a,MEMO,CORR] = cor1 ; statis[a,MEMB,CORR] = cor2 ; statis[a,MEMA,CORR] = cor3 ; allcor[a] = cor1 #0.5 * (cor2 + cor3)
#a = 3 ; (allmas[a], allalp[a], allbet[a], allsig[a], allcor[a]) = triple(cura)

          tmpstr  = "after recalibration only (using alpha and beta from the opposite set of collocations)"
MORPH && (tmpstr  = "after applying Gaussian anamorphosis, recalibration, and inverse Gaussian anamorphosis")

@printf("\nnumb = %15.0f for %s\n", linuma, ARGS[1])
@printf("cala = %15.8f %s\n", allalp[a], tmpstr)
@printf("calb = %15.8f %s\n", allbet[a], tmpstr)
@printf("mean = %15.8f %s\n", mean(allmas[a]), tmpstr)
@printf("%33s %8s %8s %8s %8s\n", " ", "allalp", "allbet", "allsig", "allcor")
@printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", allalp[a], allbet[a], allsig[a], allcor[a])
@printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", statis[a,MEMO,ALPH], statis[a,MEMO,BETA], statis[a,MEMO,SIGM], statis[a,MEMO,CORR])
@printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", statis[a,MEMB,ALPH], statis[a,MEMB,BETA], statis[a,MEMB,SIGM], statis[a,MEMB,CORR])
@printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", statis[a,MEMA,ALPH], statis[a,MEMA,BETA], statis[a,MEMA,SIGM], statis[a,MEMA,CORR])

for a = 1:linumb                                                              # apply recalibration either in Gaussian or original units
  vals = float(split(lineb[a]))                                               # using calibration parameters (and anamorphosis) from the
  EXTRA && (vals[TOTB] = (vals[TOTB] - intb) / slob ;                         # other set and report new cal/val parameters in original
            vals[TOTA] = (vals[TOTA] - inta) / sloa)                          # units
  MORPH && (vals[TOTB] = rawtonorm(cdfa, rawa, vals[TOTB]) ;
            vals[TOTA] = rawtonorm(cdfa, rawa, vals[TOTA]))
# vals[TOTB] = (vals[TOTB] - allalp[1]) / allbet[1]
# vals[TOTA] = (vals[TOTA] - allalp[1]) / allbet[1]
  vals[OCUR] = (vals[OCUR] - allalp[1]) / allbet[1]#; if vals[OCUR] < 0  vals[OCUR] = 0.0  end
##vals[TOTB] = (vals[TOTB] - statis[1,MEMB,ALPH]) / statis[1,MEMB,ALPH]
##vals[TOTA] = (vals[TOTA] - statis[1,MEMA,ALPH]) / statis[1,MEMA,BETA]
  MORPH && (vals[TOTB] = normtoraw(cdfa, rawa, vals[TOTB]) ;
            vals[TOTA] = normtoraw(cdfa, rawa, vals[TOTA]))
  curb[1,a,:] = [vals[TOTB] vals[OCUR]]
  curb[2,a,:] = [vals[TOTA] refb[a]   ]
end
a = 4 ; (mass, alp1, bet1, sig1, cor1, alp2, bet2, sig2, cor2, alp3, bet3, sig3, cor3) = triple(curb)
statis[a,MEMO,MASS] =        statis[a,MEMB,MASS] =        statis[a,MEMA,MASS] =        allmas[a] = mass
statis[a,MEMO,ALPH] = alp1 ; statis[a,MEMB,ALPH] = alp2 ; statis[a,MEMA,ALPH] = alp3 ; allalp[a] = alp1 #0.5 * (alp2 + alp3)
statis[a,MEMO,BETA] = bet1 ; statis[a,MEMB,BETA] = bet2 ; statis[a,MEMA,BETA] = bet3 ; allbet[a] = bet1 #0.5 * (bet2 + bet3)
statis[a,MEMO,SIGM] = sig1 ; statis[a,MEMB,SIGM] = sig2 ; statis[a,MEMA,SIGM] = sig3 ; allsig[a] = sig1 #0.5 * (sig2 + sig3)
statis[a,MEMO,CORR] = cor1 ; statis[a,MEMB,CORR] = cor2 ; statis[a,MEMA,CORR] = cor3 ; allcor[a] = cor1 #0.5 * (cor2 + cor3)
#a = 4 ; (allmas[a], allalp[a], allbet[a], allsig[a], allcor[a]) = triple(curb)

@printf("\nnumb = %15.0f for %s\n", linumb, ARGS222)
@printf("cala = %15.8f %s\n", allalp[a], tmpstr)
@printf("calb = %15.8f %s\n", allbet[a], tmpstr)
@printf("mean = %15.8f %s\n", mean(allmas[a]), tmpstr)
@printf("%33s %8s %8s %8s %8s\n", " ", "allalp", "allbet", "allsig", "allcor")
@printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", allalp[a], allbet[a], allsig[a], allcor[a])
@printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", statis[a,MEMO,ALPH], statis[a,MEMO,BETA], statis[a,MEMO,SIGM], statis[a,MEMO,CORR])
@printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", statis[a,MEMB,ALPH], statis[a,MEMB,BETA], statis[a,MEMB,SIGM], statis[a,MEMB,CORR])
@printf("%33s %8.3f %8.3f %8.3f %8.3f\n", " ", statis[a,MEMA,ALPH], statis[a,MEMA,BETA], statis[a,MEMA,SIGM], statis[a,MEMA,CORR])

if MORPH                                                                      # get average linear mapping metrics (as
  linrawa = (abs(rawavga11 - normtoraw(cdfb, rawb, gauavga11)) / rawvara11 +  # mean difference normalized by variance)
             abs(rawavga21 - normtoraw(cdfb, rawb, gauavga21)) / rawvara21 +
             abs(rawavga22 - normtoraw(cdfb, rawb, gauavga22)) / rawvara22) * 100 / 3
  linrawb = (abs(rawavgb11 - normtoraw(cdfb, rawb, gauavgb11)) / rawvarb11 +
             abs(rawavgb21 - normtoraw(cdfb, rawb, gauavgb21)) / rawvarb21 +
             abs(rawavgb22 - normtoraw(cdfb, rawb, gauavgb22)) / rawvarb22) * 100 / 3
  lingaua = (abs(gauavga11 - rawtonorm(cdfb, rawb, rawavga11)) / gauvara11 +
             abs(gauavga21 - rawtonorm(cdfb, rawb, rawavga21)) / gauvara21 +
             abs(gauavga22 - rawtonorm(cdfb, rawb, rawavga22)) / gauvara22) * 100 / 3
  lingaub = (abs(gauavgb11 - rawtonorm(cdfb, rawb, rawavgb11)) / gauvarb11 +
             abs(gauavgb21 - rawtonorm(cdfb, rawb, rawavgb21)) / gauvarb21 +
             abs(gauavgb22 - rawtonorm(cdfb, rawb, rawavgb22)) / gauvarb22) * 100 / 3
  linstr  = @sprintf("linrawa = %5.1f linrawb = %5.1f lingaua = %5.1f lingaub = %5.1f", linrawa, linrawb, lingaua, lingaub)
  println("\n$linstr\n")
  tmpstr *= " $linstr"
end

MORPH && (fpb = My.ouvre(ARGS[1] * ".cali.pair.morph", "a"))
MORPH || (fpb = My.ouvre(ARGS[1] * ".cali.pair",       "a"))
form = @sprintf("  mean param   MASS is %6.2f %s\n", mean(allmas[3]), tmpstr)
write(fpb, form)
form = @sprintf("  mean param   MASS is %6.2f %s\n", mean(allmas[4]), tmpstr)
write(fpb, form)
form = @sprintf("%77s %8.3f %8.3f %8.3f %8.3f\n", ARGS[1], allalp[3], allbet[3], allsig[3], allcor[3])
write(fpb, form)
form = @sprintf("%77s %8.3f %8.3f %8.3f %8.3f\n", ARGS222, allalp[4], allbet[4], allsig[4], allcor[4])
write(fpb, form)
form = @sprintf("%77s %8.3f %8.3f %8.3f %8.3f\n", "obs", statis[3,MEMO,ALPH], statis[3,MEMO,BETA], statis[3,MEMO,SIGM], statis[3,MEMO,CORR])
write(fpb, form)
form = @sprintf("%77s %8.3f %8.3f %8.3f %8.3f\n", "obs", statis[4,MEMO,ALPH], statis[4,MEMO,BETA], statis[4,MEMO,SIGM], statis[4,MEMO,CORR])
write(fpb, form)
form = @sprintf("%77s %8.3f %8.3f %8.3f %8.3f\n", "bef", statis[3,MEMB,ALPH], statis[3,MEMB,BETA], statis[3,MEMB,SIGM], statis[3,MEMB,CORR])
write(fpb, form)
form = @sprintf("%77s %8.3f %8.3f %8.3f %8.3f\n", "bef", statis[4,MEMB,ALPH], statis[4,MEMB,BETA], statis[4,MEMB,SIGM], statis[4,MEMB,CORR])
write(fpb, form)
form = @sprintf("%77s %8.3f %8.3f %8.3f %8.3f\n", "aft", statis[3,MEMA,ALPH], statis[3,MEMA,BETA], statis[3,MEMA,SIGM], statis[3,MEMA,CORR])
write(fpb, form)
form = @sprintf("%77s %8.3f %8.3f %8.3f %8.3f\n", "aft", statis[4,MEMA,ALPH], statis[4,MEMA,BETA], statis[4,MEMA,SIGM], statis[4,MEMA,CORR])
write(fpb, form)
close(fpb)
exit(0)

#=
       date    16   lat  26  lon   36  gnss  46   bef  56   now  66   aft  76  inc    87    elo   99    azo         ela         aza         gai
20150604000000   -52.500   253.750    20.250    10.005    10.442     9.874   165.8281     77.2194   -126.5719     77.6264   -128.6019     12.4070
=#
