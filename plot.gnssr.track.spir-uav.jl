#=
 = Plot an initial inspection timeseries of parameters on multiple panels
 = (and unpack the IEEC/CSIC data and convert uint variables to double to
 = allow for limited NetCDF package coverage of uint) - RD September 2016.
 =#

using My, NetCDF, Winston

if (argc = length(ARGS)) != 1
  print("\nUsage: jjj $(basename(@__FILE__)) WAV_E11_track01_L1_AB.nc\n\n")
  exit(1)
end
altfil = ARGS[1]
xyzfil = ARGS[1] * ".xyzzy"
run(`ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' $altfil $xyzfil`)

timwek = ncread(xyzfil, "gps_week")
timsow = ncread(xyzfil, "gps_sow")
timlen = length(timwek)
print("found $timlen points along the track\n")

timstr = Array(UTF8String, timlen)                                            # parse the dates as string labels
for a = 1:timlen
  timstr[a] = dategps(timwek[a], timsow[a])
end

pars = ["latitude_spec", "longitude_spec", "height_rcv", "atmospheric_delay", "SSH", "sigma_SSH", "mss_snr_max", "mss_snr_der", "elevation_txr"]
parn = length(pars)
data = Array(Float64, timlen, parn)

for (a, para) in enumerate(pars)                                              # read all data before removing the
  vara = ncread(xyzfil, para)                                                 # uint conversion file
  data[:,a] = vara
end
rm(xyzfil)

fpa = My.ouvre(altfil * ".txt", "w")                                          # echo all data to an ASCII file for
for a = 1:timlen, b = 1:parn                                                  # later geospatial plots
  form  = @sprintf("%15.8lf", data[a,b])
  form *= b == parn ? "\n" : " "
  write(fpa, form)
# print(form)
end
close(fpa)

ppp = Winston.Table(parn,1) ; setattr(ppp, "cellpadding", -0.5)               # and make timeseries plot
for (a, para) in enumerate(pars)
# print("adding $para\n")
  tpos = (a, 1)
  tmp = Winston.FramedPlot()
  ppp[tpos...] = Winston.add(tmp)
  tmp = Winston.Curve(collect(1:timlen), data[:,a])
        Winston.add(ppp[tpos...], tmp)
  tmp = Winston.PlotLabel(.53, .92, para, "size", 3.4)
        Winston.add(ppp[tpos...], tmp)
end

xyzzy = ARGS[1]*".png"
print("writing $xyzzy\n")
Winston.savefig(ppp, xyzzy, "width", 1000, "height", 1700)
exit(0)
