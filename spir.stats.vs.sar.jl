#=
 = Compute statistics (bias, stdev, correlation) between SPIR GNSS-R retrievals
 = and those of other instruments.  All tracks for each of four types of GNSS-R
 = reflections (i.e., G1 G3 E11 E19) may be taken together.  The track sections
 = falling outside the SAR grid are first omitted - RD September 2016.
 =#

using My, NetCDF, Optim

if (argc = length(ARGS)) != 2
  print("\nUsage: jjj $(basename(@__FILE__)) WAV_E11_track00_L1_ZZ.txt 2015-12-03-160500DH.00040.hdr\n\n")
  exit(1)
end
gpsfil = ARGS[1]
hdrfil = ARGS[2][1:25] * ".hdr"
netfil = ARGS[2][1:25] * ".sar.nc"
pixel_spacing = parse(Int, ARGS[2][21:25])

print("\nreading $gpsfil\n")                                                  # read the GNSS-R retrievals
data = readdlm(gpsfil)
(gpslin, gpsval) = size(data)

print("reading sice, lats, and lons from $netfil\n")                          # read the input grids
sice = ncread(netfil, "sice", start=[1,1,1], count=[-1,-1,-1])
lats = ncread(netfil, "lats", start=[1,1,1], count=[-1,-1,-1])
lons = ncread(netfil, "lons", start=[1,1,1], count=[-1,-1,-1])
(npixels, nlines) = size(sice)
ncclose(netfil)

function telescope(reflat::Float64, reflon::Float64)                          # define the telescoping to identify
  mindis = minlat = minlon = 9e9                                              # any gridpoint nearest to a reflection
  minlin = div( nlines, 2)
  minpix = div(npixels, 2)
  for a = 25:-1:0
    for b = minlin - 2^a : 2^a : minlin + 2^a
      if b >= 0 && b < nlines
        for c = minpix - 2^a : 2^a : minpix + 2^a
          if c >= 0 && c < npixels
            tmplat = lats[c+1,b+1]
            tmplon = lons[c+1,b+1]
            tmpdis = (reflat - tmplat)^2 + (reflon - tmplon)^2
            if tmpdis < mindis
              mindis = tmpdis
              minlat = tmplat
              minlon = tmplon
              minlin = b
              minpix = c
#             @printf("a b c are (%8d %8d %8d) and lat lon mindis are %lf %lf %e ***\n", a, b, c, gdlatitude, gdlongitude, mindis)
#           else
#             @printf("a b c are (%8d %8d %8d) and lat lon tmpdis are %lf %lf %e\n", a, b, c, gdlatitude, gdlongitude, tmpdis)
            end
          end
        end
      end
    end
  end
  return(mindis, minlat, minlon, minlin, minpix)
end

mask = Array(Bool, gpslin)                                                    # identify the subset of GNSS-R retrievals
cols = Array(Int64, gpslin, 2)                                                # that are collocated with the SAR grid
b = 0
for a = 1:gpslin
  mask[a] = true
  (mindis, minlat, minlon, minlin, minpix) = telescope(data[a,1], data[a,2])
  if 111000 * mindis^0.5 > pixel_spacing ; mask[a] = false ; b += 1 ; end
  cols[a,1] = minlin
  cols[a,2] = minpix
end
print("excluded $b out of $gpslin retrievals\n")
data = data[mask,:]
cols = cols[mask,:]
gpslin -= b

refa = data[:,7]                                                              # form the vectors of interest
refb = Array(Float64, gpslin)
for a = 1:gpslin
# @show lats[cols[a,2],cols[a,1]], lons[cols[a,2],cols[a,1]], data[a,1], data[a,2]
  refb[a] = sice[cols[a,2],cols[a,1]]
end
refc = data[:,8]

cora = cor(refa, refb)
corb = cor(refc, refb)

outfil = ARGS[1] * "_" * ARGS[2] * ".cor"                                     # and save the statistics to file
fpa = My.ouvre(outfil, "w")
form  = @sprintf("%8.5f %8.5f %8.5f\n", cora, corb, (cora + corb) / 2.0)
write(fpa, form)
close(fpa)
exit(0)


#=
pars = ["latitude", "longitude", "height_rcv", "atmospheric_delay",    "SSH",       "sigma_SSH",    "mss_snr_max", "mss_snr_der", "elevation_txr"]
    60.12522891     25.90870658   3018.89700000      1.64627750     15.90209092      0.00505992      0.13903438      0.12132057     55.64010066
    60.12866312     25.92446450   3017.47900000      1.64599640     16.87194781      0.00612105      0.14347620      0.12473050     55.62153628
    60.13197433     25.94014694   3021.23900000      1.64807660     15.94872116      0.00265364      0.14175899      0.12504999     55.60286170

#   print("\nreflat reflon are $(data[a,1]) $(data[a,2])\n")
#   print("minlat minlon are $minlat $minlon at line pixel $minlin $minpix\n")
#   print("allowable distance between them is about $(10 * pixel_spacing) m\n")
#   print("actual    distance between them is about $(111000 * mindis^0.5) m\n")
# end
# if 111000 * mindis^0.5 > 10 * pixel_spacing ; exit(0) ; end

@printf("nlines           %10d   npixels          %10d\n\n", nlines, npixels)
degtorad = pi / 180.0                                                         # finally get a track angle and save such

hdrfil = @sprintf("%s.%s.hdr", stem, resolution)
print("\n") ; fpa = My.ouvre(hdrfil, "w") ; print("\n")
form = @sprintf("%-23s %18.10f %18.10f %18.10f %18.10f\n", "target_lat_plane_coefs", tarreslat.minimum[1], tarreslat.minimum[2], tarreslat.minimum[3], tarreslat.minimum[4]) ; write(fpa, form)
form = @sprintf("%-23s %18.10f %18.10f %18.10f %18.10f\n", "target_lon_plane_coefs", tarreslon.minimum[1], tarreslon.minimum[2], tarreslon.minimum[3], tarreslon.minimum[4]) ; write(fpa, form)
close(fpa)
exit(0)
=#
