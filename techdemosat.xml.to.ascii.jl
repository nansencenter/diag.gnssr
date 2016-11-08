#=
 = Unpack all zip files in a directory and parse the XML files in
 = them to ASCII files in the current directory - RD October 2016.
 =#

using My, LightXML
const MISS             = -9999.0                        # generic missing value

if (argc = length(ARGS)) != 1
  print("\nUsage: jjj $(basename(@__FILE__)) /net/sverdrup-1/vol/Projects/E-GEM/L2/Tracks/RD000047\n")
  print(  "or     jjj $(basename(@__FILE__)) /home/ricani/work/workj/data_tds/L2/Tracks\n\n")
  exit(1)
end

absdir = abspath(ARGS[1])                                                     # unpack each XML file in this dir
files  = readdir(absdir)
for file in files
  full = joinpath(absdir, file)
  if isfile(full) && endswith(full, ".zip")
    run(`unzip -o -d $absdir $full`)
  end
end

allfil = Array(UTF8String, 0)                                                 # then read the list of XML files
files  = readdir(absdir)
for file in files
  full = joinpath(absdir, file)
  if isfile(full) && endswith(full, ".xml")
    push!(allfil, full)
  end
end

stem = split(dirname(allfil[1]), "/")[end]                                    # and create corresponding text files
for xmlfil in allfil                                                          # (if there is valid data to store)
  tmpfil = stem * "_" * basename(xmlfil) * ".txt"
  flag = false ; fpa = My.ouvre(tmpfil, "w")
  xmldat = find_element(root(parse_file(xmlfil)), "L2DataEpochArray")
  for xmlele in child_elements(xmldat)
    spd          = float(      content(find_element(xmlele, "WindSpeed")))
    if !isnan(spd)
      dat        =             content(find_element(xmlele, "IntegrationMidPointTime"))
      (lat, lon) = float(split(content(find_element(xmlele, "SpecularPointLatLonHeight"))))
      inc        = float(      content(find_element(xmlele, "SPIncidenceAngle")))
      elo        = float(      content(find_element(xmlele, "SPElevationORF"))) ; azo = float(content(find_element(xmlele, "SPAzimuthORF")))
      ela        = float(      content(find_element(xmlele, "SPElevationARF"))) ; aza = float(content(find_element(xmlele, "SPAzimuthARF")))
      gai        = float(      content(find_element(xmlele, "AntennaGainTowardsSpecularPoint")))
#     mss        = float(      content(find_element(xmlele, "MSS"))) ; isnan(mss) && (mss = MISS)

      dat = split(replace(replace(replace(dat, "-", ""), ":", ""), "T", ""), ".")[1]
      form = @sprintf("%s %11.4f %11.4f %11.4f %11.4f %11.4f %11.4f %11.4f %11.4f %11.4f\n", dat, lat, lon, spd, inc, elo, azo, ela, aza, gai)
      write(fpa, form)
      flag = true
    end
  end
  close(fpa)
  flag == false && rm(tmpfil)
end
exit(0)
