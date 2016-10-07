```shell
`
# get a copy
git clone git@github.com:nansencenter/diag.gnssr.git

# requirements on ubuntu 14.04
julia (http://julialang.org/)
alias  jjj 'julia $HOME/bin/\!*'
alias wrkj 'cd ~/work/workj; ls -al'
GNU parallel (http://www.gnu.org/software/parallel/)
https://github.com/rickedanielson/grads.lib (for additional GrADs scripts)

# plot overpasses of interest on the day in question
wrkj ; grads -blc "plot.gnssr.track.overpasses"

# unpack the IEEC/CSIC SPIR data and create inspection plots and text files for each type of GNSS-R reflection and track
wrkj ; cd data_skyvan
       unzip SPIR_waveformsCorr_Dec2015.zip
       parallel julia ~/bin/plot.gnssr.track.spir-uav.jl ::: WAV*nc
       cat WAV_E19*nc.txt > ../WAV_E19_track00_L1_ZZ.txt
       cat WAV_E11*nc.txt > ../WAV_E11_track00_L1_ZZ.txt
       cat WAV_G1*nc.txt > ../WAV_G1_track00_L1_ZZ.txt
       cat WAV_G3*nc.txt > ../WAV_G3_track00_L1_ZZ.txt

# following Sentinel-1 download from ESA and conversion of the full scene from SAFE to netCDF, create a subdomain
wrkj ; cd data_sentinel
       jjj sentinel.sar.cross.section.fullres.jl S1A_EW_GRDM_1SDH_20151203T160430_20151203T160530_008880_00CB22_53F6.zip 60.050 25.900 40000 45000 60
       cp 2015-12-03-160500DH.00040* ..

# mask the subdomain land and bright returns (ships) at full resolution and then reduce the subdomain resolution
wrkj ; sentinel.sar.cross.section.mask   2015-12-03-160500DH.00040.hdr -15.0 -22.0
       sentinel.sar.cross.section.reduce 2015-12-03-160500DH.00040.hdr 1
       sentinel.sar.cross.section.reduce 2015-12-03-160500DH.00040.hdr 2
       sentinel.sar.cross.section.reduce 2015-12-03-160500DH.00040.hdr 3
       sentinel.sar.cross.section.reduce 2015-12-03-160500DH.00040.hdr 4

# plot a comparison of different resolutions
wrkj ; grads -blc "plot.gnssr.track.spir.and.sar 2015-12-03-160500DH.00040.hdr"
       grads -blc "plot.gnssr.track.spir.and.sar 2015-12-03-160500DH.00080.hdr"
       grads -blc "plot.gnssr.track.spir.and.sar 2015-12-03-160500DH.00160.hdr"
       grads -blc "plot.gnssr.track.spir.and.sar 2015-12-03-160500DH.00320.hdr"
       grads -blc "plot.gnssr.track.spir.and.sar 2015-12-03-160500DH.00640.hdr"

# compute a SAR proxy for MSS following Kudryavtsev et al. (2012)
wrkj ; sentinel.sar.cross.section.mss 2015-12-03-160500DH.00040.hdr 20
       sentinel.sar.cross.section.mss 2015-12-03-160500DH.00080.hdr 20
       sentinel.sar.cross.section.mss 2015-12-03-160500DH.00160.hdr 20
       sentinel.sar.cross.section.mss 2015-12-03-160500DH.00320.hdr 20
       sentinel.sar.cross.section.mss 2015-12-03-160500DH.00640.hdr 20

# identify the highest correlation between SAR and GNSS-R MSS proxies for both reflection and track groups
wrkj ; parallel julia /home/ricani/bin/spir.stats.vs.sar.jl :::             WAV*txt ::: *hdr
       parallel julia /home/ricani/bin/spir.stats.vs.sar.jl ::: data_skyvan/WAV*txt ::: *hdr
       cat data_skyvan/*.hdr.cor > SPIR_waveformsCorr_Dec2015.cor

# plot a comparison of SPIR and Sentinel-1 MSS at each resolution and for each type of GNSS-R reflection and track
wrkj ; parallel grads -blc '"'plot.gnssr.track.spir.and.sar :::             WAV*txt ::: *hdr ::: '"'

# download JASON-2 data from AVISO and create inspection plots for the Gulf of Finland crossings
wrkj ; cd data_jason
       wget ftp://avisoftp.cnes.fr/AVISO/pub/jason-2/gdr_d/cycle_273/
       wget ftp://avisoftp.cnes.fr/AVISO/pub/jason-2/gdr_d/cycle_273/JA2_GPN_2PdP273_092_20151203_155445_20151203_165058.nc
       wget ftp://avisoftp.cnes.fr/AVISO/pub/jason-2/gdr_d/cycle_273/JA2_GPN_2PdP273_111_20151204_094249_20151204_103902.nc
       ncdump JA2_GPN_2PdP273_092_20151203_155445_20151203_165058.nc > JA2_GPN_2PdP273_092_20151203_155445_20151203_165058.ncdump
       ncdump JA2_GPN_2PdP273_111_20151204_094249_20151204_103902.nc > JA2_GPN_2PdP273_111_20151204_094249_20151204_103902.ncdump
       jjj plot.gnssr.track.altimetry.jl JA2_GPN_2PdP273_092_20151203_155445_20151203_165058.nc
       jjj plot.gnssr.track.altimetry.jl JA2_GPN_2PdP273_111_20151204_094249_20151204_103902.nc
       cp JA2*.nc.txt ..

# following ECMWF download, plot an initial view at GNSS-R collection time
wrkj ; grads -blc "plot.gnssr.track.spir.and.altim JA2_GPN_2PdP273_092_20151203_155445_20151203_165058.nc.txt WAV_G1_track00_L1_ZZ.txt windalt"
       grads -blc "plot.gnssr.track.spir.and.altim JA2_GPN_2PdP273_092_20151203_155445_20151203_165058.nc.txt WAV_G1_track00_L1_ZZ.txt wavealt"
       grads -blc "plot.gnssr.track.spir.and.altim JA2_GPN_2PdP273_092_20151203_155445_20151203_165058.nc.txt WAV_G1_track00_L1_ZZ.txt hghtalt"
       grads -blc "plot.gnssr.track.spir.and.altim JA2_GPN_2PdP273_092_20151203_155445_20151203_165058.nc.txt WAV_G1_track00_L1_ZZ.txt windgps"
       grads -blc "plot.gnssr.track.spir.and.altim JA2_GPN_2PdP273_092_20151203_155445_20151203_165058.nc.txt WAV_G1_track00_L1_ZZ.txt wavegps"
       grads -blc "plot.gnssr.track.spir.and.altim JA2_GPN_2PdP273_092_20151203_155445_20151203_165058.nc.txt WAV_G1_track00_L1_ZZ.txt hghtgps"

# plot an initial comparison of SPIR and JASON-2 SSH (omitting atmospheric corrections)
wrkj ; grads -blc "plot.gnssr.track.spir.and.altim JA2_GPN_2PdP273_092_20151203_155445_20151203_165058.nc.txt WAV_G1_track00_L1_ZZ.txt hghtall"
       grads -blc "plot.gnssr.track.spir.and.altim JA2_GPN_2PdP273_092_20151203_155445_20151203_165058.nc.txt WAV_G3_track00_L1_ZZ.txt hghtall"
       grads -blc "plot.gnssr.track.spir.and.altim JA2_GPN_2PdP273_092_20151203_155445_20151203_165058.nc.txt WAV_E11_track00_L1_ZZ.txt hghtall"
       grads -blc "plot.gnssr.track.spir.and.altim JA2_GPN_2PdP273_092_20151203_155445_20151203_165058.nc.txt WAV_E19_track00_L1_ZZ.txt hghtall"
       di plot.ECMWF.hght.JA2_GPN_2PdP273_092_20151203_155445*



# the JASON-2 waveforms and more complete data can also be obtained from AVISO
wrkj ; wget ftp://avisoftp.cnes.fr/AVISO/pub/jason-2/sgdr_d/cycle_273/
       wget ftp://avisoftp.cnes.fr/AVISO/pub/jason-2/sgdr_d/cycle_273/JA2_GPS_2PdP273_092_20151203_155445_20151203_165058.zip
       wget ftp://avisoftp.cnes.fr/AVISO/pub/jason-2/sgdr_d/cycle_273/JA2_GPS_2PdP273_111_20151204_094249_20151204_103902.zip
       wget ftp://avisoftp.cnes.fr/AVISO/pub/jason-2/gdr_d/cycle_273/
       wget ftp://avisoftp.cnes.fr/AVISO/pub/jason-2/gdr_d/cycle_273/JA2_GPN_2PdP273_092_20151203_155445_20151203_165058.nc
       wget ftp://avisoftp.cnes.fr/AVISO/pub/jason-2/gdr_d/cycle_273/JA2_GPN_2PdP273_111_20151204_094249_20151204_103902.nc
       wget ftp://avisoftp.cnes.fr/AVISO/pub/jason-2/sgdr_d/cycle_273/
       wget ftp://avisoftp.cnes.fr/AVISO/pub/jason-2/sgdr_d/cycle_273/JA2_GPS_2PdP273_092_20151203_155445_20151203_165058.zip
       wget ftp://avisoftp.cnes.fr/AVISO/pub/jason-2/sgdr_d/cycle_273/JA2_GPS_2PdP273_111_20151204_094249_20151204_103902.zip

# and for CryoSat-2, look to JPL to obtain near overpasses
wrkj ; wget ftp://podaac-ftp.jpl.nasa.gov/allData/saral/preview/L2/XOGDR-SSHA/c029/
       wget ftp://podaac-ftp.jpl.nasa.gov/allData/saral/preview/L2/XOGDR-SSHA/c029/SRL_OPRSSHA_2PTS029_0425_20151204_014859_20151204_032752.EUM.nc
       ncdump SRL_OPRSSHA_2PTS029_0425_20151204_014859_20151204_032752.EUM.nc > SRL_OPRSSHA_2PTS029_0425_20151204_014859_20151204_032752.EUM.ncdump
       unzip JA2_GPS_2PdP273_092_20151203_155445_20151203_165058.zip
       unzip JA2_GPS_2PdP273_111_20151204_094249_20151204_103902.zip
       ncdump JA2_GPS_2PdP273_092_20151203_155445_20151203_165058.nc > JA2_GPS_2PdP273_092_20151203_155445_20151203_165058.ncdump
       ncdump JA2_GPS_2PdP273_111_20151204_094249_20151204_103902.nc > JA2_GPS_2PdP273_111_20151204_094249_20151204_103902.ncdump

# for reference, the IEEC/CSIC data are
       WAV_E11_track01_L1_AB.nc WAV_E11_track02_L1_BA.nc WAV_E11_track03_L1_CD.nc WAV_E11_track04_L1_DC.nc
       WAV_E11_track05_L1_AB.nc WAV_E11_track06_L1_BA.nc WAV_E11_track07_L1_CD.nc WAV_E11_track08_L1_DC.nc
       WAV_E19_track01_L1_AB.nc WAV_E19_track02_L1_BA.nc WAV_E19_track03_L1_CD.nc WAV_E19_track04_L1_DC.nc
       WAV_E19_track05_L1_AB.nc WAV_E19_track06_L1_BA.nc WAV_E19_track07_L1_CD.nc WAV_E19_track08_L1_DC.nc
       WAV_G1_track01_L1_AB.nc WAV_G1_track02_L1_BA.nc   WAV_G1_track03_L1_CD.nc WAV_G1_track04_L1_DC.nc
       WAV_G1_track05_L1_AB.nc WAV_G1_track06_L1_BA.nc   WAV_G1_track07_L1_CD.nc WAV_G1_track08_L1_DC.nc
       WAV_G3_track01_L1_AB.nc WAV_G3_track02_L1_BA.nc   WAV_G3_track03_L1_CD.nc WAV_G3_track04_L1_DC.nc
       WAV_G3_track05_L1_AB.nc WAV_G3_track06_L1_BA.nc   WAV_G3_track07_L1_CD.nc WAV_G3_track08_L1_DC.nc

# Those starting with WAV_E# contain GALILEO reflections, those others starting with WAV_G# are for GPS reflections.  Track0# is simply the track numbering within the flight chronologically sorted. Each track corresponds to one segment, identified by its starting-ending points (e.g. AB, BC,...). Since the cross-like pattern was flew twice, tracks 1 to 4 correspond to the first round (AB, BA, CD, DC), tracks 5 to 8 to the second round (same segments again).
# The header of the netcdf should contain all the information for you to understand the data. The fields are both observables (waveforms, variable "waveform(record, range)"), intermediate products (e.g. the delay at which the maximum derivative of the waveform occurs), as well as geophysical products (altimetry given in variable "SSH(record)" (and its uncertainty ("sigma_SSH(record)").  To obtain the SSH retrievals we have corrected by tropospheric delays. The correction we have applied is under "atmospheric_delay" variable (this could be redone).
# Note that the time is given the gps-way, in gps weeks and seconds of the week. For example, the first second in the set: gps_week 1873 + gps_sow 384730 = 3 december 2015 10:52:10 UTC time (the rest of seconds are simply subsequent to this one).
# Be cautious with the different values of mean square slope (mss) as each model is not perfectly tunned to the SPIR instrument (complex and dynamic antenna beam forming, changes in attitude, etc).

# JASON-2 ssha 10
 # = altitude of satellite (alt) 9
 # - Ku band corrected altimeter range (range_ku) 7
 # - altimeter ionospheric correction on Ku band (iono_corr_alt_ku) 11
 # - model dry tropospheric correction (model_dry_tropo_corr) 12
 # - radiometer wet tropospheric correction (rad_wet_tropo_corr) 13
# - sea state bias correction in Ku band (sea_state_bias_ku) 14
# - solid earth tide height (solid_earth_tide) 15
# - geocentric ocean tide height solution 1 (ocean_tide_sol1) 16
# - geocentric pole tide height (pole_tide) 17
# - inverted barometer height correction (inv_bar_corr) 18
# - high frequency fluctuations of the sea surface topography (hf_fluctuations_corr for I/GDR off line products only) 19
# - mean sea surface (mean_sea_surface). Set to default if the altimeter echo type (alt_echo_type) is set to 1 = non ocean like, the radiometer surface type (rad_surf_type) set to 2 = land, or the rain flag (rain_flag) set to 1 = rain 20

"range_ku", 7   "range_c", 8      "alt",  9          "ssha",  10   "iono_corr_alt_ku", "model_dry_tropo_corr", "rad_wet_tropo_corr"
52378.97325     52379.01582222222 52394.173772222224 32.767        -0.0151333333333333 -2.3104277777777775     -0.0449
52393.79453     52393.89466923077 52410.02732307692  32.767        -0.0249769230769230 -2.291592307692308      -0.108
 52453.53030000  52453.64790000    52469.24250000     32.76700000   -0.02810000         -2.30980000             -0.04770000
 52340.61470000  52340.73870000    52357.37510000     32.76700000   -0.02960000         -2.29060000             -0.11050000

"sea_state_bias_ku", "solid_earth_tide", "ocean_tide_sol1", "pole_tide",   "inv_bar_corr", "hf_fluctuations_corr", "mean_sea_surface"
-0.04087222222222222 -0.00115             0.005694444444444 -0.00469444444 -0.058138888888  0.03821111111111111    17.13666111
-0.04116923076923077 -0.0560769230769230 -0.032707692307692 -0.0046         0.024846153846 -0.10313076923076923    18.25978461
 -0.02420000          -0.00200000         -0.00090000        -0.00460000    -0.05570000      0.05010000             17.59440000
 -0.03230000          -0.05620000         -0.03900000        -0.00460000     0.02940000     -0.07320000             18.71540000

"latitude_spec", "longitude_spec", "height_rcv",  "atmospheric_delay", "SSH",        "sigma_SSH",  "mss_snr_max",   "mss_snr_der", "elevation_txr"
60.03983392      25.53001704       3043.57400000  1.43702302           17.36168613   0.00220182    0.03194709       0.02887435     72.20631125
