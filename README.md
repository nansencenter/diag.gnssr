```shell
`
# unpack the IEEC/CSIC data and convert uint variables to double
wrkj ; mkdir spir ; mv SPIR_waveformsCorr_Dec2015.zip spir ; cd spir
       unzip SPIR_waveformsCorr_Dec2015.zip ; mv SPIR_waveformsCorr_Dec2015.zip ..
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_E11_track01_L1_AB.nc ../WAV_E11_track01_L1_AB.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_E11_track02_L1_BA.nc ../WAV_E11_track02_L1_BA.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_E11_track03_L1_CD.nc ../WAV_E11_track03_L1_CD.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_E11_track04_L1_DC.nc ../WAV_E11_track04_L1_DC.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_E11_track05_L1_AB.nc ../WAV_E11_track05_L1_AB.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_E11_track06_L1_BA.nc ../WAV_E11_track06_L1_BA.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_E11_track07_L1_CD.nc ../WAV_E11_track07_L1_CD.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_E11_track08_L1_DC.nc ../WAV_E11_track08_L1_DC.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_E19_track01_L1_AB.nc ../WAV_E19_track01_L1_AB.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_E19_track02_L1_BA.nc ../WAV_E19_track02_L1_BA.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_E19_track03_L1_CD.nc ../WAV_E19_track03_L1_CD.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_E19_track04_L1_DC.nc ../WAV_E19_track04_L1_DC.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_E19_track05_L1_AB.nc ../WAV_E19_track05_L1_AB.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_E19_track06_L1_BA.nc ../WAV_E19_track06_L1_BA.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_E19_track07_L1_CD.nc ../WAV_E19_track07_L1_CD.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_E19_track08_L1_DC.nc ../WAV_E19_track08_L1_DC.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_G1_track01_L1_AB.nc ../WAV_G1_track01_L1_AB.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_G1_track02_L1_BA.nc ../WAV_G1_track02_L1_BA.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_G1_track03_L1_CD.nc ../WAV_G1_track03_L1_CD.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_G1_track04_L1_DC.nc ../WAV_G1_track04_L1_DC.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_G1_track05_L1_AB.nc ../WAV_G1_track05_L1_AB.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_G1_track06_L1_BA.nc ../WAV_G1_track06_L1_BA.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_G1_track07_L1_CD.nc ../WAV_G1_track07_L1_CD.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_G1_track08_L1_DC.nc ../WAV_G1_track08_L1_DC.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_G3_track01_L1_AB.nc ../WAV_G3_track01_L1_AB.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_G3_track02_L1_BA.nc ../WAV_G3_track02_L1_BA.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_G3_track03_L1_CD.nc ../WAV_G3_track03_L1_CD.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_G3_track04_L1_DC.nc ../WAV_G3_track04_L1_DC.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_G3_track05_L1_AB.nc ../WAV_G3_track05_L1_AB.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_G3_track06_L1_BA.nc ../WAV_G3_track06_L1_BA.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_G3_track07_L1_CD.nc ../WAV_G3_track07_L1_CD.nc
       ncap2 -s 'prn=double(prn);gps_week=double(gps_week);gps_sow=double(gps_sow);msec_coh=double(msec_coh);sec_incoh=double(sec_incoh);range_samples=double(range_samples)' WAV_G3_track08_L1_DC.nc ../WAV_G3_track08_L1_DC.nc
       cd .. ; rm -r spir

# create plots and text files for each type of GNSS-R
wrkj ; parallel julia ~/bin/plot.gnssr.track.spir-uav.jl ::: WAV*nc
       cat WAV_E19*nc.txt > WAV_E19_track00_L1_ZZ.txt
       cat WAV_E11*nc.txt > WAV_E11_track00_L1_ZZ.txt
       cat WAV_G1*nc.txt > WAV_G1_track00_L1_ZZ.txt
       cat WAV_G3*nc.txt > WAV_G3_track00_L1_ZZ.txt

# download JASON-2 data from AVISO and create ncdumps
wrkj ; wget ftp://avisoftp.cnes.fr/AVISO/pub/jason-2/gdr_d/cycle_273/
       wget ftp://avisoftp.cnes.fr/AVISO/pub/jason-2/gdr_d/cycle_273/JA2_GPN_2PdP273_092_20151203_155445_20151203_165058.nc
       wget ftp://avisoftp.cnes.fr/AVISO/pub/jason-2/gdr_d/cycle_273/JA2_GPN_2PdP273_111_20151204_094249_20151204_103902.nc
       ncdump JA2_GPN_2PdP273_092_20151203_155445_20151203_165058.nc > JA2_GPN_2PdP273_092_20151203_155445_20151203_165058.ncdump
       ncdump JA2_GPN_2PdP273_111_20151204_094249_20151204_103902.nc > JA2_GPN_2PdP273_111_20151204_094249_20151204_103902.ncdump

# create data plots for Gulf of Finland crossings
wrkj ; jjj plot.gnssr.track.altimetry.jl JA2_GPN_2PdP273_092_20151203_155445_20151203_165058.nc
       jjj plot.gnssr.track.altimetry.jl JA2_GPN_2PdP273_111_20151204_094249_20151204_103902.nc







# the waveforms and more complete data can also be obtained from AVISO
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
