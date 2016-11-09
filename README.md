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

# create a mirror of the SSTL TDS-1 wind speed retrievals
       cd /net/sverdrup-1/vol/Projects/E-GEM/L2
       echo mirror "Data/L2" . | lftp ftp.merrbys.org

# convert the TDS-1 XML data files to ASCII, with use of RD33 onward safest (6-9,16-20,22-27,33,37-39,41,43-56,59-61,63-67,69,70)
# (0-5 commissioning; 12-15 orbital elements corrupt; 21 calibration; 28-32 DDM and attitude suspect; 34 unknown, 36,40,42,57,58,62,68 missing)
wrkj ; parallel julia /home/ricani/bin/techdemosat.xml.to.ascii.jl ::: /net/sverdrup-1/vol/Projects/E-GEM/L2/Tracks/*
       wc RD00000[1-5]* RD00001[2-5]* RD00002[189]* RD00003[0-2]*
       mkdir all ; mv RD* data_tds ; cat data_tds/RD0000[3-7]*txt > all/gnss.tds1.txt

# average 45 days of TDS-1 data (20150603 to 20160418) to the 6-h and 0.125-degree resolution of the ECMWF analysis (on johansen)
wrkj ; seq 335 669 | parallel -j 16 julia /Home/ricani/bin/techdemosat.wind.speed.collate.jl all/gnss.tds1.txt
       cat all/gnss.tds1.txt.2* > all/gnss.tds1.txt.all

# identify all locations with at least one daily-average observation
wrkj ; cd all ; mkdir plot.available plot.histogr plot.locate
       jjj techdemosat.wind.speed.locate.jl gnss.tds1.txt.all /home/ricani/data/ecmwf/ECMWF_Interim_invariants_0.125.nc
       xvfb-run -a grads -blc "techdemosat.wind.speed.locate gnss.tds1.txt.all.locate"
       mv plot.techdemosat.wind.speed.dots*locate*png plot.locate

# make an initial split of the daily-average observations into cal/val/extra groups (based only on insitu, and not analysis, availability)
wrkj ; cd all
       jjj techdemosat.wind.speed.collate.split.jl gnss.tds1.txt.all.locate
       jjj techdemosat.wind.speed.collate.split.jl gnss.tds1.txt.all.locate_1.0_valid
       jjj techdemosat.wind.speed.collate.split.jl gnss.tds1.txt.all.locate_1.0_valid_1.0_valid
       sort gnss.tds1.txt.all.locate_1.0_valid                     > gnss.tds1.txt.all.locate_1.0_calib_remainder
       sort gnss.tds1.txt.all.locate_1.0_valid_1.0_calib           > gnss.tds1.txt.all.locate_1.0_valid
       sort gnss.tds1.txt.all.locate_1.0_valid_1.0_valid           > gnss.tds1.txt.all.locate_1.0_valid_remainder
       sort gnss.tds1.txt.all.locate_1.0_valid_1.0_valid_1.0_calib > gnss.tds1.txt.all.locate_1.0_extra
       sort gnss.tds1.txt.all.locate_1.0_valid_1.0_valid_1.0_valid > gnss.tds1.txt.all.locate_1.0_extra_remainder
       rm gnss.tds1.txt.all.locate_1.0*1.0*
       xvfb-run -a grads -blc "techdemosat.wind.speed.locate gnss.tds1.txt.all.locate"
       xvfb-run -a grads -blc "techdemosat.wind.speed.locate gnss.tds1.txt.all.locate_1.0_calib"
       xvfb-run -a grads -blc "techdemosat.wind.speed.locate gnss.tds1.txt.all.locate_1.0_calib_remainder"
       xvfb-run -a grads -blc "techdemosat.wind.speed.locate gnss.tds1.txt.all.locate_1.0_valid"
       xvfb-run -a grads -blc "techdemosat.wind.speed.locate gnss.tds1.txt.all.locate_1.0_valid_remainder"
       xvfb-run -a grads -blc "techdemosat.wind.speed.locate gnss.tds1.txt.all.locate_1.0_extra"
       xvfb-run -a grads -blc "techdemosat.wind.speed.locate gnss.tds1.txt.all.locate_1.0_extra_remainder"
       mv plot.techdemosat.wind.speed.dots*locate*png plot.locate

# further split the in situ cal/val observations by location and store files in an insitu dir
wrkj ; mkdir insitu
       sort -k2,2 -k3,3 -k1,1 all/gnss.tds1.txt.all > gnss.tds1.txt.all.sort
       parallel julia /home/ricani/bin/techdemosat.wind.speed.collate.split.location.jl ::: all/gnss.tds1.txt.all.locate_1.0_????? ::: gnss.tds1.txt.all.sort
       cd insitu ; ls -1 ins* | grep -v GEM > z.list ; cd .. ; wc insitu/z.list
       rm gnss.tds1.txt.all.sort

# create example ncdumps of ERA Interim data (downloaded by month from http://apps.ecmwf.int/datasets/data/interim-full-daily)
wrkj ; mkdir ncdump
       ncdump /mnt/10.11.12.232/sat_auxdata/model/ecmwf/0.125-deg/interim_2015-06_oper.nc > ncdump/interim_2015-06_oper.ncdump
       ncdump /mnt/10.11.12.232/sat_auxdata/model/ecmwf/0.125-deg/interim_2015-07_oper.nc > ncdump/interim_2015-07_oper.ncdump
       ncdump /mnt/10.11.12.232/sat_auxdata/model/ecmwf/0.125-deg/interim_2015-08_oper.nc > ncdump/interim_2015-08_oper.ncdump
       ncdump /mnt/10.11.12.232/sat_auxdata/model/ecmwf/0.125-deg/interim_2015-09_oper.nc > ncdump/interim_2015-09_oper.ncdump
       ncdump /mnt/10.11.12.232/sat_auxdata/model/ecmwf/0.125-deg/interim_2015-10_oper.nc > ncdump/interim_2015-10_oper.ncdump
       ncdump /mnt/10.11.12.232/sat_auxdata/model/ecmwf/0.125-deg/interim_2015-11_oper.nc > ncdump/interim_2015-11_oper.ncdump
       ncdump /mnt/10.11.12.232/sat_auxdata/model/ecmwf/0.125-deg/interim_2015-12_oper.nc > ncdump/interim_2015-12_oper.ncdump
       ncdump /mnt/10.11.12.232/sat_auxdata/model/ecmwf/0.125-deg/interim_2016-01_oper.nc > ncdump/interim_2016-01_oper.ncdump
       ncdump /mnt/10.11.12.232/sat_auxdata/model/ecmwf/0.125-deg/interim_2016-02_oper.nc > ncdump/interim_2016-02_oper.ncdump
       ncdump /mnt/10.11.12.232/sat_auxdata/model/ecmwf/0.125-deg/interim_2016-03_oper.nc > ncdump/interim_2016-03_oper.ncdump
       ncdump /mnt/10.11.12.232/sat_auxdata/model/ecmwf/0.125-deg/interim_2016-04_oper.nc > ncdump/interim_2016-04_oper.ncdump
       ncdump /mnt/10.11.12.232/sat_auxdata/model/ecmwf/0.125-deg/interim_2016-05_oper.nc > ncdump/interim_2016-05_oper.ncdump
       ncdump /mnt/10.11.12.232/sat_auxdata/model/ecmwf/0.125-deg/interim_2016-06_oper.nc > ncdump/interim_2016-06_oper.ncdump
       ncdump /mnt/10.11.12.232/sat_auxdata/model/ecmwf/0.125-deg/interim_2016-07_oper.nc > ncdump/interim_2016-07_oper.ncdump

# return to the cal/val locations to create analysis timeseries (some locations missing too much of this timeseries might be later ignored)
wrkj ; mkdir ecmwf
       sort      all/gnss.tds1.txt.all.locate_1.0_calib    > gnss.tds1.txt.all.locate_1.0_calib.sort
       split -l 1000 gnss.tds1.txt.all.locate_1.0_calib.sort gnss.tds1.txt.all.locate_1.0_calib.sort
       sort      all/gnss.tds1.txt.all.locate_1.0_valid    > gnss.tds1.txt.all.locate_1.0_valid.sort
       split -l 1000 gnss.tds1.txt.all.locate_1.0_valid.sort gnss.tds1.txt.all.locate_1.0_valid.sort
       sort      all/gnss.tds1.txt.all.locate_1.0_extra    > gnss.tds1.txt.all.locate_1.0_extra.sort
       split -l 1000 gnss.tds1.txt.all.locate_1.0_extra.sort gnss.tds1.txt.all.locate_1.0_extra.sort
       parallel -j 16 julia /Home/ricani/bin/techdemosat.wind.speed.ecmwf.timeseries.jl ::: gnss.tds1.txt.all.locate_1.0_?????.sort?? ::: /mnt/10.11.12.232/sat_auxdata/model/ecmwf/0.125-deg
       rm gnss.tds1.txt.all.locate_1.0_?????.sor*

# verify that each subdir contains the expected number of files (e.g., 26210 + 25577 + 24982 = 76769 files with 1336 dates)
wrkj ; cd ecmwf ; ls -1 ecm* | grep -v GEM > z.list ; split -l 1000 z.list z.list ; cd .. ; wc *[a-z]/z.list

# create the forward and backward extrapolated timeseries
wrkj ; cd ecmwf ; ls z.list?? ; cd ..
       parallel -j 16 julia /Home/ricani/bin/techdemosat.wind.speed.ecmwf.timeseries.extrapolated.jl ::: ecmwf ::: z.listaa z.listab z.listac z.listad z.listae z.listaf z.listag z.listah z.listai z.listaj z.listak z.listal z.listam z.listan z.listao z.listap z.listaq z.listar z.listas z.listat z.listau z.listav z.listaw z.listax z.listay z.listaz z.listba z.listbb z.listbc z.listbd z.listbe z.listbf z.listbg z.listbh z.listbi z.listbj z.listbk z.listbl z.listbm z.listbn z.listbo z.listbp z.listbq z.listbr z.listbs z.listbt z.listbu z.listbv z.listbw z.listbx z.listby z.listbz z.listca z.listcb z.listcc z.listcd z.listce z.listcf z.listcg z.listch z.listci z.listcj z.listck z.listcl z.listcm z.listcn z.listco z.listcp z.listcq z.listcr z.listcs z.listct z.listcu z.listcv z.listcw z.listcx z.listcy

# assemble the insitu and analysis data for a triple collocation cal/val
wrkj ; parallel julia /Home/ricani/bin/techdemosat.wind.speed.assemble.insitu.jl all/gnss.tds1.txt.all ::: all/gnss.tds1.txt.all.locate_1.0_?????
       split -l 4000 all/gnss.tds1.txt.all.locate_1.0_calib_obs all/gnss.tds1.txt.all.locate_1.0_calib_obs
       split -l 4000 all/gnss.tds1.txt.all.locate_1.0_extra_obs all/gnss.tds1.txt.all.locate_1.0_extra_obs
       split -l 4000 all/gnss.tds1.txt.all.locate_1.0_valid_obs all/gnss.tds1.txt.all.locate_1.0_valid_obs
       parallel -j 16 julia /Home/ricani/bin/techdemosat.wind.speed.assemble.ecmwf.jl ::: all/gnss.tds1.txt.all.locate_1.0_?????_obs??
       cat all/gnss.tds1.txt.all.locate_1.0_calib_obs??.comb > all/gnss.tds1.txt.all.locate_1.0_calib_obs.comb
       cat all/gnss.tds1.txt.all.locate_1.0_extra_obs??.comb > all/gnss.tds1.txt.all.locate_1.0_extra_obs.comb
       cat all/gnss.tds1.txt.all.locate_1.0_valid_obs??.comb > all/gnss.tds1.txt.all.locate_1.0_valid_obs.comb
       mkdir all/zali.assemble ; mv all/*_1.0_?????_obs?? all/*_1.0_?????_obs??.comb all/zali.assemble

# perform global and local calibrations of the two extrapolations (BEF and AFT relative to NOW) using the extra collocation set
wrkj ; jjj techdemosat.wind.speed.extrapolated.histogram.jl all/gnss.tds1.txt.all.locate_1.0_extra_obs.comb uwnd
       jjj techdemosat.wind.speed.extrapolated.histogram.jl all/gnss.tds1.txt.all.locate_1.0_extra_obs.comb vwnd
       jjj techdemosat.wind.speed.extrapolated.histoplot.jl all/gnss.tds1.txt.all.locate_1.0_extra_obs.comb uwnd
       jjj techdemosat.wind.speed.extrapolated.histoplot.jl all/gnss.tds1.txt.all.locate_1.0_extra_obs.comb vwnd

# perform a paired triple collocation cal/val globally
wrkj ; jjj techdemosat.wind.speed.triple.paired.jl all/gnss.tds1.txt.all.locate_1.0_calib_obs.comb
       mkdir all/zali.recalib.paired ; mv all/*.cali.pair all/*.recalibrate all/zali.recalib.paired
       jjj analysis.evaluation.table.performance.paired.jl all/zali.recalib.paired/all.flux.daily.locate_2.0_calib.shfx*cali.pair
       jjj analysis.evaluation.table.performance.paired.jl all/zali.recalib.paired/all.flux.daily.locate_2.0_calib.lhfx*cali.pair
       jjj analysis.evaluation.table.performance.paired.jl all/zali.recalib.paired/all.flux.daily.locate_2.0_calib.wspd*cali.pair
       jjj analysis.evaluation.table.performance.paired.jl all/zali.recalib.paired/all.flux.daily.locate_2.0_calib.airt*cali.pair
       jjj analysis.evaluation.table.performance.paired.jl all/zali.recalib.paired/all.flux.daily.locate_2.0_calib.sstt*cali.pair
       jjj analysis.evaluation.table.performance.paired.jl all/zali.recalib.paired/all.flux.daily.locate_2.0_calib.shum*cali.pair
       cat all/zali.recalib.paired/*shfx*.md all/zali.recalib.paired/*lhfx*.md all/zali.recalib.paired/*wspd*.md all/zali.recalib.paired/*airt*.md all/zali.recalib.paired/*sstt*.md all/zali.recalib.paired/*shum*.md > analysis.evaluation.table.coefficients.md
       pandoc analysis.evaluation.table.coefficients.md -o analysis.evaluation.table.coefficients.html
       mv analysis.evaluation.table.coefficients.html analysis.evaluation.table.coefficients.0.00.html
       mv analysis.evaluation.table.coefficients.md   analysis.evaluation.table.coefficients.0.00.md
       mv all/zali.recalib.paired                                    all/zali.recalib.paired.0.00







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
