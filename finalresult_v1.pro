pro finalresult_v1, chains,outfile,savechains=savechains, chainsfile=chainsfile, binning=binning, trimstatement=trimstatement
;+
; NAME:
;       RECONSTRUCT_MULTI4
;       Part of MCMC_SD_V1 suite
;
; AUTHORSHIP:
;          SERGE DIETERICH, May 2018
;          Department of Terrestrial Magnetism,
;          Carnegie Institution of Washington
;
; PURPOSE:
;          1 - Plots histograms of the probability density functions for each parameter and saves files
;          2 - Computes the mean, median, and standard deviation of each parameter
;          3 - Allows trimming of samples that did not converge
;          4 - Allows for saving of the trimmed results 
;
; CALLING SEQUENCE:
;          finalresult_v1, chains, outfile[, /savechains, chainsfile="string.sav", binning=(10 element double array), trimstatement="string"] 
;
; INPUT: 
;          chains     = The array containing chains, generated by reconstruct_multi4.pro
;                             usually part of a structure: combined_results.chains
;          outfile    = String with the .ps file name to save plots: e.g. "filename.ps"
;          savechains = If set (/savechains) then the trimmed chains will be saved to a .sav file
;          chainsfile = If savechains is set then chainsfile is a string with the .sav file name
;                             where the trimmed chains will be saved: e.g. "trimmedchains.sav" 
;          binning    = A 10 element double precision vector with the histigram binning for every parameter, 
;                             in the standard index order.
;                             If not set then default values written in the code are used.
;     trimstatement   = A string to be used as the condition for a where() statement that trims bad samples.
;                             the format is, e.g, "chains0 gt 285" to trim all samples with parallax (parameter 0)
;                             greater than 285. Any conditional statement allowed in a where( ) call can be used.
;                             Use a very large or very small condition, or [set badindices = -1 in line 68 and don't
;                             set trimstatement in the call] for no trim.
;
; OUTPUT:
;         .ps file with probability density function histogram plots labeled with mean, median, and standard deviation,
;               in that order.
;          If savechains is set (/savechains) a .sav file with the trimmed chains.
;               Array format is dblarr(parameters, samples per chain, chains). Includes NaNs. 
;
;  PROCEDURES CALLED:
;          None
;
;  VERSION HISTORY:
;          CURRENT V1 - First publication version. May 2018.
;
;  BUGS:
;          None known. 
;-

sizearray = size(chains)

multiplot, /reset  
chains0 = dblarr(sizearray[2],sizearray[3])
chains1 = dblarr(sizearray[2],sizearray[3])
chains2 = dblarr(sizearray[2],sizearray[3])
chains3 = dblarr(sizearray[2],sizearray[3])
chains4 = dblarr(sizearray[2],sizearray[3])
chains5 = dblarr(sizearray[2],sizearray[3])
chains6 = dblarr(sizearray[2],sizearray[3])
chains7 = dblarr(sizearray[2],sizearray[3])
chains8 = dblarr(sizearray[2],sizearray[3])
chains9 = dblarr(sizearray[2],sizearray[3])



    chains0[*,*] = chains[0,*,*]
    chains1[*,*] = chains[1,*,*]
    chains2[*,*] = chains[2,*,*]
    chains3[*,*] = chains[3,*,*]
    chains4[*,*] = chains[4,*,*]
    chains5[*,*] = chains[5,*,*]
    chains6[*,*] = chains[6,*,*]
    chains7[*,*] = chains[7,*,*] * 180d/!pi
    chains8[*,*] = chains[8,*,*] * 180d/!pi
    chains9[*,*] = chains[9,*,*] * 180d/!pi


if ~keyword_set(trimstatement) then begin
  ;;;;;;;Enter condition for bad indices here ;;;;;;;;;;;;;;;;
  ;badindices = where(chains0 gt 285.)
  badindices = -1
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
endif else begin
  void = EXECUTE('badindices = where(' + trimstatement + ')')
endelse

help,badindices
if badindices[0] ne -1 then begin
  chains0[badindices] = !values.d_nan
  chains1[badindices] = !values.d_nan
  chains2[badindices] = !values.d_nan
  chains3[badindices] = !values.d_nan
  chains4[badindices] = !values.d_nan
  chains5[badindices] = !values.d_nan
  chains6[badindices] = !values.d_nan
  chains7[badindices] = !values.d_nan
  chains8[badindices] = !values.d_nan
  chains9[badindices] = !values.d_nan
endif

if ~keyword_set(binning) then begin
  histresult0 = histogram(chains0,binsize =  .05, locations = binstart0, /nan)
  histresult1 = histogram(chains1,binsize = .00001, locations = binstart1, /nan)
  histresult2 = histogram(chains2,binsize = .00002, locations = binstart2, /nan)
  histresult3 = histogram(chains3,binsize =   .1, locations = binstart3, /nan)
  histresult4 = histogram(chains4,binsize =  1., locations = binstart4, /nan)
  histresult5 = histogram(chains5,binsize = .0005, locations = binstart5, /nan)
  histresult6 = histogram(chains6,binsize =   1., locations = binstart6, /nan)
  histresult7 = histogram(chains7,binsize = .01, locations = binstart7, /nan)
  histresult8 = histogram(chains8,binsize = .05, locations = binstart8, /nan)
  histresult9 = histogram(chains9,binsize = .01, locations = binstart9, /nan)
endif else begin
  histresult0 = histogram(chains0,binsize = binning[0], locations = binstart0, /nan)
  histresult1 = histogram(chains1,binsize = binning[1], locations = binstart1, /nan)
  histresult2 = histogram(chains2,binsize = binning[2], locations = binstart2, /nan)
  histresult3 = histogram(chains3,binsize = binning[3], locations = binstart3, /nan)
  histresult4 = histogram(chains4,binsize = binning[4], locations = binstart4, /nan)
  histresult5 = histogram(chains5,binsize = binning[5], locations = binstart5, /nan)
  histresult6 = histogram(chains6,binsize = binning[6], locations = binstart6, /nan)
  histresult7 = histogram(chains7,binsize = binning[7], locations = binstart7, /nan)
  histresult8 = histogram(chains8,binsize = binning[8], locations = binstart8, /nan)
  histresult9 = histogram(chains9,binsize = binning[9], locations = binstart9, /nan)
endelse


set_plot, 'ps'
device,filename= outfile, /landscape

mean0 = strtrim(mean(chains0, /nan),2)
median0 = strtrim(median(chains0),2)
stdv0 = strtrim(stddev(chains0, /nan),2)
plot, binstart0,histresult0,psym=10,title="Parallax,  "+mean0+'  '+median0+'  '+stdv0,xtitle="mas",charsize=1.5,thick=2,charthick=2

mean1 = strtrim(mean(chains1, /nan),2)
median1 = strtrim(median(chains1),2)
stdv1 = strtrim(stddev(chains1, /nan),2)
plot, binstart1,histresult1,psym=10,title="RA Proper Motion, "+mean1+'  '+median1+'  '+stdv1,xtitle="mas / day",charsize=1.5,thick=2,charthick=2

mean2 = strtrim(mean(chains2, /nan),2)
median2 = strtrim(median(chains2),2)
stdv2 = strtrim(stddev(chains2, /nan),2)
plot, binstart2,histresult2,psym=10,title="DEC Proper Motion, "+mean2+'  '+median2+'  '+stdv2,xtitle="mas / day",charsize=1.5,thick=2,charthick=2

mean3 = strtrim(mean(chains3, /nan),2)
median3 = strtrim(median(chains3),2)
stdv3 = strtrim(stddev(chains3, /nan),2)
plot, binstart3,histresult3,psym=10,title="Semi-major Axis, "+mean3+'  '+median3+'  '+stdv3,xtitle="mas",charsize=1.5,thick=2,charthick=2

mean4 = strtrim(mean(chains4, /nan),2)
median4 = strtrim(median(chains4),2)
stdv4 = strtrim(stddev(chains4, /nan),2)
plot, binstart4,histresult4,psym=10,title="Orbital Period, "+mean4+'  '+median4+'  '+stdv4,xtitle="days",charsize=1.5,thick=2,charthick=2

mean5 = strtrim(mean(chains5, /nan),2)
median5 = strtrim(median(chains5),2)
stdv5 = strtrim(stddev(chains5, /nan),2)
plot, binstart5,histresult5,psym=10,title="Orbital Eccentricity, "+mean5+'  '+median5+'  '+stdv5,charsize=1.5,thick=2,charthick=2

mean6 = strtrim(mean(chains6, /nan),2)
median6 = strtrim(median(chains6),2)
stdv6 = strtrim(stddev(chains6, /nan),2)
plot, binstart6,histresult6,psym=10,title="Time of Periastron, "+mean6+'  '+median6+'  '+stdv6,xtitle="JD", charsize=1.5,thick=2,charthick=2

mean7 = strtrim(mean(chains7, /nan),2)
median7 = strtrim(median(chains7),2)
stdv7 = strtrim(stddev(chains7, /nan),2)
plot, binstart7,histresult7,psym=10,title="PA Line of Nodes, "+mean7+'  '+median7+'  '+stdv7,xtitle="Degree", charsize=1.5,thick=2,charthick=2

mean8 = strtrim(mean(chains8, /nan),2)
median8 = strtrim(median(chains8),2)
stdv8 = strtrim(stddev(chains8, /nan),2)
plot, binstart8,histresult8,psym=10,title="PA Periastron, "+mean8+'  '+median8+'  '+stdv8,xtitle="Degree", charsize=1.5,thick=2,charthick=2

mean9 = strtrim(mean(chains9, /nan),2)
median9 = strtrim(median(chains9),2)
stdv9 = strtrim(stddev(chains9, /nan),2)
plot, binstart9,histresult9,psym=10,title="Inclination, "+mean9+'  '+median9+'  '+stdv9,xtitle="Degree", charsize=1.5,thick=2,charthick=2


device, /close
set_plot, 'x'

if keyword_set(savechains) then begin
  outchains = dblarr(10,sizearray[2],sizearray[3])
  outchains[0,*,*] = chains0[*,*]
  outchains[1,*,*] = chains1[*,*]
  outchains[2,*,*] = chains2[*,*]
  outchains[3,*,*] = chains3[*,*]
  outchains[4,*,*] = chains4[*,*]
  outchains[5,*,*] = chains5[*,*]
  outchains[6,*,*] = chains6[*,*]
  outchains[7,*,*] = chains7[*,*]
  outchains[8,*,*] = chains8[*,*]
  outchains[9,*,*] = chains9[*,*]

  save,outchains,filename=chainsfile, description="Array format is dblarr(parameters, samples per chain, chains). Includes NaNs."
endif

end



