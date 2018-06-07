pro plothistograms_v1, inchains,bounds
;+
; NAME:
;       PLOTHISTOGRAMS_V1
;       Part of MCMC_SD_V1 suite
;
; AUTHORSHIP:
;          SERGE DIETERICH, May 2018
;          Department of Terrestrial Magnetism,
;          Carnegie Institution of Washington
;
; PURPOSE:
;          Plots the histograms of the probability density function
;          of each parameter after a sub-run. Meant as a graphical way
;          to judge convergence and evolution.
;
; CALLING SEQUENCE:
;          plothistograms_v1, inchains, bounds
;          All required
;          Called by mcmc_sd_wrapper_v1.pro
;          
; INPUTS:
;          inchains  = The Markov chains produced in the last sub-run.
;                          called from structure usually in form results.chains 
;          bounds    =  The X range for the histogram of each parameter. 
;                          usually set to parbounds so as to cover entire range for each parameter.
;                          see mcmc_sd_wrapper_v1.pro
;                                   
; OUTPUT:  
;          Graphic output with histogram grid.
; 
; PROCEDURES CALLED:
;          None
;
;  VERSION HISTORY:
;          CURRENT V1 - First publication version. May 2018.
;
;  BUGS:
;         Plot are rudimentary.
;         Could be much better using IDL new graphics.
;-
histresult0 = histogram(inchains[0,*,*],binsize =   .1, locations = binstart0)
histresult1 = histogram(inchains[1,*,*],binsize = .001, locations = binstart1)
histresult2 = histogram(inchains[2,*,*],binsize = .001, locations = binstart2)
histresult3 = histogram(inchains[3,*,*],binsize =   1, locations = binstart3)
histresult4 = histogram(inchains[4,*,*],binsize =  10, locations = binstart4)
histresult5 = histogram(inchains[5,*,*],binsize = .001, locations = binstart5)
histresult6 = histogram(inchains[6,*,*],binsize =   1, locations = binstart6)
histresult7 = histogram(inchains[7,*,*]*180./!pi,binsize = 1, locations = binstart7)
histresult8 = histogram(inchains[8,*,*]*180./!pi,binsize = 1, locations = binstart8)
histresult9 = histogram(inchains[9,*,*]*180./!pi,binsize = 1, locations = binstart9)

arraysize = size(inchains)
if arraysize[1] eq 11 then histresult10 = histogram(inchains[10,*,*],binsize = .01, locations = binstart10)

cgplot,  /reset
multiplot, /reset
cgerase & multiplot, [4,3]; , /square 

cgplot,binstart0,histresult0,psym=10,xrange=[bounds[0,0],bounds[1,0]]
multiplot

cgplot,binstart1,histresult1,psym=10,xrange=[bounds[0,1],bounds[1,1]]
multiplot

cgplot,binstart2,histresult2,psym=10,xrange=[bounds[0,2],bounds[1,2]]
multiplot

cgplot,binstart3,histresult3,psym=10,xrange=[bounds[0,3],bounds[1,3]]
multiplot

cgplot,binstart4,histresult4,psym=10,xrange=[bounds[0,4],bounds[1,4]]
multiplot

cgplot,binstart5,histresult5,psym=10,xrange=[bounds[0,5],bounds[1,5]]
multiplot

cgplot,binstart6,histresult6,psym=10,xrange=[bounds[0,6],bounds[1,6]]
multiplot

cgplot,binstart7,histresult7,psym=10,xrange=[bounds[0,7]*180./!pi,bounds[1,7]*180./!pi]
multiplot

cgplot,binstart8,histresult8,psym=10,xrange=[bounds[0,8]*180./!pi,bounds[1,8]*180./!pi]
multiplot

cgplot,binstart9,histresult9,psym=10,xrange=[bounds[0,9]*180./!pi,bounds[1,9]*180./!pi]
multiplot

if arraysize[1] eq 11 then begin
  cgplot,binstart10,histresult10,psym=10,xrange=[bounds[0,10],bounds[1,10]]
  multiplot
endif


end


