function MCMC_SD_WRAPPER_V1, astrofile=astrofile,n_chains=n_chains,n_iter=n_iter,n_runs=n_runs,ra_jitter=ra_jitter,dc_jitter=dc_jitter,startvals=startvals, $
                             scalingpar=scalingpar,stepscale=stepscale,stepmultiplier=stepmultiplier,tictoc=tictoc,inputbounds=inputbounds, $
                             parbounds=parbounds, commononly=commononly
;+
; NAME:
;       MCMC_SD_WRAPPER_V1
;       Part of MCMC_SD_V1 suite  
;
; AUTHORSHIP:
;          SERGE DIETERICH, May 2018
;          Department of Terrestrial Magnetism,
;          Carnegie Institution of Washington  
;  
; PURPOSE:
;          User level calling code for the MCMC_SD astrometric solver suite.
;          Establishes user set parameters, establishes common
;          variable block, runs multiple MCMC sub-runs, saves results
;          to .sav files, and outputs the last sub-run.
;
; CALLING SEQUENCE:
;          All parameters/keywords will default to those hard-wired in
;          the code below if not specified in the call sequence. The
;          call sequence can therefore contain all, any, or no
;          parameter calls:
;          Minimal calling sequence if starting from scratch:
;                                            results = mcmc_sd_wrapper_v1()
;          Calling sequence if picking up from a previous run:
;                                            results2 = mcmc_sd_wrapper_v1(startvals = results.chains)
;                                            where results.chains may be restored from a .sav file.
;
; INPUTS:
;            astrofile  =    The input data file containing the astrometric data. Format:
;                            "JD_observation   RA_displacement_mas  DEC_displacement_mas  RA_error  DEC_error	 RA_parallax_factor  DEC_parallax_factor"
;                            where RA and DEC displacement should be measured as offsets from the first epoch of observation.
;                            JD_obsewrvation must be in increasing order.
;                       
;            n_chains   =    The number of chains per session. We have been using a total of 50 chains
;                            for each result, meaning this number should be 25 if running 2 simultaneous
;                            sessions or 17 if running 3 sessions (for a total of 51 in that case).
;                            NOTE: All sessions must have the same number of chains if they are to
;                            be combined by reconstruct_multi.pro
;                            
;            n_iter     =    The number of steps per each sub-run. We have been using 10,000.
;                            n_iter x n_chains determines how much memory is needed.
;            
;            n_runs     =    The number of sub-runs in a session. The last 10 sub-runs are saved.
;                            Increasing the number of sub-runs allows one to decrease n_iter and
;                            thereby decrease memory needs.
;            
;            ra_jitter  =    Systemic RA error in milli-arcseconds to be added in quadrature to
;                            the errors in the astrometry data file
;
;            dc_jitter  =    Systemic DEC error.
;                            
;            startvals  =    If specified it will make the run start from previous results
;                            as opposed to random numbers. To start from previous
;                            results in structure previousresults do
;                            startvals=previousresults.chains .
;
;            scalingpar  =   The astrometric parameter used to set the step scale. We have been using
;                            0 (parallax), 3 (semi-major axis) also seems to work.
;            
;            stepscale   =   Sets the scale of the step size distribution for the parameter specified
;                            in scalingpar. We have been using 1, meaning steps vary by about 1 mas,
;                            but see the next parameter.
;                            
;        stepmultiplier  =   Operates in stepscale to form a distribution. We have been using
;                            stepmultiplier = 10, meaning that for every step a random
;                            number between 0 and 10 is generated and there are equal chances that
;                            stepscale will be multiplied or divided by that number  and then
;                            multiplied by either 1 or -1.
;                            
;         inputbounds    =   A 10x2 array that indicates the interval from which starting guesses
;                            should be drawn when initiating each chain. If startvals is
;                            set then inputbounds is ignored. Example:
;                            
;			               inputbounds = [[265d,285d],[10.7d,11d],[-6.8d,-6.5d],[150,350d],[3000d,7500d],[0.35d,0.6d],[2.45d6, min_jd], $
;                                           [120d*!pi/180d,180d*!pi/180d],[0d,2d*!pi],[70d*!pi/180d,!pi/2d],[150d,350d]]
;                             
;                            means that starting values for parameter 0 (parallax) will be
;                            drawn from the 265 to 285 range, and so on. Note that parameter
;                            6, the epoch of periastron, must always preceed the first epoch of
;                            observation. Because this parameter can be rather lengthy
;                            you may wish to either set it in the code or set it to an array
;                            name that can then be entered repeatedly instead of the explicit
;                            expression.
;                            
;             parbounds   =  The interval in which each parameter will be allowed to fluctuate.
;                            Same format as inputbounds. It parbounds
;                            should not be too restrictive because that takes away the MCMC's
;                            ability to find its way back to higher probabilities when
;                            a parameter happens to be off by a large amount. A good rule is to
;                            enter the plausible range given the observational setup
;			                       without taking into account specific
;			                       knowledge about the system you are
;			                       fitting. If the upper and lower limits
;			                       for a parameter are identical then that
;			                       parameter is considered a constant.
;  
;           
;                tictoc  =   If /tictoc is included in the call line then the ellapsed time for the entire
;                            run will be output by the IDL tic toc command at the end of the run.
;                            Be sure to check that your IDL implementation supports tic toc:
;			                       IDL> tic
;			                       IDL> toc
;			                       Should print the ellapsed time. 
;
;             commononly  =  If set (/commononly) then this wrapper
;                            code will establish the common variable block and exit
;                            without running any codes. This is useful
;                            in case you need to run codes individually.
;
;   
; OUTPUTS:
;        The function returns the results of the last sub-run and
;        saves the results of the last 10 sub-runs in files
;        last1_multi_results.sav ... last10_multi_results.sav
;        For a 50K iteration sub-run of 50 chains the output structure
;        is of the form:
;        ** Structure <4da368>, 3 tags, length=230000000, data length=230000000, refs=1:
;        CHAINS          DOUBLE    Array[10, 50000, 50]
;        LNPCHAINS       DOUBLE    Array[50000, 50]
;        MATCHEPOCH      LONG      Array[50000, 50]  
;        where results.chains contains the MCMC sampling results from
;        which the probability denstity fucntion will be
;        drawn. LNPCHAINS and MATCHEPOCH are mainly for
;        troubleshooting purposes.  
;
; PROCEDURES CALLED:
;        mcmcsd_spidersolver_v1(parinputs,parbounds,n_iter, n_chains,scalingpar,stepscale,stepmultiplier)  
;        plothistograms_v1,results.chains,parbounds
; 
; VERSION HISTORY:
;        CURRENT V1 - First publication version, May 2018. 
;          
;  
; BUGS:
;        Still de-bugging
;        Rather inflexible in what it saves: always the last 10 sub-runs.
;        Could use better graphical representation of sub-run histograms.
;        Erases previous results with same file names without asking.
;-   

spawn, 'rm -f last1_multi_results.sav, last2_multi_results.sav, last3_multi_results.sav, last4_multi_results.sav', garbage
spawn, 'rm -f last5_multi_results.sav, last6_multi_results.sav, last7_multi_results.sav, last8_multi_results.sav', garbage
spawn, 'rm -f last9_multi_results.sav, last10_multi_results.sav', garbage

print, "Generating mcmc_sd common variable block . . . "
  common newsolver_common_block, data_jd, data_ra, data_dc, data_raerr, data_dcerr, $
    min_jd, max_jd, pifact_ra, pifact_dc, Ec_array, sin_ec, cos_ec, n_epochs, Ec_array_lstndx, n_epochs_m1, $
    Ec_array2,sin_ec2,cos_ec2,Ec_array2_lstndx,Ec_array3,sin_ec3,cos_ec3,Ec_array3_lstndx, $
    Ec_array4,sin_ec4,cos_ec4,Ec_array4_lstndx


;;; List of parameter indices
;;;
;;; 0 - Parallax in mas
;;; 1 - RA proper motion in mas / day
;;; 2 - DEC proper motion in mas / day
;;; 3 - CTIOPI Photocentric semi-major axis in mas
;;; 4 - Period in days
;;; 5 - Eccentricity
;;; 6 - Epoch of periastron in days. Must be before first epoch and preferably within one period before it.
;;; 7 - OMEGA, position angle of the line of nodes E of N, in radians
;;; 8 - omega, position angle of periastron E of N, in radians
;;; 9 - inclination in radians
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; START USER SET PARAMETERS HERE ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;; DEFAULTS BELOW ARE SPECIFIC TO DIETERICH ET AL 2018 AND NOT MEANT TO BE GENERAL USE;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;CHANGE THESE INTO SOMETHING THAT MAKES SENSE FOR YOUR CASE;;;;;;;;;;;;;;;;;;;;;;



if ~keyword_set(astrofile) then astrofile = "ctiopi_nightly_means.txt"

  readcol,astrofile,data_jd, data_ra, data_dc, data_raerr, data_dcerr, pifact_ra, pifact_dc, $
           format =  "d,       d,       d,        d,         d,            d,         d"

 min_jd = double(min(data_jd))      
 max_jd = double(max(data_jd))

if ~keyword_set(ra_jitter) then ra_jitter = 0d  
if ~keyword_set(dc_jitter) then dc_jitter = 0d
data_raerr = (data_raerr^2d + ra_jitter^2d)^.5d
data_dcerr = (data_dcerr^2d + dc_jitter^2d)^.5d

if ~keyword_set(n_iter) then n_iter =  long(1e4)                
n_iter_m1 = n_iter - 1     

if ~keyword_set(n_chains) then n_chains = 17   ;;; 50 total is good. So 25 if running 2 simultaneous sessions, 17 if 3.
n_chains_m1 = n_chains - 1  

if ~keyword_set(scalingpar) then scalingpar = 3  ;;; See list of parameter indices above.
if ~keyword_set(stepscale) then stepscale = 1d   
if ~keyword_set(stepmultiplier) then stepmultiplier = 10d
if ~keyword_set(n_runs) then n_runs = 200  
n_runs_m1 = n_runs - 1    
n_runs_m2 = n_runs_m1 -1  

;; n_chains random starting values for each parameter will be generated within inputbounds
;; parameter                                       0           1             2           3           4            
if ~keyword_set(inputbounds) then inputbounds = [[265d,285d],[10.7d,11d],[-6.8d,-6.5d],[150,350d],[3000d,7500d], $
                                                [0.35d,0.6d],[2.45d6, min_jd],[120d*!pi/180d,180d*!pi/180d],[0d,2d*!pi],[70d*!pi/180d,!pi/2d]]
;                                                5              6                      7                    8                9      
;
;; parameters will be allowed to fluctuate within parbounds. The limit is considered a bad value. If a step touches its bound then the next step to be tested
;; will be a random value for that parameter (but not the other parameters) within parbounds.
;; parameters                                    0           1           2         3            4            5      
if ~keyword_set(parbounds) then parbounds = [[250d,300d],[10d,11.5d],[-7d,-6d],[100d,600d],[1000d,8000d],[0.2d,0.7d], $
                                             [2.448e6, min_jd],[100d*!pi/180d,200d*!pi/180d],[0,2d*!pi],[40d*!pi/180d,3d*!pi/2d]]
;                                                      6               7                          8                9     
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; END USER SET PARAMETERS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; RUN YOUR CODE! ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; Define common arrays
  
  print,"Generating common block arrays."
  Ec_array = (dindgen(1d2)+0.5d) * (2d *!pi / (1d2 + 0.5d))
  Ec_array_lstndx = n_elements(Ec_array) -1L 

  Ec_array2 = (dindgen(1d4)+0.5d) * (2d *!pi / (1d4 + 0.5d))
  Ec_array2_lstndx = n_elements(Ec_array2) -1L
  Ec_array2 = [Ec_array2[(Ec_array2_lstndx-100L):Ec_array2_lstndx],Ec_array2, Ec_array2[0:100L]]
  Ec_array2_lstndx = n_elements(Ec_array2) -1L

  Ec_array3 = (dindgen(1d6)+0.5d) * (2d *!pi / (1d6 + 0.5d))
  Ec_array3_lstndx = n_elements(Ec_array3) -1L
  Ec_array3 = [Ec_array3[(Ec_array3_lstndx-100L):Ec_array3_lstndx],Ec_array3, Ec_array3[0:100L]]
  Ec_array3_lstndx = n_elements(Ec_array3) -1L

  Ec_array4 = (dindgen(1d8)+0.5d) * (2d *!pi / (1d8 + 0.5d))
  Ec_array4_lstndx = n_elements(Ec_array4) -1L
  Ec_array4 = [Ec_array4[(Ec_array4_lstndx-100L):Ec_array4_lstndx],Ec_array4, Ec_array4[0:100L]]
  Ec_array4_lstndx = n_elements(Ec_array4) -1L

  ;; The next 2 lines add random noise ;; PROBABLY NOT NEEDED
  ;Ec_arr4_noise = (randomu(seed,n_elements(Ec_array4)) * 2d)  - 1d
  ;for i =1L, Ec_array4_lstndx do Ec_array4[i] = Ec_array4[i] + (Ec_array4[i] - Ec_array4[i-1L]) *  Ec_arr4_noise[i]

  sin_ec   = sin(Ec_array)
  cos_ec   = cos(Ec_array)
  sin_ec2   = sin(Ec_array2)
  cos_ec2   = cos(Ec_array2)
  sin_ec3   = sin(Ec_array3)
  cos_ec3   = cos(Ec_array3)
  sin_ec4   = sin(Ec_array4)
  cos_ec4   = cos(Ec_array4)

  n_epochs = n_elements(data_jd)
  n_epochs_m1 = n_epochs - 1

 
if keyword_set(startvals) then begin
  parinputs = dblarr(10,n_chains)
  for i=0,9 do begin
    for j=0,n_chains_m1 do parinputs[i,j] = startvals[i,-1,j]
  endfor
endif else begin
  parinputs = randomu(seed,10,n_chains)
  for i =0, 9 do parinputs[i,*] = parinputs[i,*] * (inputbounds[1,i] - inputbounds[0,i]) + inputbounds[0,i]
endelse

if keyword_set(commononly) then begin
   print, "COMMON VARIABLE BLOCK ESTABLISHED."
   print, "RETURNING NaN AND EXITING."
   return, !values.d_nan
endif

if keyword_set(tictoc) then tic

print,"Starting run 1"

results_firstiter =  mcmcsd_spidersolver_v1(parinputs,parbounds,n_iter, n_chains,scalingpar,stepscale,stepmultiplier)

;;output = create_struct('chains',parvalues, 'lnPchains',lnPValues,  'RAmodel',modRAvalues, $
;  'DCmodel',modDCvalues, 'partials',partialD,'Matchepoch',matchedepochs) ;,'LastChange',lastchange,'UsedSteps',usedsteps)

print,"Mean results after this sub-run: "
print,mean(results_firstiter.chains[0,*,*]),mean(results_firstiter.chains[1,*,*]),mean(results_firstiter.chains[2,*,*]),mean(results_firstiter.chains[3,*,*]), $
  mean(results_firstiter.chains[4,*,*]),mean(results_firstiter.chains[5,*,*]),mean(results_firstiter.chains[6,*,*]),mean(results_firstiter.chains[7,*,*]), $
  mean(results_firstiter.chains[8,*,*]),mean(results_firstiter.chains[9,*,*])
  print,"Median results after this sub-run: "
  print,median(results_firstiter.chains[0,*,*]),median(results_firstiter.chains[1,*,*]),median(results_firstiter.chains[2,*,*]),median(results_firstiter.chains[3,*,*]), $
    median(results_firstiter.chains[4,*,*]),median(results_firstiter.chains[5,*,*]),median(results_firstiter.chains[6,*,*]),median(results_firstiter.chains[7,*,*]), $
    median(results_firstiter.chains[8,*,*]),median(results_firstiter.chains[9,*,*])
  print,"Standard deviation of results after this sub-run: "
  print,stddev(results_firstiter.chains[0,*,*]),stddev(results_firstiter.chains[1,*,*]),stddev(results_firstiter.chains[2,*,*]),stddev(results_firstiter.chains[3,*,*]), $
    stddev(results_firstiter.chains[4,*,*]),stddev(results_firstiter.chains[5,*,*]),stddev(results_firstiter.chains[6,*,*]),stddev(results_firstiter.chains[7,*,*]), $
    stddev(results_firstiter.chains[8,*,*]),stddev(results_firstiter.chains[9,*,*])

print,"Matches to epoch error bar for each chain: "
matchvector = lonarr(n_chains)
for c = 0, n_chains_m1 do matchvector[c] = max(results_firstiter.matchepoch[*,c])
print,matchvector

plothistograms_v1,results_firstiter.chains,parbounds

save,results_firstiter, filename="results_firstiter.sav",description="First sub-run"

for par = 0, 9 do begin
  for v =0, n_chains_m1 do parinputs[par,v] = results_firstiter.chains[par,-1,v]
endfor
delvar,results_firstiter

for i=1,n_runs_m1 do begin
  Print,"Starting run ", long(i+1)

  results = mcmcsd_spidersolver_v1(parinputs,parbounds,n_iter,n_chains,scalingpar,stepscale,stepmultiplier)

    for par = 0, 9 do begin
      for v =0, n_chains_m1 do parinputs[par,v] = results.chains[par,-1,v]
    endfor
    print,"Mean results after this sub-run: "
    print,mean(results.chains[0,*,*]),mean(results.chains[1,*,*]),mean(results.chains[2,*,*]),mean(results.chains[3,*,*]), $
      mean(results.chains[4,*,*]),mean(results.chains[5,*,*]),mean(results.chains[6,*,*]),mean(results.chains[7,*,*]), $
      mean(results.chains[8,*,*]),mean(results.chains[9,*,*])
    print,"Median results after this sub-run: "
    print,median(results.chains[0,*,*]),median(results.chains[1,*,*]),median(results.chains[2,*,*]),median(results.chains[3,*,*]), $
      median(results.chains[4,*,*]),median(results.chains[5,*,*]),median(results.chains[6,*,*]),median(results.chains[7,*,*]), $
      median(results.chains[8,*,*]),median(results.chains[9,*,*])
    print,"Standard deviation of results after this sub-run: "
    print,stddev(results.chains[0,*,*]),stddev(results.chains[1,*,*]),stddev(results.chains[2,*,*]),stddev(results.chains[3,*,*]), $
      stddev(results.chains[4,*,*]),stddev(results.chains[5,*,*]),stddev(results.chains[6,*,*]),stddev(results.chains[7,*,*]), $
      stddev(results.chains[8,*,*]),stddev(results.chains[9,*,*])

    print,"Matches to epoch error bars for each chain: "
    matchvector = lonarr(n_chains)
    for c = 0, n_chains_m1 do matchvector[c] = max(results.matchepoch[*,c])
    print,matchvector
     plothistograms_v1,results.chains,parbounds

    spawn, 'mv -f last9_multi_results.sav last10_multi_results.sav', garbage     ;; keeps the last 10 runs
    spawn, 'mv -f last8_multi_results.sav last9_multi_results.sav', garbage      
    spawn, 'mv -f last7_multi_results.sav last8_multi_results.sav', garbage
    spawn, 'mv -f last6_multi_results.sav last7_multi_results.sav', garbage
    spawn, 'mv -f last5_multi_results.sav last6_multi_results.sav', garbage
    spawn, 'mv -f last4_multi_results.sav last5_multi_results.sav', garbage
    spawn, 'mv -f last3_multi_results.sav last4_multi_results.sav', garbage
    spawn, 'mv -f last2_multi_results.sav last3_multi_results.sav', garbage
    spawn, 'mv -f last1_multi_results.sav last2_multi_results.sav', garbage
    save,results,filename="last1_multi_results.sav"
    delvar,results  
endfor

print,"End of all runs. Use reconstruct_multi.pro to combine results from the last 10 sub-runs."

if keyword_set(tictoc) then toc

return,results

end
