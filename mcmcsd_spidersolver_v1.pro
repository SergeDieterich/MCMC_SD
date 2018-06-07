function mcmcsd_spidersolver_v1, parinputs,parbounds,n_iter, n_chains,scalingpar,stepscale,stepmultiplier
;+
; NAME:
;       MCMCSD_SPIDERSOLVER_V1
;       Part of MCMC_SD_V1 suite  
;
; AUTHORSHIP:
;          SERGE DIETERICH, May 2018
;          Department of Terrestrial Magnetism,
;          Carnegie Institution of Washington  
;  
; PURPOSE:
;          The MCMC iterator for the MCMC_SD suite. Its functions are:
;          1 - Produces the MCMC steps.
;          2 - Runs the astrometric model with new values based on
;              the latest step.  
;          3 - Decides if a location in parameter space is accepted
;              based on the results returned by the astrometric model.
;          4 - Returns the Markov chains.  
;
; CALLING SEQUENCE:
;          result = mcmcsd_spidersolver_v1(parinputs,parbounds,n_iter, n_chains,scalingpar,stepscale,stepmultiplier)
;          All parameters required.   
;          Meant to be called by MCMC_SD_WRAPPER_V1, which also
;          establishes the necessary common variable block.  
;         
; INPUTS:
;            parinputs   =  Starting location for each parameter for each
;                           chain, calculated by MCMC_SD_WRAPPER.
;;            n_chains   =  The number of chains per session. We have been using a total of 50 chains
;                            for each result, meaning this number should be 25 if running 2 simultaneous
;                            sessions or 17 if running 3 sessions (for a total of 51 in that case).
;                            NOTE: All sessions must have the same number of chains if they are to
;                            be combined by reconstruct_multi.pro
;
;            n_iter     =    The number of steps per each sub-run. We have been using 10,000.
;                            n_iter x n_chains determines how much memory is needed.
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
;             parbounds   =  The interval in which each parameter will be allowed to fluctuate.
;                            Same format as inputbounds (see MCMC_SD_WRAPPER_V1.pro). Parbounds
;                            should not be too restrictive because that takes away the MCMC's
;                            ability to find its way back to higher probabilities when
;                            a parameter happens to be off by a large amount. A good rule is to
;                            enter the plausible range given the observational setup
;                            without taking into account specific
;                            knowledge about the system you are
;                            fitting. If the upper and lower limits
;                            for a parameter are identical then that
;                            parameter is considered a constant.
;    
;    
;   
; OUTPUTS:
;        The function returns the results of the last sub-run
;        to the wrapper function, which then saves it.
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
;        astromodel.pro
;        partials.pro  
;        MCMC_SD_V1 common block must already be established by mcmc_sd_wrapper_v1.pro
;        
; VERSION HISTORY:
;        CURRENT V1 - First publication version. May 2018.
;  
; BUGS:
;        Eccentricity (parameter index 5) still not well scaled, but
;        seems to work.
;        Does not establish common block automatically if called without wrapper.
;            Use result = mcmc_sd_wrapper_v1(/commononly)
;-   

n_chains_m1 = n_chains - 1
n_iter_m1 = n_iter - 1

  common newsolver_common_block, data_jd, data_ra, data_dc, data_raerr, data_dcerr, $
    min_jd, max_jd, pifact_ra, pifact_dc, Ec_array, sin_ec, cos_ec, n_epochs, Ec_array_lstndx, n_epochs_m1, $
    Ec_array2,sin_ec2,cos_ec2,Ec_array2_lstndx,Ec_array3,sin_ec3,cos_ec3,Ec_array3_lstndx, $
    Ec_array4,sin_ec4,cos_ec4,Ec_array4_lstndx

fixedpars =  where(parbounds[0,*] eq parbounds[1,*])  

;;;;;;; Establishes step array  ;;;;;;;;;;;;;;
steparray = randomu(seed,10,n_iter,n_chains)*stepmultiplier

stepexp =   randomu(seed,10,n_iter,n_chains)*2 -1
stepexp = stepexp/abs(stepexp)

n_iter_m5 = n_iter - 5
for x = 0, n_iter_m5, 4 do stepexp[*,x,*] = -1d    ;; forces all steps to be < 1 at every 4th step.

steparray = steparray^stepexp

delvar,stepexp           ;; free up memory after operations are done.

stepsign = randomu(seed,10,n_iter,n_chains)*2 -1
stepsign = stepsign/abs(stepsign)

steparray = steparray * stepsign  ;; Equal chances of divide or multiply, + or -

delvar,stepsign

if fixedpars[0] ne -1 then  steparray[fixedpars,*,*] = 0d    ;; If the parameter is fixed then all steps for that parameter are zero

ln_chainadvance = alog(randomu(seed,n_iter,n_chains))    ;; random number to be compared to ratio of probabilities

print, "Iterating MCMC."

;; For every iteration the following must be recorded:
;; The input parameters
;; The ln(P)
;; The RA and DC models outputs

parvalues    =dblarr(10,n_iter,n_chains)
lnPvalues    =dblarr(n_iter,n_chains)
modRAvalues  =dblarr(n_epochs,n_iter,n_chains)
modDCvalues  =dblarr(n_epochs,n_iter,n_chains)
partialD     =dblarr(10,n_iter,n_chains)
matchedepochs = lonarr(n_iter,n_chains)
spidercounter = intarr(n_chains)
spiderleg     = intarr(n_chains)
;; Do the first iteration, which does not require partial derivative scaling.

for j = 0, n_chains_m1 do begin 
  if fixedpars[0] ne -1 then parinputs[fixedpars,j] = parbounds[0,fixedpars]   ;; Fixes all fixed parameters to their values.
   solverresults = astromodel(parinputs[*,j])
  ;; solverresults form: ('p_model_ra', model_ra, 'p_model_dc', model_dc, 'ln_prob', -ln_probability)

  parvalues[*,0,j]   = parinputs[*,j]
  lnPValues[0,j]     = solverresults.ln_prob
  modRAvalues[*,0,j] = solverresults.p_model_ra
  modDCvalues[*,0,j] = solverresults.p_model_dc
   dlnPdPar          = partials(parvalues[*,0,j],modRAvalues[*,0,j],modDCvalues[*,0,j])
  partialD[*,0,j]    = dlnPdPar

  parvalues[*,1,j]   = parvalues[*,0,j]    ;; So that the big loop can start at 2 and check if the values have changed in the previous 2 iterations
  lnPValues[1,j]     = lnPValues[0,j] 
  modRAvalues[*,1,j] = modRAvalues[*,0,j]
  modDCvalues[*,1,j] = modDCvalues[*,0,j]
  partialD[*,1,j]    = partialD[*,0,j] 
endfor

for i = 2, n_iter_m1 do begin   ;; This loop over iterations

 if i mod 200 eq 0 then begin 
   spiderstep =1               ;; Indicates that this will be a spider step for all chains.
;   spidercounter = spidercounter * 0    ;; zeroes the spider step counter.
 endif  else spiderstep = 0

 for j =0, n_chains_m1 do begin  ;; This loop is for chains
    if spidercounter[j] eq 100 then begin ;; Evaluates if the last spider step put the chain in an overall better vicinity.
       spidercounter[j] = 0
      if mean(lnPValues[(i-199):(i-101),j]) gt mean(lnPValues[(i-100):(i-1),j])  then begin 
        lnPValues[i-1,j] = lnPValues[i-52,j]   ;; go back to before spider jump if 50 step average got worse
                 print, '                              SPIDER RAN AWAY!',i,j,spiderleg[j];,mean(parvalues[spiderleg[j],(i-100):(i-1),j])
      endif else print,'              SPIDER FOUND NEW HOME!',i,j,spiderleg[j],mean(parvalues[spiderleg[j],(i-199):(i-101),j]),' ----> ',mean(parvalues[spiderleg[j],(i-100):(i-1),j])  
 ;;     stop
    endif
   jolt = 0  ;; if things were jolted changes to 1 and then step always needs to be improvement to accept
    if lnPValues[i-1,j] eq lnPValues[i-2,j]  or spiderstep eq 1 then begin  ;; compute new partials only if chain has advanced and it's not a spider. If not copy previous values.
      dlnPdPar = partialD[*,i-1,j]
    endif else begin
      dlnPdPar = partials(parvalues[*,i-1,j], modRAvalues[*,i-1,j], modDCvalues[*,i-1,j])       ;; get partial derivatives
    endelse
    partialD[*,i,j] = dlnPdPar
    partialscale = (dlnPdPar[scalingpar] /dlnPdPar)
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; partialscale[5] = 0.005       ;;;;; Because eccentricity does not scale nicely.
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    steps = steparray[*,i,j] * partialscale * stepscale ;; uncomment to use partials

    steps[5] = steps[5]/100d
    ;;      0   1      2   3  4    5   6   7        8        9
    ; stepscale = [ 1d,.001d,.001d,1d,1d,.001d,1d, 180d/!pi,180d/!pi,180d/!pi]

    ;; steps = steparray[*,i,j] *stepscale  ;; uncomment to not use partials

 ; endelse


 if spiderstep eq 1 then begin
  newparvalues = parvalues[*,i-1,j]
 spiderleg[j] = fix(randomu(seed)*11 - 1)
  spiderlegfixed = where(spiderleg[j] eq fixedpars)
  while  spiderlegfixed[0] ne -1 do begin            ;; This while loop makes sure a fixed parameter is not chosen as the spider parameter.
    spiderleg[j] = fix(randomu(seed)*11 - 1)
    spiderlegfixed = where(spiderleg[j] eq fixedpars)
  endwhile
                              
  newparvalues[spiderleg[j]] = (parbounds[1,spiderleg[j]] - parbounds[0,spiderleg[j]]) * randomu(seed)   + parbounds[0,spiderleg[j]] ;;; SPIDERSTEP!!!!    Good for fixed parameters too.
 endif else begin  ;; advance by regular step mechanism if it is not a spider step.
    newparvalues = (parvalues[*,i-1,j] + steps)    ;; STEP!
;;;;;;    if fixedpars[0] ne -1 then newparvalues[fixedpars] = parbounds[0,fixedpars]   ;; Fixes all fixed parameters to their values.   This line is not needed. Steps are zero for fixed pars.
    for k=0, 9 do begin
      while newparvalues[k] lt parbounds[0,k] or newparvalues[k] gt parbounds[1,k] do begin  ;; While loop changes out of bound parameter values
        jolt =1
        steps[k] = (parbounds[1,k] - parbounds[0,k]) * (randomu(seed)*2d - 1d) ;; step becomes +- up to 1 bound interval
        newparvalues[k] = parvalues[k,i-1,j] + steps[k]
      endwhile
    endfor
  endelse

  solverresults = astromodel(newparvalues)    ;; Run the model
  matchedepochs[i,j] = solverresults.goodmatches
  if matchedepochs[i,j] ge n_epochs * 0.68d and matchedepochs[i-1,j] lt n_epochs * 0.68d then print, "Chain ", j, "converged at iteration ",i

   ;;;;;;; Probability of advancing a step ;;;;;;;;;;;;;;;;;
     probratio = solverresults.ln_prob - lnPValues[i-1,j]
     if probratio ge ln_chainadvance[i,j] then begin ;; advance chain if ratio of probabilities >= uniform random number
       parvalues[*,i,j] = newparvalues
       lnPValues[i,j]  = solverresults.ln_prob
       modRAvalues[*,i,j] = solverresults.p_model_ra
       modDCvalues[*,i,j] = solverresults.p_model_dc
       if spiderstep eq 1 then begin
        spidercounter[j] = 1      ;; Start counting steps towards the spider average only if the spider step was accepted.
        print, 'SPIDER ON THE MOVE!',i,j,spiderleg[j] ;; ,parvalues[spiderleg[j],i,j]
       endif
     endif else begin            ;; go back to previous value if chain did not advance.
       parvalues[*,i,j] = parvalues[*,i-1,j]
       lnPValues[i,j]  =  lnPValues[i-1,j]
       modRAvalues[*,i,j] = modRAvalues[*,i-1,j]
       modDCvalues[*,i,j] = modDCvalues[*,i-1,j]
      endelse   ;; This completes 1 iteration for 1 chain. 
   if spidercounter[j] ne 0 then spidercounter[j] ++     ;; Keeps track of how many steps taken ever since last spider step. 
  endfor   ;; Ends the loop over chains
endfor   ;;  Ends the loop over iterations


delvar, steparray, ln_chainadvance, modRAvalues, modDCvalues,partialD       ;; free up space to create return structure

;;print, "Minimum and Maximum epochs matched ",min(matchedepochs),max(matchedepochs)
;output = create_struct('chains',parvalues, 'lnPchains',lnPValues,  'RAmodel',modRAvalues, $ 
;                      'DCmodel',modDCvalues, 'partials',partialD,'Matchepoch',matchedepochs) ;,'LastChange',lastchange,'UsedSteps',usedsteps)

output = create_struct('chains',parvalues, 'lnPchains',lnPValues,'Matchepoch',matchedepochs)

return, output

End
