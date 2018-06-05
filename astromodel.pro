function astromodel, par_a
;+
; NAME:
;       ASTROMODEL
;       Part of MCMC_SD_V1 suite
;
; AUTHORSHIP:
;          SERGE DIETERICH, May 2018
;          Department of Terrestrial Magnetism,
;          Carnegie Institution of Washington
;
; PURPOSE:
;          The 10 parameter astrometric model. Calculates the probability that
;          a given set of values matches the data given the uncertainty on the data.
;          Good for a single data set where all dispacements are relative to the first
;          epoch of observation. See instructions.txt and/or mcmc_sd_wrapper_v1.pro for
;          data format and general considerations.
;
; CALLING SEQUENCE:
;           result = astromodel(par_a)
;           Called by mcmcsd_spidersolver_v1.pro
;           
; INPUTS:
;        par_a
;              parameter values to be tested against the data.
;              10 element double precision vector. 
;              Indices are as follows (units of mas, radians, days):
;          0: parallax
;;         1: ra proper motion
;;         2: dec proper motion
;;         3: semi-major axis
;;         4: period
;;         5: eccentricity
;;         6: time of periastron passage (JD before first epoch)
;;         7: big omega (position angle of the line of nodes)
;;         8: little omega (position angle of periastron)
;;         9: inclination (i)
;
;  OUTPUT:
;        outputs a structure of the form
;               create_struct('p_model_ra', model_ra, 'p_model_dc', model_dc, 'ln_prob', ln_probability, $
;                             'goodmatches',goodmatches)
;              p_model_ra:  computed RA displacement for each epoch
;              p_model_dc:  computed DEC displacement for each epoch
;              ln_prob:     natural log of the probability that these parameters fit the data
;              goodmatches: number of epochs within error bars for given epoch
;
;  PROCEDURES CALLED:
;           ecca.pro
;        MCMC_SD_V1 common block must already be established by mcmc_sd_wrapper_v1.pro
;
;  VERSION HISTORY:
;          CURRENT V1 - First publication version. May 2018.
;
;  BUGS:   
;          Epoch of periastron passage, parameter 6, must be before first epoch.
;           
;-

common newsolver_common_block, data_jd, data_ra, data_dc, data_raerr, data_dcerr, $
  min_jd, max_jd, pifact_ra, pifact_dc, Ec_array, sin_ec, cos_ec,  n_epochs, Ec_array_lstndx, n_epochs_m1, $
  Ec_array2,sin_ec2,cos_ec2,Ec_array2_lstndx,Ec_array3,sin_ec3,cos_ec3,Ec_array3_lstndx, $
  Ec_array4,sin_ec4,cos_ec4,Ec_array4_lstndx


;; Now compute the fit model -- Hlimit2.0 notebook, page 108
;;  Define Thiele-Innes constants

Athi = par_a[3] * ( cos(par_a[7]) * cos(par_a[8]) - sin(par_a[7]) * sin(par_a[8]) * cos(par_a[9]))
Bthi = par_a[3] * ( sin(par_a[7]) * cos(par_a[8]) + cos(par_a[7]) * sin(par_a[8]) * cos(par_a[9]))
Fthi = par_a[3] * (-cos(par_a[7]) * sin(par_a[8]) - sin(par_a[7]) * cos(par_a[8]) * cos(par_a[9]))
Gthi = par_a[3] * (-sin(par_a[7]) * sin(par_a[8]) + cos(par_a[7]) * cos(par_a[8]) * cos(par_a[9]))

model_ra = dblarr(n_epochs)
model_dc = dblarr(n_epochs)

;model_ra = dblarr(n_epochs + caps_epochs, /nozero)
;model_dc = dblarr(n_epochs + caps_epochs, /nozero)


;;;;;;;;;;; Kepler's equation ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  Solved now by function Ecca
;; (E - e sin E) * (P/2pi) = (littlet - bigT) ;;  0 <= E < 2pi    BUT ----->   ;;
;;                                                                            ;;
;; (E - e sin E) / 2pi =  ((littlet - bigT) mod P) / P ;; if not mod then E unbound ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; Solve first epoch separately in order to compute barycentric correction

Ecc_a_index4 = Ecca(par_a[5],par_a[4],data_jd[0],par_a[6])

XXfirst = cos_ec4[Ecc_a_index4] - par_a[5]                       ;; Define XX and YY eliptical rectangular coordinates
YYfirst = (1d - par_a[5]^2d)^(.5d) * sin_ec4[Ecc_a_index4]       ;; Hilditch, 50

;; Barycentric correction for first epoch
BC_ra = pifact_ra[0] * par_a[0] + (Bthi * XXfirst + Gthi * YYfirst) ; + par_a[10]
BC_dc = pifact_dc[0] * par_a[0] + (Athi * XXfirst + Fthi * YYfirst) ; + par_a[11]

for epochs = 0, n_epochs_m1 do begin   ;; Solve Kepler's eqn and compute model for all other epochs

  Ecc_a_index4 = Ecca(par_a[5],par_a[4],data_jd[epochs],par_a[6])
  XX = cos_ec4[Ecc_a_index4] - par_a[5]                       ;; Define XX and YY eliptical rectangular coordinates   
  YY = (1d - par_a[5]^2d)^(.5d) * sin_ec4[Ecc_a_index4]       ;; Hilditch, 50
  
    ;;Solutions computed here, notebook p. 108.
  model_ra[epochs] = par_a[1] * (data_jd[epochs] - data_jd[0]) + pifact_ra[epochs] * par_a[0] + (Bthi * XX + Gthi * YY) - BC_ra
  model_dc[epochs] = par_a[2] * (data_jd[epochs] - data_jd[0]) + pifact_dc[epochs] * par_a[0] + (Athi * XX + Fthi * YY) - BC_dc
endfor

chi2_ra = total((data_ra - model_ra)^2d / (2d * data_raerr^2d))
chi2_dc = total((data_dc - model_dc)^2d / (2d * data_dcerr^2d))

ln_likelihood = -chi2_ra/2d - chi2_dc/2d
 
ln_probability = ln_likelihood  ;+ total(ln_prior)


goodmatchesarray = where(abs(data_ra - model_ra) le data_raerr and abs(data_dc - model_dc) le data_dcerr)

goodmatches = 0
if goodmatchesarray[0] ne -1 then goodmatches = n_elements(goodmatchesarray)

returnstruct = create_struct('p_model_ra', model_ra, 'p_model_dc', model_dc, 'ln_prob', ln_probability, $ 
                              'goodmatches',goodmatches)
;return, -ln_probability
return, returnstruct

End
