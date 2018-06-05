function partials, pval, p_model_ra, p_model_dc
;+
; NAME:
;       PARTIALS
;       Part of MCMC_SD_V1 suite
;
; AUTHORSHIP:
;          SERGE DIETERICH, May 2018
;          Department of Terrestrial Magnetism,
;          Carnegie Institution of Washington
;
; PURPOSE:
;         Approximates the partial derivatives of the astrometric model
;         with respect to each parameter. Uses d(total displacement)/d(parameter)
;         as a rough proxy for d(ln(P))/d(parameter). See equations A2 and A3
;         of Dieterich et al. 2018.
;         Used for step scaling by mcmcsd_spidersolver_v1.pro
;
; CALING SEQUENCE:         
;         result = partials(pval,p_model_ra,p_model_dc)
;         All parameters required
;         Called by mcmc_spidersolver_v1.pro
;
; INPUTS:        
;         pval        = vector with values for each astrometric parameter for the previous step,
;                       following the standard index convention 
;         p_model_ra  = vector with previous step's computed RA for each epoch
;         p_model_dc  = vector with previous step's computed DEC for each epoch
;          
;  OUTPUT:
;         sum_partials = 10 element double precision vector containing the derivative of total displacement
;                        (RA + DEC) with respect to each astrometric parameter.
;  PROCEDURES CALLED:
;        ecca.pro
;        MCMC_SD_V1 common block must already be established by mcmc_sd_wrapper_v1.pro
;
;  VERSION HISTORY:
;          CURRENT V1 - First publication version. May 2018.
;
;  BUGS:
;          Not a true approximation of d(ln(P))/d(parameter)
;          Relies on previous step even if a jump makes that step far
;                away from current location in probability space.
;          Problems with the bound nature of eccentricity (parameter 5),
;          which for now must be divided by a large number by mcmcsd_spidersolver_v1.pro       
;-

common newsolver_common_block, data_jd, data_ra, data_dc, data_raerr, data_dcerr, $
  min_jd, max_jd, pifact_ra, pifact_dc, Ec_array, sin_ec, cos_ec, n_epochs, Ec_array_lstndx, n_epochs_m1, $
  Ec_array2,sin_ec2,cos_ec2,Ec_array2_lstndx,Ec_array3,sin_ec3,cos_ec3,Ec_array3_lstndx, $
  Ec_array4,sin_ec4,cos_ec4,Ec_array4_lstndx

; Define Thiele-Innes elements

Athi = pval[3] * (cos(pval[7])*cos(pval[8]) - sin(pval[7])*sin(pval[8])*cos(pval[9]))
Bthi = pval[3] * (sin(pval[7])*cos(pval[8]) + cos(pval[7])*sin(pval[8])*cos(pval[9]))
Fthi = pval[3] * (-cos(pval[7])*sin(pval[8]) - sin(pval[7])*cos(pval[8])*cos(pval[9]))
Gthi = pval[3] * (-sin(pval[7])*sin(pval[8]) + cos(pval[7])*cos(pval[8])*cos(pval[9]))

; Define partials for Thiele-Innes elements, notebook p. 113

dAda  = cos(pval[7])*cos(pval[8]) - sin(pval[7])*sin(pval[8])*cos(pval[9])
dAdOM = pval[3] * (-sin(pval[7])*cos(pval[8]) - cos(pval[7])*sin(pval[8])*cos(pval[9]))
dAdo  = pval[3] * (-cos(pval[7])*sin(pval[8]) - sin(pval[7])*cos(pval[8])*cos(pval[9]))
dAdi  = pval[3] * (sin(pval[7])*sin(pval[8])*sin(pval[9]))

dBda  = sin(pval[7])*cos(pval[8]) + cos(pval[7])*sin(pval[8])*cos(pval[9])
dBdOM = pval[3] * (cos(pval[7])*cos(pval[8]) - sin(pval[7])*sin(pval[8])*cos(pval[9]))
dBdo  = pval[3] * (-sin(pval[7])*sin(pval[8]) + cos(pval[7])*cos(pval[8])*cos(pval[9]))
dBdi  = -pval[3] * (cos(pval[7])*sin(pval[8])*sin(pval[9]))

dFda  = -(cos(pval[7])*sin(pval[8]) + sin(pval[7])*cos(pval[8])*cos(pval[9]))
dFdOM = pval[3] * (sin(pval[7])*sin(pval[8]) - cos(pval[7])*cos(pval[8])*cos(pval[9]))
dFdo  = pval[3] * (-cos(pval[7])*cos(pval[8]) + sin(pval[7])*sin(pval[8])*cos(pval[9]))
dFdi  = pval[3] * (sin(pval[7])*cos(pval[8])*sin(pval[9]))

dGda  = (-sin(pval[7])*sin(pval[8]) + cos(pval[7])*cos(pval[8])*cos(pval[9]))
dGdOM = -pval[3] * (cos(pval[7])*sin(pval[8]) + sin(pval[7])*cos(pval[8])*cos(pval[9]))
dGdo  = -pval[3] * (sin(pval[7])*cos(pval[8]) + cos(pval[7])*sin(pval[8])*cos(pval[9]))
dGdi  = -pval[3] * (cos(pval[7])*cos(pval[8])*sin(pval[9]))

;; Now compute the partial derivatives of the models, which must be done for each epoch. 

partials_ra = dblarr(10,n_epochs)    ; 10 parameters, n epochs, each parameter is one column (1st index)
partials_dc = dblarr(10,n_epochs)

for epochs = 0, n_epochs_m1 do begin
   Ec_index = Ecca(pval[5],pval[4],data_jd[epochs],pval[6])    ;; Eccentric anomaly
   
   dEEdP = (Ec_array4[Ecca(pval[5],(pval[4]+100d),data_jd[epochs],pval[6])] - Ec_array4[Ecca(pval[5],(pval[4]-100d),data_jd[epochs],pval[6])])/200d      ;; Numerical aprox of partials of eccntric anomaly
   dEEdT = (Ec_array4[Ecca(pval[5],pval[4],data_jd[epochs],(pval[6]+100d))] - Ec_array4[Ecca(pval[5],pval[4],data_jd[epochs],(pval[6]-100d))])/200d

   dEEde = (Ec_array4[Ecca((pval[5]+0.01d < 0.9d ),pval[4],data_jd[epochs],pval[6])] - Ec_array4[Ecca((pval[5]-0.01d > .05d),pval[4],data_jd[epochs],pval[6])]) $
            /((pval[5]+0.1d < 0.9d ) - (pval[5]-0.1d > .05d))

   X = cos_ec4[Ec_index] - pval[5]
   Y = (1d - pval[5]^2d)^.5d * sin_ec4[Ec_index]

   dXdP = -sin_ec4[Ec_index]*dEEdP
   dXdT = -sin_ec4[Ec_index]*dEEdT
   dXde = -sin_ec4[Ec_index]*dEEde -1d
   
   dYdEE = -(1d - pval[5]^2d)^.5d * cos_ec4[Ec_index]
   dYdP  = dYdEE * dEEdP
   dYdT  = dYdEE * dEEdT
   dYde  = -pval[5]*((1d - pval[5]^2d)^(-.5d)*sin_ec4[Ec_index] - ((1d - pval[5]^2d)^.5d *cos_ec4[Ec_index] * dEEde))

;; Now enter the partial derivatives. Notebook page 117.
   partials_ra[0,epochs] = pifact_ra[epochs]                         ;; parallax     
   partials_ra[1,epochs] = data_jd[epochs] - min_jd                  ;; RA proper motion
   partials_ra[2,epochs] = 0d                                        ;; DC proper motion
   partials_ra[3,epochs] = dBda*X + dGda*Y                           ;; CTIOPI semi-major axis
   partials_ra[4,epochs] = Bthi*dXdP + Gthi*dYdP                     ;; period
   partials_ra[5,epochs] = Bthi*dXde + Gthi*dYde                     ;; eccentricity
   partials_ra[6,epochs] = Bthi*dXdT + Gthi*dYdT                     ;; time of periastron passage
   partials_ra[7,epochs] = dBdOM*X + dGdOM*Y                         ;; OMEGA, PA of line of nodes
   partials_ra[8,epochs] = dBdo*X + dGdo*Y                           ;; omega, PA of periastron
   partials_ra[9,epochs] = dBdi*X + dGdi*Y                           ;; inclination
   
   partials_dc[0,epochs] = pifact_dc[epochs]                         ;; parallax
   partials_dc[1,epochs] = 0d
   partials_dc[2,epochs] = data_jd[epochs] - min_jd 
   partials_dc[3,epochs] = dAda*X + dFda*Y 
   partials_dc[4,epochs] = Athi*dXdP + Fthi*dYdP   
   partials_dc[5,epochs] = Athi*dXde + Fthi*dYde   
   partials_dc[6,epochs] = Athi*dXdT + Fthi*dYdT  
   partials_dc[7,epochs] = dAdOM*X + dFdOM*Y   
   partials_dc[8,epochs] = dAdo*X + dFdo*Y    
   partials_dc[9,epochs] = dAdi*X + dFdi*Y   
     
endfor  
  
;; Then the CAPSCam stuff comes here, but let's check with just CTIOPI first.

;; Getting the final derivatives. Bottom of page 118

;; 2 ways to measure the influence of changing a parameters. First way is to measure how much the RA and DEC
;; will change as a result of changing that parameter. Second is to carry the derivaties through the chi square
;; and compute d(ln(P))/d(parameter). Second way is in theory best, since d(ln(P))/d(parameter) is what the MCMC
;; ultimately cares about, but it can cause extreme slopes. First way may be sufficient as a way of scaling how much
;; model predictions move as parameters change.

;; FIRST WAY
sum_partials_ra = dblarr(10)
sum_partials_dc = dblarr(10)
for parameters = 0, 9 do begin
  sum_partials_ra[parameters] = total(abs(partials_ra[parameters,*]))
  sum_partials_dc[parameters] = total(abs(partials_dc[parameters,*]))
endfor

sum_partials = sum_partials_ra + sum_partials_dc
return, sum_partials 
 
;; SECOND WAY
;; Chi_ra_partials = dblarr(10)
;; Chi_dc_partials = dblarr(10)   
;;  
;; for parameters = 0, 9 do begin
;;  Chi_ra_partials[parameters] = total((-(data_ra - p_model_ra) * partials_ra[parameters,*])/(2d* data_raerr^2d))
;;  Chi_dc_partials[parameters] = total((-(data_dc - p_model_dc) * partials_dc[parameters,*])/(2d* data_dcerr^2d))
;; endfor

;;dlnPdPar = dblarr(10)

;; chi2_ra = total((data_ra - p_model_ra)^2d / (2d * data_raerr)^2d)
;; chi2_dc = total((data_dc - p_model_dc)^2d / (2d * data_dcerr)^2d)

;; ln_likelihood = -chi2_ra/2d - chi2_dc/2d

;; ;;Prob = exp(-chi2_ra/2d - chi2_dc/2d)

;; ;dProbdPar = -.5d*(Chi_ra_partials + Chi_dc_partials)*exp(-chi2_ra/2d - chi2_dc/2d)

;; dlnPdPar = -0.5d*(Chi_ra_partials +  Chi_dc_partials)

;; return, dlnPdPar
  
 end