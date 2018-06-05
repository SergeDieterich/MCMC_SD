function Ecca, e, P, littlet, BigT
;+
; NAME:
;       ECCA
;       Part of MCMC_SD_V1 suite
;
; AUTHORSHIP:
;          SERGE DIETERICH, May 2018
;          Department of Terrestrial Magnetism,
;          Carnegie Institution of Washington
;
; PURPOSE:
;          Computes the eccentric anomaly by approximating the solution to Kepler's equation
;          Returns its index in MCMC_SD_V1 common arrays
;
;          ;;;;;;;;;; Kepler's equation ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;          ; (E - e sin E) * (P/2pi) = (littlet - bigT) ;;  0 <= E < 2pi    BUT ----->   ;;
;          ;                                                                            ;;
;          ; (E - e sin E) / 2pi =  ((littlet - bigT) mod P) / P ;; if not mod then E unbound ;;
;          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; CALLING SEQUENCE:
;          result = ecca(e,P,littlet,BigT)
;          all required
;          Called by astromodel.pro and partials.pro
;
; INPUTS:
;          e   =  eccentricity
;          P   =  orbital period
;      littlet =  epoch at which eccentric anomaly is to be computed
;         BigT =  epoch of periastron passage. Must be before littlet.
;
; OUTPUT:
;         The index of common arrays  Ec_array4, sin_ec4, and cos_ec4
;         that contains the corresponding value of interest.
;               Returning the index rather than the eccentric anomaly itself is more
;               efficient.
;
; PROCEDURES CALLED:
;           none.
;           MCMC_SD_V1 common block must already be established by mcmc_sd_wrapper_v1.pro
;           
; VERSION HISTORY:
;          CURRENT V1 - First publication version. May 2018.
;
; BUGS:
;          Probably solves things to unecessarily high precision at great computational cost.
;               Whether lower precision works needs to be tested.
;-
common newsolver_common_block, data_jd, data_ra, data_dc, data_raerr, data_dcerr, $
  min_jd, max_jd, pifact_ra, pifact_dc, Ec_array, sin_ec, cos_ec, n_epochs, Ec_array_lstndx, n_epochs_m1, $
  Ec_array2,sin_ec2,cos_ec2,Ec_array2_lstndx,Ec_array3,sin_ec3,cos_ec3,Ec_array3_lstndx, $
  Ec_array4,sin_ec4,cos_ec4,Ec_array4_lstndx


leftside = (Ec_array - (e * sin_ec)) / (2d *!pi)
;;
;;*******************CTIOPI*****************
;;
;;;;  To make thing faster go down to eccentric anomaly increments of 1d-6 in 3 steps of 1d-2.

rightside = ((littlet - BigT) mod P) / P      ;;This and the next 2 lines solve Kepler's equation
minresult = min(rightside - leftside,Ecc_a_index, /absolute)

central2 = ((Ecc_a_index)*1d2)+100L
minindex2 = ((central2 - 100L) > 0L)
maxindex2 = ((central2 + 100L) < Ec_array2_lstndx)
Ec_array22 =Ec_array2[minindex2:maxindex2]
sin_ec22 = sin_ec2[minindex2:maxindex2]
leftside2 = (Ec_array22 - (e * sin_ec22)) / (2d *!pi)
minresult2 = min(rightside - leftside2,Ecc_a_index2, /absolute)
Ecc_a_index2 = Ecc_a_index2 + minindex2

central3 = (((Ecc_a_index2 - 100L) > 0L)*1d2)+100L
minindex3 = ((central3 - 100L) > 0L)
maxindex3 = ((central3 + 100L) < Ec_array3_lstndx)
Ec_array33 =Ec_array3[minindex3:maxindex3]
sin_ec33 = sin_ec3[minindex3:maxindex3]
leftside3 = (Ec_array33 - (e * sin_ec33)) / (2d *!pi)
minresult3 = min(rightside - leftside3,Ecc_a_index3, /absolute)
Ecc_a_index3 = Ecc_a_index3 + minindex3

central4 = (((Ecc_a_index3 - 100L) > 0L)*1d2)+100L
minindex4 = ((central4 - 100L) > 0L)
maxindex4 = ((central4 + 100L) < Ec_array4_lstndx)
Ec_array44 =Ec_array4[minindex4:maxindex4]
sin_ec44 = sin_ec4[minindex4:maxindex4]
leftside4 = (Ec_array44 - (e * sin_ec44)) / (2d *!pi)
minresult4 = min(rightside - leftside4,Ecc_a_index4, /absolute)
Ecc_a_index4 = Ecc_a_index4 + minindex4

;;return, Ec_array4[Ecc_a_index4] * 180d/!pi

return, Ecc_a_index4

end