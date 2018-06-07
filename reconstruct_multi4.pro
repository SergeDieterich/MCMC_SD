function reconstruct_multi4,nsessions
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
;          Combines the last 10 sub-runs of an MCMC_SD_V1 run into
;          a single structure. If the MCMC was run in simultaneous
;          IDL sessions then combine those sessions as well into
;          a master structure containing the run results.
;          
; CALLING SEQUENCE:          
;          combined_results = reconstruct_multi4(nsessions)
;              nsessions required
;          
; INPUT:
;          nsessions =  1, 2, 3, or 4.
;                       The number of IDL sessions that were run in parallel to speed up
;                       the computations. If nsessions > 1 then the other sessions must 
;                       have IDENTICAL calling sequences for mcmc_sd_wrapper_v1.pro
;                       and must be run from subdirectories ./session2, ./session3, and ./session4
; 
; OUTPUT:
;         A structure of the form create_struct('chains',chains, 'lnPchains',lnpchains,'Matchepoch',matchepoch)
;                        where the number of chains has been combined from the multiple IDL sessions (if more than 1)
;                        and the length of the chains has been combined from the last 10 subsessions.
;
; PROCEDURES CALLED:
;          None
;
;  VERSION HISTORY:
;          CURRENT V1 - First publication version. May 2018.
;
;  BUGS:
;         Currently limited to 4 simultaneous IDL sessions. Each session generally
;              used 100% of a processor core.
;         Combines only the last 10 subsessions.
;         Resulting structure must be saved manually into .sav file.
;-

;; first (main session)
restore,'last1_multi_results.sav'

sizearray = size(results.chains)
npars = sizearray[1]
niters = sizearray[2]
nchains = sizearray[3]
niters_m1 = niters - 1
nchains_m1 = nchains - 1

results1 = results
restore,'last2_multi_results.sav'
results2 = results
restore,'last3_multi_results.sav'
results3 = results
restore,'last4_multi_results.sav'
results4 = results
restore,'last5_multi_results.sav'
results5 = results
restore,'last6_multi_results.sav'
results6 = results
restore,'last7_multi_results.sav'
results7 = results
restore,'last8_multi_results.sav'
results8 = results
restore,'last9_multi_results.sav'
results9 = results
restore,'last10_multi_results.sav'
results10 = results

chains1 = dblarr(npars,niters * 10 ,nchains)
lnpchains1 = dblarr(niters*10,nchains)
matchepoch1 = dblarr(niters*10,nchains)

for i =0, niters_m1 do begin
   chains1[*,i,*] = results1.chains[*,i,*]
   chains1[*,i + niters,*] = results2.chains[*,i,*]
   chains1[*,i + 2*niters,*] = results3.chains[*,i,*]
   chains1[*,i + 3*niters,*] = results4.chains[*,i,*]
   chains1[*,i + 4*niters,*] = results5.chains[*,i,*]
   chains1[*,i + 5*niters,*] = results6.chains[*,i,*]
   chains1[*,i + 6*niters,*] = results7.chains[*,i,*]
   chains1[*,i + 7*niters,*] = results8.chains[*,i,*]
   chains1[*,i + 8*niters,*] = results9.chains[*,i,*]
   chains1[*,i + 9*niters,*] = results10.chains[*,i,*]

   lnpchains1[i,*] = results1.lnpchains[i,*]
   lnpchains1[i + niters,*] = results2.lnpchains[i,*]
   lnpchains1[i + 2*niters,*] = results3.lnpchains[i,*]
   lnpchains1[i + 3*niters,*] = results4.lnpchains[i,*]
   lnpchains1[i + 4*niters,*] = results5.lnpchains[i,*]
   lnpchains1[i + 5*niters,*] = results6.lnpchains[i,*]
   lnpchains1[i + 6*niters,*] = results7.lnpchains[i,*]
   lnpchains1[i + 7*niters,*] = results8.lnpchains[i,*]
   lnpchains1[i + 8*niters,*] = results9.lnpchains[i,*]
   lnpchains1[i + 9*niters,*] = results10.lnpchains[i,*]

   matchepoch1[i,*] = results1.matchepoch[i,*]
   matchepoch1[i + niters,*] = results2.matchepoch[i,*]
   matchepoch1[i + 2*niters,*] = results3.matchepoch[i,*]
   matchepoch1[i + 3*niters,*] = results4.matchepoch[i,*]
   matchepoch1[i + 4*niters,*] = results5.matchepoch[i,*]
   matchepoch1[i + 5*niters,*] = results6.matchepoch[i,*]
   matchepoch1[i + 6*niters,*] = results7.matchepoch[i,*]
   matchepoch1[i + 7*niters,*] = results8.matchepoch[i,*]
   matchepoch1[i + 8*niters,*] = results9.matchepoch[i,*]
   matchepoch1[i + 9*niters,*] = results10.matchepoch[i,*]
endfor

chains = chains1
lnPchains = lnpchains1
Matchepoch =matchepoch1
;;  combined_results = create_struct('chains',chains1, 'lnPchains',lnpchains1,'Matchepoch',matchepoch1)

if nsessions gt 1 then begin
  ;; session 2

  restore,'session2/last1_multi_results.sav'
  results1 = results
  restore,'session2/last2_multi_results.sav'
  results2 = results
  restore,'session2/last3_multi_results.sav'
  results3 = results
  restore,'session2/last4_multi_results.sav'
  results4 = results
  restore,'session2/last5_multi_results.sav'
  results5 = results
  restore,'session2/last6_multi_results.sav'
  results6 = results
  restore,'session2/last7_multi_results.sav'
  results7 = results
  restore,'session2/last8_multi_results.sav'
  results8 = results
  restore,'session2/last9_multi_results.sav'
  results9 = results
  restore,'session2/last10_multi_results.sav'
  results10 = results

  chains2 = dblarr(npars,niters * 10 ,nchains)
  lnpchains2 = dblarr(niters*10,nchains)
  matchepoch2 = dblarr(niters*10,nchains)

  for i =0, niters_m1 do begin
    chains2[*,i,*] = results1.chains[*,i,*]
    chains2[*,i + niters,*] = results2.chains[*,i,*]
    chains2[*,i + 2*niters,*] = results3.chains[*,i,*]
    chains2[*,i + 3*niters,*] = results4.chains[*,i,*]
    chains2[*,i + 4*niters,*] = results5.chains[*,i,*]
    chains2[*,i + 5*niters,*] = results6.chains[*,i,*]
    chains2[*,i + 6*niters,*] = results7.chains[*,i,*]
    chains2[*,i + 7*niters,*] = results8.chains[*,i,*]
    chains2[*,i + 8*niters,*] = results9.chains[*,i,*]
    chains2[*,i + 9*niters,*] = results10.chains[*,i,*]

    lnpchains2[i,*] = results1.lnpchains[i,*]
    lnpchains2[i + niters,*] = results2.lnpchains[i,*]
    lnpchains2[i + 2*niters,*] = results3.lnpchains[i,*]
    lnpchains2[i + 3*niters,*] = results4.lnpchains[i,*]
    lnpchains2[i + 4*niters,*] = results5.lnpchains[i,*]
    lnpchains2[i + 5*niters,*] = results6.lnpchains[i,*]
    lnpchains2[i + 6*niters,*] = results7.lnpchains[i,*]
    lnpchains2[i + 7*niters,*] = results8.lnpchains[i,*]
    lnpchains2[i + 8*niters,*] = results9.lnpchains[i,*]
    lnpchains2[i + 9*niters,*] = results10.lnpchains[i,*]

    matchepoch2[i,*] = results1.matchepoch[i,*]
    matchepoch2[i + niters,*] = results2.matchepoch[i,*]
    matchepoch2[i + 2*niters,*] = results3.matchepoch[i,*]
    matchepoch2[i + 3*niters,*] = results4.matchepoch[i,*]
    matchepoch2[i + 4*niters,*] = results5.matchepoch[i,*]
    matchepoch2[i + 5*niters,*] = results6.matchepoch[i,*]
    matchepoch2[i + 6*niters,*] = results7.matchepoch[i,*]
    matchepoch2[i + 7*niters,*] = results8.matchepoch[i,*]
    matchepoch2[i + 8*niters,*] = results9.matchepoch[i,*]
    matchepoch2[i + 9*niters,*] = results10.matchepoch[i,*]
  endfor
endif
  
if nsessions gt 2 then begin
    ;; session 3

    restore,'session3/last1_multi_results.sav'
    results1 = results
    restore,'session3/last2_multi_results.sav'
    results2 = results
    restore,'session3/last3_multi_results.sav'
    results3 = results
    restore,'session3/last4_multi_results.sav'
    results4 = results
    restore,'session3/last5_multi_results.sav'
    results5 = results
    restore,'session3/last6_multi_results.sav'
    results6 = results
    restore,'session3/last7_multi_results.sav'
    results7 = results
    restore,'session3/last8_multi_results.sav'
    results8 = results
    restore,'session3/last9_multi_results.sav'
    results9 = results
    restore,'session3/last10_multi_results.sav'
    results10 = results

    chains3 = dblarr(npars,niters * 10 ,nchains)
    lnpchains3 = dblarr(niters*10,nchains)
    matchepoch3 = dblarr(niters*10,nchains)

    for i =0, niters_m1 do begin
      chains3[*,i,*] = results1.chains[*,i,*]
      chains3[*,i + niters,*] = results2.chains[*,i,*]
      chains3[*,i + 2*niters,*] = results3.chains[*,i,*]
      chains3[*,i + 3*niters,*] = results4.chains[*,i,*]
      chains3[*,i + 4*niters,*] = results5.chains[*,i,*]
      chains3[*,i + 5*niters,*] = results6.chains[*,i,*]
      chains3[*,i + 6*niters,*] = results7.chains[*,i,*]
      chains3[*,i + 7*niters,*] = results8.chains[*,i,*]
      chains3[*,i + 8*niters,*] = results9.chains[*,i,*]
      chains3[*,i + 9*niters,*] = results10.chains[*,i,*]

      lnpchains3[i,*] = results1.lnpchains[i,*]
      lnpchains3[i + niters,*] = results2.lnpchains[i,*]
      lnpchains3[i + 2*niters,*] = results3.lnpchains[i,*]
      lnpchains3[i + 3*niters,*] = results4.lnpchains[i,*]
      lnpchains3[i + 4*niters,*] = results5.lnpchains[i,*]
      lnpchains3[i + 5*niters,*] = results6.lnpchains[i,*]
      lnpchains3[i + 6*niters,*] = results7.lnpchains[i,*]
      lnpchains3[i + 7*niters,*] = results8.lnpchains[i,*]
      lnpchains3[i + 8*niters,*] = results9.lnpchains[i,*]
      lnpchains3[i + 9*niters,*] = results10.lnpchains[i,*]

      matchepoch3[i,*] = results1.matchepoch[i,*]
      matchepoch3[i + niters,*] = results2.matchepoch[i,*]
      matchepoch3[i + 2*niters,*] = results3.matchepoch[i,*]
      matchepoch3[i + 3*niters,*] = results4.matchepoch[i,*]
      matchepoch3[i + 4*niters,*] = results5.matchepoch[i,*]
      matchepoch3[i + 5*niters,*] = results6.matchepoch[i,*]
      matchepoch3[i + 6*niters,*] = results7.matchepoch[i,*]
      matchepoch3[i + 7*niters,*] = results8.matchepoch[i,*]
      matchepoch3[i + 8*niters,*] = results9.matchepoch[i,*]
      matchepoch3[i + 9*niters,*] = results10.matchepoch[i,*]
    endfor

endif

    if nsessions eq 4 then begin
      ;; session 4

      restore,'session4/last1_multi_results.sav'
      results1 = results
      restore,'session4/last2_multi_results.sav'
      results2 = results
      restore,'session4/last3_multi_results.sav'
      results3 = results
      restore,'session4/last4_multi_results.sav'
      results4 = results
      restore,'session4/last5_multi_results.sav'
      results5 = results
      restore,'session4/last6_multi_results.sav'
      results6 = results
      restore,'session4/last7_multi_results.sav'
      results7 = results
      restore,'session4/last8_multi_results.sav'
      results8 = results
      restore,'session4/last9_multi_results.sav'
      results9 = results
      restore,'session4/last10_multi_results.sav'
      results10 = results

      chains4 = dblarr(npars,niters * 10 ,nchains)
      lnpchains4 = dblarr(niters*10,nchains)
      matchepoch4 = dblarr(niters*10,nchains)

      for i =0, niters_m1 do begin
        chains4[*,i,*] = results1.chains[*,i,*]
        chains4[*,i + niters,*] = results2.chains[*,i,*]
        chains4[*,i + 2*niters,*] = results3.chains[*,i,*]
        chains4[*,i + 3*niters,*] = results4.chains[*,i,*]
        chains4[*,i + 4*niters,*] = results5.chains[*,i,*]
        chains4[*,i + 5*niters,*] = results6.chains[*,i,*]
        chains4[*,i + 6*niters,*] = results7.chains[*,i,*]
        chains4[*,i + 7*niters,*] = results8.chains[*,i,*]
        chains4[*,i + 8*niters,*] = results9.chains[*,i,*]
        chains4[*,i + 9*niters,*] = results10.chains[*,i,*]

        lnpchains4[i,*] = results1.lnpchains[i,*]
        lnpchains4[i + niters,*] = results2.lnpchains[i,*]
        lnpchains4[i + 2*niters,*] = results3.lnpchains[i,*]
        lnpchains4[i + 3*niters,*] = results4.lnpchains[i,*]
        lnpchains4[i + 4*niters,*] = results5.lnpchains[i,*]
        lnpchains4[i + 5*niters,*] = results6.lnpchains[i,*]
        lnpchains4[i + 6*niters,*] = results7.lnpchains[i,*]
        lnpchains4[i + 7*niters,*] = results8.lnpchains[i,*]
        lnpchains4[i + 8*niters,*] = results9.lnpchains[i,*]
        lnpchains4[i + 9*niters,*] = results10.lnpchains[i,*]

        matchepoch4[i,*] = results1.matchepoch[i,*]
        matchepoch4[i + niters,*] = results2.matchepoch[i,*]
        matchepoch4[i + 2*niters,*] = results3.matchepoch[i,*]
        matchepoch4[i + 3*niters,*] = results4.matchepoch[i,*]
        matchepoch4[i + 4*niters,*] = results5.matchepoch[i,*]
        matchepoch4[i + 5*niters,*] = results6.matchepoch[i,*]
        matchepoch4[i + 6*niters,*] = results7.matchepoch[i,*]
        matchepoch4[i + 7*niters,*] = results8.matchepoch[i,*]
        matchepoch4[i + 8*niters,*] = results9.matchepoch[i,*]
        matchepoch4[i + 9*niters,*] = results10.matchepoch[i,*]
      endfor

    ;; combine chains in the case of 4 sessions

    chains = dblarr(npars,niters * 10, nchains * 4)
    lnpchains = dblarr(niters * 10, nchains * 4)
    matchepoch = dblarr(niters * 10, nchains * 4)

    for i=0, nchains_m1 do begin
      chains[*,*,i] = chains1[*,*,i]
      chains[*,*,i+nchains] = chains2[*,*,i]
      chains[*,*,i+ 2d*nchains] = chains3[*,*,i]
      chains[*,*,i+ 3d*nchains] = chains4[*,*,i]
      
      lnpchains[*,i] = lnpchains1[*,i]
      lnpchains[*,i+nchains] = lnpchains2[*,i]
      lnpchains[*,i+ 2d*nchains] = lnpchains3[*,i]
      lnpchains[*,i+ 3d*nchains] = lnpchains4[*,i]
      
      matchepoch[*,i] = matchepoch1[*,i]
      matchepoch[*,i+nchains] = matchepoch2[*,i]
      matchepoch[*,i+ 2d*nchains] = matchepoch3[*,i]
      matchepoch[*,i+ 3d*nchains] = matchepoch4[*,i]
      
    endfor
  endif

if nsessions eq 2 then begin
  chains = dblarr(npars,niters * 10, nchains * 2)
  lnpchains = dblarr(niters * 10, nchains * 2)
  matchepoch = dblarr(niters * 10, nchains * 2)

  for i =0, nchains_m1 do begin
    chains[*,*,i] = chains1[*,*,i]
    chains[*,*,i+nchains] = chains2[*,*,i]
    lnpchains[*,i] = lnpchains1[*,i]
    lnpchains[*,i+nchains] = lnpchains2[*,i]
    matchepoch[*,i] = matchepoch1[*,i]
    matchepoch[*,i+nchains] = matchepoch2[*,i]
  endfor
endif
if nsessions eq 3 then begin
  
  chains = dblarr(npars,niters * 10, nchains * 3)
  lnpchains = dblarr(niters * 10, nchains * 3)
  matchepoch = dblarr(niters * 10, nchains * 3)

  for i=0, nchains_m1 do begin
    chains[*,*,i] = chains1[*,*,i]
    chains[*,*,i+nchains] = chains2[*,*,i]
    chains[*,*,i+ 2d*nchains] = chains3[*,*,i]
    lnpchains[*,i] = lnpchains1[*,i]
    lnpchains[*,i+nchains] = lnpchains2[*,i]
    lnpchains[*,i+ 2d*nchains] = lnpchains3[*,i]
    matchepoch[*,i] = matchepoch1[*,i]
    matchepoch[*,i+nchains] = matchepoch2[*,i]
    matchepoch[*,i+ 2d*nchains] = matchepoch3[*,i]
  endfor
endif
  combined_results = create_struct('chains',chains, 'lnPchains',lnpchains,'Matchepoch',matchepoch)
  
  return, combined_results
end
