; calculate the time lag between L-VOD and LAI seaosnality

; inputs: L-VOD and LAI seasonality resulted from 'concatenating_multi_year_data_into_mean_seasonality.pro'
;         mask of the pixels which meet the cretiria for time lag analysis
        
; outputs: time lag (number of days)
;          highest Spearman correlation coefficient and p-value

; by Feng Tian, 2017 November 

pro time_lag_analysis

  compile_opt idl2
  envi, /restore_base_save_files
  envi_batch_init

  file_validmask=' ' ; mask of the pxiels included in the analysis (value=1: included; value=0:not included)
  file_vi=' ' ; LAI seasonality data path (365 bands)
  file_smosvod=' ' ; L-VOD seasonality data path (365 bands)
  outfile=' '  ; output data path

  envi_open_file, file_vi, r_fid=fid_vi, /NO_REALIZE
  envi_open_file, file_smosvod, r_fid=fid_smosvod, /NO_REALIZE
  envi_open_file, file_validmask, r_fid=fid_mask, /NO_REALIZE
  envi_file_query, fid_vi, dims=dims, nb=nb, ns=ns, nl=nl

  data_vi=make_array(ns,nl,nb,/float,value='nan')
  data_smosvod=make_array(ns,nl,nb,/float,value='nan')

  data_mask = envi_get_data(fid=fid_mask, dims=dims, pos=0)
  idx=where(data_mask eq 1)
  
  for pp=0,nb-1 do begin  ; VOD and LAI will have the same data gaps after this process
    
    tmp=envi_get_data(fid=fid_vi, dims=dims, pos=pp)
    tmp1=envi_get_data(fid=fid_smosvod, dims=dims, pos=pp)
    idx_nan=where(finite(tmp+tmp1) eq 0)
    tmp[idx_nan]='nan'
    tmp1[idx_nan]='nan'
    data_vi[*,*,pp]=tmp
    data_smosvod[*,*,pp]=tmp1
    
  endfor
  
  ;extent the L-VOD time seires for moving
  tmp=make_array(ns,nl,nb+30+183,/float,value='nan') 
  tmp[*,*,[0:29]]=data_smosvod[*,*,[(nb-30):(nb-1)]]
  tmp[*,*,[30:(nb+30-1)]]=data_smosvod
  tmp[*,*,[(nb+30):(nb+30+183-1)]]=data_smosvod[*,*,[0:182]]
  data_smosvod=tmp
  tmp=fltarr(1)  ;reduce the memory
  
  r=make_array(ns,nl,/float,value='nan')
  r_p=make_array(ns,nl,/float,value='nan')
  lag=make_array(ns,nl,/float,value='nan')

  envi_report_init, ['Processing...'], title='Wait a moment...', base=base ;add processing toolbar

  nb1=nb+30+183 ; extended number of bands
  l=1 ;for the processing bar
  
  foreach ele, idx do begin  ;use foreach loop to perform the analysis on selected pixrels (based on the mask)
    
    envi_report_stat, base, l, n_elements(idx)
    envi_report_inc, base,float(n_elements(idx))/(l)  ; handel the processing toolbar
    s=ele mod ns ;column number
    t=ele/ns     ; row number

    vi=reform(data_vi[s,t,*])
    smosvod=reform(data_smosvod[s,t,*])
    results=fltarr(2,nb1-nb+1)-2   ;forwards one month 30 days
    
    idx_vi=where(finite(vi) eq 1, Countvi)

    if Countvi eq nb then begin ;pixels with complete seasonality over the course of the year
      
      ; moving time lag based on the maximum correlation coefficient
      for p=0,nb1-nb do results[*,p] = r_correlate(vi,smosvod[p:p+nb-1])
      r[s,t]=max(results[0,*],Max_subscript)
      r_p[s,t]=results[1,Max_subscript]
      lag[s,t]=Max_subscript-30
      l=l+1
      
    endif else begin ; pixels with data gaps during frozen soil conditions
      
      ; moving time lag based on the maximum correlation coefficient
      maxvi=max(vi,Max_Subscript,/nan)
      
      for p=0,nb1-nb do begin
        
        idx_vod=where(finite(smosvod[p:(nb-1+p)]) eq 1)
        match, idx_vi, idx_vod, suba, subb, COUNT = countcor
        
        if countcor gt 60 and idx_vi[suba[-1]]-Max_Subscript gt 30 then begin
          idx_cor=idx_vi[suba]
          results[*,p] = r_correlate(vi[idx_cor],smosvod[idx_cor+p])
        endif else continue
        
      endfor
      
      r[s,t]=max(results[0,*],Max_subscript)
      r_p[s,t]=results[1,Max_subscript]
      lag[s,t]=Max_subscript-30
      l=l+1
      
    endelse
    
  endforeach

  outfile=outfile+'_r'    ; highest correlation coefficent
  outfile1=outfile+'_p'   ; p-value
  outfile2=outfile+'_lag' ; time lag

  openw, unit, outfile, /Get_LUN
  writeu,unit,r
  free_lun,unit

  openw, unit, outfile1, /Get_LUN
  writeu,unit,r_p
  free_lun,unit

  openw, unit, outfile2, /Get_LUN
  writeu,unit,lag
  free_lun,unit

  inherit = envi_set_inheritance(fid_vi, dims, /SPATIAL)
  envi_setup_head, fname=outfile, $
    ns=ns, nl=nl, nb=1, $
    data_type=4,INTERLEAVE =0, $
    offset=0,inherit=inherit, /write

  envi_setup_head, fname=outfile1, $
    ns=ns, nl=nl, nb=1, $
    data_type=4,INTERLEAVE =0, $
    offset=0,inherit=inherit, /write

  envi_setup_head, fname=outfile2, $
    ns=ns, nl=nl, nb=1, $
    data_type=4,INTERLEAVE =0, $
    offset=0,inherit=inherit, /write

  envi_report_init, base=base, /finish        ;destroy the processing toolbar

end

