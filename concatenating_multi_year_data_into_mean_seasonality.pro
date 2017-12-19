; purposes: Concatenating the six-year (2011-2016) data into a mean yearly daily time series
;           for each pixel to obtain a clear and reliable seasonality

; input: layer stacked satellite image series (6 years)

; output: smoothed mean seasonality (365 bands) of the input data

; by Feng Tian; 2017 November

pro concatenating_multi_year_data_into_mean_seasonality

  compile_opt idl2
  envi, /restore_base_save_files
  envi_batch_init

  file_in=' ' ;stacked file of raw daily time series during 2011-2016

  envi_open_file, file_in, r_fid=fid, /NO_REALIZE
  envi_file_query, fid, dims=dims, nb=nb, ns=ns, nl=nl

  frequency=365 ; the number of observations in each year (exclude the last day in leap years)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; concatenating raw time series into mean seasonality
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  outfile_concat=file_dirname(file_in, /mark_directory)+file_basename(file_in)+'_concatenated'

  envi_report_init, ['step1: seasonality'], title='Wait a moment...', base=base ;add processing toolbar

  openw, unit, outfile_concat, /Get_LUN

  for j=0,frequency-1 do begin

    envi_report_stat, base, j+1, frequency
    envi_report_inc, base,frequency/(j+1)  ; handel the processing toolbar

    data_y=make_array(ns,nl, nb/frequency, /float, VALUE = 'NaN')

    for i=0, nb-frequency, frequency do data_y[*,*,i/frequency] = float(envi_get_data(fid=fid, dims=dims, pos=i+j))

    writeu,unit,median(data_y, DIMENSION=3, /EVEN)

  endfor

  free_lun,unit

  ; witre header for the ENVI file
  inherit = envi_set_inheritance(fid, dims, /SPATIAL)
  envi_setup_head, fname=outfile_concat, $
    ns=ns, nl=nl, nb=frequency, $
    data_type=4,INTERLEAVE =0, $
    offset=0,inherit=inherit, /write

  envi_report_init, base=base, /finish        ;destroy the processing toolbar

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Denoise the raw seasonaltiy and then calculate a smoothed seasonality
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  window_width=30  ; set the window size to 30 days
  
  file_in=outfile_concat

  outfile_denoise=file_in+'_denoised'  ; used for determing where we have a clear seasonaltiy

  outfile_season=file_in+'_smooth_seasonality'

  envi_open_file, file_in, r_fid=fid, /NO_REALIZE

  envi_file_query, fid, dims=dims, nb=nb, ns=ns, nl=nl

  data=make_array(ns,nl,nb, /float,value='nan')

  data_denoised=make_array(ns,nl,nb, /float,value='nan')

  for p=0,nb-1 do data[*,*,p]=envi_get_data(fid=fid,dims=dims,pos=p)
  Zeros=total(data,3, /nan)

  envi_report_init, ['Step 2: Denoise and smooth seasonality'], title='Wait a moment...', base=base ;add processing toolbar

  for i=0,nl-1 do begin

    envi_report_stat, base, i, nl
    envi_report_inc, base,float(nl)/(i+1)  ; handel the processing toolbar

    for j=0, ns-1 do begin

      ; delete the pixels with all the values are 0 (background pixels)
      if Zeros[j,i] eq 0 then data[j,i,*]='nan'

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;Denoise begin
      idx=where(finite(data[j,i,*]) eq 1, Count) ; test the number of valid data points

      if Count gt 60 then begin ; go for the analysis for pixels with enough observations (Note there are gaps in the data)

        ; extended data for moving median, etc
        ts=reform(data[j,i,*])
        ts=[ts[(nb-2*window_width):(nb-1)],ts,ts[0:(2*window_width-1)]]

        running_median=my_median(ts,window_width,1)

        anom=(ts-running_median) ; calculate the anomaly

        std=stddev(anom,/nan) ;standard deviation

        idx_noise=where(abs(anom) gt std)  ; throw out abnormal values
        if idx_noise[0] ne (-1) then ts[idx_noise]='nan'

        data[j,i,*]=ts[(2*window_width):(nb-1+2*window_width)]

        data_denoised[j,i,*]=data[j,i,*]
        ;Denoise end

        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;smooth begin
        sg_width=20
        savgolFilter = savgol(sg_width, sg_width, 0, 1)

        idx60=where(finite(data[j,i,0:59]) eq 1, Count1) ; number of valid obs in January and Febuary

        ;;;;;;; only for pixels in the northern hemesphere (no data during frozen soil conditions)
        ; and non-background (valid observation > 60) ;;;;;;;
        if Count1 lt 20 then begin

          ;running average begin
          ts=reform(data[j,i,*])
          ; extended data for moving average, etc
          ts_in=[ts[(nb-sg_width):(nb-1)],ts,ts[0:(window_width-1)]]

          out=my_average(ts_in,window_width,1.3)

          data[j,i,*]=out[window_width:(nb-1+window_width)]

          ;running average end and SG begin
          ts=reform(data[j,i,*])

          ; extended data for moving average, etc
          ts_in=[ts[(nb-sg_width):(nb-1)],ts,ts[0:(sg_width-1)]]

          SG=make_array(nb+2*sg_width, /float, VALUE = 'NaN')

          idx = fix(finite(ts_in))

          idx_end=where(ts_diff(idx,1) eq 1, Count3)

          idx_start=where(ts_diff(idx,1) eq (-1), Count4)+1

          if Count3 gt Count4 then idx_start=[0,idx_start]

          if Count3 lt Count4 then idx_end=[idx_end,nb-1]

          for tt=0,n_elements(idx_start)-1 do begin
            if idx_end[tt]-idx_start[tt] ge 60 then SG[idx_start[tt]:idx_end[tt]]=convol(ts_in[idx_start[tt]:idx_end[tt]], savgolFilter,/EDGE_TRUNCATE)
          endfor

          data[j,i,*]=SG[sg_width:(nb-1+sg_width)]
          ;SG end

        endif else begin
          ;;;;;;; For pixels in the tropics where we have valid retrivals over the course of the year ;;;;;;;

          ;running average begin
          ts=reform(data[j,i,*])

          ; extended data for moving median, etc
          ts_in=[ts[(nb-window_width):(nb-1)],ts,ts[0:(window_width-1)]]

          out=my_average(ts_in,window_width,1.7)

          data[j,i,*]=out[window_width:(nb-1+window_width)]

          ;running median end and SG begin
          ts=reform(data[j,i,*])

          ; extended data for moving median, etc
          ts_in=[ts[(nb-sg_width):(nb-1)],ts,ts[0:(sg_width-1)]]

          SG = convol(ts_in, savgolFilter,/EDGE_TRUNCATE)

          data[j,i,*]=SG[sg_width:(nb-1+sg_width)]
          ; SG end
        endelse

      endif else begin
        data[j,i,*]='nan'
      endelse
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;endif
    endfor
  endfor

  openw, unit, outfile_denoise, /Get_LUN
  writeu,unit,data_denoised
  free_lun,unit

  openw, unit, outfile_season, /Get_LUN
  writeu,unit,data
  free_lun,unit

  inherit = envi_set_inheritance(fid, dims, /SPATIAL)
  envi_setup_head, fname=outfile_denoise, $
    ns=ns, nl=nl, nb=nb, $
    data_type=4,INTERLEAVE =0, $
    offset=0,inherit=inherit, /write

  envi_setup_head, fname=outfile_season, $
    ns=ns, nl=nl, nb=nb, $
    data_type=4,INTERLEAVE =0, $
    offset=0,inherit=inherit, /write

  envi_report_init, base=base, /finish        ;destroy the processing toolbar


end



