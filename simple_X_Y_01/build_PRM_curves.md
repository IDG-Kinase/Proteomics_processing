Process Simple For-Rev
================
Matthew Berginski
March 11, 2019

Read In Data
------------

``` r
process_excel_sheet <- function(excel_data_file, sheet = 1) {

  excel_data = read_excel(excel_data_file,sheet = sheet)

  split_header = str_split(names(excel_data)[1],"-")
  kinase = split_header[[1]][1]

  split_peptide_type = str_split(split_header[[1]][2]," ")

  peptide_start = split_peptide_type[[1]][1]
  curve_type = split_peptide_type[[1]][2]

  excel_data = read_excel(excel_data_file,sheet = sheet,skip=1) %>%
    mutate(kinase=kinase,curve_type=curve_type,peptide_start=peptide_start)
  return(excel_data)
}

deep_dive_data = process_excel_sheet(here('simple_X_Y/Deep Dive kinome plot data 022019-2.xlsx'),sheet = 1)
```

    ## New names:
    ## * `` -> `..2`

``` r
for (i in 2:length(excel_sheets(here('simple_X_Y/Deep Dive kinome plot data 022019-2.xlsx')))) {
  deep_dive_data = rbind(deep_dive_data,
                         process_excel_sheet(here('simple_X_Y/Deep Dive kinome plot data 022019-2.xlsx'),sheet = i))
}
```

    ## New names:
    ## * `` -> `..2`
    ## New names:
    ## * `` -> `..2`
    ## New names:
    ## * `` -> `..2`
    ## New names:
    ## * `` -> `..2`
    ## New names:
    ## * `` -> `..2`
    ## New names:
    ## * `` -> `..2`
    ## New names:
    ## * `` -> `..2`
    ## New names:
    ## * `` -> `..2`

    ## New names:
    ## * `` -> `..2`
    ## * `` -> `..3`

    ## New names:
    ## * `` -> `..2`
    ## New names:
    ## * `` -> `..2`
    ## New names:
    ## * `` -> `..2`
    ## New names:
    ## * `` -> `..2`
    ## New names:
    ## * `` -> `..2`
    ## New names:
    ## * `` -> `..2`
    ## New names:
    ## * `` -> `..2`
    ## New names:
    ## * `` -> `..2`
    ## New names:
    ## * `` -> `..2`

``` r
curve_summary = read_excel(here('simple_X_Y/Deep Dive kinome Peptide parameters 022019.xlsx'),skip=1) %>%
  clean_names() %>%
  mutate(peptide_start = str_extract(peptide,'.{3}')) %>%
  filter(!is.na(protein_name)) %>%
  mutate(correlation_plot_slope = NA,
         correlation_plot_intercept = NA,
         response_curve_slope = NA,
         response_curve_intercept = NA)

for (row_num in 1:dim(curve_summary)[1]) {
  this_left_hand_side = str_split(curve_summary$correlation_plot_equation[row_num],'=')[[1]][2]
  slope_intercept = str_split(this_left_hand_side,'x')[[1]]
  curve_summary$correlation_plot_slope[row_num] = as.numeric(slope_intercept[1])
  curve_summary$correlation_plot_intercept[row_num] = as.numeric(slope_intercept[2])
  
  this_left_hand_side = str_split(curve_summary$response_curve_equation[row_num],'=')[[1]][2]
  slope_intercept = str_split(this_left_hand_side,'x')[[1]]
  curve_summary$response_curve_slope[row_num] = as.numeric(slope_intercept[1])
  curve_summary$response_curve_intercept[row_num] = as.numeric(slope_intercept[2])
}

deep_dive_data = deep_dive_data %>% left_join(curve_summary)
```

    ## Joining, by = "peptide_start"

Plotting
========

``` r
deep_dive_reverse = deep_dive_data %>%
  filter(curve_type == "Reverse")

reverse_curves = list()

build_peptide_reverse_curve <- function(this_peptide_data) {
  this_peptide_low_vals = this_peptide_data %>%
    filter(Theor_conc < 20)
  
  full_range_plot <- ggplot(this_peptide_data,aes(x=Theor_conc,y=Meas_conc)) +
    geom_point() +
    theme_berginski() +
    labs(x="Theoretical Concentration (fmol/\U00B5L)",y="Measured Concentration (fmol/\U00B5L)") +
    geom_abline(slope = this_peptide_data$correlation_plot_slope[1],
                intercept = this_peptide_data$correlation_plot_intercept[1],
                color='blue') +
    geom_abline(slope = this_peptide_data$response_curve_slope[1],
                intercept = this_peptide_data$response_curve_intercept[1]) +
    ggtitle(paste0(this_peptide_data$kinase," - ", this_peptide_data$peptide)) +
    geom_text(aes(x = max(this_peptide_data$Theor_conc),
                  y = 1,
                  label = paste0('y=',this_peptide_data$correlation_plot_slope[1],'x',this_peptide_data$correlation_plot_intercept[1])),
              hjust="inward",vjust="inward",size=8) +
    NULL
  
  zoom_plot <- ggplot(this_peptide_low_vals,aes(x=Theor_conc,y=Meas_conc)) +
    geom_vline(xintercept = this_peptide_low_vals$lod[1], color='red') +
    geom_text(aes(x=this_peptide_low_vals$lod[1],
                  y=max(this_peptide_low_vals$Meas_conc),
                  label = paste0('LOD=',this_peptide_low_vals$lod[1])),
              color='red',hjust="inward",vjust="inward",size=6) +
    geom_point() +
    theme_berginski() +
    labs(x="",y="") +
    geom_abline(slope = this_peptide_low_vals$correlation_plot_slope[1],
                intercept = this_peptide_low_vals$correlation_plot_intercept[1],
                color='blue') +
    geom_abline(slope = this_peptide_low_vals$response_curve_slope[1],
                intercept = this_peptide_low_vals$response_curve_intercept[1])
  
  combined_plot <- ggdraw() +
    draw_plot(full_range_plot) +
    draw_plot(zoom_plot,x = 0.075,y = 0.5,width = 0.4,height = 0.4)
  return(combined_plot)
}

dir.create(here('simple_X_Y','figures'),showWarnings = F)
for (this_peptide in unique(deep_dive_reverse$peptide)) {
  this_peptide_data = deep_dive_reverse %>% filter(peptide == this_peptide)
  
  this_kinase = this_peptide_data$kinase[1]

  reverse_curves[[this_kinase]][[this_peptide]] <- build_peptide_reverse_curve(this_peptide_data)
  reverse_curves[[this_kinase]][[this_peptide]]
  ggsave(here('simple_X_Y','figures',paste0(this_kinase,'-',this_peptide,'-Reverse.svg')), 
         reverse_curves[[this_kinase]][[this_peptide]])
}
```

    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image

``` r
# deep_dive_reverse = deep_dive_data %>%
#   filter(curve_type == "Reverse")
```
