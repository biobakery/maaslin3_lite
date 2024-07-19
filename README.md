## Usage: 

#### STEP 1: Generate lefse like results using maaslin3
```
Rscript lefse_testing.R --input hmp_aerobiosis_small.txt
```
#### Output 
```
lefse_style_results_abundance.res
lefse_style_results_prevalence.res
```

#### STEP2:  Generating lefse like plots 
```
lefse_plot_res.py lefse_style_results_prevalence.res hmp_prevalance_aerobiosis_small.png
```
![hmp_prevalance_aerobiosis_small](./output/hmp_prevalance_aerobiosis_small.png)

```
lefse_plot_res.py lefse_style_results_abundance.res hmp_abundance_aerobiosis_small.png
```
![hmp_abundance_aerobiosis_small](./output/hmp_abundance_aerobiosis_small.png)

#### STEP3:  Generating lefse like cladogram 
```
lefse_plot_cladogram.py lefse_style_results_prevalence.res hmp_prevalence_aerobiosis_small.cladogram.png --format png
```
![hmp_prevalence_aerobiosis_small.cladogram](./output/hmp_prevalence_aerobiosis_small.cladogram.png)

```
lefse_plot_cladogram.py lefse_style_results_abundance.res hmp_abundance_aerobiosis_small.cladogram.png --format png
```
![hmp_abundance_aerobiosis_small.cladogram](./output/hmp_abundance_aerobiosis_small.cladogram.png)
