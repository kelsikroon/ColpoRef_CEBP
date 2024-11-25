# ColpoRef_CEBP
Code for analysis and results in paper "Colposcopy Referral and CIN3+ Risk of Human Papillomavirus Genotyping Strategies in Cervical Cancer Screening" (Kroon et al. 2024, CEBP)

- `code/`: codes used in this analysis
- `plots/`: plots used in this study
- `results/`: output csv files
- `manuscript/`: pre-peer reviewed version of the manuscript 

## Codes
- data_prep.R
	- Data cleaning
	  - fixing dates 
		- column selection
		- wide -> long
		- inclusion/exclusion criteria
- analysis.R
  - main analysis function that runs each separate script
- referral_rates.R
  - code used to calculate referral rates for each of the 14 strategies and produce Figure 2 and Table S1
- PPV_analysis.R
  - code used to calculate positive predictive values (PPV) for each of the 14 strategies for endpoint CIN3+ (and CIN2+ in suppl.), Table 2 and Table S2
- NPV_analysis.R
  - code used to calculate negative predictive values (NPV) for each of the 14 strategies for endpoint CIN3+ (and CIN2+ in suppl.), Table 3 and Table S3
- INNR_analysis.R
  - code used for incremental analysis and to produce the efficient frontier in Figure 3

## Plots
- Figure_1.png
- Figure_3.png
- Figure_3.png
- Figure_S1.png

## Results 
- Table1_pop.csv
- Table2_PPVresults.csv
- Table3_NPVresults.csv
- TableS1_refRates.csv
- TableS2_PPVresults_CIN2plus.csv
- TableS3_NPVresults_CIN2plus.csv
