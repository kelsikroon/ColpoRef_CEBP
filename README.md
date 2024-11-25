# ColpoRef_CEBP
Code for analysis and results in paper "Colposcopy Referral and CIN3+ Risk of Human Papillomavirus Genotyping Strategies in Cervical Cancer Screening" (Kroon et al. 2024, CEBP)

- `code/`: codes used in this analysis
- `plots/`: plots used in this study
- `results/`: output csv files
- `manuscript/`: pre-peer reviewed version of the manuscript 

## Codes
- `data_prep.R`
	- Data cleaning
	  - fixing dates 
		- column selection
		- wide -> long
		- inclusion/exclusion criteria
- `analysis.R`
  - main analysis function that runs each separate script
- `referral_rates.R`
  - code used to calculate referral rates for each of the 14 strategies and produce Figure 2 and Table S1
- `PPV_analysis.R`
  - code used to calculate positive predictive values (PPV) for each of the 14 strategies for endpoint CIN3+ (and CIN2+ in suppl.), Table 2 and Table S2
- `NPV_analysis.R`
  - code used to calculate negative predictive values (NPV) for each of the 14 strategies for endpoint CIN3+ (and CIN2+ in suppl.), Table 3 and Table S3
- `INNR_analysis.R`
  - code used for incremental analysis and to produce the efficient frontier in Figure 3

## Plots
- `Figure_1.png`: Flowchart of the POBASCAM women included in this analysis. “7 types” refers to the seven hrHPV types from the nonavalent HPV vaccine (HPV16/18/31/33/45/52/58).
- `Figure_2.png`: Cumulative referral rate of the 14 genotyping strategies. Referral rates are shown over three screening time points (baseline, first repeat test, and second round) for (left) more aggressive strategies than the reference and (right) more conservative strategies than the reference.
- `Figure_3.png`: Efficient frontier of the fourteen referral strategies. The number of first round referrals are plotted against the number of first round CIN3+ cases detected for each of the 14 strategies. The efficient frontier is shown in red, with labels specifying the efficient strategies and their respective marginal PPV (mPPV).
- `Figure_S1.png`: Example referral rate calculation 

## Results 
- `Table1_pop.csv`: Population characteristics of the intervention and control arm of the POBASCAM study.
- `Table2_PPVresults.csv`: PPV for CIN3+ of 14 different strategies in hrHPV-positive women stratified by baseline cytology.
- `Table3_NPVresults.csv`: NPV for CIN3+ of 14 different strategies in hrHPV-positive women.
- `TableS1_refRates.csv`: Colposcopy referral rates of the 14 different genotyping strategies
- `TableS2_PPVresults_CIN2plus.csv`: PPV for endpoint CIN2+ of 14 different strategies in hrHPV-positive women 
- `TableS3_NPVresults_CIN2plus.csv`: NPV for endpoint CIN2+ of 14 different strategies in hrHPV-positive women 
