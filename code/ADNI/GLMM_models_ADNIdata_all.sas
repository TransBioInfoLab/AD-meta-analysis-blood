
/**************************************************************
* Project: Cross-tissue meta-analysis of blood and brain 
*  epigenome-wide association studies in Alzheimer’s disease
*
* Author: Lily Wang, PhD
* Date: 7-2021
***************************************************************/


%let dir = C:\Users\lxw391\TBL Dropbox\AD-meta-analysis-blood-samples\LW\ADNI_GLMM\AD_CN\DATA_for_SAS_models; 

%let res = C:\Users\lxw391\TBL Dropbox\AD-meta-analysis-blood-samples\LW\ADNI_GLMM\AD_CN\results; 

%macro p; 

%do i = 1 %to 729; 

PROC IMPORT OUT=ad_cn DATAFILE= "&dir\ADNI_1000_cpgs_&i..csv" REPLACE; GETNAMES=yes; RUN;


proc sort data = ad_cn; by cpg rid PlateNumber; run;

ods exclude all; ods noresults;

proc glimmix data = ad_cn; 
by cpg; 
	class RID PlateNumber PTGENDER; 

	model DIAGNOSIS (descending)  = beta age_at_visit  PTGENDER  PlateNumber B NK CD4T CD8T Mono Neutro / dist = binary s; 

	random int /subject = RID; 

	*ods output Tests3 = res_beta; 

	ods output ParameterEstimates = res_beta; 
run; 
ods exclude none; ods results;

data res_beta; set res_beta; if Effect = "beta"; format Probt e10.; drop PTGENDER PlateNumber; run;

proc export data=res_beta dbms=csv outfile="&res\glmm_using_beta_&i..csv" replace; run;

%end; 

%mend; 

%p; 
