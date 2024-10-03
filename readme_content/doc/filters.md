# Filters selection

<p align="justify"> The thresholds and everything you want to include or exclude for filtering are present in the conf.ini file. Below is what our configuration includes and excludes.

| Filters | Type | Include
|-------------------------|----------------| :---:|
|CLIN_SIG <br> | <p align="justify">"risk_factor","pathogenic","likely_pathogenic","drug_response"| Yes
|CONSEQUENCES <br> | <p align="justify">"splice_region_variant","frameshift_variant","splice_donor_variant","stop_gained","splice_acceptor_variant","start_lost","inframe_insertion","inframe_deletion"| Yes
|ONCOKB_FILTER <br> | <p align="justify">"Likely Oncogenic","Oncogenic"| Yes
|BENIGN <br> | <p align="justify">benign\likely_benign| No
|POLYPHEN <br> | <p align="justify">"benign" | No
|IMPACT <br> |<p align="justify">   "LOW"| No
|SIFT <br> | <p align="justify">"tolerated"| No
 
 ⚠️ *There are other thresholds that involve a range of inclusion and exclusion, and they are:t_VAF_min=0.02;t_VAF_min_novel=0.05; t_VAF_max=0.98;AF=<0.0003*
