Frozen Protocol remains v1.1. This is an IMPLEMENTATION patch only.

What changed:
1) If second-stage external support has exactly ONE variable, glmnet previously failed
   with: "x should be a matrix with 2 or more columns".
2) The code now appends an all-zero dummy column ONLY for glmnet/cv.glmnet and
   explicitly excludes that dummy predictor. The refined debiasing step still uses the
   original one-variable design. Therefore the statistical target is unchanged.
3) Two diagnostic columns are recorded:
   singleton_glmnet_patch_site1
   singleton_glmnet_patch_site2
4) Empty external support is NOT repaired. It remains a genuine finite-sample failure
   and is recorded as before.
5) Added 03B_RUN_PILOT_R500.R.

Recommended next run:
source("03B_RUN_PILOT_R500.R")
