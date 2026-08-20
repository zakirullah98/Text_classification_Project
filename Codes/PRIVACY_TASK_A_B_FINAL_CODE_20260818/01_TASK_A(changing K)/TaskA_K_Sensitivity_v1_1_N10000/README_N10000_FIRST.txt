TASK A -- K SENSITIVITY -- FORMAL N=10,000 PER K

Run order:
1. The K=2 backward-compatibility test and Quick/Pilot were already passed in v1.0.
2. Run: 05_RUN_FORMAL_KGRID_N10000.R
3. The formal run attempts exactly N=10,000 repetitions for EACH K = 2,3,5,10,20,50,100.
4. Results are saved under: results_formal_Kgrid_N10000/
5. After ALL K values finish, run: 06_SUMMARIZE_TASKA_N10000.R
6. Teacher-facing outputs will appear under:
   results_formal_Kgrid_N10000/teacher_summary/

Main teacher outputs:
- TABLE1_TEACHER_MAIN_R2_R3_XY.csv
- TABLE2_R2_R3_XY_UNCERTAINTY.csv
- FIG1_MEMBER_ACCURACY_R2_R3_VS_RANDOM.png
- FIG2_MEMBER_MINUS_RANDOM_R2_R3.png
- FIG3_MEMBER_MINUS_FRESH_R2_R3.png
- FIG4_R2XY_MEMBER_FRESH_RANDOM.png
- FIG5_EXTERNAL_SUPPORT_SIZE.png
- FIG6_HESSIAN_REGULARIZATION_RATE.png

The run uses batch checkpoints. If RStudio stops, rerun the formal script; completed batches are skipped.
