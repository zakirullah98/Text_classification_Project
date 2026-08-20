TASK B -- POPULATION ROBUSTNESS OF FED-FDR SOURCE-SITE ATTACK
Version: v1.0

THE CORE RULE
-------------
ATTACK FIXED, POPULATION CHANGED.
The Fed-FDR implementation and R0--R3 attack functions are imported directly
from the frozen K=2 source-site Protocol v1.1.1. Only X generation changes.

MAIN SETTINGS
-------------
K=2, n1=n2=50, p=50, true support 1:5.
Formal: N=10,000 attempted repetitions PER MAIN POPULATION.
Main populations P0--P4 are defined in TaskB_FROZEN_PROTOCOL_v1_0.md.

RUN ORDER
---------
1) 02_TEST_P0_BACKWARD_COMPATIBILITY.R
   Must print Overall: PASS.
   This verifies the new P0 branch reproduces Frozen v1.1.1 exactly.

2) 03_RUN_QUICK_MAIN_POPULATIONS.R
   N=5/population, oracle diagnostic M=1000.
   Purpose: code paths + DGP sanity only. DO NOT interpret attack accuracy.

3) 04_RUN_PILOT_MAIN_POPULATIONS.R
   N=100/population, oracle diagnostic M=2000, 8 workers.
   Inspect oracle signal diagnostics, P3 covariance conditioning, failures,
   and runtime. Do not tune based on whether an attack accuracy looks exciting.

4) 05_RUN_FORMAL_MAIN_POPULATIONS_N10000.R
   Only after Pilot is judged non-degenerate.
   N=10,000 attempted repetitions per P0--P4 = 50,000 attempts total.
   Checkpointed and resumable.

5) 06_SUMMARIZE_TASKB_FORMAL.R
   Produces teacher-ready formal tables and figures under:
   results_taskB_formal_N10000/teacher_summary/

OPTIONAL ONLY
-------------
7) 07_RUN_SENSITIVITY_TEMPLATE.R
   P1_IID_UNIFORM and P4_POISSON_RAW.
   Do not run unless the main P1/P4 result needs interpretation.

IMPORTANT P2 DETAIL
-------------------
P2 uses the SAME fixed block placement for both sites and all repetitions.
Variables 1:5 are deliberately dispersed into five different covariance blocks.
The mapping is exported as P2_FIXED_PERMUTATION.csv.

IMPORTANT P4 DETAIL
-------------------
The main P4 is a grid-count representation of a homogeneous spatial Poisson
process. lambda_k is expected count PER CELL. With 50 equal cells in [0,1]^2,
process intensity is 50*lambda_k per unit area. Main P4 uses theoretical
standardization (C-lambda_k)/sqrt(lambda_k). lambda_k is generator-only and is
NOT released to the attacker.

FORMAL VS PILOT DIAGNOSTICS
---------------------------
The large oracle diagnostic sample is used only in Quick/Pilot. It is OFF in
formal N=10,000 because it is not part of the attack estimand and would waste
substantial computation.
