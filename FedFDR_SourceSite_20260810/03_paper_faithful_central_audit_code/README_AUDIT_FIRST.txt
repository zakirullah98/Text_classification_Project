Fed-FDR Paper-Faithful Central-Step Audit v1
============================================

WHAT THIS AUDIT DOES
--------------------
This is NOT a new privacy attack experiment and does NOT change Frozen Protocol v1.1.
It keeps all upstream simulation/fitting steps identical and re-audits only the central
Fed-FDR selection step.

Three central variants are computed on the SAME repetition:

A: CODE_THETA
   theta_mirror = raw refined beta / sqrt(SE)
   + the public-repository-compatible central logic used in the formal privacy run.

B: CODE_RAWBETA
   raw zero-padded refined beta
   + the same public-repository-compatible central logic.
   Comparing A vs B isolates the mirror-input scaling difference.

C: PAPER_RAWBETA
   raw zero-padded refined beta
   + Algorithm 1 as printed in the JRSSB paper.
   Comparing B vs C isolates differences in threshold/aggregation implementation.

WHY BOTH PAIRWISE AND FINAL RESULTS ARE SAVED
---------------------------------------------
K=2 has only one site pair. Every feature selected by that pair receives the same inclusion
rate 1/|S^(12)|. Therefore the final inclusion-rate step in the printed Algorithm 1 has severe
ties and can behave degenerately. The audit reports:
  * pairwise mirror support / FDP / Power
  * final inclusion-aggregated support / FDP / Power
separately.

RUN ORDER
---------
1) Open this folder as the RStudio working directory.
2) Run:
      source("08_RUN_AUDIT_QUICK_R20.R")
3) Send/check the R20 results if desired.
4) If the quick audit is clean, run:
      source("09_RUN_AUDIT_FORMAL_R10000.R")
5) Then run:
      source("10_ANALYZE_CENTRAL_AUDIT.R")

COMPUTING SETTINGS
------------------
In 09_RUN_AUDIT_FORMAL_R10000.R you may change only:
  N_WORKERS <- 8L
  BATCH_SIZE <- 100L
These affect runtime only, not the statistical audit.

IMPORTANT
---------
The main privacy attack quantities (R1/R2/R3 source-site attack scores) are NOT recalculated
or changed here. This audit only asks how final Fed-FDR support/FDP/Power change when the
central mirror step follows the paper's raw refined beta input and Algorithm 1 text.

v1.0.1 implementation patch:
- Fixed QUICK R20 cluster cleanup. The original quick script could try to stop the same PSOCK cluster twice / invoke on.exit at script scope, causing: defaultCluster(cl): no cluster supplied or registered.
- No statistical method, seed, Fed-FDR calculation, or audit definition changed.


VERSION 1.0.2 NOTE
------------------
Before the formal R=10000 audit, the paper-threshold boundary convention was
made explicit. Because Algorithm 1 writes min over tau>0 with strict M<-tau
and M>tau, a strictly positive minimum need not exist when the FDP estimate is
already <= alpha for arbitrarily small positive tau. Version 1.0.2 records the
boundary-limit tau=0 whenever #{M<0}/(#{M>0} v 1) <= alpha, not only when
there are no negative mirror statistics. No upstream fitting or privacy-attack
definition is changed.
