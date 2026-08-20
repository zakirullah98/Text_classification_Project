# Task A Frozen Attack Definition v1.0

## Target

For a known training member \(\widetilde Z=(\widetilde X,\widetilde Y)\), infer
\(J\in\{1,\ldots,K\}\), with \(P(J=k)=1/K\).

For every regime, construct site scores \(Q_1,\ldots,Q_K\) and predict

\[
\widehat J = \arg\max_k Q_k.
\]

If several sites attain the maximum, use expected uniform tie credit.

## R0

\[
Q_{0,k}=1/K.
\]

## R1: support access

Let

\[
s_{kj}=\mathbf 1\{j\in\widehat{\mathcal S}^{(k)}\}.
\]

Define the site-specific selection weight

\[
w_{kj}=s_{kj}\left[1-\frac{1}{K-1}\sum_{\ell\ne k}s_{\ell j}\right].
\]

Then

\[
u_{S,k}(x)=
\begin{cases}
\dfrac{\sum_{j=1}^p w_{kj}x_j^2}{\sum_{j=1}^p w_{kj}},&\sum_jw_{kj}>0,\\
0,&\sum_jw_{kj}=0,
\end{cases}
\]

and

\[
Q_{S,k}(x)=\frac{u_{S,k}(x)}{\sum_{\ell=1}^K u_{S,\ell}(x)+\varepsilon_0}.
\]

R1 is identical for X-only and (X,Y).

## R2: refined-coefficient access

### X-only

\[
r_k(x)=x^\top\widehat\beta^{(k)},\qquad
u_{\beta,X,k}=|r_k(x)|,
\]

\[
Q_{\beta,X,k}=\frac{u_{\beta,X,k}}{\sum_{\ell=1}^K u_{\beta,X,\ell}+\varepsilon_0}.
\]

### (X,Y)

Let \(y^*=2y-1\). Define

\[
m_k(x,y)=y^*x^\top\widehat\beta^{(k)},
\]

\[
u_{\beta,XY,k}=\operatorname{expit}(m_k),
\]

\[
Q_{\beta,XY,k}=\frac{u_{\beta,XY,k}}{\sum_{\ell=1}^K u_{\beta,XY,\ell}+\varepsilon_0}.
\]

Because expit is strictly increasing, this has the same decision as
\(\arg\max_k m_k\), equivalently the smallest slope-only logistic loss.

## R3: support + coefficient

\[
Q_{\mathrm{Full},X,k}=Q_{S,k}+Q_{\beta,X,k},
\]

\[
Q_{\mathrm{Full},XY,k}=Q_{S,k}+Q_{\beta,XY,k}.
\]

No weight is tuned.

## K=2 backward compatibility

When \(K=2\),

\[
w_{1j}=s_{1j}(1-s_{2j}),\qquad
w_{2j}=s_{2j}(1-s_{1j}),
\]

so R1 reduces exactly to the old asymmetric-support definition. Moreover,
for every normalized component,

\[
Q_1-Q_2=D,
\]

and therefore the K-general R3 decision reduces exactly to the old
\(D_S+D_\beta\) rule.

## Main metrics

\[
\operatorname{Acc}_{\rm mem}(K),\qquad
\operatorname{Acc}_{\rm mem}(K)-1/K,
\]

\[
\operatorname{Acc}_{\rm fresh}(K),\qquad
\operatorname{Acc}_{\rm mem}(K)-\operatorname{Acc}_{\rm fresh}(K).
\]
