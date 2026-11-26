# Writer's helper (loop 0): assembles main.tex from per-section strings; not part of the paper.
import re, sys
S = {}

S['abstract'] = r'''\begin{abstract}
Terahertz time-domain spectroscopy can in principle return the thickness and complex refractive index of every layer of a stratified sample from one transmission trace, but the inverse problem is unstable once the sample is thin or has several layers. We fit a differentiable transfer-matrix model to the measured pulse with the gradient optimisers Adam and L-BFGS and extract $n$, $\kap$ and $d$ of every layer. Simulations with a measured reference pulse locate the film thickness below which a thin film can no longer be separated from its substrate, and a polyimide film on silicon at that threshold is recovered in agreement with profilometry. On a quartz-on-silicon bilayer the sensitivity matrix of the forward model constrains the total optical path three orders of magnitude better than the extinction coefficients and separates $n$ from $d$ only through the Fabry--P\'erot echoes; the optimiser trajectory, projected onto the principal planes of its own path, shows Adam creeping along a flat $\kap$ valley for most of its budget, whereas L-BFGS completes the well-conditioned part of the descent in a few hundred evaluations. Cram\'er--Rao bounds show that the negative silicon $\kap$ of the frequency-independent fit is a model artefact, removed by a six-knot dispersive model of a quartz--air--silicon stack.
\end{abstract}
'''

# ---- target 190 words
S['intro'] = r'''
% =====================================================================================
\section{Introduction}
Terahertz time-domain spectroscopy (THz-TDS) yields the complex refractive index of a sample from the transmitted field of a picosecond pulse \cite{Grischkowsky1990,Duvillaret1996,Naftaly2007}. For an optically thick slab the classical frequency-by-frequency inversion \cite{Duvillaret1996} is stable and noise-limited \cite{Withayachumnankul2008}; for thin samples or several layers internal reflections overlap, the transfer function admits several solutions and the objective becomes non-convex \cite{Scheller2009}, so thickness and index can no longer be separated, precisely in the thickness-gauging regime where THz-TDS is most in demand \cite{Yasui2005}. Fitting the whole time trace \cite{Peretti2019} replaces the inversion by a non-convex optimisation whose behaviour is rarely examined.

Here the extraction is treated as an inverse problem \cite{Aster2018} with the transfer-matrix method (TMM) \cite{BornWolf1999} as forward model, implemented with automatic differentiation \cite{Paszke2019} so that $d$, $n$ and $\kap$ of every layer are fitted simultaneously with Adam \cite{Kingma2017} or L-BFGS \cite{Nocedal2006}. Simulations with a measured reference pulse locate where a thin film on a substrate ceases to be extractable, and a polyimide film on silicon confirms the prediction. The loss landscape of an experimental quartz-on-silicon fit is then analysed through the sensitivity matrix of the forward model and, borrowing from deep learning \cite{Li2018}, through projections of the optimiser trajectory onto the planes in which it actually moves; a six-knot dispersive model of a quartz--air--silicon stack closes the paper.
'''

# ---- target 80 words + eqs
S['tmm'] = r'''
% =====================================================================================
\section{Forward model: transfer-matrix method}\label{sec:tmm}
For a pulse entering layer $j$ of complex index $\tilde n_j = n_j - i\kap_j$ and thickness $d_j$ from a medium of index $n_0$, with $D_j=\omega d_j/c$, continuity of the tangential fields gives \cite{BornWolf1999}
\begin{align}
  r_j &= \frac{i}{2}\,\frac{\left(\frac{n_0}{\tilde n_j}-\frac{\tilde n_j}{n_0}\right)\sin(\tilde n_j D_j)}{\cos(\tilde n_j D_j)+\frac{i}{2}\left(\frac{n_0}{\tilde n_j}+\frac{\tilde n_j}{n_0}\right)\sin(\tilde n_j D_j)}, \label{eq:r}\\
  t_j &= \frac{1}{\cos(\tilde n_j D_j)+\frac{i}{2}\left(\frac{n_0}{\tilde n_j}+\frac{\tilde n_j}{n_0}\right)\sin(\tilde n_j D_j)}, \label{eq:t}
\end{align}
combined recursively over the layers into the transmission $T(\omega)$ of the stack. The sample pulse is the measured reference pulse $x(t)$ filtered by $T$,
\begin{equation}
  y(t)=\mathcal{F}^{-1}\!\left[T(\omega)\,\mathcal{F}[x(t)]\right],
  \label{eq:forward}
\end{equation}
on a zero-padded grid ($4L$ points for $L$ samples). Every operation is differentiable.
'''

# ---- target 250 words + eqs
S['inverse'] = r'''
% =====================================================================================
\section{Inverse problem: parameters, loss and optimisers}\label{sec:inverse}
Until Section~\ref{sec:disp} each layer has a frequency-independent index, so $N$ layers give $3N$ parameters $\{n_j,\kap_j,d_j\}$, optimised as dimensionless factors $s_i$ with $p_i = s_i\,p_i^{(0)}$ about the initial guess $p^{(0)}$, so that a step is the same fractional change for every parameter; the landscape analysis uses the same scaling, $x_i=(p_i-p_i^{\mathrm{opt}})/p_i^{(0)}$. The prediction $\hat y(t;\theta)$ is compared with the measured pulse by the root-mean-square error for the simulations and for the experiments by
\begin{equation}
  \Lt_{\mathrm{exp}}(\theta)=\frac{1}{L}\sum_i\left|\hat y_i-y_i\right| + \frac{1}{L}\sum_i \ln\cosh\!\left(\hat y_i-y_i\right),
  \label{eq:lossexp}
\end{equation}
which tolerates the few large residuals near the main pulse that a frequency-independent model cannot remove and is V-shaped around the optimum. Throughout, the noise floor $\sigma$ is the standard deviation of the first \SI{10}{\pico\second} of the sample trace, before the pulse (of the measured reference for the simulated sweep).

Gradients come from automatic differentiation through the forward model. Adam \cite{Kingma2017} keeps moving averages $\hat m_t$, $\hat v_t$ of the gradient and of its square and updates
\begin{equation}
  \theta_{t+1}=\theta_t-\eta\,\hat m_t/(\sqrt{\hat v_t}+\epsilon),
  \label{eq:adam}
\end{equation}
so every parameter moves at roughly the same rate $\eta$ however weakly the loss depends on it, and in a flat direction Adam keeps walking as long as the gradient keeps its sign. L-BFGS \cite{Nocedal2006} builds a curvature estimate from the gradient history and sets the step by a strong-Wolfe line search. Both are local and the basin of a THz fit is narrow, so each fit is run from several starts around the nominal guess and the lowest loss kept; Bayesian optimisation of the start, used in earlier work, was abandoned because in fifteen dimensions it moved the start into the wrong echo-assignment basin at a cost of more evaluations than it saved.
'''

# ---- target 50 words
S['sim5'] = r'''
% =====================================================================================
\section{Simulated extractions and a thin film at the limit}\label{sec:sim}
A noise-free five-layer stack with fifteen free parameters ($n=3.4/2.5/5.4/4.2/1.7$, $\kap=0.006/0.004/0.005/0.001/0.01$, $d=70/200/100/200/\SI{50}{\micro\metre}$), fitted by multi-start L-BFGS from a hand guess within $n\pm0.3$, $\kap\pm0.005$, $d\pm\SI{30}{\percent}$, is recovered to \SI{0.001}{\percent} in $n$, \SI{0.002}{\percent} in $d$ and \SI{0.3}{\percent} in $\kap$: the method handles the stack, and the extinction coefficients are two orders of magnitude less accurate than the indices.
'''

# ---- target 230 words
S['limits'] = r'''
For the most demanding geometry, a thin film on a thick substrate, we scan the film thickness across the point at which the film's internal reflection separates from the main pulse, using the measured reference of Section~\ref{sec:exp} and white noise at its $\sigma$ (three draws, six free parameters, multi-start L-BFGS). Figure~\ref{fig:limits}(a) shows the film's errors: $n$ and $d$ degrade gracefully as the film thins and fail steeply below a round trip of about \SI{0.24}{\pico\second} (\SI{10}{\micro\metre}), once the film echo no longer leaves the main pulse; $\kap$ is one to two orders of magnitude worse at every thickness. The converged loss is $(0.98\text{--}1.00)\,\sigma$ throughout, so the errors are set by the noise, and the echoes \cite{Scheller2009,Yasui2005} are what separates $n$ from $d$.

\begin{figure}[tb]\centering
  \includegraphics[width=0.72\linewidth]{fig2_thin_film}
  \caption{Thin film on a substrate. (a) Simulated sweep with the measured reference (film $n=3.6$, $\kap=0.0203$; substrate $n=2.213$, $\kap=0.0125$, $d=\SI{496}{\micro\metre}$; $\sigma=\num{1.58e-3}$, power SNR \SI{44}{dB}): median (symbols) and range (bands) of the film's errors over three noise draws; top axis, round trip $2n_{\mathrm f}d/c$. Medians at 5, 10, \SI{20}{\micro\metre}: $n$ 1.9, 0.29, \SI{0.03}{\percent}; $d$ 4.8, 0.99, \SI{0.08}{\percent}; $\kap$ 25, 8.3, \SI{6.3}{\percent}. (b) Polyimide on silicon, measured and fitted pulse with residual.}
  \label{fig:limits}
\end{figure}

The experimental counterpart is a spin-coated polyimide film (HD Microsystems PI-2525) on a high-resistivity silicon wafer \todo{[sample details]}, its thickness known from profilometry and its visible index from white-light interferometry (Table~\ref{tab:pi}). Its round trip of \SI{0.30}{\pico\second} sits at the threshold of Fig.~\ref{fig:limits}(a), so the extraction should be marginal in $n$ and $d$ separately and good in $n\,d$, and so it is: six of eight starts agree to within \SI{1}{\percent} at a loss of $2.7\sigma$, but the split puts $d_{\mathrm{PI}}$ \SI{4.6}{\percent} above the profilometer value and $n_{\mathrm{PI}}$ \SI{0.8}{\percent} below $n_\parallel(\SI{633}{\nano\metre})$, with $\rho(n_{\mathrm{PI}},d_{\mathrm{PI}})=-0.93$ and a discrepancy eighteen times the Cram\'er--Rao bound. Two starts end in a $\kap$-driven local minimum (Table~\ref{tab:pi}), and a plain full-band fit from \SI{16}{\micro\metre} collapses to $d_{\mathrm{PI}}=\SI{22.3}{\micro\metre}$, $n_{\mathrm{PI}}=1.93$ at \SI{48}{\percent} higher loss: the simulated failure mode in the flesh.

\begin{table}[tb]\centering\footnotesize\setlength{\tabcolsep}{3.5pt}
  \caption{Polyimide on silicon ($\sigma=\num{1.35e-3}$; thicknesses in \si{\micro\metre}): best of eight L-BFGS starts, Cram\'er--Rao bound, the local minimum of two starts, and the independent values ($n_\parallel$ at \SI{633}{\nano\metre}, profilometer, silicon index from Ref.~\cite{Dai2004}).}
  \label{tab:pi}
  \begin{tabular}{l cccccc c c}\toprule
    & $n_{\mathrm{PI}}$ & $\kap_{\mathrm{PI}}$ & $d_{\mathrm{PI}}$ & $n_{\mathrm{Si}}$ & $\kap_{\mathrm{Si}}$ & $d_{\mathrm{Si}}$ & $n_{\mathrm{PI}}d_{\mathrm{PI}}$ & rMSE$/\sigma$\\\midrule
    THz fit (best) & 1.712 & $<10^{-6}$ & 25.84 & 3.259 & 0.0006 & 400.2 & 44.26 & 2.69\\
    Cram\'er--Rao $1\sigma$ & 0.0024 & 0.0018 & 0.062 & 0.0002 & 0.0001 & 0.030 & 0.052 & \\
    local minimum & 1.611 & 0.0107 & 28.47 & 3.255 & $<10^{-6}$ & 401.3 & 45.87 & 3.34\\
\midrule
    independent & 1.726 & & 24.7 & 3.4175 & & & 42.6 & \\
    \bottomrule
  \end{tabular}
\end{table}
'''

# ---- target 150 words
S['exp'] = r'''
% =====================================================================================
\section{Experimental extraction: quartz on silicon}\label{sec:exp}
A quartz plate bonded to a silicon wafer \todo{[bonding method, nominal thicknesses and supplier to be added]} was measured with a Menlo Systems TeraSmart spectrometer under nitrogen purge (1000 waveforms over \SI{15}{\minute}, 2425 points at \SI{33}{\femto\second}; $\sigma=\num{2.9e-4}$, peak-to-$\sigma$ ratio $7.9\times10^3$) and fitted with Eq.~(\ref{eq:lossexp}) and Adam at $\eta=10^{-4}$ for $2\times10^4$ iterations from three starts bracketing the expected quartz index (Fig.~\ref{fig:expfit}, Table~\ref{tab:exp}). The starts converge to the same point to four figures in $n$ and $d$; the quartz index is that of fused silica \cite{Grischkowsky1990,Naftaly2007} rather than crystalline quartz \todo{[confirm sample material]}, the silicon index is within a percent of the literature value \cite{Dai2004}, and the silicon extinction coefficient is negative, which is unphysical. The loss falls by an order of magnitude in the first \SI{8}{\percent} of the iterations and only slowly thereafter, and the mean absolute residual at the optimum is $15\sigma$ (Fig.~\ref{fig:expfit}(b)).

\begin{figure}[tb]\centering
  \includegraphics[width=0.72\linewidth]{fig3_exp_fit}
  \caption{Quartz/Si extraction. (a) Measured reference and sample pulses with the converged fit (start 1); inset, the echoes. (b) Residual against $\pm\sigma$, concentrated at the main pulse and the first echo. (c) Amplitude spectra; above \SI{1.5}{THz} the fit lies above the measurement, the signature of absorption increasing with frequency.}
  \label{fig:expfit}
\end{figure}

\begin{table}[tb]\centering\footnotesize\setlength{\tabcolsep}{3.5pt}
  \caption{Quartz/Si fits ($2\times10^4$ Adam iterations) from three starts ($n_1=1.982/1.900/2.010$, $d_1=1000/900/\SI{900}{\micro\metre}$, $\kap_1=0.0071$; silicon $3.40/0.0054/\SI{400}{\micro\metre}$), from start 1 with the softplus constraint $\kap\ge0$, and the Cram\'er--Rao bounds (Section~\ref{sec:cr}). Fits 2 and 3 lie within $2\times10^{-4}$ ($n_1$), \SI{0.15}{\micro\metre} ($d_1$), $1.5\times10^{-3}$ ($n_2$) and \SI{0.16}{\micro\metre} ($d_2$) of fit 1.}
  \label{tab:exp}
  \begin{tabular}{l cccccc c}\toprule
    & $n_1$ & $\kap_1$ & $d_1$ (\si{\micro\metre}) & $n_2$ & $\kap_2$ & $d_2$ (\si{\micro\metre}) & loss\\\midrule
    fit 1 & 1.9639 & 0.00541 & 994.28 & 3.4478 & $-0.00144$ & 381.41 & $4.4522\times10^{-3}$\\
    fit 1, $\kap\ge0$ (softplus) & 1.9595 & 0.00478 & 996.60 & 3.4681 & 0.00011 & 379.11 & $4.4646\times10^{-3}$\\\midrule
    Cram\'er--Rao $1\sigma$ & $6\times10^{-5}$ & $9\times10^{-6}$ & 0.032 & $3\times10^{-4}$ & $1.9\times10^{-5}$ & 0.031 & \\
    \bottomrule
  \end{tabular}
\end{table}
'''

# ---- target 190 words
S['sens'] = r'''
% =====================================================================================
\section{The loss landscape of the experimental fit}\label{sec:landscape}
\subsection{Sensitivity matrix and the two-phase descent}\label{sec:sens}\label{sec:pca}
In the scaled coordinates $x_i$ about the optimum of Table~\ref{tab:exp}, the geometry of the problem is described, independently of the cone-shaped loss, by the sensitivity matrix of the forward model \cite{Kay1993,Aster2018},
\begin{equation}
  G=\frac{1}{L}J^\top J,\qquad J_{ti}=\frac{\partial \hat y(t_t)}{\partial x_i},
  \label{eq:G}
\end{equation}
whose eigenvectors are the parameter combinations along which the pulse changes most and least and whose singular values, 8.66, 1.54, 0.69, 0.073, 0.019 and 0.0016, span a factor of 5500. The best-constrained direction is the total optical path, fixed by the pulse delay; after two directions that redistribute optical path comes a gap of an order of magnitude to the anti-phase combination ($n_2$ up, $d_2$ down) that leaves $n_2d_2$ unchanged, the optical-thickness degeneracy ($\rho(n_1,d_1)=0.97$, $\rho(n_2,d_2)=0.99$), lifted only by the echoes, whose Fresnel coefficients depend on $n$ alone. Another order of magnitude down are the two extinction coefficients: per unit fractional change, $\kap_1$ is 280 times and $\kap_2$ 890 times less visible than $n_1$.
'''

# ---- target 190 words
S['pca'] = r'''
How an optimiser moves through this landscape is best seen in the plane the trajectory itself chooses \cite{Li2018}: the two leading principal components of the centred iterates span it, the loss is evaluated on a grid in it and the path projected onto it. Over the whole run almost all of the path variance lies along $\kap_2$, yet the loss barely depends on that motion, so we split the run where every start has completed \SI{95}{\percent} of its loss decrease, iteration 1615 (Fig.~\ref{fig:twophase}), with one plane per start in curvature-normalised coordinates (axes scaled by $\sqrt{G_{ii}}$) for Phase~I and a shared plane for Phase~II. In Phase~I the loss drops by \SIrange{81}{92}{\percent} as each start descends a bowl whose axes are the best-constrained combinations of $n$ and $d$; in Phase~II, \SI{92}{\percent} of the budget, the starts drift along the same $\kap$ valley at Adam's full speed for a further \SIrange{22}{29}{\percent}.

\begin{figure}[tb]\centering
  \includegraphics[width=0.72\linewidth]{fig4_adam_two_phase}
  \caption{The Adam descent of the quartz/Si fit in the principal planes of its own path ($+$ optimum, $\bullet$ start, $\star$ end of phase). (a)--(c) Phase~I, iterations 1--1615, one panel per start in curvature-normalised coordinates. (d) Phase~II, shared plane (\SI{99.7}{\percent} of the path variance): the merged starts creep along the valley $0.97\kap_2+0.24\kap_1$, taking $\kap_2$ through zero. (e) Loss versus iteration.}
  \label{fig:twophase}
\end{figure}
'''

# ---- target 250 words
S['cr'] = r'''
\subsection{Noise-limited uncertainties, and why $\kap_{\mathrm{Si}}<0$ is a model artefact}\label{sec:cr}\label{sec:constrained}
Whether the negative silicon extinction coefficient is a statistical fluctuation or a systematic effect the sensitivity matrix decides, step by step. (i) If the model were exact and the residual white noise of standard deviation $\sigma$, the covariance of any unbiased estimator would obey the Cram\'er--Rao bound $\mathrm{cov}(\hat x)\ge\sigma^2 (J^\top J)^{-1}$ \cite{Kay1993}, a property of the measurement, not of the optimiser; (ii) $J$ is Eq.~(\ref{eq:G}) at the optimum. (iii) The square roots of the diagonal are the bounds in Table~\ref{tab:exp}: the indices at the $10^{-4}$ level, the thicknesses at \SI{31}{\nano\metre}, $\kap_2$ at $1.9\times10^{-5}$. (iv) The fitted $\kap_2=-1.5\times10^{-3}$ is therefore $77\sigma$ from zero: noise cannot have produced the sign. (v) Forcing $\kap_2\ge0$ and re-optimising everything else raises the loss by only \SI{0.29}{\percent}: the negative $\kap_2$ buys almost nothing. (vi) The sign is therefore the frequency-independent model absorbing the high-frequency amplitude excess of Fig.~\ref{fig:expfit}(c) through the one parameter the data leave nearly free: a model artefact, which the dispersive fit of Section~\ref{sec:disp} removes. The constraint in (v) is the map $\kap_j=\kap_j^{(0)}\,\mathrm{softplus}(\beta\theta_j)/\mathrm{softplus}(\beta)$, $\beta=5$, which keeps a gradient everywhere, whereas a rectifier freezes the parameter the moment it reaches zero; the softplus fit (Table~\ref{tab:exp}) drives $\kap_{\mathrm{Si}}$ to zero, as expected for high-resistivity silicon \cite{Dai2004}, and moves $n$ and $d$ of each layer in opposite directions by up to \SI{0.6}{\percent} while $n_1d_1$ and $n_2d_2$ change by \SI{0.015}{\percent}: along the valley, not out of it.
'''

# ---- target 90 words
S['constrained'] = r'''
'''

# ---- target 160 words
S['lbfgs'] = r'''
% =====================================================================================
\section{L-BFGS versus Adam}\label{sec:lbfgs}
Adam creeps along the $\kap$ valley at a fixed normalised speed, whereas a quasi-Newton method learns the curvature. We therefore repeated the fit with L-BFGS (strong-Wolfe line search, history 20) from the same three starts with the same scaled parameters and loss (Fig.~\ref{fig:lbfgs}). L-BFGS cuts across the Phase-I bowl instead of tracing Adam's L-shaped route and reaches, in 250--300 evaluations, the loss that Adam reaches after 5000--5900 iterations; then it stops. The stop is a kink, not a barrier: the line search finds no decrease even for steps of $10^{-8}$, because the mean-absolute term of Eq.~(\ref{eq:lossexp}) is not differentiable where residual samples change sign, yet the loss decreases monotonically along the straight line to Adam's end point, \SI{3.7}{\percent} lower. The two end points differ by \SI{5.7}{\micro\metre} in $d_1$ and $d_2$, with $\kap_{\mathrm{Si}}$ still positive at the kink, while $n_1d_1$ and $n_2d_2$ agree to \SI{0.012}{\percent} and \SI{0.034}{\percent}: what L-BFGS reaches twenty times cheaper is the part of the fit the data constrain; what it misses is the flat $\kap$ valley of Section~\ref{sec:landscape}, where the end point is decided by how the loss weights the model error. A smooth loss is the right companion for a quasi-Newton method: with the rMSE, or the log-cosh term alone, L-BFGS converges without a stall to yet another point of the $n\,d$ valley. A combined time- and frequency-domain loss adds no information, its Cram\'er--Rao bounds being those of the time trace.

\begin{figure}[tb]\centering
  \includegraphics[width=0.72\linewidth]{fig5_lbfgs_trajectory}
  \caption{L-BFGS on the quartz/Si fit in the planes of Fig.~\ref{fig:twophase}, every iterate marked, Adam's path drawn faintly. (a)--(c) Phase-I planes; (d) Phase-II plane, where L-BFGS stops at the kink before the valley; (e) loss versus function evaluations. End points (start 1; evaluations over the three starts): Adam, $2\times10^4$ iterations, $d_1=994.3$, $d_2=\SI{381.4}{\micro\metre}$, $\kap_2=-0.0014$, $n_1d_1=1.9526$, $n_2d_2=\SI{1.3150}{\milli\metre}$; L-BFGS on the same loss, 250--299 evaluations, 999.6, \SI{376.1}{\micro\metre}, $+0.0044$, 1.9528, \SI{1.3146}{\milli\metre}; L-BFGS on the rMSE (181--314 evaluations) or log-cosh (178--484) loss, 977.7, \SI{398.1}{\micro\metre}, $-0.0073$, 1.9523, \SI{1.3153}{\milli\metre}.}
  \label{fig:lbfgs}
\end{figure}

'''

# ---- target 330 words
S['disp'] = r"""
% =====================================================================================
\section{Letting the model be right: quartz--air--silicon with six frequency knots}\label{sec:disp}
Every diagnostic so far points to the missing frequency dependence of the optical constants, which we add in an interpretable, well-scaled form: $n_j(f)$ and $\kap_j(f)$ of each layer are the natural cubic interpolants of their values at six knots spanning the band in which the sample spectrum stays above the noise floor, \SIrange{0.2}{3.0}{THz}, held at the band-edge values outside it; each knot value is $\theta^{\mathrm{init}}+s\,\vartheta$ with a trainable $\vartheta$ and a fixed scale $s$, and $\kap\ge0$ is imposed by a hinge penalty $\lambda_\kap\langle\max(0,-\kap_j(f))\rangle$, $\lambda_\kap=1$. Adam runs in two stages, the knots of each curve tied in stage~1 (the frequency-independent problem) and freed in stage~2. The sample is the quartz plate and silicon wafer of Section~\ref{sec:exp} re-mounted with an air gap, a third layer with $n=1$, $\kap=0$ frozen and only its thickness free (Figs.~\ref{fig:qasfit}, \ref{fig:qasnk}, Table~\ref{tab:qas}).

\begin{figure}[tb]\centering
  \includegraphics[width=0.62\linewidth]{fig6_qas_fit}
  \caption{Quartz--air--silicon stack, six-knot model. (a) Measured and fitted pulse (stage~2); (b) residual after stage~1 (knots tied) and stage~2 (knots free) against $\pm\sigma$; (c) amplitude spectra, fit band shaded. Freeing the knots lowers the loss from \num{3.62e-3} to \num{2.93e-3} (MAE from $4.0\sigma$ to $3.2\sigma$) and removes the pulse-shape residual; what remains carries the Fabry--P\'erot periodicity of the echoes.}
  \label{fig:qasfit}
\end{figure}

\begin{figure}[tb]\centering
  \includegraphics[width=0.57\linewidth]{fig7_qas_nk}
  \caption{Extracted $n(f)$ (top) and $\kap(f)$ (bottom) of quartz (left) and silicon (right); knots marked, stage~1 (tied) dotted, fit band shaded. Silicon $\kap$ stays below $1.2\times10^{-3}$ except at the \SI{3}{THz} knot, where the spectrum reaches the noise floor.}
  \label{fig:qasnk}
\end{figure}

\begin{table}[tb]\centering\footnotesize\setlength{\tabcolsep}{3.5pt}
  \caption{Quartz--air--silicon fit: knot values after stage~2, stage-1 (tied) values, thicknesses, optical thicknesses at \SI{1}{THz} and losses ($\sigma=\num{9.0e-4}$).}
  \label{tab:qas}
  \begin{tabular}{l l cccccc c}\toprule
    & & \multicolumn{6}{c}{knot frequency (THz)} & stage 1\\\cmidrule(lr){3-8}
    & & 0.20 & 0.76 & 1.32 & 1.88 & 2.44 & 3.00 & (tied)\\\midrule
    quartz & $n$    & 1.9652 & 1.9712 & 1.9727 & 1.9745 & 1.9780 & 1.9917 & 1.9753\\
           & $\kap$ & 0.0006 & 0.0034 & 0.0068 & 0.0089 & 0.0110 & 0.0114 & 0.0050\\
    Si     & $n$    & 3.3963 & 3.3926 & 3.3928 & 3.3940 & 3.3916 & 3.3818 & 3.3903\\
           & $\kap$ & 0.0012 & 0.0000 & 0.0009 & 0.0002 & 0.0007 & 0.0048 & 0.0002\\\midrule
    \multicolumn{2}{l}{$d_{\mathrm{q}}$, $d_{\mathrm{air}}$, $d_{\mathrm{Si}}$ (\si{\micro\metre})} & \multicolumn{6}{c}{976.6, 125.9, 388.9} & 974.5, 126.8, 389.1\\
    \multicolumn{2}{l}{$n\,d$ quartz, Si (\si{\micro\metre})} & \multicolumn{6}{c}{1926.0, 1319.3} & 1924.9, 1319.3\\
    \multicolumn{2}{l}{loss, MAE$/\sigma$} & \multicolumn{6}{c}{\num{2.93e-3}, 3.2} & \num{3.62e-3}, 4.0\\
    \bottomrule
  \end{tabular}
\end{table}

Three points follow. First, the air gap is far better conditioned than the per-layer $n$--$d$ split, because it enters through the timing of two interfaces, not through an optical path that an index could compensate: a scan of the loss against $d_{\mathrm{air}}$ shows a single sharp minimum at \SI{126}{\micro\metre}, 5.5 times lower than at the nominal \SI{50}{\micro\metre}. Second, the quartz extinction rises almost linearly with frequency, as for fused silica rather than crystalline quartz \cite{Grischkowsky1990,Naftaly2007} \todo{[confirm sample material]}; this rising $\kap$ is what drove $\kap_{\mathrm{Si}}$ negative in Section~\ref{sec:cr}, and with it in the model the silicon is lossless. Third, the silicon index remains \SI{0.7}{\percent} below the literature value 3.4175 \cite{Dai2004}: its optical thickness is identical to \SI{0.01}{\percent} between stages, but the $n$--$d$ split rests on the echo amplitudes alone, and a loss scan against $n_{\mathrm{Si}}$ with $d_{\mathrm{Si}}$ re-optimised puts the valley floor near 3.37 (literature value: $+\SI{4}{\percent}$ in loss).
"""

S['disc'] = r'''
% =====================================================================================
\section{Discussion and conclusions}\label{sec:conclusions}
A differentiable TMM fitted by gradient optimisers degrades gracefully in $n$ and $d$ and erratically in $\kap$ as a film thins, failing once the film echo no longer leaves the main pulse; on the quartz/Si bilayer it converges reproducibly but slowly. All of this follows from the sensitivity spectrum of the measurement: the total optical path is constrained thousands of times better than the extinction coefficients, the per-layer $n$--$d$ split only by the echoes. L-BFGS reaches the data-constrained part of the optimum twenty times cheaper than Adam and, with a smooth loss, is the better default. Uncertainties from $(J^\top J)^{-1}$ \cite{Kay1993,Aster2018,Withayachumnankul2008} must be read alongside the ratio of residual to noise floor: when it is far above unity, the honest uncertainty on $d$ is the spread between fits that weight the model error differently, here \SI{2}{\percent}. \todo{[mechanically measured thickness to be included?]}

'''

ORDER = ['abstract','intro','tmm','inverse','sim5','limits','exp','sens','pca','cr','constrained','lbfgs','disp','disc']

def words(t):
    t = re.sub(r'\\begin\{(figure|table|align|equation)\}.*?\\end\{\1\}', '', t, flags=re.S)
    t = re.sub(r'\\cite\{[^}]*\}', '', t)
    t = re.sub(r'\\[a-zA-Z]+', '', t)
    return len(re.sub(r'[{}$^_\\]', ' ', t).split())

if __name__ == '__main__':
    tot = 0
    for k in ORDER:
        w = words(S[k]); tot += w
        caps = re.findall(r'\\caption\{(.*?)\}\n', S[k], flags=re.S)
        print(f'{w:5d}  {k:12s}  captions {[len(c.split()) for c in caps]}')
    print('prose total', tot)
    if '--write' in sys.argv:
        src = open('/home/claude/work/paper_v3/main.tex').read()
        head = src[:src.index('\n\\begin{abstract}')+1]
        tail = src[src.index(r'\begin{backmatter}'):]
        open('/home/claude/work/paper_v3/main.tex', 'w').write(head + ''.join(S[k] for k in ORDER) + tail)
        print('written')
