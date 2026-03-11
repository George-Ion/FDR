\documentclass[12pt]{article}
\usepackage[utf8]{inputenc}
\usepackage{amsmath, amssymb, amsthm}
\usepackage{geometry}
\usepackage{listings}
\usepackage{xcolor}
\usepackage{hyperref}

\geometry{a4paper, margin=1in}

\definecolor{codegreen}{rgb}{0,0.6,0}
\definecolor{codegray}{rgb}{0.5,0.5,0.5}
\definecolor{codepurple}{rgb}{0.58,0,0.82}

\lstset{
    language=R,
    basicstyle=\ttfamily\small,
    keywordstyle=\color{blue},
    stringstyle=\color{codepurple},
    commentstyle=\color{codegreen},
    breaklines=true,
    showstringspaces=false
}

\title{Computational Statistics: Repository Analysis and Notes}
\author{Documentation for GitHub: George-Ion/FDR}
\date{March 2026}

\begin{document}

\maketitle

\section{Overview}
This repository serves as a computational laboratory for high-dimensional statistical inference and stochastic simulation. The primary focus is on managing the \textbf{False Discovery Rate (FDR)} and implementing iterative algorithms for parameter estimation.

\section{Core Concepts and Algorithms}

\subsection{False Discovery Rate (FDR)}
In multiple hypothesis testing, we control the expected proportion of false rejections. The Benjamini-Hochberg (BH) procedure ensures:
\begin{equation}
E\left[ \frac{V}{R \vee 1} \right] \le \frac{m_0}{m} \alpha \le \alpha
\end{equation}
Where $V$ is the number of false positives and $R$ is the total number of rejections.



\subsection{Expectation-Maximization (EM) Algorithm}
Used for Maximum Likelihood Estimation in the presence of latent variables $\mathbf{Z}$. The algorithm iterates between:
\begin{itemize}
    \item \textbf{E-Step:} Calculate the expected log-likelihood $Q(\theta|\theta^{(t)}) = E_{\mathbf{Z}|\mathbf{X},\theta^{(t)}}[\log L(\theta; \mathbf{X}, \mathbf{Z})]$.
    \item \textbf{M-Step:} Update the parameter $\theta^{(t+1)} = \arg \max_{\theta} Q(\theta|\theta^{(t)})$.
\end{itemize}



\subsection{Markov Chain Monte Carlo (MCMC)}
The \texttt{metropolis\_hastings.R} script implements sampling from a target distribution $P(x)$. The acceptance probability $\alpha$ for a move from $x$ to $x'$ is given by:
\begin{equation}
\alpha(x, x') = \min \left( 1, \frac{P(x') q(x|x')}{P(x) q(x'|x)} \right)
\end{equation}



\section{Implementation Details}
The following R scripts are provided in the root directory:
\begin{description}
    \item[\texttt{EM\_algorithm.R}] Handles Gaussian Mixture Models (GMM) and missing data convergence.
    \item[\texttt{metropolis\_hastings.R}] Generates posterior samples for non-conjugate priors.
    \item[\texttt{bootstrap\_jackknife.R}] Computes non-parametric confidence intervals and bias correction.
\end{description}

\section{Summary of Computational Trade-offs}
\begin{center}
\begin{tabular}{|l|l|l|}
\hline
\textbf{Method} & \textbf{Advantage} & \textbf{Complexity} \\ \hline
MLE (Direct) & Efficient for simple models & Often non-analytical \\ \hline
EM Algorithm & Guaranteed likelihood increase & Sensitive to initialization \\ \hline
MCMC & Samples from any distribution & High autocorrelation / Burn-in \\ \hline
Bootstrap & No distributional assumptions & $O(B \cdot n)$ iterations \\ \hline
\end{tabular}
\end{center}

\end{document}
