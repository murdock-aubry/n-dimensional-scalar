# Table of Contents

- $\texttt{scalar_n.f90}$: Fortran file which contains the main routines used in the algorithm for the general $n$ dimensional case.

- $\texttt{experiment\_i.f90}$: Testing the the algorithm against the lower dimensional code, assuring accuracy.

- $\texttt{odesolve\_n.tex}$: A brute force ODE solver for $n$ dimensional equations used as a step in this algorithm.

- $\texttt{partitions.nb}$: A mathematical file with contains for for computing all partitions of positive integers. This data is then compiled in an convenient way and used in a step of the algorith proper.


# Introduction

A popular and well-understood class of ordinary differential equations, commonly refered to as homogenous scalar equations, take on the form
$$y^{(n)}(t) + q_{n-1}(t)y^{(n-1)} + \cdots + q_{1}(t)y'(t) + q_0(t)y(t) = 0.$$

The cost of numerically representing solutions to this class of equations using standard methods increases with the magnitude of the coefficient functions $\{q_i\}_{i=0}^{n-1}$. However, the phase functions of these equations are solvable, independent of the magnitude of the coefficient functions. More information on the theory, algorithm description, implementation and testing can be found in this [paper](https://arxiv.org/abs/2308.03288), authored by James Bremer and myself. 

As you may notice, the implementation and description of the algorithm in the paper only considers for the case where $n=2, 3$, or $4$. Utitlizing the code for the lower dimensional cases, I independently generalized the algorithm to the case where the value of $n$ is unknown at compile time. Further, I conducted several experiments, both comparisons to the lower dimensions to assure accuracy, and to the general cases for large values of $n$. The results are conclusive to run-time independence of the magnitude of the coefficient functions, and quadratic proportionality to the value of $n$.

# Algorithm Description

We now describe the global Levin method; the first of two algorithms considered in this paper. This algorithm computes the $n$ non-oscillatory phase functions $\{r_i(t)\}_{i=1}^n$ of the $n^\text{th}$ order Riccati equation, provided their values at some point $c$ in the domain of interest $[a, b]$. These functions can then be used to find the solutions of (\ref{introduction:scalarode}) via the substitution (\ref{introduction:phaserep2}). It operates by in two steps. first by computing initial estimates for the basis functions of the solution space, and secondly by using this initial estimate to initialize an adaptivereccursion utilizing the Levin method.

## The global Levin Method

Before describing the algorithm proper, we first discuss a subprocedure which estimates the solution of the Riccati equation on a subinterval $[a_0, b_0]\subset [a, b]$. This procedure is used as a step in the iterative process of the algorithm proper, and operates by slightly perturbing a current estimate of the solution through the addition of some small function $\delta(t)$:
$$
    \begin{equation}
        r_{m+1}(t) = r_m(t) + \delta(t)
    \end{equation}
$$
Inserting $(2)$ into the $n^{\text{th}}$ order Riccati equation, we obtain a discretizable equation. For example, in the second order case, we obtain from $(1)$
$$
    \begin{equation}
        \delta'(t) + (2r_m(t) + q_1(t))\delta = - r_m'(t) - r_m^2(t) - q_1r_m(t) - q_0(t)
    \end{equation}
$$
A general expression for the discretizable Riccati equation is derived bellow, and is a pivotal step in the completetion of this project. This subprocedure takes as input the interval $[a_0, b_0]$, an integer $\ell$ controlling the number of Chebyshev nodes used in the discretization process, estimates of the value of $r$ at the $\ell-$point Chebyshev quadrature, and an external subroutine for the evaluation of the coefficient functions $q_0, q_1, \ldots, q_{n-1}$. In return, the algorithm provides estimates for the values of the solution to the initial value problem at the Chebyshev quadrature. The subroutine proceeds as follows.

1. Construct an initial $k$ point Chebyshev extremal quadrature $\{t_{i, \ell}\}_{i=1}^\ell$ on the interval $[a, b]$, and use the provided external subroutine to evaluate the $n$ coefficient functions $q_0,q_1, \ldots, q_{n-1}$ on said quadrature.

2. Estimate the vector values of the first $(n-1)$ derivatives of $r$ at the Chebyshev nodes through repeated application the spectral differentiation matrix $\mathcal{D}_\ell$ to the vector values of $r$. In particular,
$$
    \begin{equation}
        \begin{pmatrix} r^{(i)}(t_{i, \ell}) \\ \vdots \\ r^{(i)}(t_{\ell, \ell}) \end{pmatrix} = 
        \mathcal{D}_k^{i}\begin{pmatrix} r(t_{i, \ell}) \\ \vdots \\ r(t_{\ell, \ell}) \end{pmatrix}
    \end{equation}
$$

3. Construct the $\ell \times\ell$ matrix $A$ and vector $[\tilde{r}]$ which discretizes the perturbed $n^{th}$ order Riccati equation. For example, in the case where $n = 2$, we obtain from $(3)$
$
    \begin{equation}
        \begin{aligned}
            \mathcal{A} &= \mathcal{D}_\ell + \text{diag}[2r + q_1] \\
            [\tilde{r}] &= -([r'] + [r^2] + [q_1r] + [q_0])
        \end{aligned}
    \end{equation}
$
which discretizes the left and right hand sides of $(3)$;
$
    \begin{equation}
        \mathcal{A}[\delta] = [\tilde{r}]
    \end{equation}
$    
        
4. Solve the system (\ref{levin:disc4}) for $[\delta]$ via a QR decomposition.
        
5. Return the values $[r] + [\delta]$ as the new estimate for the solution on the interval $[a_0, b_0]$.
    
The algorithm proper takes as input the domain of interest $[a, b]$, the desired values of the basis functions at some point $c\in [a, b]$, as well as an external subroutine used to evaluate the coefficient functions. The algorithm then returns the values of the non-oscillatory phase functions on a piecewise Chebyshev extremal grid. It maintains two lists of intervals; one storing the intervals which are yet to be utilized, and one storing the list of accepted intervals. Initially, the first list of intervals only contains the original interval $[a, b]$. As long as the first list is non-empty the following steps are repeated:
    
1. Remove an entry $[a_0, b_0]$ from the first list, and construct an initial $\ell-$ point grid Chebyshev extremal quadrature $\{t_{i, \ell}\}_{i=1}^\ell$ on the entire interval
$[a, b]$. Use the provided external subroutine to evaluate the $n$ coefficient functions $q_0,q_1, \ldots, q_{n-1}$ at on said quadrature. 

2. If the list of accepted interval is empty, then for each quadrature point $t_{i, k}$ use the previously mentioned root solver to compute the $n$ roots of the $n^{\text{th}}$ order complex polynomial
$$
    \begin{equation}
        p_i(z) = z^n + q_{n-1}(t_{i,\ell}) + \cdots + q_1(t_{i,\ell})z + q_0(t_{i, \ell}).
    \end{equation}
$$
Ordering the roots in such a way as to minimize the proximity of each root to the roots of $p_{i=1}(z)$, we obtain the $\ell$ values of the $n$ initial guess functions $[r_{10}], \ldots [r_{1n}]$. If the list of accepted intervals is non-empty, the utilize the values of the solutions at the endpoint of the adjacent accepted interval as the constant initial guess functions

1. Using these initial estimates, new estimates $[r_{11}], \ldots [r_{n1}]$ are computed using the subprocedure described above.

2. Compute the $n$ convergence paramaters
$
    \begin{equation}
        \xi_i = \frac{\|[\delta]\|_{L^2}}{\|[r_{0i}]\|_{L^2}} \ \ \ \ i = 1, \ldots, n.
    \end{equation}
$
If $\text{max}_{i=1, \ldots, n}\xi_i < \epsilon_0$, then exit the alogrithm and return the current estimate.

3. Otherwise, subdivide the interval $[a, b]$ into the equal sectors $[a, (a+b)/2], [(a+b)/2, b]$ and exit this loop.

4. Otherwise, we have found an accepted collection of values $[r_1], \ldots, [r_n]$ the some subinterval $[a_0, b_0]$. Add the subinterval $[a_0, b_0]$ to the list of accepted intervals. Steps 3 and 4 are now repeated with the interval directly adjecent to the one just accepted, with the initial guess being given by the constant functions of value equal to $r_1(b_0), \ldots, r_n(b_0)$ which are taken from the value which were just accepted.

In the end, we obtain accurate estimates of the values of the $n$ non-oscillatory phase functions on a peicewise Chebyshev structure, which can then be used to interpolate these functons
anywhere on the entire interval $[a, b]$.



## Local Levin method

In this section we describe the local Levin method, the second of the two alogorthms described and utilized in this paper. The algorithm operates in three steps. First, reasonable initial estimates are made for the values of the derivatives of the $n$ phase functions $r_1,\ldots, r_n$ on a Chebyshev quadrature, denoted $[r_1], \ldots, [r_n]$. Secondly, these initial estimates are iterated upon via the subprocedure described in the previous section to obtain an accurate estimate of the values of the $n$ solutions at some point $c$ in the domain of interest $[a, b]$. Finally, the general ODE solver described in Appendix (1) of the aforemention paper, is used, in conjuction with the estimated initial values, to solve for the derivatives of all $n$ phase functions on the entire interval. The algorithm takes as input the entire interval of interest $[a, b]$, an integer $\ell$ controlling the number of Chebyshev nodes used in the discretization process, a small subinterval $[a_0, b_0]\subset [a, b]$ containing the point $c_0$ at which to compute estimates of the first $(n-1)$ derivatives of the solutions $r_1, \ldots, r_n$, and an external subroutine used to compute the values of the coefficient functions $q_0, q_1, \ldots, q_{n-1}$. The algorithm outputs the values of the first $(n-1)$ derivatives of the solutions $r_1, \ldots, r_n$ on a piecewise Chebyshev structure over the interval $[a, b]$. The first step of the algorithm proceeds as follows.

1. Construct an $\ell-$ point Chebyshev extremal quadrature $\{t_{i, \ell}\}_{i=1}^\ell$ on the subinterval $[c_0, b_0]$, and use the provided external subroutine to evaluate the $n$ coefficient functions $q_0,q_1, \ldots, q_{n-1}$ on said quadrature. The algorithm can be easily modified to instead utilize the subinterval $[a_0, c_0]$.

2. For each quadrature point $t_{i, k}$ use the previously mentioned root solver to compute the $n$ roots of the $n^{\text{th}}$ order complex polynomial
$
    \begin{equation}
        p_i(z) = z^n + q_{n-1}(t_{i,\ell})z^{n-1} + \cdots + q_1(t_{i,\ell})z + q_0(t_{i, \ell}).
    \end{equation}
$
Ordering the roots in such a way as to minimize the proximity of each root to the roots of $p_{i-1}(z)$, we obtain the $\ell$ values of the $n$ initial guess functions $[r_{1,0}], \ldots [r_{n, 0}]$.

The second step of the algorithm maintains lists of estimates for the values of the derivatives of each of the $n$ phase functions on the $\ell-$ point quadrature. Initially, these lists only contain the initial estimates made above. As long as $\max_i \xi_i \geq \epsilon_0$, the following steps are repeated.
    
1. Use the previous iterations' estimates $[r_1], \ldots,[r_n]$ and the subprocedure described in the previous section to compute improved estimates $[r_{1, \text{new}}], \ldots, [r_{n, \text{new}}]$ of the solutions at the quadrature points.

2. Compute the $n$ convergence paramaters
$
\begin{equation}
    \xi_i = \frac{\|[\delta]\|_{L^2}}{\|[r_{i}]\|_{L^2}} \ \ \ \ i = 1, \ldots, n.
\end{equation}
$
and append the new values $[r_{i, \text{new}}]$ to the list of estimates.

The above iterative process provides accurate estimates $[r_1], \ldots, [r_n]$ of the values of the derivatives of the $n$ phase functions of $(1)$ at the constructed quadrature. The values of the  first $(n-1)$ derivatives of each solution can again be computed through repeated application of the spectral differentiation matrix $\mathscr{D}_\ell$ (see (4)). The first element of each of the vectors of values provides accurate estimates of the values of the solutions at the endpoint $a_0$. The algorithm is completed by applying the ODE solver of Appendix (1) of the aforementioned paper to the Riccati equation with the initial values $r_1(a_0), \ldots, r_n(a_0)$ to obtain the values of the global solutions on a peicewise Chebyshev structure.


## Generalization to $n^{\text{th}}$ degree

In this section, the derivation of the analag of $(3)$ for the arbitrary $n$ is shown.







