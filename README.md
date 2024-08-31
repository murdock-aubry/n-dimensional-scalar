# Table of Contents

- ``scalar-n.f90`` Fortran file which contains the main routines used in the algorithm for the general $n$ dimensional case.

- ``experiment-i.f90``: Testing the algorithm against the lower dimensional code, assuring accuracy.

- ``odesolve-n.tex``: A brute force ODE solver for $n$ dimensional equations used as a step in this algorithm.

- ``partitions.nb``: A mathematical file with contains for for computing all partitions of positive integers. This data is then compiled in an convenient way and used in a step of the algorith proper.

Please note that the files ``chebyshev.f90``, ``chebpw.f90``, ``legendre.f90``, ``linalg0.f90`` and ``utils.f90`` are written by my supervisor James Bremer, and are contain pivotal routines in the functionality of my algorithm.


# Introduction

A popular and well-understood class of ordinary differential equations, commonly refered to as homogenous scalar equations, take on the form
$$y^{(n)}(t) + q_{n-1}(t)y^{(n-1)} + \cdots + q_{1}(t)y'(t) + q_0(t)y(t) = 0. \tag{1}$$

The cost of numerically representing solutions to this class of equations using standard methods increases with the magnitude of the coefficient functions $\lbrace q_i\rbrace_{i=0}^{n-1}$. However, the phase functions of these equations are solvable, independent of the magnitude of the coefficient functions. More information on the theory, algorithm description, implementation and testing can be found in this [paper](https://arxiv.org/abs/2311.08578), authored by James Bremer and myself. 

As you may notice, the implementation and description of the algorithm in the paper only considers for the case where $n=2, 3$, or $4$. Utitlizing the code for the lower dimensional cases, I independently generalized the algorithm to the case where the value of $n$ is unknown at compile time. Further, I conducted several experiments, both comparisons to the lower dimensions to assure accuracy, and to the general cases for large values of $n$. The results are conclusive to run-time independence of the magnitude of the coefficient functions, and quadratic proportionality to the value of $n$.


## Generalization to $n^{\text{th}}$ degree

In this section, the derivation of the analag of $(3)$ for the arbitrary $n$ is shown.

The $n^{\text{th}}$ derivative chain rule of the composition of functions $f\circ g$ is provided by the following formula:

$$
    \frac{d^n}{dx^n}f(g(x)) = \sum_S A_{m_1, \ldots, m_n} f^{(m_1 + \cdots + m_n)}(g(x)) \prod_{j=1}^n(g^{(j)}(x))^{m_j}, \tag{12}
$$

where the sum is taken over the collection $S$ of all sequences $(m_1, \ldots, m_n)$ such that

$$
    m_1 + 2m_2 + \ldots + nm_n = \sum_{i=1}^n im_i = n \tag{13}
$$

and

$$
    A_{m_1, \ldots, m_n} = \frac{n!}{m_1!(1!)^{m_1}\cdot m_2!(2!)^{m_2} \cdots m_n!(n!)^{m_n}}\tag{14}
$$

Substituting $f(t) = \exp(t)$ and $g(t) = \int r(t)dt$, we obtain

$$
    \frac{d^n}{dx^n}\exp\left(\int r(t)dt \right) = \exp\left(\int r(t)dt \right) \sum_S A_{m_1, \ldots, m_n} \prod_{j=1}^n(r^{(j-1)}(x))^{m_j}\tag{15}
$$

Therefore, the $n^{\text{th}}$ order Riccati equation for the $n^{\text{th}}$ order equation

$$
    y^{(n)}(t) + k^nq(t)y(t) = 0 \tag{16}
$$

is given by

$$
    \sum_S A_{m_1, \ldots, m_n} \prod_{j=1}^n(r^{(j-1)})^{m_j} + k^nq = 0 \tag{17}
$$    

Deploying Newton's method via perturbation by a function of small magnitude $\delta$;

$$
    r_{n+1}(t) = r_n(t) + \delta(t) \tag{18}
$$

we wish to find (\ref{chain n}) in the case where $g(t) = \int (r(t) + \delta(t))dt$. A direct substitution of (18) into (12) yields

$$
    \frac{d^n}{dx^n}\exp\left(\int (r(t) + \delta(t))dt \right) = \exp\left(\int (r(t) + \delta(t))dt \right) \sum_S A_{m_1, \ldots, m_n} \prod_{j=1}^n(r^{(j-1)} + \delta^{(j-1)})^{m_j}\tag{19}
$$

and so the linearized Riccati equation becomes

$$
    \sum_S A_{m_1, \ldots, m_n} \prod_{j=1}^n(r^{(j-1)} + \delta^{(j-1)})^{m_j} + k^nq = 0 \tag{20}
$$

The binomial series provides

$$
    (r^{(j-1)} + \delta^{(j-1)})^{m_j} = \sum_{k=0}^{m_j} \binom{m_j}k(r^{(j-1)})^{m_j - k}(\delta^{(j-1)})^k. \tag{21}
$$

For all $i_1, i_2, j_1, j_2 \geq 1$, we make the assumption that $(\delta^{(i_1)})^{j_1}(\delta^{(i_2)})^{j_2} \approx 0$ and thus we only consider the terms $k = 0, 1$ in (21), yielding the following

$$
    (r^{(j-1)} + \delta^{(j-1)})^{m_j} \approx (r^{(j-1)})^{m_j} + m_j(r^{(j-1)})^{m_j - 1}\delta^{(j-1)} \tag{22}
$$

Further, under the assumption $b_ib_j \approx 0$ for all $i, j$, we have

$$
    \prod_{i=1}^n(a_i + b_i) \approx \prod_{i=1}^n a_i + \sum_j b_j\left(\prod_{i=1, i\neq j}^n a_i \right) \tag{23}
$$

Therefore, utilizing (22) and (23), we obtain

$$
    \prod_{j=1}^n(r^{(j-1)} + \delta^{(j-1)})^{m_j} \approx \prod_{j=1}^n (r^{(j-1)})^{m_j} + \sum_{j=1}^nm_j(r^{(j-1)})^{m_j - 1}\delta^{(j-1)}\left[\prod_{i=1, i\neq j}^n(r^{(i-1)})^{m_i} \right] \tag{24}
$$

Plugging this into the linearized Riccati equation (9), we obtain

$$
        \sum_S A_{m_1, \ldots, m_n} \left(\prod_{j=1}^n (r^{(j-1)})^{m_j} + \sum_{j=1}^nm_j(r^{(j-1)})^{m_j - 1}\delta^{(j-1)}\left[\prod_{i=1, i\neq j}^n(r^{(i-1)})^{m_i} \right]\right) + k^nq = 0 \tag{25}
$$

$$
    \implies \boxed{\sum_{j=1}^n\left(\sum_Sm_j A_{m_1, \ldots, m_n}(r^{(j-1)})^{m_j - 1}\left[\prod_{i=1, i\neq j}^n(r^{(i-1)})^{m_i}\right]\right)\delta^{(j-1)} = - \sum_S A_{m_1, \ldots, m_n}\prod_{j=1}^n(r^{(j-1)})^{m_j} - k^nq} \tag{26}
$$

The above equation is discretizable, just as in the lower dimensional cases. The discretization code can be found in the interative steps of the algorithm contained in the subroutine ``scalar-levin-n`` in the file  ``scalar-n.f90
``.







