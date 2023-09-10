# n-dimensional-scalar


## Table of Contents

- $\texttt{scalar\_n.f90}$: Hello

- $\texttt{experiment\_n.nb}$

- $\texttt{partitions.nb}$: A mathematical file with contains for for computing all partitions of positive integers. This data is then compiled in an convenient way and used in a step of the algorith proper.

## Introduction

A popular and well-understood class of ordinary differential equations, commonly refered to as homogenous scalar equations, take on the form
$$y^{(n)}(t) + q_{n-1}(t)y^{(n-1)} + \cdots + q_{1}(t)y'(t) + q_0(t)y(t) = 0.$$


The cost of numerically representing solutions to this class of equations using standard methods increases with the magnitude of the coefficient functions $\{q_i\}_{i=0}^{n-1}$. However, the phase functions of these equations are solvable, independent of the magnitude of the coefficient functions. More information on the theory, algorithm description, implementation and testing can be found in this [paper](https://arxiv.org/abs/2308.03288), authored by James Bremer and myself. 

As you may notice, the implementation and description of the algorithm in the paper only considers for the case where $n=2, 3$, or $4$. Utitlizing the code for the lower dimensional cases, I independently generalized the algorithm to the case where the value of $n$ is unknown at compile time. Further, I conducted several experiments, both comparisons to the lower dimensions to assure accuracy, and to the general cases for large values of $n$. The results are conclusive to run-time independence of the magnitude of the coefficient functions, and quadratic proportionality to the value of $n$.

Included in this repository 


Fortran algorithm which computes a basis of the solution space of n^th order scalar equations, time-independent of the magnitude and complexity of the coefficient functions. 
