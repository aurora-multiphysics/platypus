# Darcy Problem

## Summary

Solves a simple 2D/3D mixed Darcy problem corresponding to the saddle point system.
We discretize with Raviart-Thomas finite elements
(velocity $ \mathbf{u} $) and piecewise discontinuous polynomials (pressure $ p $).
This example is based on
[MFEM Example 5](https://mfem.org/examples/).

## Description

This problem solves a mixed Darcy problem with the following form:

\begin{align}
\mathbf{K} \mathbf{u} + \nabla p &= \mathbf{f} \\
\nabla \cdot \mathbf{u} &= g
\end{align}

with natural boundary condition $ -p = \text{"given pressure"} $.
Here we use a given exact solution $ (\mathbf{u}, p) $ and compute the
corresponding right-hand side $ (\mathbf{f}, g) $.
We discretize with Raviart-Thomas finite elements (velocity $ \mathbf{u} $) and
piecewise discontinuous polynomials (pressure $ p $).

The example demonstrates the use of transpose in the input file.

In this example, we solve this using the weak form

\begin{equation}
(\mathbf{K} \mathbf{u}, \mathbf{v})_\Omega + (\nabla p, \mathbf{v})_\Omega - (\nabla \cdot \mathbf{u}, q)_\Omega
= (\mathbf{f}, \mathbf{v})_\Omega
\quad \forall \mathbf{v} \in V, \forall q \in W
\end{equation}

where

\begin{equation}
\begin{split}
\mathbf{u} \in H(\mathrm{div})(\Omega) &: \mathbf{u} \cdot \hat n = \mathbf{g} \quad \text{on} \quad \partial \Omega \\
p \in L^2(\Omega) &: p = p_0 \quad \text{on} \quad \Gamma_\mathrm{D}
\end{split}
\end{equation}

### Block Darcy Operator

Assemble the finite element matrices for the Darcy operator:

\begin{equation}
D = \begin{bmatrix}
M & B^T \\
B & 0
\end{bmatrix}
\end{equation}

where:

\begin{equation}
M = \int_\Omega k \mathbf{u}_h \cdot \mathbf{v}_h \, d\Omega \quad \mathbf{u}_h, \mathbf{v}_h \in R_h
\end{equation}

\begin{equation}
B = -\int_\Omega \nabla \cdot \mathbf{u}_h \, q_h \, d\Omega \quad \mathbf{u}_h \in R_h, \quad q_h \in W_h
\end{equation}

The transpose of the matrix $ B $ is used in the weak form to ensure the correct assembly of the system.

## Example File

!listing kernels/darcy.i
