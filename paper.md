---
title: 'iceFEM: A FreeFem package for wave induced ice-shelf vibrations'
tags:
  - FreeFem
  - C++
  - Ice-Shelves
  - Fluid-Structure Interaction
authors:
  - name: Balaje Kalyanaraman
    orcid: 0000-0001-6181-1876
    affiliation: 1
  - name: Michael H. Meylan
    affiliation: 1
  - name: Bishnu P. Lamichhane
    affiliation: 1
  - name: Luke G. Bennetts
    affiliation: 2
affiliations:
  - name: School of Mathematical and Physical Sciences, University of Newcastle, Callaghan, NSW 2308, Australia
    index: 1
  - name: School of Mathematical Sciences, University of Adelaide, Adelaide, SA 5005, Australia
    index: 2

date: 18 December 2020
bibliography: paper.bib
---

# Summary

Vibrations of ice-shelves in response to ocean waves were first
investigated by [@holdsworth1978iceberg] who proposed that resonant
vibrations lead to icebergs calving from the shelf front. Since then
seismometric measurements on the Ross ice-shelf, the largest Antarctic
ice-shelf, confirmed the presence of this ocean wave-induced ice-shelf
vibration [@bromirski2015ross; @Massom2018]. The period of vibration
ranged from the long infragravity/tsunami waves to shorter, swell
waves. More recently, [@brunt2011antarctic] presented the first
observational evidence that a Northern Hemisphere tsunami triggered
calving on the Sulzberger ice-shelf. Mathematical models based on
linear wave theory have been proposed to study these ocean-wave
induced ice-shelf vibrations.

![Geometry and the Governing
equations. \label{fig:0}](Images/PaperImages/iceGeo.png)

Consider the geometry and the schematic shown in \autoref{fig:0}. The
motion of the ice-shelf (solid) is governed by the elastodynamic
equations which are coupled with the fluid motion, governed using the
linear potential flow theory. The fluid is modelled as a semi-infinite
rectangular region of uniform depth, whereas the sub-shelf cavity
region is assumed to be non-uniform. The time-domain problem is
converted to the frequency domain by applying a transformation
\begin{equation}
  \Phi(x,z,t) = \text{Re}\left\{ \phi(x,z)e^{-i \omega t}
  \right\},\quad
  \mathbf{u}(x,z,t) = \text{Re}\left\{ \eta(x,z)e^{-i \omega t} \right\}
\end{equation}
at a prescribed frequency $\omega$.

The semi-infinite region is truncated by constructing analytic
expressions and deriving a non-local boundary condition of the form
\begin{equation}
  \partial_x \phi = \mathbf{Q}\phi + \chi \quad \text{on} \quad \Gamma_f^{(4)}.
\end{equation}
The computational domain is then restricted to the sub-shelf cavity
region where the finite element method is used to solve the resulting
governing equations. The finite dimensional weak formulation of
the coupled problem is to find $(\phi_h,\mathbf{w}_h) \in V_h\times
W_h$ such that
\begin{eqnarray}
\left(\nabla\phi_h,\nabla\psi\right)_{\Omega_f}
    &=-i\omega\left\langle\mathbf{w}_h,\psi\right\rangle_{\Gamma_f^{(3)}} +
    \left\langle\mathbf{Q}\phi_h,\psi\right\rangle_{\Gamma_f^{(4)}} +
    \left\langle\chi,\psi\right\rangle_{\Gamma_f^{(4)}} \label{eq:1new}\\
    (\sigma(\mathbf{w}_h)\,:\,\epsilon(\mathbf{v}))_{\Omega_s}
    &=\rho_s\omega^2\left(\mathbf{w}_h,\mathbf{v}\right)_{\Omega_s} +
    \left\langle\mathbf{w}_h,\mathbf{v}\cdot\mathbf{n}\right\rangle_{\Gamma_s^{(1)}}
    -i\omega\left\langle\phi_h,\mathbf{v}\right\rangle_{\Gamma_f^{(3)}}
\end{eqnarray}
where $V_h$ and $W_h$ are appropriate finite element spaces.

## Numerical Method

The numerical method is based on the modal expansion technique. The
displacement and the velocity potential can be written as
\begin{equation}\label{eq:2new}
 \phi_h(x,z) = \phi_0(x,z) +
  \sum_{j=1}^{M}\lambda_j\phi_j(x,z),\quad \mathbf{w}_h(x,z) = \sum_{j=1}^{M}
  \lambda_j\eta_j(x,z)
\end{equation}
where $\lambda_j$'s being the unknown dofs. Substituting equation
\eqref{eq:2new} into the weak formulation of the linear elasticity
equations \eqref{eq:1new}, we obtain
\begin{align*}
  \sum_{j=1}^{M}\lambda_j\bigg[
  (\sigma(\eta_j)\,:\,\epsilon(\mathbf{w}_h))_{\Omega_s}
  &- \rho_s\omega^2(\eta_j,\mathbf{w}_h)_{\Omega_s}\\ -
  \left\langle\eta_j,\mathbf{w}_h\cdot\mathbf{n}\right\rangle_{\Gamma_s^{(1)}}
  &+ i\omega\left\langle\phi_j,\mathbf{w}_h\right\rangle_{\Gamma_f^{(3)}}
  \bigg] =
  - i\omega\left\langle\phi_0,\mathbf{w}_h\right\rangle_{\Gamma_f^{(3)}}
\end{align*}
which corresponds to the (reduced) linear system
\begin{equation}\label{eq:1}
[\mathbf{K}-\omega^2\mathbf{M}+i\omega\mathbf{B}]\{\lambda\} = \{\mathbf{f}\}
\end{equation}
where the entries of the matrices are given by
\begin{align*}
  \mathbf{K}_{jk} =
  \left(\sigma(\eta_j)\,:\,\epsilon(\eta_k)\right)_{\Omega_s},&\quad
  \mathbf{M}_{jk} = \rho_s(\eta_j,\eta_k)_{\Omega_s}, \\
  \mathbf{C}_{jk} =
  \left\langle\eta_j,\eta_k\cdot\mathbf{n}\right\rangle_{\Gamma_s^{(1)}},&\quad
  \mathbf{B}_{jk} = \left\langle\phi_j,\eta_k\right\rangle_{\Gamma_f^{(3)}}.
\end{align*}

The functions $\eta_j \in W_h$ are the in--vacuo vibration
modes of the ice--shelf which corresponds to solving the eigenvalue problem
\begin{equation*}
  \left(\sigma(\eta):\epsilon(\mathbf{v})\right)_{\Omega_s} =
  \rho_s\,\beta^2\left(\eta,\mathbf{v}\right)_{\Omega_s}
\end{equation*}
for all $\mathbf{v} \in W_h$. The diffraction potential $\phi_0 \in
V_h$ and the radiation potential $\phi_j \in V_h$ corresponding to the
vibration mode $\eta_j$ can be obtained by solving:
\begin{align*}
  (\nabla\phi_0,\nabla\psi)_{\Omega_f} &= \left\langle\mathbf{Q}\phi_0,\psi\right\rangle_{\Gamma_f^{(4)}}
  + \left\langle\chi,\psi\right\rangle_{\Gamma_f^{(4)}}\\
  (\nabla\phi_j,\nabla\psi)_{\Omega_f} &= \left\langle\mathbf{Q}\phi_j,\psi\right\rangle_{\Gamma_f^{(4)}}
  - i\omega \left\langle\eta_j,\psi\right\rangle_{\Gamma_f^{(3)}}
\end{align*}
The diffraction and radiation potentials can be computed in parallel
once the in-vacuo modes are obtained. The entries in the linear system
(\autoref{eq:1}) are analytic functions of the incident frequency
$\omega$. The dimension of the linear system is much
smaller than the finite element degrees of freedom and can be solved
efficiently. Using the analyticity of the resulting system, a large
number of frequency domain solutions can be obtained by interpolating
the linear system without having to solve the much larger finite
element problem on a finer frequency grid.

## Example

For the iceberg motion example in the program `iceberg.edp`, the
following bash script solves a small set of problems using the finite
element method for different values of incident wave frequencies:

``` bash
#!/bin/bash
# Generate the working directory
./genDir.sh 1_ICEBERG/;

# Solve the finite element-frequency domain problems.
for i in $(seq 15 2.5 80)
do
    mpirun -np 2 FreeFem++-mpi -v 0 iceberg.edp -N1 20 -N2 30 -Tr $i
    -L 3000 -H 2000 -h 200 -nev 8 -iter $(echo $i/2.5-5 | bc) >
    /dev/null;

    echo "Done $i";
done
```

![Figure showing (Top) the value of the reflection coefficients on a
coarse $\omega$-space (blue,+) and on a fine $\omega$-space
(red,solid) obtained after solving the interpolated system. (Middle)
The value of the reflection and transmission coefficients as a
function of the incident frequency. The energy conservation result
$|T|^2+|R|^2=1$ is also verified. (Bottom) Modal contribution,
$|\lambda_j|$ of the various in-vacuo modes as a function of
frequency. \label{fig:1}](Images/PaperImages/coeffs.png)

This computes the solution on a coarse $\omega$-space, ($\omega_c$)
containing 27 points. The program then writes a set of files containing the
real and imaginary part of the LHS and RHS of the linear system
(\ref{eq:1}) inside the working directory along with the
diffraction and radiation reflection coefficients. Once the reduced
system is obtained, the entries inside the files can be interpolated
on a finer $\omega$-space ($\omega_f$) to obtain the solution. Once
the solutions for $\lambda$ on the finer grid is obtained,
quantities like the reflection coefficients can be computed using the
new solution and the diffraction and radiation reflection
coefficients (\autoref{fig:1}). Several `MATLAB` routines are
available in the `modules` folder to compute the coefficients and the
reflection coefficients in real and complex $\omega$ space. For
example, to perform interpolation on the real $\omega$-space and
obtain the solution $\lambda_j$, the function `interpolateFreq()`
could be used. The [PDF
manual](https://github.com/Balaje/iceFem/blob/ParIceFem/manual.pdf)
contains more details on the `MATLAB` interface along with tutorials
on the `macros` available within the package. More examples can be
found in the [`README.md`](https://github.com/Balaje/iceFem#icefem).

# Statement of need
FreeFem [@ffpp] is an open-source domain specific language to implement
finite element methods and is based on C++. FreeFem is an excellent choice
for studying ice-shelf vibrations problems due to its flexibility and
ease of implementation. Similar packages include deal.II [@dealII92],
FEniCS project [@LoggMardalEtAl2012a],
which are examples of automated differential equation solvers and more modern
packages like Gridap [@Badia2020], based on Julia. `iceFEM` is a
FreeFem package for simulating ice-shelf vibrations, heavily inspired
by the `ffddm` module available in FreeFem for implementing the domain
decomposition methods. Numerical Methods have been proposed to study
the vibrations of these ice-shelves, predominantly based on the
thickness averaged thin-plate model for the ice-shelf and the depth
averaged shallow-water models for the fluid flow in the sub--shelf
cavity region [@Sergienko2013; @Meylan2017]. More recently, finite
element methods have been proposed to model non-uniform shelf/cavity
regions [@Ilyas2018; @kalyanaraman2020coupled]. The method also
involves implementing a general non-local boundary condition based on
analytic expressions to handle semi-infinite regions. The numerical
methods can be extended to complex valued inputs for the incident
frequency and hence require flexible and versatile finite element
solvers while being easy to implement. `iceFEM`
implements the method to solve the vibration problem for complex
valued incident frequencies which are useful to study resonances in
the complex plane [@kalyanaraman2020coupled], shown in \autoref{fig:4}.

![Reflection and Transmission Coefficient as a function of complex
frequencies. The white dots indicate the poles which are points where
there is an abrupt change in color. The poles correspond to the
resonance frequencies of the shelf/cavity system. The color denotes
the phase angle and the brightness, the magnitude of the complex number
[@wegert2012visual]. \label{fig:4}](Images/PaperImages/Figure_1.png)

While *a priori* knowledge of finite elements are useful to write
scripts using `iceFEM`, some macros which yield the essential
stiffness matrix and load vector are available. More work is being
done to improve the user-friendliness of `iceFEM` in future
releases. Real-life examples using the BEDMAP2 dataset can also be
imported and solved using `iceFEM`. The finite element algorithms in
the `iceFEM` package was validated using the thin-plate solutions
obtained using the eigenfunction matching methods in
[@KALYANARAMAN2019].

## References
