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
  - name: University of Newcastle, Australia
    index: 1
  - name: University of Adelaide, Australia
    index: 2

date: 18 December 2020
bibliography: paper.bib
---

# Summary

Vibrations of ice-shelves in response to ocean waves were first investigated by [@holdsworth1978iceberg] who proposed that resonant vibrations lead to icebergs calving from the shelf front. Since then seismometric measurements on the Ross ice-shelf, the largest Antarctic ice-shelf, confirmed the presence of this ocean wave-induced ice-shelf vibration [@bromirski2015ross; @Massom2018]. The period of vibration ranged from the long infragravity/tsunami waves to shorter, swell waves. More recently, [@brunt2011antarctic] presented the first observational evidence that a Northern Hemisphere tsunami triggered calving on the Sulzberger ice-shelf. Mathematical models based on linear wave theory have been proposed to study these ocean-wave induced ice-shelf vibrations.

The numerical method is based on finite element method and the modal expansion technique. This results in a system of equations
$$[\mathbf{K}-\omega^2\mathbf{M}+i\omega\mathbf{B}]\{\lambda\} = \{\mathbf{f}\}$$
whose entries are analytic functions of the incident frequency $\omega$ which has a lot of advantages over domain decomposition methods. The dimension of the above linear system is much smaller than the finite element degrees of freedom and can be solved efficiently. Further, using the analyticity of the resulting system, a large number of frequency domain solutions can be obtained by interpolating the linear system without having to solve the much larger finite element problem on a finer frequency grid as shown in Figure \autoref{fig:1}. Further use of the package can be found in the ``README`` file located on the GitHub page.

![Figure showing (Top) the value of the reflection coefficients on a coarse $\omega$-space (blue,+) and on a fine $\omega$-space (red,solid) obtained after solving the interpolated system. (Middle) The value of the reflection and transmission coefficients as a function of the incident frequency. (Bottom) Modal contribution, $\lambda_j$ of the various in-vacuo modes as a function of frequency. \label{fig:1}](Images/PaperImages/coeffs.png)


# Statement of need
FreeFem [@ffpp] is an open-source programming language to implement finite element methods and is based on C++. FreeFem is an ideal choice for studying ice-shelf vibrations problems due to its flexibility and ease of implementation. `iceFEM` is a FreeFem package for simulating ice-shelf vibrations, heavily inspired by the `ffddm` module available in FreeFem. Numerical Methods have been proposed to study the vibrations of these ice-shelves, predominantly based on the thickness averaged thin-plate model for the ice-shelf and the depth averaged shallow-water models for the fluid flow in the sub--shelf cavity region [@Sergienko2013; @Meylan2017]. More recently, finite element methods have been proposed to model non-uniform shelf/cavity regions [@Ilyas2018; @kalyanaraman2020coupled]. The method also involves implementing a general non-local boundary condition based on analytic expressions to handle semi-infinite regions. `iceFEM` also implements the method to solve the vibration problem for complex valued incident frequencies which are useful to study resonances in the complex plane [@kalyanaraman2020coupled]. Real-life examples using the BEDMAP2 dataset can also be imported and solved using `iceFEM`. The `iceFEM` package was validated using the thin-plate solutions obtained using the eigenfunction matching methods in [@KALYANARAMAN2019].

## References
