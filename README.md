# iceFEM++

To use iceFEM++, you must have FreeFem++ installed. Currently, the
program is written to model vibrations of ice-shelves (1D, 2D and 3D), an
example of a fluid-structure interaction problem. The schematic and
the governing equations are shown in the Figure below.

<p style='text-align: center;'>
<img width="760" height="345" src="./Images/iceGeo.png" border="0">
</p>


In future, the
program will be extended to solve more problems in fluid-structure
interaction.

There are different ice-shelf examples that can be solved.

1. `iceshelf_submerged_moving.edp` assumes that the ice--shelf is a
1D thin-plate and the vibrations are modelled using the
Euler-Bernoulli beam theory. The vibration of the ice-shelf and the
velocity potential in the cavity region for an
incident wave-forcing of 200 s (in the frequency domain) is shown below. Figure on the left shows the solution on a non uniform cavity and on the right, the solution on a uniform cavity. In this case, the problem is solved using the modal expansion technique, which is used for solving hydro-elasticity problems of large container ships.

| ![Non-Uniform cavity](./Images/femEB1.png) | ![Uniform Cavity](./Images/femEB2.png) |
| ---------------------------------- | ------------------------------ |


2. `iceSpline.edp` uses the 2D linear elasticity equations under plane strain
conditions for the ice-shelf. Figure on the left shows the finite element meshes used for the cavity and the ice-shelves (both non-uniform).
The governing equations are solved using the combined approach of modal expansion and the finitie element method.

| ![Meshes](./Images/femLEmesh.png) | ![Solution](./Images/femLE.png) |
| ---------------------------------- | ------------------------------ |

The solution to the linear elasticity problem agrees with the thin-plate solution when the ice-shelf is uniform and thin!

| ![Meshes](./Images/femLEvsEB3.png) | ![Solution](./Images/femLEvsEB4.png) |
| ---------------------------------- | ------------------------------ |


MATLAB routines
```matlab
leSolu.m
femEBvsFull.m
thinVsFull.m
```
assist with the visualization of the solution. A fancier example is to find solution to the linear elasticity problem in the frequency domain for complex incident frequencies. Run `refCoeff_cplx.m` to obtain the Figure below, which shows the reflection coefficient as a function of the incident frequency in the complex plane. The white dots indicate the resonance frequencies.

<p style='text-align: center;'>
<img width="530" height="345" src="/Images/resonance3.png" border="0">
</p>

*Do this in your personal laptop at your own risk! The code takes a long time to run on a personal laptop. My laptop almost died; but, fortunately she hung on.*

More examples coming soon.
