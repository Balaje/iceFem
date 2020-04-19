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

There are three different ice--shelf examples that can be solved.

1. `iceshelf_submerged_moving.edp` assumes that the ice--shelf is a
1D thin-plate and the vibrations are modelled using the
Euler-Bernoulli beam theory. The vibration of the ice-shelf and the
velocity potential in the cavity region for an
incident wave-forcing of 200 s is shown below.

<p style='text-align: center;'>
<img width="525" height="388" src="./Images/femEB.png" border="0">
</p>


2. `iceSpline.edp` uses the 2D linear elasticity equations under plane strain
conditions for the ice-shelf. Figure on 

| ![Image 1](./Images/femLEmesh.png) | ![Image 1](./Images/femLE.png) |
| ---------------------------------- | ------------------------------ |
