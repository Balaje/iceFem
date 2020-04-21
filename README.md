# iceFEM++

To use iceFEM++, you must have FreeFem++ installed. Currently, the
program is written to model vibrations of ice-shelves (1D, 2D and 3D), an
example of a fluid-structure interaction problem. The schematic and
the governing equations are shown in the Figure below.

<p style='text-align: center;'>
<img width="760" height="345" src="./Images/iceGeo.png" border="0">
</p>

To setup and use the package, cd into the package directory, open up terminal and type

```shell
chmod +x genDir.sh
./genDir.sh 1_Forced
./genDir.sh 1_ThinPlate
./genDir.sh 1_DirectFEM
export FF_INCLUDEPATH="./include" #Add the include/ folder into the FreeFem++ environment variable. Needs to be run for the first time.
FreeFem++ -ne iceSpline.edp
```

`genDir.sh` is a shell script to generate the directories used by the FreeFem++ routine to write the appropriate files. The last command computes the solution for a default uniform ice-shelf and cavity with

* Length L=20 km
* Thickness h=200 m
* Cavity Depth H=800 m
* Incident Wave period T=200 s

To specify more inputs, one can use the following command

```shell
FreeFem++ -ne -v 0 iceSpline.edp -L [LENGTH] -H [DEPTH OF CAVITY] -h [THICKNESS OF ICE]
                                 -N [MESH PARAM]
                                 -Tr [REAL(PERIOD)] -Ti [IMAGE(period)]
                                 -iter [SOL. INDEX]
                                 -isUniIce [ON/OFF UNIFORM/NON UNIFORM ICE]
                                 -isUniCav [ON/OFF UNIFORM/NON UNIFORM CAVITY]
                                 -isForced [ON/OFF SHELF-FRONT FORCES]
```

**Example: Run**

```shell
FreeFem++ -ne -v 0 iceSpline.edp -L 10000 -H 800 -h 200 -N 4 -Tr 100 -Ti 0 -iter 0 -isUniIce 1 -isUniCav 1 -isForced 0    

FreeFem++ -ne -v 0 iceSpline.edp -L 15000 -H 800 -h 200 -N 4 -Tr 200 -Ti 0 -iter 0 -isUniIce 1 -isUniCav 0 -isForced 0
```

MATLAB can be used to visualize the solution generated using FreeFem++. A FreeFem++ to MATLAB converter is available online. If using this package, the required files are located in the `modules/` directory. To generate good PDF plots, it is recommended to use `export_fig` MATLAB package. To visualize the solution obtained by FreeFem++, we use this MATLAB code:

```Matlab
% From the FreeFem to MATLAB converter.
[pts1,seg1,tri1] = importfilemesh('1_Forced/2_Deformation/movedIce0.msh');
[pts2,seg2,tri2] = importfilemesh('1_Forced/2_Deformation/cavityMesh.msh');
PHILE = importfiledata('1_Forced/2_Potential/potential0.bb');
% Plot using PDEPLOT
figure(1)
subplot(2,1,1);
pdeplot(pts1,seg1,tri1)
axis equal tight
grid on
subplot(2,1,2)
pdeplot(pts2,seg2,tri2,'XYData',real(PHILE)','colormap','jet');
```
The following plots are generated.

  | ![Example 1](./Images/eg1.png) | ![Example 2](./Images/eg2.png) |
  | ---------------------------------- | ------------------------------ |

The Reflection Coefficients are tabulated below.

  | Reflection Coefficient, R | abs(R)  |
  | ---------------------------------- | ------------------------------ |
  | (0.8507259058288464,0.525609582438974) | 0.9999999999999919 |
  | (-0.3166231272563836,0.9485514194213012) | 0.9999999999998886 |

There are four demo MATLAB scripts that are written to generate the plots shown below.
```Matlab
thinVsFull.m
femEBvsFull.m
leSolu.m
refCoeff_cplx.m
```

**Note: Before using the MATLAB scripts, it is essential to run the `genDir.sh` script to generate the following directories.**

```shell
./genDir.sh 1_Forced
./genDir.sh 1_ThinPlate
./genDir.sh 1_DirectFEM
```

**NOTE:** `refCoeff_cplx.m` **computes the Analytic extension of the solution for complex frequencies. The routine computes the solution for a large number of incident frequencies and hence takes a long time to run on a personal laptop.**

There are different ice-shelf examples that can be solved.

1. `iceshelf_submerged_moving.edp` assumes that the ice--shelf is a
  1D thin-plate and the vibrations are modelled using the
  Euler-Bernoulli beam theory. The vibration of the ice-shelf and the
  velocity potential in the cavity region for an
  incident wave-forcing of 200 s (in the frequency domain) is shown below. Figure on the left shows the solution on a non uniform cavity and on the right, the solution on a uniform cavity. In this case, the problem is solved using the modal expansion technique, which is used for solving hydro-elasticity problems of large container ships. To use this module, generate the appropriate directory by running

  ```shell
  mkdir -v 1_ThinPlate
  ```

  | ![Non-Uniform cavity](./Images/femEB1.png) | ![Uniform Cavity](./Images/femEB2.png) |
  | ---------------------------------- | ------------------------------ |


2. `iceSpline.edp` uses the 2D linear elasticity equations under plane strain
  conditions for the ice-shelf. Figure on the left shows the finite element meshes used for the cavity and the ice-shelves (both non-uniform).
  The governing equations are solved using the combined approach of modal expansion and the finite element method. Run

  ```shell
  ./genDir.sh 1_Forced
  ```

  to generate the solution directories before running the solver.

  | ![Meshes](./Images/femLEmesh.png) | ![Solution](./Images/femLE.png) |
  | ---------------------------------- | ------------------------------ |  

3. `icefem.edp` is used to solve the coupled problem on uniform meshes. Since the mesh points match exactly at the shelf/cavity interface, the coupled problem is solved directly using the finite element method without modal expansions.  To generate the appropriate directories, run

```shell
./genDir.sh 1_DirectFEM
```

The solution to the linear elasticity problem agrees with the thin-plate solution when the ice-shelf is uniform and thin!

| ![Meshes](./Images/femLEvsEB3.png) | ![Solution](./Images/femLEvsEB4.png) |
| ---------------------------------- | ------------------------------ |

Will support 3D models in the future. Some examples currently in the works are shown below.

 | ![3D_1](./Images/Displacement.png) | ![3D_2](./Images/Velocity.png) |
 | ---------------------------------- | ------------------------------ |

More coming soon.

Contact: Balaje K,  Email: [balaje6@gmail.com](mailto:balaje6@gmail.com)
