#!/bin/bash

dirName=$1

mkdir -v $dirName
cd $dirName
mkdir -v 2_Deformation 2_Modes 2_ModesMatrix 2_Potential 2_RefCoeff 2_Stresses
cd 2_ModesMatrix
mkdir -v Interpolated_F Interpolated_H Interpolated_L
cd ..
cd 2_RefCoeff
mkdir -v Interpolated_R RefCoeff_Dif RefCoeff_Rad

mkdir -v Meshes/ Meshes/BEDMAP2
