#!/bin/bash

dirName=$1


mkdir -p $dirName
cd $dirName
mkdir -p 2_Deformation 2_Modes 2_ModesMatrix 2_Potential 2_RefCoeff 2_Stresses
cd 2_ModesMatrix
mkdir -p Interpolated_F Interpolated_H Interpolated_L
cd ..
cd 2_RefCoeff
mkdir -p Interpolated_R RefCoeff_Dif RefCoeff_Rad
