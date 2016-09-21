# L-BFGS-Cutoff-sd

## Introduction 
This is the final version for the cutoff implementation of L-BFGS algorithm for minimizing energy on s^d sphere. It uses angle coordinate
and is being developed and tested right now.

The program is compiled using cmake. The code was mainly completed with Alex Vlasiuk in 2016 fall.

## Algorithm
The code is implemented mainly based on the L-BFGS algorithm implemented by PatWie from https://github.com/PatWie/CppNumericalSolvers.git.
His implementation of the algorithm provides general line search and the l-bfgs solver class. We modify the solver and the line search 
based on our needs. Both classes are revised based on **Numerical Optimization** by Norcedal and Wright. In addition the solver class is 
 revised according to the stopping condition (convergence) in our case. We also added output functionality in the solver class.


## Compile
To compile, just open a terminal and run   
    `cmake . && make`  
A executable binary file called L-BFGS will be genearated.

## Run
To run the program, you need to provide a input file and a control file. The format of the input file should be .txt and it should contains
the initial coordinates of the point configuration on 3d sphere. The name of it is not required. For the control file, the name and format
are strictly required. The name has to be control.inp and an example of the format can be found in the folder. It is recommaned that you 
download the control.inp and modify it based on your input.

The exact command to run the program is the same as any binary executable file:   
    `./L-BFGS`   
By default, the output file generated will be named as output.txt. This implementation also allows to generate output files containing the coordinates. The number of files generated will be exactly the partitionnum. To obtain the outputs, create a folder name output under the L-BFGS-Cutoff-s2 directory. 

Note: the outputs will be overwritten after each execution.

## Input file or random configuration
In the control.inp file, there is an option of infile. 0 indicates no input file and 1 otherwise. When 0 is put, the inputfile should be exactly the inputfile name. When 1 is selected, the program will just ignore the inputfile line and generate a random configuration on 3d sphere.
