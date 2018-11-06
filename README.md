# G17-ASU-MathLab

## Introduction ##

This is a college project which aims to make a C++ console app that works as MATLAB does.

It's a mathematical software similar to MATLAB, it mainly performs most of operations.
It's develepod in two phases.


## Specifications ##

- It works on linux.
- It has the same interface as MATLAB.
- It has a fast response time for large inputs.
- Provide a result identical to MATLAB.

## Phase 1

### Features ###

- Creating matrices of any size.
- Support mathematical core operations (addition, subtraction, multiplication, transpose and division).

## Phase 2

### Features ###

- Support mathematical functions (Trignometric, logarithmic, roots and power).
- Support element wise operations.
- Support matrix parsing. 
- Support error handling.

## Input example

A = [1 2 3; 4 5.7 6; 3 4 2];

B = [1.5 4.1 5.4; 3.1 4.2 1.2; 3.2 4.3 2.2];

C = A + B

D = 5.5 + 12 * sin(0.4) + 2.2^4;

E = [1.2 2.3 D;[1.3 2.4;4.6 1.3],[3.2;7.8]];

F = (C^3 * sqrt(1./D))^(0.1)

G = A * B

H = A’


## How to use:

- First option: The user writes the commands in the console and gets the results directly.
- Second option: The user makes a file containing the commands and the program reads it and produces the results.   
