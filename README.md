# Tarea_4_ANPI
Luis Fernando Murillo 
Jorge Aguero
*****************************************************************************************************************************
**********************************************  Tarea 4 **************************************************
*****************************************************************************************************************************

README

In this task we implemented three methods for matrix decomposition.
1.LU Crout
2.LU Doolitlle
3 QR
Which are used to solve the equations Ax = b

CONTACT

If you have problems or comments with this program you
can contact please contact us.

This project can also be found at GitHub in:
https://github.com/luisf1997cr/Tarea_4_ANPI

-----------------------------------------------------------------------------------------------------------------------
------------------------------------------- RUN INSTRUCTIONS ----------------------------------------------
-----------------------------------------------------------------------------------------------------------

Unzip the project, open a terminal and change your working directory to the unzipped folder

Create a directory build:

> mkdir build;

Go into that directory

> cd build;

You can choose to build a release version with:

> cmake ../ -DCMAKE_BUILD_TYPE=Release

or a debug version with

> cmake ../ -DCMAKE_BUILD_TYPE=Debug

And build everything with

> make

The executables will be stored at build/bin.

To execute the program go to the /build/bin directory

> cd bin

If you can, running all test of task you can use

> ./tester 

Or you can see all unit tests available using

> ./tester --list_content

Examples of unit test

> ./tester -t QR

> ./tester -t Matrix

> ./tester -t LU/Crout

You can run the benchmarks used to select the LU algorithm by running

>./benchmarks -t LU

________________________________________________________________________________________________________________________
____________________________________________________ Dependencies _______________________________________
_________________________________________________________________________________________________________________________

You need CMAKE and Boost to build the program

> sudo apt-get install libboost-all-dev
> sudo apt-get -y install cmake


