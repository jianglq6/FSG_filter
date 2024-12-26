#******************************************************************************#
#*                                  FSG_filter                                *#
#*                           =========================                        *#
#*                                                                            *#
#*  Fully Staggered Grid Finite-Difference 2D Elastic Aniso Wave Modeling     *#
#*                                                                            *#
#*                                                                            *#
#*     A signed non-commercial agreement is required to use this program.     *#
#*             Free for non-commercial academic research ONLY.                *#
#*       This program is distributed WITHOUT ANY WARRANTY whatsoever.         *#
#*       Do not redistribute this program without written permission.         *#
#*                                                                            *#
#*   Written by:                                                              *#
#*     Luqian Jiang, Jun Tian, and Wei Zhang                                  *#
#******************************************************************************#

1 Introduction

  This program simulates elastic aniso wave propagation in 2D complex media by using  
  fully staggered grid finite-difference method. 

2 Direcotry layout

    Makefile 
       makefile for this package

    example/
       run script to demonstrate the usage

    forward/
       source codes

    mfile/
       matlab scriptes to show the result

    lib/
       lib source codes

    media/
       codes for discreting model, written by


3 Installation

    3.1 set compiler env, install NetCDF

    3.2 Compiling:
        
        edit Makefile

        make

4 Usage
  
  1) cd example/
  2) modify the parameters in run_*.sh
  3) run run_*.sh
  4) use mfiles/.m to show results

# vim:ft=conf:ts=4:sw=4:nu:et:ai:
