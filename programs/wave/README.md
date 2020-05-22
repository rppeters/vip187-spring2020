#Gaussian Pulse Wave Propogation

A simple 2-d simulation of a gaussian pulse in the water modeled in a gif file

##Dependencies
Sequential:
```c 
    #include <stdlib.h>
    #include <stdio.h>
    #include <math.h>
    #include <gd.h>
    #include <assert.h>
    #include <time.h>
```

Parallel: same as sequential but also the ability to run Cuda code.

##Compiling
Using the Makefile provided

Sequential:
    For normal build:
    ```
    make wave
    ```

    For debugging build:
    ```
    make wave_debug
    ```

    Running:
    ```
    make wave_small
    make wave_large
    ```

Modifications will most likely need to be made for parallel/Cuda.

#Gif
Excercise caution opening larger gif files 