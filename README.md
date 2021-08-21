## simrx2 in a nutshell

**simrx** is an HPC (High Performance Computing) X-ray projection simulator based on **PENELOPE**, a particle simulator developed at the University of Barcelona. It was originally implemented in Fortran, and I extended it with a detector to capture images, i.e., image formation, and a C++ module using MPI in order to execute it in a cluster of CPUs.

The first version was intended to produce only one projection, just one image. Now, it is able to produce multiple images at different angles so that they can be used later to feed them into an algorithm for volumetric reconstruction of the sample, also know as tomography or CT-scans.

Below there is an example of an experiment with 4 projections of a cube filled with water, and a silicon cylinder inside. The simulation took 10<sup>8</sup> X-ray photons (monochromatic photons) with 50 KeV of energy. The geometry of the X-ray beam is conical.

![](https://github.com/marselan/simrx2/blob/refactor/result_first_refactor/png/image1.png)
![](https://github.com/marselan/simrx2/blob/refactor/result_first_refactor/png/image2.png)
![](https://github.com/marselan/simrx2/blob/refactor/result_first_refactor/png/image7.png)
![](https://github.com/marselan/simrx2/blob/refactor/result_first_refactor/png/image8.png)

The values for the attenuation coefficient measured in the images vs theoretical:

||Measured|Theoretical|
|-|---|---|
|Water|0.2214 cm<sup>-1</sup>|[0.2269 cm<sup>-1</sup>](https://physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/water.html)|
|Silicon|0.9903 cm<sup>-1</sup>|[1.0217 cm<sup>-1</sup>](https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z14.html)|

[![simrx2 sample simulation](https://github.com/marselan/simrx2/blob/master/misc/cube_video.jpg)](https://youtu.be/ZdYvHYo7Ff4)

## Installation and configuration instructions

Please refer to the [wiki here.](https://github.com/marselan/simrx2/wiki)

