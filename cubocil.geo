0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   1)
INDICES=( 0, 0, 0, 0, 0)
    AXX=( 1.000000000000000E-00,   0)
     A0=(-1.000000000000000E+02,   0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   2)
INDICES=( 0, 0, 0, 0, 0)
    AYY=( 1.000000000000000E-00,   0)
     A0=(-1.000000000000000E+02,   0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   3)
INDICES=( 0, 0, 1, 0,-1)
Z-SHIFT=( 2.000000000000000E+01,   0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   4)
INDICES=( 0, 1, 1, 0,-1)
Y-SHIFT=( 1.500000000000000E+00,   0)
Z-SHIFT=( 1.500000000000000E+00,   0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   5)
INDICES=( 0, 0, 0, 0, 0)
    AXX=( 1.000000000000000E-00,   0)
     A0=(-1.000000000000000E+00,   0)
1111111111111111111111111111111111111111111111111111111111111111
X-SHIFT=( 1.500000000000000E+00,   0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   6)
INDICES=( 0, 0, 0, 0, 0)
    AXX=( 1.000000000000000E-00,   0)
     A0=(-9.000000000000000E-00,   0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   7)
INDICES=( 0, 0, 0, 0, 0)
    AYY=( 1.000000000000000E-00,   0)
     A0=(-9.000000000000000E-00,   0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   8)
INDICES=( 0, 0, 0, 0, 0)
    AZZ=( 1.000000000000000E-00,   0)
     A0=(-9.000000000000000E-00,   0)
0000000000000000000000000000000000000000000000000000000000000000
BODY    (   1)  DETECTOR
MATERIAL(   1)
SURFACE (   1), SIDE POINTER=(-1)
SURFACE (   2), SIDE POINTER=(-1)
SURFACE (   3), SIDE POINTER=(-1)
0000000000000000000000000000000000000000000000000000000000000000
BODY    (   2)  CYLINDER
MATERIAL(   1)
SURFACE (   4), SIDE POINTER=(-1)
SURFACE (   5), SIDE POINTER=(-1)
0000000000000000000000000000000000000000000000000000000000000000
MODULE  (   3)  CUBE
MATERIAL(   2)
SURFACE (   6), SIDE POINTER=(-1)
SURFACE (   7), SIDE POINTER=(-1)
SURFACE (   8), SIDE POINTER=(-1)
BODY    (   2)
0000000000000000000000000000000000000000000000000000000000000000
END      0000000000000000000000000000000000000000000000000000000

