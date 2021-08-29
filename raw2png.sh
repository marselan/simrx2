#!/bin/bash
rm -f results/image*
rm -f results/png/*.png
cp image* results
octave to_png.m
