#!/bin/bash

scp -r dhruvsrikanth@midway3.rcc.uchicago.edu:/home/dhruvsrikanth/Motif_4/output.zip ./output

unzip ./output.zip -d ./

rm ./output.zip