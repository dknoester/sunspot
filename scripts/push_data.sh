#!/bin/bash

rsync -avz -e ssh ./data/ dk@hpc.msu.edu:data/sunspot/
