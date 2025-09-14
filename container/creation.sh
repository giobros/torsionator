#!/bin/bash

# Step 1: Build the Docker image
docker build -t ubuntu22_cuda11.2 .

# Step 2: Save the image as a tar file
docker save -o cuda11.2.tar ubuntu22_cuda11.2

# Step 3: Convert the image to a Singularity (Apptainer) image
apptainer build image.sif docker-archive://./cuda11.2.tar

