This code module applies a "simple Spectral Band Adjustment Factor (SBAF)" on LISS III and AWiFS to calibrate them with the help
of a reference image like Landsat 8 and Sentinel 2.

Inputs:

inpf_liss = path to folder LISS III or AWiFS bands
inpf_ref = path to folder containing reference image bands

Output:
Final output folder: Georeferenced

Four temporary images are also generated which are automatically deleted after the execution. These are two composite images,
one resampled image, and one clipped image

Note: Please make sure that both input and reference image have the same projection system. 
