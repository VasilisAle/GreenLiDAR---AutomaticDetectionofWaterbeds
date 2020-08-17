## Introduction
The code in this repository was developed during the graduation project of Vasileios Alexandridis for the MSc Geomatics at TU Delft, the Netherlands (the thesis is available [here: to be added]())

Bathymetric Airborne LiDAR technology is used to map the depth of the water bodies. The green light sensor is able to penetrate the water surface and reach the bottom of the interesting water areas. 
Factors such as the water clarity, the water turbidity (waves) and the vegetation are some of the water conditions that affect the capability of the green laser penetration. 

A workflow presents the required procedures to process raw green LiDAR point cloud of water bodies, then classify them into **THREE** classes: **_water surface_**, **_underwater_** and **_bottom_** points. **Pulse-** and **Voxel-based** methods were implemented to perform a classification process with high level of automation. Points' characteristics such as *intensity*, *return number*, *number of returns*, *GPS time* were analysed per pulse. Voxel The spatial distribution of points into voxels - water columns was examined based on factors: *elevation*, *density*, *distance*, *intensity*.

## Dependencies
The automatic detection of waterbeds programm is written in ``` julia - 1.3.1 ```. The implementation is dependent on the following libraries:

[LazIO v0.2.1](https://github.com/evetion/LazIO.jl)

[LasIO v0.3.6](https://github.com/visr/LasIO.jl)

[StaticArrays v0.12.4](https://github.com/JuliaArrays/StaticArrays.jl)

[JSON v0.21.0](https://github.com/JuliaIO/JSON.jl)

[FileIO v1.4.0](https://github.com/JuliaIO/FileIO.jl)

[Plots v1.5.8](https://github.com/JuliaPlots/Plots.jl)

Some parts are written in ``` Python 3.8 ```. The python libraries are: 

[rasterio](https://rasterio.readthedocs.io/en/latest/intro.html)

[startin 0.6.0](https://github.com/hugoledoux/startin)

[numpy 1.18.1](https://pypi.org/project/numpy/1.18.4/)

## Structure
The repository is structured as follows: 

- **/code**: the Julia scripts that contain the voxel and pulse-based implementations and some plots (eg. histograms)

## Usage
- Pulse-based method:  ``` pulse_approach.jl ``` is used to categorize the points into laser pulses based on their z value, intensity, return number (rn), number of returns (nr) and GPS time. Specifically, it reads  the input dataset from the ``` pulse.json ``` file and groups the points into triplets or doubles. Then, it checks if the last point's intensity is bigger than the middle's one in order to mainly identify the bottom points.
- Voxel-based method:  ``` voxelization.jl ``` is used to divide the 3D space into voxels with specific x, y and z dimensions, which are specified into the ``` voxelization.json ``` file. It reads the boundaries of the input dataset and then identifies the points that fall into each voxel. For each voxel a histogram with specific bin size is created and the bottom points are identified. Based on the density, distance and intensity values, the confidence values are computed for each potential bottom point. These confidence values are stored as classification code in the LAS file. 

## Parameters
The ``` voxelization.json ``` and ``` pulse.json ``` files contain the parameters that can be set by the user for both implemenentations. 
The options for the ``` voxelization.json ``` each are the following:

- ``` dataname ```: the name of the LAS file
- ``` vx ```: the voxelsize in x dimension
- ``` vy ```: the voxelsize in y dimension
- ``` output_name ```: the name of the classified point cloud
- ``` raster_output ```: output raster with 5 bands: z_value, density, distance, intensity, confidence value

The options for the ``` pulse.json ``` each are the following:

- ``` dataname ```: the name of the LAS file
- ``` output_name ```: the name of the classified point cloud

## License
[MIT](https://choosealicense.com/licenses/mit/)
