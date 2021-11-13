Analysis glioma
===============

## Step 1: tracking

The tracking of the cells has been done using [ImageJ][1] with the plugin [TrackMate][2]. We use the DoG detection (Difference of Gaussian) with a radius of 20μm and a threshold of 1. As a result, we obtain two CSV files (see an example and an illustration in the folder `data_tracking/step1_data_trackMate/NPA_stich_3`) containing the positions $(x,y)$ of the cells (see the file `*spots.csv`) along with the trajectories (`*tracks.csv`).

## Step 2: smoothing trajectories

We then smooth the trajectories obtained with *TrackMate*. Smoothing is necessarily to make an estimate of  the velocity of each cell at any time. Without smoothing, the trajectories are too noisy. We use a Gaussian kernel as a filter with standard deviation $σ=2$ and a stencil of $9$ points. The implementation has been done in [Julia](https://julialang.org/), the script is in the folder `src`:
```julia
	> include("step2_filter_data.jl")
```
The velocity of the cells  ${\bf v}_i(t)$ at a given time $t$ is estimated using the classical formula:

```math
$${\bf v}_i(t) = \frac{{\bf x}_i(t+Δt)-{\bf x}_i(t-Δt)}{2Δt}$$
```

Thus, we have now an estimation of the positions ${\bf x}_i(t)$ and  velocities ${\bf v}_i(t)$ of the cells $i$. From the velocity ${\bf v}_i(t)$, we also deduce the velocity direction $θ_i(t)$. All the trajectories are saved in `jld2` file in the folder `data_tracking/step2_data_filtered`.

## Step 3: statistics in zone

We now perform statistical analysis in specific zones of the experiment. These zones are specified manually (see an illustration `data_tracking/step3_data_zone/NPA_stich_3/NPA_stich_3_zones.pdf`) and store in a *json* file called `coordinates.json`. This file is used in the julia following script:
```julia
	> include("step3_zone_create_df.jl")
```
which creates dataframes containing the analysis on each zone. This script has to be run twice to get both the velocities analysis and correlation functions (see the line 13 in the script).

The results of the analysis can be visualized using the two scripts:
```julia
	> include("step3_zone_velocity_plot.jl")
	> include("step3_zone_correlation_plot.jl")
```
## Step 4: classification flock/stream/swarm

To classify a zone as a flock, a stream or a swarm, we take the collection of angle $\theta_i$ inside a specific zone and compare the distribution to three distributions. The three distributions are a (wrapped) Gaussian for a flock, a symmetrize Gaussian for a stream and the constant function for a swarm. We use maximum likelihood and Akaike weight to select which distribution matches the best the distribution of $\theta_i$. From a practical point of view, it suffices to run the following [R](https://www.r-project.org/) script:
```R
	> source("step4_zone_testing_distribution.r")
```

## Additional remarks

The Julia scripts have been tested using Julia version 1.6.2 (stable) with the following packages:
- computing: `Revise`, `ProgressMeter`, `CSV`, `JSON`, `DataFrames`, `FileIO`
- visualization: `PyPlot`, `LaTeXStrings`
- math: `LinearAlgebra`, `Statistics`, `DSP`

The R script requires the library `circular`.


## References

1. Schindelin, Johannes, et al. *Fiji: an open-source platform for biological-image analysis.* _Nature methods_ 9.7 (2012): 676-682.
2. Tinevez, Jean-Yves, et al. *TrackMate: An open and extensible platform for single-particle tracking.* _Methods_ 115 (2017): 80-90.

[1]: https://imagej.net/software/fiji/ "Fiji: an open-source platform for biological-image analysis."
[2]: https://imagej.net/plugins/trackmate/ "TrackMate: An open and extensible platform for single-particle tracking"

## Legal notice information

The programs are distributed under the GNU GPL license version 2. More information are available in the file COPYING.txt. For any information or bugs, please contact me at: `smotsch[at]asu.edu`.


