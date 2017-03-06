# Extended-Cepstral-Distance
A time series distance measure for efficient clustering of input/output signals by their underlying dynamics
## Summary
The code in this repository represents a minimal working example in Matlab for the simulations in "A time series distance measure for efficient clustering of input/output signals by their underlying dynamics, a summary of which can be found on ARXIV LINK, by Oliver Lauwers and Bart De Moor, which deals with the question of finding a distance measure for time series that captures aspects of the underlying dynamics of the time series.

The simulations consist of the following steps:
 - Constructing two electrical circuits (generateelectriccircuits.m)
 - Generating input signals for these circuits, applying them to the circuits, and capturing the output signals (generateelectriccircuits.m)
 - Constructing distance matrices out of the input/output data using the Euclidean distance, the Keogh Lower Bound to Dynamic Time Warping (dLBKeogh.m), the Martin cepstral distance, an extension of the Martin cepstral distance, a distance based on the H2-norm and a distance based on the H-infinity-norm. (clusterperformance.m)
 - Clustering the data using hierarchical clustering. (clusterperformance.m)
 - Comparing the clustering results with the ground truth through use of the Adjusted Rand Index (adjrandindex.m)

The experiments were run by the scrip "runExperiments.m". The different functions are then nested in each other.

Bear in mind that this code is not meant as a fully working software package, but serves merely as an illustration accompanying the manuscript mentioned earlier.

##Reference
When using this code or discussing results of the extended cepstral distance measure, please refer to ARXIV LINK