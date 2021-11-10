# Topographic-Horizons

This set of functions, based on the paper in IEEE Geoscience and Remote Sensing Letters [1],  computes the angles to the horizons from an elevation grid, measured in degrees upward from the horizontal. The one-dimensional problem uses an order N algorithm [2]. Horizons for arbitrary azimuths are derived by rotating the elevation grid, then calculating the horizons along the columns of the rotated grid and re-rotating the grid back to its original orientation. The code supports grids in either projected or geographic format; it also calculates the distances to the horizon.

The examples in the demo folder reproduce Figures 1 through 4 in the associated manuscript [1]. The demos require some functions from my Sun Position repository that are included in the demo folder, but the horizon functions themselves do not. To reduce the size of the DEM (digital elevation model), run the demos using MainDemo, which uses just a 0.25°×0.25° section of that paper’s topographic 1°×1° tile. Alternatively, run any of the Demo_...m functions in the demo folder with any geographic elevation grid Z and raster reference R.

[1]	J. Dozier, Revisiting the topographic horizon problem in the era of big data and parallel computing, IEEE Geosci. Remote Sens. Lett., 2021, doi: 10.1109/LGRS.2021.3125278.

[2]	J. Dozier, J. Bruno, and P. Downey, A faster solution to the horizon problem, Comp. Geosci., vol. 7, pp. 145-151, 1981, doi: 10.1016/0098-3004(81)90026-1.
