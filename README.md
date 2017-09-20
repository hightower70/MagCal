# MagCal
Magnetic sensor calibration in c# using Q. Li ellipsoid fitting. The original algorithm has been published by Q. LI here: https://www.mathworks.com/matlabcentral/fileexchange/23377-ellipsoid-fitting

This implementation is a straight conversion of the algorithm written in C and published here: https://sites.google.com/site/sailboatinstruments1/step-1

The calibrated magnetic sensor data can be obtained by the following formula:

![calibration formula](http://latex.codecogs.com/gif.latex?h_%7Bcal%7D%20%3D%20A%5E%7B-1%7D%20*%20%28h%20-%20b%29)
