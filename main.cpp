// Given is the following code structure/
// Data Files

// measurement.txt: The measurements from the sonar rangefinder sensors attached to the robot at each time stamp recorded over a period of 413 seconds. (timestamp, measurement 1:8).
// poses.txt: The exact robot poses at each timestamp recorded over a period of 413 seconds. (timestamp, x, y, ϴ).
// Global Functions

// inverseSensorModel(): You'll code this function as part of your second quiz after learning the inverse sensor model for sonar rangefinder sensors.
// occupancyGridMapping(): You'll code this function as part of your first quiz.
// Main Function

// File Scan: Scanning both the measurement and poses files to retrieve the values. At each time stamp, the values are passed to the occupancy grid mapping function.
// Display Map: After processing all the measurements and poses, the map is displayed.

#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include "src/matplotlibcpp.h" //Graph Library

using namespace std;
namespace plt = matplotlibcpp;

// Sensor characteristic: Min and Max ranges of the beams
double Zmax = 5000, Zmin = 170;
// Defining free cells(lfree), occupied cells(locc), unknown cells(l0) log odds values
double l0 = 0, locc = 0.4, lfree = -0.4;
// Grid dimensions
double gridWidth = 100, gridHeight = 100;
// Map dimensions
double mapWidth = 30000, mapHeight = 15000;
// Robot size with respect to the map 
double robotXOffset = mapWidth / 5, robotYOffset = mapHeight / 3;
// Defining an l vector to store the log odds values of each cell
vector< vector<double> > l(mapWidth/gridWidth, vector<double>(mapHeight/gridHeight));

double inverseSensorModel(double x, double y, double theta, double xi, double yi, double sensorData[])
{
     //******************Code the Inverse Sensor Model Algorithm**********************//
    // Defining Sensor Characteristics
    double Zk, thetaK, sensorTheta;
    double minDelta = -1;
    double alpha = 200, beta = 20;

    //******************Compute r and phi**********************//
    double r = sqrt(pow(xi - x, 2) + pow(yi - y, 2));
    double phi = atan2(yi - y, xi - x) - theta;

    //Scaling Measurement to [-90 -37.5 -22.5 -7.5 7.5 22.5 37.5 90]
    for (int i = 0; i < 8; i++) {
        if (i == 0) {
            sensorTheta = -90 * (M_PI / 180);
        }
        else if (i == 1) {
            sensorTheta = -37.5 * (M_PI / 180);
        }
        else if (i == 6) {
            sensorTheta = 37.5 * (M_PI / 180);
        }
        else if (i == 7) {
            sensorTheta = 90 * (M_PI / 180);
        }
        else {
            sensorTheta = (-37.5 + (i - 1) * 15) * (M_PI / 180);
        }

        if (fabs(phi - sensorTheta) < minDelta || minDelta == -1) {
            Zk = sensorData[i];
            thetaK = sensorTheta;
            minDelta = fabs(phi - sensorTheta);
        }
    }

    //******************Evaluate the three cases**********************//
    // Also consider the cells with Zk > Zmax or Zk < Zmin as unkown states
    if ((r > std::min(Zmax, (Zk + (alpha/2)))) || (fabs(phi - thetaK) > (beta / 2)) || ( Zk > Zmax) || (Zk < Zmin)) {
        return l0;
    }
    else if ((r < Zmax) && (fabs(r - Zk) < (alpha/2))){
        return locc;
    }
    else if (r <= Zk){
        return lfree;
    }
    return 0.4;
}

void occupancyGridMapping(double Robotx, double Roboty, double Robottheta, double sensorData[])
{
    int t = 1;
    for (int i=0; i<300;i++)
    {
        for(int j=0;j<150;j++)
        {
            double xi = i * gridWidth + gridWidth / 2 - robotXOffset;
            double yi = -(j * gridHeight + gridHeight / 2) +robotYOffset;
            double xd = Robotx - xi;
            double yd = Roboty - yi;
            if (( xd*xd + yd*yd ) <= (Zmax*Zmax))
            {
                l[i][j] = l[i][j] + inverseSensorModel( Robotx, Roboty, Robottheta, xi, yi, sensorData ) - l0;
            }
        }
    }
    
}

void visualization()
{
    //Initialize a plot named Map of size 300x150
		//namespace plt = matplotlibcpp;
		plt::title("MAP");
		//plt::figure_size(300, 150);
		plt::xlim(0, 300);
		plt::ylim(0, 150);
    
    //Loop over the log odds values of the cells and plot each cell state. 
    //Unkown state: green color, occupied state: black color, and free state: red color
    //Help: https://matplotlib.org/2.1.1/api/_as_gen/matplotlib.pyplot.plot.html 
		for(int x = 0; x < 300; x++){
			for(int y = 0;y < 150; y++){
				if (l[x][y] == 0){				
					plt::plot({x},{y},"gs"); // Unknown state
				}
				else if (l[x][y] > 0){				
					plt::plot({x},{y},"ks"); // Occupied state
				}
				else if (l[x][y] < 0){				
					plt::plot({x},{y},"rs");	// Free state	
				}
			}
		}
    
    //Save the image and close the plot 
		//plt::show();
		plt::save("./Images/Map.png");
		plt::clf();
}

int main()
{
    double timeStamp;
    double measurementData[8];
    double robotX, robotY, robotTheta;

    FILE* posesFile = fopen("poses.txt", "r");
    FILE* measurementFile = fopen("measurement.txt", "r");

    // Scanning the files and retrieving measurement and poses at each timestamp
    while (fscanf(posesFile, "%lf %lf %lf %lf", &timeStamp, &robotX, &robotY, &robotTheta) != EOF) {
        fscanf(measurementFile, "%lf", &timeStamp);
        for (int i = 0; i < 8; i++) {
            fscanf(measurementFile, "%lf", &measurementData[i]);
        }
        occupancyGridMapping(robotX, robotY, (robotTheta / 10) * (M_PI / 180), measurementData);
    }
    
    // Displaying the map
    for (int x = 0; x < mapWidth / gridWidth; x++) {
        for (int y = 0; y < mapHeight / gridHeight; y++) {
            cout << l[x][y] << " ";
        }
    }
    
    return 0;
}

/*
* void occupancyGridMapping(double Robotx, double Roboty, double Robottheta, double sensorData[])
* {
*    //******************Code the Occupancy Grid Mapping Algorithm**********************//
*   for (int x = 0; x < mapWidth / gridWidth; x++) {
*        for (int y = 0; y < mapHeight / gridHeight; y++) {
*            double xi = x * gridWidth + gridWidth / 2 - robotXOffset;
*            double yi = -(y * gridHeight + gridHeight / 2) + robotYOffset;
*            if (sqrt(pow(xi - Robotx, 2) + pow(yi - Roboty, 2)) <= Zmax) {
*                l[x][y] = l[x][y] + inverseSensorModel(Robotx, Roboty, Robottheta, xi, yi, sensorData) - l0;
*            }
*        }
*    }
* }
*/
    
    
    
    
