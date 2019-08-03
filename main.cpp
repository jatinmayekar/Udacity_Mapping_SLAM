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
using namespace std;

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
    // You will be coding this section in the upcoming concept! 
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
void occupancyGridMapping(double Robotx, double Roboty, double Robottheta, double sensorData[])
{
    //******************Code the Occupancy Grid Mapping Algorithm**********************//
    for (int x = 0; x < mapWidth / gridWidth; x++) {
        for (int y = 0; y < mapHeight / gridHeight; y++) {
            double xi = x * gridWidth + gridWidth / 2 - robotXOffset;
            double yi = -(y * gridHeight + gridHeight / 2) + robotYOffset;
            if (sqrt(pow(xi - Robotx, 2) + pow(yi - Roboty, 2)) <= Zmax) {
                l[x][y] = l[x][y] + inverseSensorModel(Robotx, Roboty, Robottheta, xi, yi, sensorData) - l0;
            }
        }
    }
}
*/
