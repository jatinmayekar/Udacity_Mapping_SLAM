#include <iostream>
#include <math.h>
using namespace std;

const int mapWidth =  2;
const int mapHeight = 2;

void sensorFusion(double m1[][mapWidth], double m2[][mapWidth])
{
    // Fuse the measurments of the two maps using DeMorgan's Law and print the resulting 
    //map in a matrix form:
    //a  b
    //c  d
    double m3[mapHeight][mapWidth] = { {0, 0}, {0, 0} };
    for (int i = 0; i < mapWidth; i++){
        for (int j = 0; j < mapHeight; j++){
            m3[i][j] = 1 - ((1 - m1[i][j]) * (1 - m2[i][j]));
            
            if ( j == mapWidth - 1){ 
            cout << m3[i][j] << " \n";
            }
            else {
            cout << m3[i][j] << " ";    
            }
        }
    }
}

int main()
{

    double m1[mapHeight][mapWidth] = { { 0.9, 0.6 }, { 0.1, 0.5 } };
    double m2[mapHeight][mapWidth] = { { 0.3, 0.4 }, { 0.4, 0.3 } };
    sensorFusion(m1, m2);

    return 0;
}
