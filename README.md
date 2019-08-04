# Udacity_Mapping_SLAM

## SLAM 
Robot must construct map of the environment while simulataneoulsy localizing

### Categories of SLAM Algorithms:
1. Extended Kalman Filter SLAM(EKF)
2. Sparse Extended Information Filter(SEIF)
3. Extended Information Filter(EIF)
4. FastSLAM - Particle filter with low dimesnional EKF
5. GraphSLAM - Uses constraints to relate between robot poses and the environment,
and then tries to resolve all this contsraints to create the most likely map given the data. Implementation of Graph-SLAM - RTABmap (Real Time Appearance Mapping)

## Ocupancy Grid Mapping

- Localization
  * Assumption - Known Map
  * Estimation - Robot's Trajectory
  
- Mapping
  * Assumption - Robot's Trajectory
  * Estimation - Known Map
  
 - Mapping is important in dynamic environment and as well as in static                environments because sesnor measurements are always noisy and accumulate errors,    and may at times produce a low accuracy map. It's always better to map a static    environment and aim to correctthe noise to obtain a high accuracy map.
 
 - Challenges & Difficulties of Mapping
    1. Unkown map
    2. Huge Hypothesis space - Continous spread of map
    3. Size of environment
    4. Noise in sensor measurements
    5. Perceptual Ambiguity - Similar objects in environment
    6. Cycles - Traversing the same path twice back & forth
    
 - Mapping with known poses 
   * The problem of generating a mapunder the assumption that the robot poses are
   known and non-noisy is referred to as mapping with known poses. In SLAM,          mapping is done with unknown poses and then occupancy grid mapping uses these    exact same posesand noisy measurements to generate a (posterior)map for path      plannig and navigation.
   
 - Posterior Probaility 
   * Mapping wth known poses can be represented as P(m | z<sub>1:t</sub>,x<sub>1:t</sub>) function
     x - pose nd z - measurements, so P is the posterior probablility over the        map given all the measurements(z) up to time t and all the poses(x) up to        time t represented by the robot's trajectory
   * For 2D maps - use laser rangefinder - capture slice of 3D world - merge at        instant - partitioned into grid cells - estimate posterior through occupancy
     grid mapping algorithm
   * For 3D maps -  can be done by occupancy grid mapping - but - much higher          computational memory - large no of noisy measurements - need filtering
   
 -  |     Robotics Challenge      |         Probability Equations       | 
    | --------------------------- | ----------------------------------- |
    |        Localization         |     P( x<sub>1:t</sub> \| z<sub>1:t</sub>, m, u<sub>1:t</sub> )|
    |        Mapping              |     P( m \| z<sub>1:t</sub>, x<sub>1:t</sub> )        |
    |        SLAM                 |     P( x<sub>1:t</sub>, m \| z</sub>1:t</sub>, u<sub>1:t</sub> )|
    
 - Grid Cells
    * Each grid - holds one value - either 0 or 1 - binary
    * Free or unoccupied space - 0
    * Obstacle or occupied space - 1
    * No of possible maps for a grid m\*n (rows \* columns) = 2 ^ (m\*n)
    * No of grid cells = m \* n
    
 - Three approcahes to calculate the posterior
    * P( m \| z<sub>1:t</sub>, x<sub>1:t</sub> ) - Very high computational power
    * P( m<sub>i</sub> \| z<sub>1:t</sub>, x<sub>1:t</sub> ) - Compute probablity of each cell independently
    * \|\| <sub>i</sub> P( m<sub>i</sub> \| z<sub>1:t</sub>, x<sub>1:t</sub> ) - Best approach - Relates cells & overcomes the huge computational memory to estimate the map with the product of marginals or factorization.
    
 - Measurement Model Selection
    * Due to factorization, a binary estimation problem has to be solved to identify grid cells holding static state that does not change over time. By static, it means that the state of the system does not change during sensing. 
    * Binary Bayes filter solves this problem. It uses the log odds ratio of the belief. With static state, the belief is now a function of measurements only.
    * bel<sub>t</sub> (x) = p(x | z<sub>1:t</sub>) -  Inverse measurement model 
    * bel - Represents binay state of the model w.r.t. measurements 
    * Depending upon the measuremnet value -  state of the grid is updated
    * Log odds ratio form = l<sub>t</sub> = p(x | z<sub>t</sub>) / (1 - p( x | z<sub>t</sub>) )
    * Advantage of using log odds ratio is to avoid probablility intabilities near 0 or 1. Another advantage relates to system speed, accuracy, and simplicity
    * Forward vs. Inverse Measurement Model
    * Forward Measurement Model - P(z<sub>1:t</sub> | x): Estimating a posterior over the measurement given the system state.
    * Inverse Measurement Model - P(x | z<sub>1:t</sub>): Estimating a posterior over the system state given the measurement.
    * The inverse measurement model is generally used when measurements are more complex than the system's state.
    
 - Binary Bayes Filter Algorithm
    * l<sub>t</sub> = l<sub>t-1</sub> + log( p(x | z<sub>1:t</sub>) / (1 - p(x | z<sub>1:t</sub>) ) - log( p(x) / 1 - p(x) )
    * New belief = l<sub>t</sub>
    * Previous belief = l<sub>t-1</sub>
    * Log odd of the inverse measurements model = log( p(x | z<sub>1:t</sub>) / (1 - p(x | z<sub>1:t</sub>) )
    * Initial belief = log( p(x) / 1 - p(x) ) - Represents the initial state of the system before taking any sensor measurements into consideration
    
 - Occupancy Grid Mapping Algorithm
    * for all cells m<sub>i</sub> do
    * if m<sub>i</sub> in perceptual filed of z<sub>t</sub> then
    * l<sub>t,i</sub> = l<sub>t-1,i</sub> + log( p( m<sub>i</sub> | z<sub>t</sub>, x<sub>t</sub>) / (1 - p( m<sub>i</sub> | z<sub>t</sub>, x<sub>t</sub> ) - log( p(m<sub>i</sub>) / 1 - p(m<sub>i</sub>)
    * else 
    * l<sub>t,i</sub> = l<sub>t-1,i</sub>
    * end if
    * endfor
    
  - Viusalization of Map created using Occupancy grid mapping algorithm
  * Legend
       * Black - Occupied state
       * Red -  Free state
       * Green - Unknown state
       
  ![Map](https://github.com/gonfreces/Udacity_Mapping_SLAM/blob/master/Map.png)
  
  Ignore this warning when compiling
  ![Warining](https://github.com/gonfreces/Udacity_Mapping_SLAM/blob/master/ogm_warning.png)
  
  
    
    
    
    
    
    
    
