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

## Ocupancy Grip Mapping

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
   * Mapping wth known poses can be represented as P(m | z(1:t),x(1:t)) function
     x - pose nd z - measurements, so P is the posterior probablility over the        map given all the measurements(z) up to time t and all the poses(x) up to        time t represented by the robot's trajectory.
   * For 2D maps - use laser rangefinder - capture slice of 3D world - merge at        instant - partitioned into grid cells - estimate posterior through occupancy
     grid mapping algorithm.
   * For 3D maps -  can be done by occupancy grid mapping - but - much higher          computational memory - large no of noisy measurements - need filtering
   
 -  |     Robotics Challenge      |         Probability Equations       | 
    | --------------------------- | ----------------------------------- |
    |        Localization         |     P( x(1:t) \| z(1:t), m, u(1:t) )|
    |        Mapping              |     P( m \| z(1:t), x(1:t) )        |
    |        SLAM                 |     P( x(1:t), m \| z(1:t), u(1:t) )|
    
 - Grid Cells
    * Each grid - hold one value - either 0 or 1 - binary values.
    * Free or unoccupied space - 0
    * Obstacle or occupied space - 1
    * No of possible maps for a grid m\*n (rows \* columns) = 2 ^ (m\*n)
    * No of grid cells = m \* n
    
 - Three approcahes to calculate the posterior
    * P( m \| z(1:t), x(1:t) ) - Very high computational power
    * P( m(i) \| z(1:t), x(1:t) ) - Compute probablity of each cell independently
    * \|\| (i) P( m(i) \| z(1:t), x(1:t) ) - Relates cells & overcomes the huge computational memory to estimate the map with the product of marginals or factorization.
    
    
    
    
