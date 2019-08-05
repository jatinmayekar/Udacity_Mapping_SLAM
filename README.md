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
  
 - Mapping is important in dynamic environment and as well as in static environments because sesnor measurements are always noisy and accumulate errors, and may at times produce a low accuracy map. It's always better to map a static environment and aim to correct the noise to obtain a high accuracy map.
 
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
  
 - MultiSensor Fusion
    * LIDAR sensor and RGBD sensor
    * Intuitive way - Build the occupancy grid mapping algorithm and solve for each sensor model(integrate in Zk) but this will fail as 
       * Each sensor has different characteristics
       * Different senstivity wrt to obstacles
    * Best approach - build seperate maps for each sensor model independent of each model and then integrate them
    * Use DeMorgans's Law(best) or take the maximum probability or perform a null operation
    * DeMorgan's Law: p(m<sub>i</sub>) = 1 - ||<sub>k</sub> (1 - p(m<sub>i</sub><sup>k</sup>))
    * where k is the no of sensors hence the no of maps
 
 - 3D Mapping
    * 2D and 3D maps(more) - computationally expensive to build and maintain
    * Collect 3D data using various tech:
       * 3D Lidar - single sensor array of beams stacked horizontally 
       * 2D Lidar - tilted or rotated 360 deg to obtain 3D coverage
       * RGBD camera  - single visual camera + laser rangefinder or infrared depth sensor - allows for  determination of depth of image - hence distance of object
       * Stereo camera - A  pair of offset cameras - used to directly infer distance of close objects - same as human eyes
       * Single camera - Cheap - Small - Complex monocular SLAM - Depth cannot be inferred directly -  calculated by analysing data from a sequence of frames in a video
     
  - 3D data representations
     * Desired charactersitics:
        * Probabilistic -  accomodate sensor noise and dynamic environments
        * Distinguish data between free and unknown space  - enable robot to plan an unobstructed path and build a complete map
        * Memory efficient - memory on mobile robot is limited - therefore map should be accesible to robot's main memory while mapping over a large are for a large period of time -  to accoplish this - need a memory data representation that is compact and allows for efficient updates and queries
        
     * Types of data representations:
         * Point clouds -  set of data points corresponding to range measurement at defined x,y,z coordinates - flaw - exists only where things are in the world - does not distinguish between unoccupied and unknown states - stores a large amount of measurement points - each scan requires more memory - memory inefficient
         * Voxels - Volumetric data representation using a grid of cubic volumes of equal size - probabilistic memory representation - so can distinguish between free and unoccupied states - flaw - needs info on size of area known or approximated before measurement - seconly, complete map must be allocated in memory so memory requirement is high
         * Octrees - memory efficient tree based data representation - can be dynamically expanded to different resoltuions and different areas - where every voxel can be subdivided into 8 voxels recursively - because map volumes are not initialzed until you need to add new measurement - have been used to adapt occupancy grid mapping from 2D to 3D - introducing probabilistic represenation of occupied vs free space  
         * 2.5D maps
         * Elevation maps
         * Extended elevation maps
         * Multi-level surface (MLS)
         
    * OctoMap
       * OctoMap framework is a open-source C++ library and ROS package based on Octrees
       * can be used to generate volumetric 3D models
       * Not a 3D SLiM solution - mapping framework and needs a pose estimate
       * It converts and integrates point clouds into 3D occupancy maps
       * Uses a probabilistic occupancy estimation modelled as a recursive binary bayes filter - static filter which assumes the environment does not change 
       * Efficient updates are achieved using the log odds notation
       * Represented volumeterically with modeling of free, occupied and unmapped areas
       * Upper and lower limits are placed on the log odds notation to limit the no. of updates required to change the state of voxel
       * Supports multi-resolution map queries - min. voxel size - determines the resolution
       * Tree pruining is also used to remove the redudancy between discrete occupancy states - achieved by defining a threshold probability for a free or occupied voxel - childrren that are identical to the parent can be pruned - 
       * Memory efficiency is achieved by using a compression methods that produces compact map files - coherent vlumes are locally combined including both mapped free areas and occupied space 
       
   - SLAM 
     * Challenging as both poses and map are unknown
     * Also called as "CLAM - Concurrent Localization and Mapping"
     * Input = Measurment +  Controls
     * Output = Map + Trajectory
     * Two forms of SLAM 
        * Online SLAM
           * At time t, the robot will estimate its new pose xt and the map m given only its current measurements zt and controls ut
           * It solves instantaneous poses using current measurements and controls independently from previous measurements ad controls - estimate variables that occur at time t only
           * Problem can be modeled by probability equation: p ( x<sub>t</sub>, m : z<sub>1:t</sub>, u<sub>1:t</sub>) where the posterior is solved by instantaneous pose and the map given the current measurements and controls
           * Posterior over the current pose
           ![Online SLAM](https://github.com/gonfreces/Udacity_Mapping_SLAM/blob/master/Online%20SLAM_1.png)
        * Full SLAM
           * Estimate entire path upto time t instead of a instantaneous pose given all measurements and controls
           * Problem can be modeled by probability equation: p ( x<sub>1:t</sub>, m : z<sub>1:t</sub>, u<sub>1:t</sub>)
           * Posterior over the entire path 
           ![Offline SLAM](https://github.com/gonfreces/Udacity_Mapping_SLAM/blob/master/Full_Offline%20SLAM.png)
           
        * Realtion between Online and Full SLAM
        ![Offline SLAM](https://github.com/gonfreces/Udacity_Mapping_SLAM/blob/master/RelationSLAM.png)
        
   * Nature of SLAM
      * Continous -  Object locations(sense environment for landmarks or objects) + Poses(odometry) 
      * Discrete - SLAM algorithms have to identify if a relation exists between the newly detected object and the previously detected object - helps robot identify if it has been in this location before - Answer to this binary -yes or no - that's what makes it discreet - Correspondence
      
   * Correspondence
      * Online SLAM - p (x<sub>t</sub>, m, c<sub>t</sub> : z<sub>1:t</sub>, u<sub>1:t</sub>)
      * Full SLAM - p ( x<sub>1:t</sub>, m, c<sub>1:t</sub> : z<sub>1:t</sub>, u<sub>1:t</sub>)
      * Reason to consider it - To help robot better understand its position by establishing relation between objects
      ![Correspondence](https://github.com/gonfreces/Udacity_Mapping_SLAM/blob/master/RelationSLAM_C.png)
         
   * SLAM Challenges
      * Continuous - encounter many objects - keep track of them - no. of variables increases - makes the problem highly dimensional and challenging to compute the posterior
      * Discreet - Large correspondence values  - increase exponentially 
      * Therfore SLAM algorithms will have to rely on approximation while estimating a posterior in order to conserve computational memory
      
   * FastSLAM
      * Adding map as a new dimension to particle in MCL fails - as map is modeled with any variables - high dimensionality - particle filter approach to SLAM in this current form will scale exponentially and is doomed to fail
      * It uses a custom particle filter approach to solve full SLAM problem with known correspondences      
      * Uses a particles to estimate posterior over the robot path along the map
      * Each of this particle holds the robot trajectory which will give an advantage to SLAM to solve the problem of mapping with known poses
      * Each of the particle also holds a map and each feature of the map is represented by a local Gaussian
      * The main problem is now divide into seperate independent problem, each of which aims to solve the problem of estimating the features of the map
      * To solve this independent mini problem, fastSLAM uses low dimensional extended kalman filter
      * While map features are expressed independently, dependency exists only between the robot pose uncertainity
      * This custom approach of representing posterior with particle filter and Gaussian is called Rao-Blackwellized Particle Filter One
      * Estimating the Trajectory: FastSLAM estimates a posterior over the trajectory using a particle filter approach. This will give an advantage to SLAM to solve the problem of mapping with known poses
      * Estimating the Map: FastSLAM uses a low dimensional Extended Kalman Filter to solve independent features of the map which are modeled with local Gaussian
      * Capable of solving both the Full SLAM and Online SLAM problems
           * FastSLAM estimates the full robot path, and hence it solves the Full SLAM problem
           * On the other hand, each particle in FastSLAM estimates instantaneous poses, and thus FastSLAM also solves the Online SLAM problem
      * Three diferent instances of Fast SLAM:
          * FastSLAM 1.0
               - Simple and easy to implement
               - Inefficient as particle filters generate sample inefficiency
               - Landmark based algo - therefore not able to model an arbitrary environment
          * FastSLAM 2.0
               - overcomes this inefficiency by implementing different distribution which generates low no of particles
               - Keep in mind that both of the FastSLAM 1.0 and 2.0 algorithms use a low dimensional Extended Kalman filter to estimate the posterior over the map features
               - Landmark based algo - therefore not able to model an arbitrary environment
          * Grid-based FastSLAM
               - Adpats FastSLAM to grid maps - extends to occupancy grid mapping
               - Non-landmark based algo - can model and solve in an arbitrary environment 
               
   * Grid-based FastSLAM
      * Probability equation: p(x<sub>0:t</sub>, m : z<sub>1:t</sub>, u<sub>1:t</sub>) = p(x<sub>0:t</sub>: z<sub>1:t</sub>, u<sub>1:t</sub>) * p (m : x<sub>1:t</sub>, z<sub>1:t</sub>) [Robot trajectory * map]
      * The grid-based FastSLAM algorithm will update each particle by solving the mapping with known poses problem using the occupancy grid mapping algorithm
      * Three different techniques represented by three different probability functions to adapt fast SLAM to grid mapping:
          * Sampling motion - p(x<sub>t</sub>: x<sub>t</sub><sup>k</sup>, u<sub>t</sub>) Estimates the current pose given the k-th particle previous pose and the current controls u - solved by MCL
          * Map estimation - p(m<sub>t</sub>: z<sub>t</sub>, x<sub>t</sub><sup>k</sup>, m<sub>t-1</sub><sup>k</sup>) - Estimates the current map given the current measurements, the current k-th particle pose, and the previous k-th particle map - solved by Occupancy Grid Mapping 
          * Importance weight - p(z<sub>t</sub>: x<sub>t</sub><sup>k</sup>, m<sub>t</sub><sup>k</sup>)Estimates the current likelihood of the measurement given the current k-th particle pose and the current k-th particle map - solved by MCL
          
      * gmapping ROS package
         * Provies laser based SLAM - feed its node with the robot laser measurements and odometry values and expect it to provide you with a 2D occupancy grid map of the environment. The map will be updated as the robot moves and collect sensory information using its laser range finder sensor.


       
    
    
    
    
    
    
    
