# Run Away Robot with Unscented Kalman Filter

Self-Driving Car Engineer Nanodegree Program (Bonus Challenge)

---

<img src="./assets/video.gif?raw=true" width="500">




Making use of the previously implemented [UKF](https://github.com/camigord/Unscented-Kalman-Filter) to catch an escaped car driving in a circular path.

_The run away car is sensed by a stationary sensor, that is able to measure both noisy LiDAR and radar data. The capture vehicle will need to use these measurements to close in on the run away car. To capture the run away car the capture vehicle needs to come within 0.1 unit distance of its position. However the capture car and the run away car have the same max velocity, so if the capture vehicle wants to catch the car, it will need to predict where the car will be ahead of time._

### Strategy

We will make use of the UKF to estimate the current location of the run-away-car and to predict its future location. We start by setting an estimated _"catching time"_ equal to 3 seconds, which is equivalent to the time when we expect to intercept the run-away-car. We then employ the UKF to predict the future location of the car after _catching time_ seconds and aim our _"hunter car"_ in that direction. As new measurements are received, we update our expected _catching time_ and repeat the prediction and aiming steps. As the expected interception time decreases, our prediction error should decrease, thus giving us higher chances of catching the run-away-car.

### Running the Code

This project involves the Term 2 Simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases)

This repository includes two files that can be used to set up and intall uWebSocketIO for either Linux or Mac systems. For windows you can use either Docker, VMware, or even Windows 10 Bash on Ubuntu to install uWebSocketIO.

Once the install for uWebSocketIO is complete, the main program can be built and ran by doing the following from the project top directory.

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF`

### Dependencies

* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4
* uWebSocketIO
