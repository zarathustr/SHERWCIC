# Simultaneous Hand-eye/Robot-World/Camera-IMU Calibration
This repository specifies the data acquisition and computation for the work "Simultaneous Hand-eye/Robot-World/Camera-IMU Calibration". The target platform we use for calibration is the Universal Robot 5 (UR5) industrial robot. First, clone the repository to the ROS workspace (after cloning, please use
``git submodule update --init --recursive`` 
to update the repo). Then, enter the librealsense and make install the repository, which includes the driver of Realsense camera D435i. After that, compile the ROS workspace and run the UR5 starting script from the launch file. A planning UI will be displayed and the user can plan the trajectory of the robot. Launch the Realsense data gathering program via the roslaunch file contained in the librealsense-ros repository. Finally, log all required data via rosbag record. To compute the calibration parameters, the logged rosbag should be converted into .csv and .png files. This repository consists of these converters which can be launched via rosrun or roslaunch easily.


# MATLAB Script Usage
Sample Data:
data1.mat and data2.mat contain synchronized data sequences of camera, robot and IMU.

Usage:
1. Make sure that the MATLAB version is beyond R2014b.
2. Setup the MATLAB parallel computation preferences. Set the thread numbers for the paralellization (e.g. 12 concurrent threads).
3. Run the file 'test_continuous_hec_rw_cam_imu_optXW.m'.


# Reference
Wu, J., Wang, M., Jiang, Y., Yi, B., Liu, M. (2020) Simultaneous Hand-eye/Robot-World/Camera-IMU Calibration. Submitted to IEEE/ASME Transactions on Mechatronics.
