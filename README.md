# Synthesis of Control Barrier Functions Using a Supervised Machine Learning Approach

This repository contains code to simulate results of our 2020 IROS submission titled, "Synthesis of Control Barrier Functions Using a Machine Learning Approach". The simulator used is called "Simpled Two Dimensional Robot (STDR) Simulator". 

## Getting Started
Please follow the instructions below to instantiate the environment for the simulation:

* Download the STDR simulator available at 
```
http://stdr-simulator-ros-pkg.github.io/
```

* Create a robot according to your own specifications by creating a .yaml file and saving it in 
```
/opt/ros/kinetic/share/stdr_resources/resources/robots/
```

* Create your own world. This requires an image of .jpg or .png extension. Resize the PNG based on your needs. Create a .yaml file with the same name as the image and mention the resolution in the image as needed. Changing the resolution in the .yaml file and the overall size of the image are the 2 ways in which the length & width of the domain in STDR can be controlled.

* Save both the image and the .yaml file in 
```
/opt/ros/kinetic/share/stdr_resources/maps/
```

* Run the command
```
/opt/ros/kinetic/setup.bash 
```
and 
```
source devel/setup.bash
```

Below we provide an example which illustrates the steps required to load a map and robot along with desired sensors in the STDR simulator

# Example
The following is an example to show how one can load maps and a robot with a LiDAR sensor in the STDR simulator:

* Run the following command in a command line/terminal window
```
roslaunch stdr_launchers server_no_map.launch
```

* In a new command line/terminal window, load the five obstacle map, by running the following command
```
rosrun stdr_server load_map src/stdr_simulator/stdr_resources/maps/FiveObs.yaml
```

* You should receive the following message once the previous command has been executed
```
[ INFO] [1581962946.363200566]: Loading map from image "src/stdr_simulator/stdr_resources/maps/FiveObs.png"
[ INFO] [1581962946.370371252]: Read a 832 X 522 map @ 0.004 m/cell
[ INFO] [1581962946.400851447]: Map successfully loaded
```

* In the same command line/terminal window, run the following command to place the robot at a specific initial condition
```
rosrun stdr_robot robot_handler add src/stdr_simulator/stdr_resources/resources/robots/omni_stdr.xml x y theta
```
where x, y are the position coordinates of the robot and theta is its orientation. You should receive the following message upon successful spawning of the robot
```
[ INFO] [1581963139.789220573]: New robot spawned successfully, with name robot0.
```

* In a new command line/terminal window, run the following command to launch the map with the robot at the specified initial condition
```
roslaunch stdr_gui stdr_gui.launch
```
