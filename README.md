#  Optimal Trajectory Planning for a Tractor Maneuver in a Field

This project focuses on optimizing the trajectory of a tractor maneuvering in a headland to improve the efficiency and accuracy of agricultural tasks. The optimization ensures that the tractor transitions between rows with minimal time and within constraints like headland boundaries, speed, and steering limits, aiming to reduce time spent on idle maneuvers while maintaining row alignment. The project has been developed as part of *CONSTRAINED AND NUMERICAL OPTIMIZATION FOR ESTIMATION AND CONTROL* exam held by Professor Fagiano at *Politecnico di Milano* during the 2023/24 academic year. See  [[1]](#1) as reference. 

[Read the full report](Optimal%20Trajectory%20Planning%20for%20a%20Tractor%20Maneuver%20in%20a%20Field_v1.1.pdf)
## Code
It's possibile to run effectively the optimization routine 

## Problem Statement
This work focuses on a specific maneuver that a tractor needs to perform multiple times. Since efficient use of the terrain is crucial in this field, maintaining rows parallel and equally spaced according to the row width is one of the most important features. To achieve this the tractor needs to arrive at the final point, the one from which the row is then worked, with a specific orientation. 
The operational area where the tractor must perform the maneuver before entering the next row, is called _headland_. The final pose of the tractor is heading parallel with respect starting point but in opposite direction. A schematic representation is provided in images below, depicting an ideal case without considering the real dynamics of the tractor.  
Manuever schematic            |  Headland particular
:-------------------------:|:-------------------------:
![schematic](Images\problem_description.png) | ![schematic](Images\Headland_particular.png) 

Here it's shown a realistic example of a headland, taken from a satellite image of a field in southern Milan, Italy. The green area defines the operational space, typically constrained by a prohibited area and the actual field, both of which must be avoided during the maneuvering phase and are here colored in red. The objective of this work is to develop an algorithm to compute the optimal trajectory for executing the described maneuver, minimizing the execution time.

## Possible Manuvers
To perform the manuever there are main possibilities, for two consecutive rows we're considering only the **bulb**, characterized by a great exploitation of the headland but with low execution time, and the **fishtail**, that exploits the reverse gear and can be considered as dual case.

Bulb            |  Fishtail
:-------------------------:|:-------------------------:
![schematic](Images\bulb.svg) | ![Fishtail](Images\fishtail.svg)

To deal with the selection of the proper according to the headland configuration a double optimization routine has been developed to better adapt according to different headland configuration.

## Results
The cost function that has been used to balances two competing objectives:

- Minimizing the time taken to complete the maneuver.
- Minimizing the error in the tractor’s final state (position and orientation).

The results demonstrate successful optimization of tractor maneuvers under various headland conditions, including the trailer:

Bulb (with trailer)            |  Fishtail (with trailer)
:-------------------------:|:-------------------------:
![schematic](Images\01____m00__q16\Trajectory.svg) | ![Fishtail](Images\01____m05__q8\Trajectory.svg)


## References
<a id="1">[1]</a> 
[ Xuyong Tu, Lie Tang] (2019). 
 Headland turning optimisation for agricultural vehicles and those with towed implements, [link](https://www.sciencedirect.com/science/article/pii/S2666154319300092).
