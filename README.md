# Node_method
 This repository provides a *physical-point-of-view implementation* of the famous *node method* in thermal dynamic calculation.

Typically, people inclines to code node method with complex index operations to obtain the average temperature of the grey-shadowed mass flow blocks in the following Figure as suggested by the [paper](https://ieeexplore.ieee.org/document/7243359). And the method applies only to the scenarios where mass flow rates are fixed.

<img src="https://user-images.githubusercontent.com/93114918/194892740-efc07bed-56df-4c32-92ca-f86b1c60790b.jpg" width="600">

From the viewpoint of data-structure, the pipe in node method can be viewed as a simple `queue`, that is, each time we push a mass flow block into the inlet and some mass flow blocks pop out of the outlet. Actually we only have to calculate the average temperatures of those mass flow blocks popped out. Here, we store the temperatures and mass flows of those blocks in the attribute `NMpipe.m_buffer` and `NMpipe.T_buffer`. The function `NMpipe.generate_Tout()` sums the multiplications of temperature and mass flow buffers, and then outputs the outlet temperature by dividing the total mass flows in `NMpipe.m_buffer`. 

This kind of relization of node method has distinct physical meanings and can easily cope with scenarios of variable mass flow rates. The `NM_test.m` file gives a simple test case of [a real DHS located in Shijiazhuang, China](https://ieeexplore.ieee.org/document/9904489) where both inlet temperatures and mass flow rates are variables.

The code in this repository can deal with the cases where there are more than one pipes. In other words, one can use the code to perform the simulation of a district heating systems. But numerical simulations suggest that node method has [bad convergence performance](https://ieeexplore.ieee.org/document/9904489). 

Moreover, it should be pointed out that the object-oriented programming in MATLAB has a lot of limitations. It is a better choice to complement the code with the aid of Python, a more friendly language in terms of OOP. 

If you want to contribute to this repository, please contact me!
