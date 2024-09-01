# Obstacle Problem Using Random Minimizing Movement

This repository contains the code used for the numerical simulations presented in our paper [*Mini-Batch Descent in Semiflows*](https://arxiv.org/abs/2407.07556). The code implements the **Random Minimizing Movement** scheme applied to solving the obstacle problem.

## Problem Description

In our paper, we explore the application of **Random Minimizing Movement**, a method where in each iteration of the minimizing movement, a randomly selected sub-functional of the total functional is minimized rather than following the full functional. This approach is illustrated using the obstacle problem, comparing the obtained solutions with real solutions.


<p align="center">
    <img src="https://github.com/Martinshs/mini-batch_obtacle_problem/blob/main/comparation_solutions.gif" alt="" /></a>
</p>



## Repository Structure

- `utils.py`: Contains several functions. Among them are:
  * smoothmax : Used to approximate the max function
  * two_mountains : 2d function considered as the obstacle. It has the shape of two mountains
  * obs_1d : 1d function considered as the obstacle for the 1d case.
  * make_gif :  Generate gifs
  * batches_gen : Generate a set that contains sets randomly generated with a certain probability. 
  * unity_partition_1d : Partition of the unity 1d, defined in the interval [-1,1]
  * unity_partition_2d : Partition of the unity 2d, defined in the interval [-1,1]x[-1,1]
              
  
- `solvers.py`: Contain the solvers. Namely,
  * solver_op :  Obstacle problem solver
  * sol_one_realization :  Solve one realization of the random minimizing movement
  * solver_avg_mb_op : Compute the average of a given number of realizations of the random minimizing movement
- `mini_batch_obstacle_problem.ipynb`: Jupyter Notebooks with examples demonstrating the code.

## References

For more details on the theory and numerical results, please refer to Section 4.3 of our paper on arXiv: [Mini-Batch Descent in Semiflows](https://arxiv.org/abs/2407.07556).
