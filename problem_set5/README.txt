For all the following questions, I took the cross section of a cylinder (2D circle).

Question 1:

Solving the potential for all space:

The plot of the potential is saved as "potential1.png".

The charge density of the numerical solution is saved as "density1.png".
It ressembles a line charge that is at the edges of the circle (2D cylinder).


As mentionned in the script, to find the true potential, the following steps were done:

we know that at the edge of the box, the potential is 0 and 
that in the circle the potential is 1 (or at the edge of the circle)

the radius of the box is just half the size of a side of the box

we can calculate the difference in potential and find the constant that is
multiplying the expression (lambda). We can then isolate for the constant 
to allow ourselves to have a potential of 0 on the edges of the box and a
potential of 1 in the circle.

The true potential was plotted and shown in graph : "potential_true1.png"
The true density was shown in "density_true1.png".
We can see that we don't have an exact analytic solution for this problem
since we are imposing that the values at the edges of the box are 0. 

(Would work if box >> than circle)

The solution is still far from the analytic answer, the potential effects 
are still close to the center.

For the error analysis presented in question 2, 
The script for question1 prints out : 
"on iteration 9774 residual is 0.009999580841263362"


+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Question 2:

I have set a limit so that the shown residual is less than 0.01.
The residuals are represented by r=b-Ax (and by taking the sum of the square elements
of that matrix)

For the conjugate gradient, the script prints out:

on iteration 75 residual is 0.009916414688149178

For a tolerance of 0.01, 
For the conjugate gradient, we get 75 steps for convergence 
For the method used in question 1, we get 9774 steps for convergence.

We can see that the conjugate gradient is much faster. 

Plot of potential is saved as "potential2.png"

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Question 3:

For the resolution, i started at 128 and went up by a factor of 2 every time 
I wanted to upgrade the resolution. 

The code printed out the following for the convergence of each resolution:

on iteration 42 residual is 0.009885762865525453
on iteration 1 residual is 0.007147441416496474
on iteration 1 residual is 0.006798765112495565
on iteration 1 residual is 0.007731419343935307
on iteration 1 residual is 0.007952709797850383
on iteration 1 residual is 0.008486286385191135


The higher the resolution, the less steps it takes to converge.

All plots were saved as "res_ 128 X 128 _potential3.png"
(changing the 128 by the actual resolution number)


++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Question 4: 

The solver from 3 was repeated with similar results for the potential and the iteration numbers.

The electric field was plotted on top of the potential. The image is saved as :"elec_fied_and_potential4.png"
The field lines and amplitudes are represented by arrows. 

I also plotted a 1D crossection of the 2D electric field (all values of the electric field located at the middle of the box, from top to bottom)
This is also where the bump is situated. 

We can see in the crossection plotted that there is a spike in the electric field at the position of the bump.
This plot is saved as "elec_bump.png". This is expected since the gradient of the bump is obviously larger when the radius is smaller. 
Having bumps increases the electric field strength at that point and loses more energy. No bumps are favorable.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Question 5:

For this question, since it is clearly asking a plot from the center of the heated side to the center of the opposite side, 
we can solve the 1D heat equation just for that line. (all other dimension contibutions cancel out since we are taking the center)
We are also not fixing our box to be 0 on the edges and hence we are not taking dirichlet boundary conditions. We are taking both
edges to be free. 

I am setting a constant k to be of the power of 10e-4.
If this constant is of order 10e-3 or higher, the heat function diverges and we get weird results.

The plot for the evolution of the heat equation in time from the center of the box (with linearly increasing potential)
to the center of the other side is saved as "heat_equation.png" . As time increases, the value of the heat on the left side increases linearly.
This plot shows different colors which show the evolution in time (superposition of plots in time). As time increases, the amplitude increases.