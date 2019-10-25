Question 1

Chisquare was found to be : 
1588.4366720631901

which is printed in the script. Using the definition of chisquare

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Question 2:

The new parameters, keeping a fixed tau are : 

the arays will be presented in the following order: 
H0,ombh2,omch2,As,ns 

array([6.31721130e+01, 2.20404629e-02, 1.22012412e-01, 2.08134880e-09,
       9.49420085e-01])

with errors: (taken from the covariant matrix)

array([1.95280321e+00 4.92017979e-04 5.03091301e-03 3.90389095e-11
 1.21777174e-02])

for fixed tau of 0.5

with a chisquare of :1245.8089154720751

We know the derivatives are good since it lowered the initial chisquare.
The derivative with respect to the first parameter (H0) was taken and shown in (plot_deriv_Q2)

I have also shown the derivative with respect to tau(plot_dervi_tau_free.PNG)


-----------------------------------------------------
Now taking tau to be the a non-fixed parameter, we get: 

CAMBError: Error in Fortran called from calc_transfer:
Reionization did not converge to optical depth
tau does not converge

first guess:(we would expect the errors to be higher since there are more free parameters,
 however this is not necessarily true)

from question 3, we can see that tau is correlated with As, hence making tau a
non-fixed parameter, would increase the error on As but not on the other parameters

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Question 3 

*****In the plots, the parameters are shown in order H0,ombh2,omch2,tau,As,ns *********
*****for questions 3 and 4*******

new params foro 5000 steps and scale factor of 1/2:
array([7.45879105e+01, 2.36651494e-02, 1.04413600e-01, 1.96868785e-01,
       2.68738349e-09, 1.01251400e+00])

take parameters by taking the average of the chain at half way (better estimate)

with error: (found from covariance matrix)
array([3.51642250e+00, 6.83128239e-04, 6.41326186e-03, 4.25936278e-02,
       1.92723132e-10, 1.93471883e-02])

We can see that the chains converge since they start looking like random noise, with
fluctuations that aren't very big. (Can view chains in plot_Q3.PNG)
We can also see that the fast fourrier transforms converge (become constant a later times)
view (FFT_Q3.PNG)

We can also see the correlation between variables by observing the corner plot (corner_plots_Q3.PNG)
we can see that there is almost no correlation except between tau and As 
(which eventually tends to a gaussian as well after the chain had been run for a long time)


++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Question 4 

new params for 5000 steps and scale factor of 1/2

array([7.08997410e+01, 2.27668159e-02, 1.09789379e-01, 5.42885685e-02,
       2.04512148e-09, 9.79046880e-01])

with error (constraints):
The error was taken from the covariance matrix

array([1.90589237e+00, 4.56490793e-04, 3.96390631e-03, 3.95030844e-03,
       3.13964475e-11, 1.23960567e-02])

----------------------------------------
Part 2 of question 4 

The parameter tau was weighted using a gaussian function. 
The mean was taken to be the given value of 0.0544 with sigma = 0.0073

I took each data point from my chain found in Question 3 for tau
to find the gaussian weight on each data point

Then i applied the gaussian weight on the chains and calculated a weighted average. 

The new parameters are now: 
[6.89953393e+01 2.27358714e-02 1.14533841e-01 8.56586935e-02
 2.22119300e-09 9.73844503e-01]

With error: 
[3.65705113e+00 5.48308342e-04 4.82568067e-03 1.70288607e-02
 9.00456316e-11 1.67407239e-02]


All the new parameters found are within 2 sigma of the parameters for question 3. 
They seem to be more precise, since we have weighed them knowing that one of the parameters behaves a certain way

