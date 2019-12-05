Please refer to "part3_periodic.py" script as the main script. (this part is the best commented as well)
All the other scripts are copied from that one with tweaks to answer each question.


PART1

Script for this part is shown in "part1.py"

Here are the parameters that were set to get the observed results:

m=20		mass
soft=1		softening
dt=1		time step
grid_size=30	size of the grid
n=1		number of particles

The potential was plotted and shown in:
"potential_part1.png" This is directly over the particle position of the particle.

To show that the particle isn't moving, I plotted the kinetic energy as a function 
of time steps, shown in : "Kinetic_energy_part1.png". It is always 0 (no change in velocity)
It means it feels no exterior force.

Substracting new particle position to initial position after 500 time steps:

print(x_new-part.x)
print(y_new-part.y)

both returned: [0.]
===================================================================
PART2

For 2 particles view script "part2.py"

The initial position and velocities were tweaked for the particle to start in an orbit.
(this can be seen in the code)
initial conditions: 
    part.x[0]=grid_size//2-15
    part.x[1]=grid_size//2+1
    
    part.y[0]=grid_size//2-15
    part.y[1]=grid_size//2+1
    
    part.vx[0]=0
    part.vx[1]=0
    
    part.vy[0]=1
    part.vy[1]=-1

Here are the parameters that were set to get the observed results:

m=20		mass
soft=1		softening
dt=1		time step
grid_size=150	size of the grid
n=2		number of particles

The animation is stored under: 

"2particles.gif"


===================================================================
PART3

Periodic Boundary conditions:

To view the results for this part, you can view the script 
"part3_periodic.py" and run it.

Here are the parameters that were set to get the observed results:

m=1		mass
soft=1		softening
dt=0.1		time step
grid_size= 500	size of the grid
n=100000	number of particles


for the particles with periodic boundary condition,
the file was too large so I had to upload it to the drive to show the animation
(if you download the link it will show a beautiful simulation), the file is 313MB
https://drive.google.com/open?id=1jqHlg60DRxmTnsOgUdTMLXfaSTSEjDgS

here is a description of the animation: many small galaxy like
objects are formed (clusters) since the density is not equal everywhere. Theses galaxies
approach each other and form bigger clusters. Eventually they all gather in a huge cluster
and explode (spreading everywhere, with big amounts of kinetic energy), from then on, the 
kinetic energy is much higher than the potential and particles are flying everywhere with no
specific pattern

A plot of the energy vs time is showed in "energy_periodic.png".
We can see that the energy is somewhat constant and then explodes releasing a lot 
of kinetic energy (huge peak in plot) almost doubles the initial energy. Then it stabilizes
at a higher constant than the initial energy.
It is observed that energy isn't conserved very well using this scheme. 


Non-periodic Boundary conditions:

The script: "part3_non_periodic.py"

We obtain a very similar scenario to when we have periodic boundary conditions, 
however we can see that the particles are all heading towards the center and are converging
there faster. When in the center they explode once again and a similar effect happens 
then in the periodic case.

m=1		mass
soft=1		softening
dt=1		time step
grid_size= 500	size of the grid
n=100000	number of particles


For the video file size to be reasonable, i had to increase the time steps:

shown in : "nparticles_nonperiodic.mp4"

A plot of the energy vs time is showed in "energy_nonperiodic.png".
We can see that the energy is somewhat constant and then explodes releasing a lot 
of kinetic energy (huge peak in plot) almost doubles the initial energy. Then it stabilizes
at a lower constant than the initial energy. (the particles get dispersed and get out of the grid)
It is observed that energy isn't conserved very well using this scheme. 


===================================================================
PART4