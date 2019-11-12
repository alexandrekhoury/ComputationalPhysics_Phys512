README

Assignment 4 phys 512
Alexandre Khoury

Question 1

a) The noise model was a power spectrum that was windowed and that was smoothed. 

The window was applied before fourier transforming the data to get the power spectrum. 
The window is a slightly modified cos function (with edges that go to 0).
To apply the window, just multiply the cos function mentionned above with the data then fourier transform to analyse further results.

The power spectrum was found by taking the fourier transform of the strain data and taking the absolute value squared

To smooth the power spectrum, I assumed that the noise was uncorrelated (gaussian). I then used a built in python function to 
apply a gaussian filter to the data points. The plots of the smoothed power spectrum (semi log) can be shown under Power_spectrum.PNG
Ideally it would be good to filter out the high and low range frequencies. However it didn't have much effect on our data. 

b)The matched filters are shown in Matched_filter.PNG
First column represents the different Hanford events and second column represent the different Livingston events.

c)An estimate of the noise of each even can be taken from the power spectrum.
The Signal to noise ratio was found using the matched filter and that estimated noise and is shown in SNR_plots.PNG
First column is Hanford, second column is Livingston and third column is a combination of both (SNR_H^2+SNR_L^2)^(1/2)

d)The signal to noise ratio was plotted (analytically) in: SNR_analytic.PNG
We get more signal in analytic model (without the matched filter). 
We have made a noise model using the match filter who's goal is to get as close as possible to the analytic solution.
The closer the matched filter analysis is to the actual noise, the close both SNR will be.

e)

For 1e, the frequencies are for Hanford:[ 92.03125 101.78125  78.78125  77.09375]
For 1e, the frequencies are for Livingston:[ 76.   115.90625  96.8125  107.21875]

These represents the frequencies of the events in order in hertz: 
1: 'H-H1_LOSC_4_V1-1167559920-32.hdf5','L-L1_LOSC_4_V1-1167559920-32.hdf5'
2: 'H-H1_LOSC_4_V2-1126259446-32.hdf5','L-L1_LOSC_4_V2-1126259446-32.hdf5'
3: 'H-H1_LOSC_4_V2-1128678884-32.hdf5','L-L1_LOSC_4_V2-1128678884-32.hdf5'
4: 'H-H1_LOSC_4_V2-1135136334-32.hdf5','L-L1_LOSC_4_V2-1135136334-32.hdf5'

f)We can see that for the first event the uncertainty in the matched filter is of order 10^-3 seconds
The analysis can be repeated for the different events and would yield error of same order. 

We can see by the code that I wrote that the difference in time detected by Hanford and Livingston is of the order 10^-3 seconds.
If we multiply that number by the speed of light, we get something of order 10^2. Now since sky coordinates are represented by angles, 
we can see that the gravitational wave hit closer to the Livingston detector then hit the Hanford detector but not in a straight line trajectory.
There is some angle associated to the position of the grav wave with deltay=10^2 and deltax=10^3. We can get an apporximate location in the sky
taking the arctangent. The error associated would be the error found previously for the matched filters times the speed of light
(order a little less than 10^2).
