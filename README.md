This is a code announcement board for Stat. Mech's umbrella sampling code project.

**Ver. 3.0.0**------
Molecular Dynamic Simulation algorithm version 3.0.0
Bug fixed in Ver 3.0.1 - Temperature was evaluated incorrectly. Now should be good. 
- Biggest update so far.
- Many Bugs have been fixed:
    - Non-uniformity of unit have been address. See my code and later posted pdf for more info.
    - Box_size is increased to better address the number of particles inside the box
    - New PBC condition was added to correctly placed particles (Thanks to Lauren who pointed out the functionality of rint)
- Ver 3.1.0 (expected Dec. 10th or 11th) would include the new bias potential and reactive look up. Currently these two functions are temporaly removed from Ver 3.0.0.
- Ver 3.1.0 also would have functions to visualize particle movements via animation.
- Ver 4.0.0 (expected Thursday or Friday) would include functions to combine the bias probability together.
- I am tired when typing this so mistakes may occurred. 

**Ver. 2.0.0**------ 
Molecular Dynamic Simulation algorithm version 2.0.0
- I am happy to annouce that umbrella sampling is finally, and successfully implemented!
- My test is conducted with the file fixed_test.csv which you can find this in my folder as well.
- If you want to see what it looks like, feel free to see the first_attempt_..._ figure. (with Ver.2.0.0)
- I would love to know the feedback with this code to see if there is any incompatibility issue.
- Ver.1.0.0 is still a great starting point if you want to make your own umbrella sampling algorithm. It contained most of the necessary details without the annoying bias potential.
- I would now head to work on grading lab reports and Dynamic's paper. Let me know if there is any feedback with the code to perfect it. 
  
- Ver.1.2.0 First attempt of implementing simple harmonic bias potential
  - A new set of compute_force_US and main_US is modified and optimized for umbrella sampling. 
  - Currently both are only available in the ver2.0 python file.
  - bonding distance plot is added to visualize the bond distance after bias is added. Check the png file to take a peek on my result!

- Ver.1.1.0 New functionality was added in for both python codes.
  - You can now monitor the distance between particles of solutes. (It currently only program to differentiate water from other particles)
  - reaction coordinates can now be monitored in a step-based processed.
  - Plotting is now a built in function in md_run_2.0 (so more controllability I hope)

**Ver.1.0.0**------
Molecular Dynamic Simulation algorithm version 1.0.0
- Sample input file can reference test_2 and test_14.
- Version 1.0.0 stands for (major update; minor update; bug fixes)
- Update Roadmap:
  - Conduct a simulation with two neutral solutes and track the distance between two solute molecules (Ver. 1.1.0)
  - Fix one of the solute as surface (Ver. 1.2.0)
  - Umbrella sampling implementation would be updated in version 2.0.0
