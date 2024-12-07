This is a code announcement board for Stat. Mech's umbrella sampling code project.

Ver. 1.2.0------ This is a non-official, still tuning version

- Ver.1.2.0 First attempt of implementing simple harmonic bias potential
  - A new set of compute_force_US and main_US is modified for specific usage. (Differentiate them)
  - Currently both are only available in the ver2.0 python file.
  - bonding distance plot is added to visualize the bond distance after bias is added. Check the png file to take a peek on my result!

- Ver.1.1.0 New functionality was added in for both python codes.
  - You can now monitor the distance between particles of solutes. (It currently only program to recognize water as the only type of solvent)
  - reaction coordinates can now be monitored in a step-based processed.
  - Plotting is now a built in function in md_run_2.0 (so more controllability I hope)

Ver.1.0.0------
Molecular Dynamic Simulation algorithm version 1.0.0
- Sample input file can reference test_2 and test_14.
- Version 1.0.0 stands for (major update; minor update; bug fixes)
- Update Roadmap:
  - Conduct a simulation with two neutral solutes and track the distance between two solute molecules (Ver. 1.1.0)
  - Fix one of the solute as surface (Ver. 1.2.0)
  - Umbrella sampling implementation would be updated in version 2.0.0

* Update Dec. 6th (Ver. 1.0.0)
  + Congrats to myself! MD Code version 1.0.0 is now published!
    - Data can now be saved as output.xyz, at each dump frequency the energy, temperature and position data would be recorded.
    - I also added a matplotlib python code where you can visualize the energy change over steps. 
  + Two test files were provided, test_2 contain one solute and one solvent particle so you can play around with ease; test_14 contained one solute and 13 solvent particles.
  + Random seed modification was added to F24_... python file so you can controlled the random (Thanks to Nik)
  + Currently epsilon is in kJ. Multiply this value by Avogadro's number if you want per mole basis.

* Update Dec. 5th (11:53 am) :
   + Complete the set up before MD running
   + Some test statements are added to double check if import process is successfully used
   + Some new functions are imported in the F24_MD_code_Tim_Tan file, you can check to see the difference.
       - normalization and shift of com are now conducted there.
       - velocity, acceleration and dimension are now stored within the class for convenience accessibility.
       - md_run is still in progress unfortunately.

* Update Dec. 4th (5:31 pm) : 
  + Class Particle can be used to store many useful information without worrying array index issue too much.
  + Overload the minus operator to handle the combination and difference in position between a particle and the rest of the particle in the list.
  + Compute forces has completed.
  - md_run is still in progress.

