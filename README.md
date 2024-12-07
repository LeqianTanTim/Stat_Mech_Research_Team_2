This is a code announcement board for Stat. Mech's umbrella sampling code project.

Ver. 1.2.0------ This is a non-official, still tuning version
* Warning!!
* Just a quick update on the ver 1.2.0 code if you decided to play around with it. There is currently a fatal flaw where the bias harmonic force is not correctly assigned (both solute had been assigned with the same force that basically just dragged them around in the same distance) and the set r0 value is not normalized by the box size as well, causing serious issues when computing. I would try to fix that in the meantime. That explains why this mess occurred.
* If you do encounter problems let me know as well.

- Ver.1.2.0 First attempt of implementing simple harmonic bias potential
  - A new set of compute_force_US and main_US is modified and optimized for umbrella sampling. 
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
