This is a code update for Stat. Mech's umbrella sampling code. 

* Update Dec. 6th
  + Incomplete the process of saving data. If you know how to modify it then this code is basically complete.
  + Two test files were provided, test_2 contain one solute and one solvent particle so you can play around with ease; test_14 is a much bigger file, recommend to use when you are much more confident with the current code implementation.
  + Random seed modification was added to F24_... python file so you can controlled the random (Thanks to Nik)
  + Currently epsilon is in kJ/mol. Divide this value by Avogadro's number if you want per particle. 

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

