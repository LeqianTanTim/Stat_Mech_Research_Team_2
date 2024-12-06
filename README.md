This is a code update for Stat. Mech's umbrella sampling code. 

* Update Dec. 4th (5:31 pm) : 
  + Class Particle can be used to store many useful information without worrying array index issue too much.
  + Overload the minus operator to handle the combination and difference in position between a particle and the rest of the particle in the list.
  + Compute forces has completed.
  - md_run is still in progress.

 * Update Dec. 5th (11:53 am) :
   + Complete the set up before MD running
   + Some test statements are added to double check if import process is successfully used
   + Some new functions are imported in the F24_MD_code_Tim_Tan file, you can check to see the difference.
       - normalization and shift of com are now conducted there.
       - velocity, acceleration and dimension are now stored within the class for convenience accessibility.
       - md_run is still in progress unfortunately.
  * Update Dec. 6th
    + Incomplete the process of saving data. If you know how to modify it then this code is basically complete for its capability of md_run (# I hope)
