An implementation of the Cornell Off-Lattice Particle-based Simulator (COPS)

Points to remember when editing .txt files:

particlenos.txt
-used to specify number of each type of species
-make sure these numbers are valid with specification in particleconfig.txt. Any errors will result in segmentation fault

reactions.txt
-reactions must be in the format (ax+by=cz) (where a,b,c are integers and x,y,z are lower case alphabets)
-there must not be any spaces in the reaction

diffusion.txt
-holds the diffusion constants and distance constants for diffusion (in that order)

particleconfig.txt
-contains initial setup of particles (particle type: fixed/random (if fixed)number)
-do not forget newline after last line

rateconstants.txt
-the rate constants for unimolecular reactions must come first, followed by bimolecular reactions

fixedparticles.txt
-contains floating point values for the particles of type which have been specified as fixed.
-must contain exact number of total fixed particles (as specified in particleconfig.txt)
-must have 3 space separated columns of floats