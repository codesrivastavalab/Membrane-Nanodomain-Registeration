# Contents 

This code will remove water molecules from within the lipid bilayer chains between a defined Z coordinate per leaflet headgroup region (an approximation). It is extremely useful when constructing new bilayer Molecular Dynamics system. The script has been designed for **GROMACS** and with water molecules defined as `SOL` and with an `OW` atom type present.  Providing you’ve built your bilayer normal to the z-axis, after solvating your system you can determine the two z points between which water should be removed by using the other code which extracts the highest and the lowest Z coordinates of the P headgroup atom. 

Please note that this code assumes you are have **MDAnalysis** installed.
