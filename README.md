# Change of variables in 2D and 3D
This project introduces a method to establish a mapping between u-space and Cartesian space. Electronic structure calculations can first be performed in u-space using the plane-wave method. The developed transformation is then used to interpret results in real-space coordinates.
# Files:
1.Benzene: Change of variables in 2D, currently tested with benzene.
2.Single Atom: Change of variables in 3D, currently tested with single atom.
# Scripts:
1. Benzen.py: Tests the coordinate transformation. Once executed, it generates a plot showing the mapping from (v, w) to (x, y). (Noted that it performs better when the atoms are aligned along the y-axis, and less effectively when they are arranged along the x-axis.)
2. Jacobian_determinate.py: Converts calculations from u-space to Cartesian space using the Jacobian determinant. (Functionality not yet fully tested.)
3. Single_Atom.py: Tests the coordinate transformation. Once executed, it generates a plot showing the mapping from (u, v, w) to (x, y, z). (Noted that it performs better when the atoms are aligned along the y-axis and z-axis, and less effectively when they are arranged along the x-axis.)
4. Other scripts: Contain supporting functions used by Benzene.py and Jacobian_determinant.py.
