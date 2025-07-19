# Change of variables in 2D, currently only benzene available
This project introduces a method to establish a mapping between u-space and Cartesian space. Electronic structure calculations can first be performed in u-space using the plane-wave method. The developed transformation is then used to interpret results in real-space coordinates.
# Scripts:
1. Benzen.py: Tests the coordinate transformation. Once executed, it generates a plot showing the mapping from (v, w) to (x, y).
2. Jacobian_determinate.py: Converts calculations from u-space to Cartesian space using the Jacobian determinant. (Functionality not yet fully tested.)
3. Other scripts: Contain supporting functions used by Benzene.py and Jacobian_determinant.py.
