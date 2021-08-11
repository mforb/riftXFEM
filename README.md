This code is a modification of an XFEM implementation made for academic purposes and distributed at: https://sourceforge.net/projects/xfem/
The description for the the original code is as follows:


The extended finite element method (XFEM) classified, one of the partition of unity method (PUM), allows discontinuities to be simulated independently of the mesh. This is possible by adding appropriate functions to the FE approximation basis, for example, the Heaviside function. The discontinuities can evolve in time, without a need for a conforming mesh. A MATLAB implementation of the XFEM written by VP Nguyen, is given here. The interaction of cracks and crack-inclusion interaction is modelled with XFEM framework. The elements intersected by discontinuity surface are sub-divided into quadrature subcells aligned with the discontinuity and higher order quadrature is adopted.

The implementation is described in the following article:

Meshless methods: a review and computer implementation aspects
VP Nguyen, T Rabczuk, S Bordas, M Duflot, Mathematics and computers in simulation 79 (3), 763-813.

%--------------------------------------------------------------------
If using these codes for research or industrial purposes, please cite:
Title: An extended finite element library
Author(s): BORDAS, S; NGUYEN, PV; DUNANT, C; et al.
Source: INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN ENGINEERING
Volume: 71 Issue: 6 Pages: 703-732 Year: AUG 6 2007 DOI: 10.1002/nme.1966

Title: Numerical integration over arbitrary polygonal domains based on
Schwarz-Christofel conformal maping.
Author(s): S Natarajan; S Bordas; D.R. Mahapatra
Souce: INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN ENGINEERING
DOI: 10.1002/nme.2589

Author(s): Sundar, Stephane
Other contributor(s): Stefano, Phu, David
Modified: Jan 2011
%--------------------------------------------------------------------

Modifications of this code have been made by Martin Forbes (PhD Candidate at the university of Otago)
contact: martin.forbes@postgrad.otago.ac.nz
