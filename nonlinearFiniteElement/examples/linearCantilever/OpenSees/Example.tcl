set ElementNumberX 2900
set ElementNumberY 290
###########################################################################
# OpenSees Example: Cantilever beam modeled with 2D quadrilateral elements
###########################################################################

#-----------------------------
# Define the model
# ----------------------------

# Create ModelBuilder with 2 dimensions and 2 DOF/node
#---------------------
model BasicBuilder -ndm 2 -ndf 2

# create the material
#---------------------
nDMaterial ElasticIsotropic 1 100000 0.25 0.0

# dimensions
#---------------------
set dimensionX 4000.0
set dimensionY 400.0

set incrementX [expr $dimensionX/$ElementNumberX];
set incrementY [expr $dimensionY/$ElementNumberY];

set l1 [expr ($ElementNumberX + 1)*($ElementNumberY + 1)]
set l2 [expr ($ElementNumberX + 1)]

set counter 1

# define nodes
#---------------------
for {set j 0} {$j <= $ElementNumberY} {incr j 1} {
	for {set i 0} {$i <= $ElementNumberX} {incr i 1} {
		node $counter [expr $incrementX*$i] [expr $incrementY*$j]
		set counter [expr $counter + 1 ]
	}
}

# define boundary conditions
#---------------------
for {set i 0} {$i <= $ElementNumberY} {incr i 1} {
	fix [expr $i*($ElementNumberX+1)+1] 1 1
}

# define element (quad)
#---------------------
set counter 1
set xElement [expr $ElementNumberX + 1]
for {set i 0} {$i < $ElementNumberY} {incr i 1} {
	for {set j 0} {$j < $ElementNumberX} {incr j 1} {
		element quad $counter [expr $i*$xElement+$j + 1] [expr $i*$xElement+1+$j + 1] [expr ($i+1)*$xElement+1+$j + 1] [expr ($i+1)*$xElement+$j + 1] 300 "PlaneStress2D" 1
		set counter [expr $counter + 1 ]
	}
}

# define load pattern
#---------------------
set loadVal [expr 1000.0/($ElementNumberY*2)]
pattern Plain 1 Linear {
		load 1 0.0 $loadVal
		load $l1 0.0 $loadVal
	for {set i 2} {$i < [expr $ElementNumberY+1]} {incr i 1} {
		load [expr $i*($ElementNumberX+1)] 0.0 [expr 2*$loadVal]
	}
}

# define the recorder
#---------------------
recorder Node -file Data/Node.out -node $l1 -dof 2 disp


# --------------------------------------------------------------------
# Start of static analysis (creation of the analysis & analysis itself)
# --------------------------------------------------------------------

# Load control with variable load steps
#                      init Jd min max
integrator LoadControl 1.0 1 1.0 10.0

# Convergence test
#              tolerance maxIter displayCode
test EnergyIncr 1.0e-12    10         0

# Solution algorithm
algorithm Newton

# DOF numberer
numberer RCM

# Cosntraint handler
constraints Plain

# System of equations solver
system ProfileSPD

# Type of analysis analysis
analysis Static

# Perform the analysis
analyze 1

# --------------------------
# End of static analysis
# --------------------------
