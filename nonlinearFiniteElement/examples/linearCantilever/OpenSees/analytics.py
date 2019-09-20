# Units 
mm = 1.0
MPa = 1.0
N = 1.0
# parameters
length = 4000.0*mm
hight = 400.0*mm
width = 300*mm
elasticModulus = 100000.0*MPa
load = 1000.0*N
# Analyse
moment_of_inertia = width*(hight**3)/12.0
tip_displacement = load * (length ** 3 ) / (3 * elasticModulus * moment_of_inertia)
print("Tip disp:", round(tip_displacement,5), "mm")

