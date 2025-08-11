import os; # set up paths
if not ("msRoot" in os.environ):
  print("try again after running:\nsource .../src/bashrc"); exit(-1);
from msrc import *;


runSh('.', "rm -rf dropInCylinder*")

with open("dropInCylinder.mhd", 'w') as f1:
	f1.write("""
DimSize = 30 30 30
Offset =      0    0    0
ElementSize = 1   1   1
Unit =      1
replaceRange 0 255 1
Paint   cylinder 0  15 15   20 15 15  10 0; //< point1(on-axis)  point2  radius value=0 water
PaintAdd  sphere 15 15 15   12  2; //< x_centre y_centre z_centre  radius value=2(oil)
replaceRange    3 255  1; //< range toValue, set outside of cylinder as rock=1
reset s 1e-6
write dump.tif
""");#ElementDataFile = NO_READ

runSh('.', "AllRunContAngle")
# FIXME: Number 55.0 is not analytical
errc =  fileFloatDiffersFrom("AllRunContAngle.log","ContactAngMax", 55.0, 0.2)
errc += fileFloatDiffersFrom("AllRunContAngle.log","ContactAngMin", 55.0, 0.2)
errc += fileFloatDiffersFrom("AllRunContAngle.log","ContactAngAvg", 55.0, 0.1)
exit(errc)
