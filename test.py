import os

if not ("msRoot" in os.environ):
  print("try again after running:\nsource .../src/bashrc"); exit(-1);
from pnmkit.runtime import *;
import image3kit as ik

run_sh('.', "rm -rf dropInCylinder*")

img = ik.VxlImgU8((30, 30,30), 1)
img.replace_range(0, 255, 1)
img.paint(ik.cylinder((0,15,15), (20,15,15), 10, 0))       # point1(on-axis), point2, radius, value=0 water
img.paint_add(ik.sphere((15,15,15), 12, 2))                # x_centre, y_centre, z_centre, radius, value=2(oil)
img.replace_range(3, 255, 1)                            # set outside of cylinder as rock=1
img.spacing = ik.dbl3(1, 1, 1)
img.write("dropInCylinder.mhd")

run_sh('.', "rm -rf  dropInCylinder")
run_sh('.', "AllRunContAngle")
# reference: results of the original code (porescale b37dfe1), not analytical
errc =  read_deviates("AllRunContAngle.log","ContactAngMax", 62.7336, 0.01)
errc += read_deviates("AllRunContAngle.log","ContactAngMin", 46.3434, 0.01)
errc += read_deviates("AllRunContAngle.log","ContactAngAvg", 55.0465, 0.01)
exit(errc)
