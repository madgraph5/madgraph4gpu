# This file was automatically created by FeynRules 2.3.35
# Mathematica version: 12.1.0 for Linux x86 (64-bit) (March 18, 2020)
# Date: Thu 7 Jan 2021 13:56:43


from object_library import all_vertices, Vertex
import particles as P
import couplings as C
import lorentz as L


V_1 = Vertex(name = 'V_1',
             particles = [ P.a, P.W__minus__, P.W__plus__ ],
             color = [ '1' ],
             lorentz = [ L.VVV1, L.VVV3, L.VVV5, L.VVV6, L.VVV7, L.VVV9 ],
             couplings = {(0,1):C.GC_471,(0,0):C.GC_533,(0,5):C.GC_283,(0,4):C.GC_282,(0,2):C.GC_3,(0,3):C.GC_532})

V_2 = Vertex(name = 'V_2',
             particles = [ P.a, P.W__minus__, P.W__plus__ ],
             color = [ '1' ],
             lorentz = [ L.VVV5 ],
             couplings = {(0,0):C.GC_465})

V_3 = Vertex(name = 'V_3',
             particles = [ P.a, P.W__minus__, P.W__plus__ ],
             color = [ '1' ],
             lorentz = [ L.VVV5 ],
             couplings = {(0,0):C.GC_474})

V_4 = Vertex(name = 'V_4',
             particles = [ P.a, P.W__minus__, P.W__plus__ ],
             color = [ '1' ],
             lorentz = [ L.VVV5 ],
             couplings = {(0,0):C.GC_486})

V_5 = Vertex(name = 'V_5',
             particles = [ P.W__minus__, P.W__plus__, P.Z ],
             color = [ '1' ],
             lorentz = [ L.VVV2, L.VVV3, L.VVV4, L.VVV5, L.VVV7, L.VVV9 ],
             couplings = {(0,1):C.GC_534,(0,0):C.GC_470,(0,5):C.GC_128,(0,4):C.GC_127,(0,3):C.GC_205,(0,2):C.GC_468})

V_6 = Vertex(name = 'V_6',
             particles = [ P.W__minus__, P.W__plus__, P.Z ],
             color = [ '1' ],
             lorentz = [ L.VVV5 ],
             couplings = {(0,0):C.GC_528})

V_7 = Vertex(name = 'V_7',
             particles = [ P.W__minus__, P.W__plus__, P.Z ],
             color = [ '1' ],
             lorentz = [ L.VVV5 ],
             couplings = {(0,0):C.GC_529})

V_8 = Vertex(name = 'V_8',
             particles = [ P.W__minus__, P.W__plus__, P.Z ],
             color = [ '1' ],
             lorentz = [ L.VVV5 ],
             couplings = {(0,0):C.GC_535})

V_9 = Vertex(name = 'V_9',
             particles = [ P.g, P.g, P.g ],
             color = [ 'f(1,2,3)' ],
             lorentz = [ L.VVV3, L.VVV5, L.VVV7, L.VVV8 ],
             couplings = {(0,0):C.GC_478,(0,3):C.GC_32,(0,2):C.GC_31,(0,1):C.GC_7})

V_10 = Vertex(name = 'V_10',
              particles = [ P.g, P.g, P.g, P.g ],
              color = [ 'f(-1,1,2)*f(3,4,-1)', 'f(-1,1,3)*f(2,4,-1)', 'f(-1,1,4)*f(2,3,-1)' ],
              lorentz = [ L.VVVV1, L.VVVV10, L.VVVV12, L.VVVV13, L.VVVV2, L.VVVV3, L.VVVV4, L.VVVV6, L.VVVV9 ],
              couplings = {(0,7):C.GC_147,(1,6):C.GC_147,(2,5):C.GC_147,(0,4):C.GC_146,(1,3):C.GC_146,(2,2):C.GC_146,(1,8):C.GC_8,(0,0):C.GC_8,(2,1):C.GC_8})

V_11 = Vertex(name = 'V_11',
              particles = [ P.g, P.g, P.g, P.g, P.g ],
              color = [ 'f(-2,1,2)*f(-1,-2,3)*f(4,5,-1)', 'f(-2,1,2)*f(-1,-2,4)*f(3,5,-1)', 'f(-2,1,2)*f(-1,-2,5)*f(3,4,-1)', 'f(-2,1,3)*f(-1,-2,2)*f(4,5,-1)', 'f(-2,1,3)*f(-1,-2,4)*f(2,5,-1)', 'f(-2,1,3)*f(-1,-2,5)*f(2,4,-1)', 'f(-2,1,4)*f(-1,-2,2)*f(3,5,-1)', 'f(-2,1,4)*f(-1,-2,3)*f(2,5,-1)', 'f(-2,1,4)*f(-1,-2,5)*f(2,3,-1)', 'f(-2,1,5)*f(-1,-2,2)*f(3,4,-1)', 'f(-2,1,5)*f(-1,-2,3)*f(2,4,-1)', 'f(-2,1,5)*f(-1,-2,4)*f(2,3,-1)', 'f(-2,2,3)*f(-1,-2,1)*f(4,5,-1)', 'f(-2,2,3)*f(-1,-2,4)*f(1,5,-1)', 'f(-2,2,3)*f(-1,-2,5)*f(1,4,-1)', 'f(-2,2,4)*f(-1,-2,1)*f(3,5,-1)', 'f(-2,2,4)*f(-1,-2,3)*f(1,5,-1)', 'f(-2,2,4)*f(-1,-2,5)*f(1,3,-1)', 'f(-2,2,5)*f(-1,-2,1)*f(3,4,-1)', 'f(-2,2,5)*f(-1,-2,3)*f(1,4,-1)', 'f(-2,2,5)*f(-1,-2,4)*f(1,3,-1)', 'f(-2,3,4)*f(-1,-2,1)*f(2,5,-1)', 'f(-2,3,4)*f(-1,-2,2)*f(1,5,-1)', 'f(-2,3,4)*f(-1,-2,5)*f(1,2,-1)', 'f(-2,3,5)*f(-1,-2,1)*f(2,4,-1)', 'f(-2,3,5)*f(-1,-2,2)*f(1,4,-1)', 'f(-2,3,5)*f(-1,-2,4)*f(1,2,-1)', 'f(-2,4,5)*f(-1,-2,1)*f(2,3,-1)', 'f(-2,4,5)*f(-1,-2,2)*f(1,3,-1)', 'f(-2,4,5)*f(-1,-2,3)*f(1,2,-1)' ],
              lorentz = [ L.VVVVV1, L.VVVVV10, L.VVVVV11, L.VVVVV12, L.VVVVV13, L.VVVVV14, L.VVVVV15, L.VVVVV17, L.VVVVV18, L.VVVVV19, L.VVVVV2, L.VVVVV20, L.VVVVV21, L.VVVVV22, L.VVVVV23, L.VVVVV24, L.VVVVV25, L.VVVVV28, L.VVVVV29, L.VVVVV3, L.VVVVV30, L.VVVVV31, L.VVVVV33, L.VVVVV34, L.VVVVV35, L.VVVVV36, L.VVVVV37, L.VVVVV4, L.VVVVV40, L.VVVVV41, L.VVVVV42, L.VVVVV43, L.VVVVV44, L.VVVVV46, L.VVVVV47, L.VVVVV48, L.VVVVV49, L.VVVVV5, L.VVVVV50, L.VVVVV51, L.VVVVV53, L.VVVVV54, L.VVVVV6, L.VVVVV7, L.VVVVV9 ],
              couplings = {(27,37):C.GC_154,(24,8):C.GC_155,(21,12):C.GC_154,(18,11):C.GC_155,(15,9):C.GC_154,(12,27):C.GC_154,(28,42):C.GC_154,(25,15):C.GC_155,(22,14):C.GC_154,(9,16):C.GC_154,(6,13):C.GC_155,(3,43):C.GC_155,(29,0):C.GC_155,(19,20):C.GC_154,(16,18):C.GC_155,(10,17):C.GC_154,(7,21):C.GC_155,(0,44):C.GC_154,(26,10):C.GC_154,(20,1):C.GC_155,(13,24):C.GC_154,(11,2):C.GC_155,(4,22):C.GC_154,(1,23):C.GC_154,(23,19):C.GC_155,(17,4):C.GC_154,(14,25):C.GC_155,(8,3):C.GC_154,(5,28):C.GC_155,(2,26):C.GC_155,(24,29):C.GC_153,(21,30):C.GC_152,(18,30):C.GC_153,(15,29):C.GC_152,(28,6):C.GC_153,(22,34):C.GC_153,(9,34):C.GC_152,(3,6):C.GC_152,(29,7):C.GC_153,(16,35):C.GC_153,(10,35):C.GC_152,(0,7):C.GC_152,(26,39):C.GC_152,(20,38):C.GC_152,(4,38):C.GC_153,(1,39):C.GC_153,(25,33):C.GC_153,(6,33):C.GC_152,(19,36):C.GC_153,(7,36):C.GC_152,(23,41):C.GC_152,(17,40):C.GC_152,(5,40):C.GC_153,(2,41):C.GC_153,(27,5):C.GC_153,(12,5):C.GC_152,(13,31):C.GC_153,(11,31):C.GC_152,(14,32):C.GC_152,(8,32):C.GC_153})

V_12 = Vertex(name = 'V_12',
              particles = [ P.g, P.g, P.g, P.g, P.g, P.g ],
              color = [ 'f(-3,1,2)*f(-2,3,4)*f(-1,-2,-3)*f(5,6,-1)', 'f(-3,1,2)*f(-2,3,5)*f(-1,-2,-3)*f(4,6,-1)', 'f(-3,1,2)*f(-2,3,6)*f(-1,-2,-3)*f(4,5,-1)', 'f(-3,1,2)*f(-2,4,5)*f(-1,-2,-3)*f(3,6,-1)', 'f(-3,1,2)*f(-2,4,6)*f(-1,-2,-3)*f(3,5,-1)', 'f(-3,1,2)*f(-2,5,6)*f(-1,-2,-3)*f(3,4,-1)', 'f(-3,1,3)*f(-2,2,4)*f(-1,-2,-3)*f(5,6,-1)', 'f(-3,1,3)*f(-2,2,5)*f(-1,-2,-3)*f(4,6,-1)', 'f(-3,1,3)*f(-2,2,6)*f(-1,-2,-3)*f(4,5,-1)', 'f(-3,1,3)*f(-2,4,5)*f(-1,-2,-3)*f(2,6,-1)', 'f(-3,1,3)*f(-2,4,6)*f(-1,-2,-3)*f(2,5,-1)', 'f(-3,1,3)*f(-2,5,6)*f(-1,-2,-3)*f(2,4,-1)', 'f(-3,1,4)*f(-2,2,3)*f(-1,-2,-3)*f(5,6,-1)', 'f(-3,1,4)*f(-2,2,5)*f(-1,-2,-3)*f(3,6,-1)', 'f(-3,1,4)*f(-2,2,6)*f(-1,-2,-3)*f(3,5,-1)', 'f(-3,1,4)*f(-2,3,5)*f(-1,-2,-3)*f(2,6,-1)', 'f(-3,1,4)*f(-2,3,6)*f(-1,-2,-3)*f(2,5,-1)', 'f(-3,1,4)*f(-2,5,6)*f(-1,-2,-3)*f(2,3,-1)', 'f(-3,1,5)*f(-2,2,3)*f(-1,-2,-3)*f(4,6,-1)', 'f(-3,1,5)*f(-2,2,4)*f(-1,-2,-3)*f(3,6,-1)', 'f(-3,1,5)*f(-2,2,6)*f(-1,-2,-3)*f(3,4,-1)', 'f(-3,1,5)*f(-2,3,4)*f(-1,-2,-3)*f(2,6,-1)', 'f(-3,1,5)*f(-2,3,6)*f(-1,-2,-3)*f(2,4,-1)', 'f(-3,1,5)*f(-2,4,6)*f(-1,-2,-3)*f(2,3,-1)', 'f(-3,1,6)*f(-2,2,3)*f(-1,-2,-3)*f(4,5,-1)', 'f(-3,1,6)*f(-2,2,4)*f(-1,-2,-3)*f(3,5,-1)', 'f(-3,1,6)*f(-2,2,5)*f(-1,-2,-3)*f(3,4,-1)', 'f(-3,1,6)*f(-2,3,4)*f(-1,-2,-3)*f(2,5,-1)', 'f(-3,1,6)*f(-2,3,5)*f(-1,-2,-3)*f(2,4,-1)', 'f(-3,1,6)*f(-2,4,5)*f(-1,-2,-3)*f(2,3,-1)', 'f(-3,2,3)*f(-2,1,4)*f(-1,-2,-3)*f(5,6,-1)', 'f(-3,2,3)*f(-2,1,5)*f(-1,-2,-3)*f(4,6,-1)', 'f(-3,2,3)*f(-2,1,6)*f(-1,-2,-3)*f(4,5,-1)', 'f(-3,2,3)*f(-2,4,5)*f(-1,-2,-3)*f(1,6,-1)', 'f(-3,2,3)*f(-2,4,6)*f(-1,-2,-3)*f(1,5,-1)', 'f(-3,2,3)*f(-2,5,6)*f(-1,-2,-3)*f(1,4,-1)', 'f(-3,2,4)*f(-2,1,3)*f(-1,-2,-3)*f(5,6,-1)', 'f(-3,2,4)*f(-2,1,5)*f(-1,-2,-3)*f(3,6,-1)', 'f(-3,2,4)*f(-2,1,6)*f(-1,-2,-3)*f(3,5,-1)', 'f(-3,2,4)*f(-2,3,5)*f(-1,-2,-3)*f(1,6,-1)', 'f(-3,2,4)*f(-2,3,6)*f(-1,-2,-3)*f(1,5,-1)', 'f(-3,2,4)*f(-2,5,6)*f(-1,-2,-3)*f(1,3,-1)', 'f(-3,2,5)*f(-2,1,3)*f(-1,-2,-3)*f(4,6,-1)', 'f(-3,2,5)*f(-2,1,4)*f(-1,-2,-3)*f(3,6,-1)', 'f(-3,2,5)*f(-2,1,6)*f(-1,-2,-3)*f(3,4,-1)', 'f(-3,2,5)*f(-2,3,4)*f(-1,-2,-3)*f(1,6,-1)', 'f(-3,2,5)*f(-2,3,6)*f(-1,-2,-3)*f(1,4,-1)', 'f(-3,2,5)*f(-2,4,6)*f(-1,-2,-3)*f(1,3,-1)', 'f(-3,2,6)*f(-2,1,3)*f(-1,-2,-3)*f(4,5,-1)', 'f(-3,2,6)*f(-2,1,4)*f(-1,-2,-3)*f(3,5,-1)', 'f(-3,2,6)*f(-2,1,5)*f(-1,-2,-3)*f(3,4,-1)', 'f(-3,2,6)*f(-2,3,4)*f(-1,-2,-3)*f(1,5,-1)', 'f(-3,2,6)*f(-2,3,5)*f(-1,-2,-3)*f(1,4,-1)', 'f(-3,2,6)*f(-2,4,5)*f(-1,-2,-3)*f(1,3,-1)', 'f(-3,3,4)*f(-2,1,2)*f(-1,-2,-3)*f(5,6,-1)', 'f(-3,3,4)*f(-2,1,5)*f(-1,-2,-3)*f(2,6,-1)', 'f(-3,3,4)*f(-2,1,6)*f(-1,-2,-3)*f(2,5,-1)', 'f(-3,3,4)*f(-2,2,5)*f(-1,-2,-3)*f(1,6,-1)', 'f(-3,3,4)*f(-2,2,6)*f(-1,-2,-3)*f(1,5,-1)', 'f(-3,3,4)*f(-2,5,6)*f(-1,-2,-3)*f(1,2,-1)', 'f(-3,3,5)*f(-2,1,2)*f(-1,-2,-3)*f(4,6,-1)', 'f(-3,3,5)*f(-2,1,4)*f(-1,-2,-3)*f(2,6,-1)', 'f(-3,3,5)*f(-2,1,6)*f(-1,-2,-3)*f(2,4,-1)', 'f(-3,3,5)*f(-2,2,4)*f(-1,-2,-3)*f(1,6,-1)', 'f(-3,3,5)*f(-2,2,6)*f(-1,-2,-3)*f(1,4,-1)', 'f(-3,3,5)*f(-2,4,6)*f(-1,-2,-3)*f(1,2,-1)', 'f(-3,3,6)*f(-2,1,2)*f(-1,-2,-3)*f(4,5,-1)', 'f(-3,3,6)*f(-2,1,4)*f(-1,-2,-3)*f(2,5,-1)', 'f(-3,3,6)*f(-2,1,5)*f(-1,-2,-3)*f(2,4,-1)', 'f(-3,3,6)*f(-2,2,4)*f(-1,-2,-3)*f(1,5,-1)', 'f(-3,3,6)*f(-2,2,5)*f(-1,-2,-3)*f(1,4,-1)', 'f(-3,3,6)*f(-2,4,5)*f(-1,-2,-3)*f(1,2,-1)', 'f(-3,4,5)*f(-2,1,2)*f(-1,-2,-3)*f(3,6,-1)', 'f(-3,4,5)*f(-2,1,3)*f(-1,-2,-3)*f(2,6,-1)', 'f(-3,4,5)*f(-2,1,6)*f(-1,-2,-3)*f(2,3,-1)', 'f(-3,4,5)*f(-2,2,3)*f(-1,-2,-3)*f(1,6,-1)', 'f(-3,4,5)*f(-2,2,6)*f(-1,-2,-3)*f(1,3,-1)', 'f(-3,4,5)*f(-2,3,6)*f(-1,-2,-3)*f(1,2,-1)', 'f(-3,4,6)*f(-2,1,2)*f(-1,-2,-3)*f(3,5,-1)', 'f(-3,4,6)*f(-2,1,3)*f(-1,-2,-3)*f(2,5,-1)', 'f(-3,4,6)*f(-2,1,5)*f(-1,-2,-3)*f(2,3,-1)', 'f(-3,4,6)*f(-2,2,3)*f(-1,-2,-3)*f(1,5,-1)', 'f(-3,4,6)*f(-2,2,5)*f(-1,-2,-3)*f(1,3,-1)', 'f(-3,4,6)*f(-2,3,5)*f(-1,-2,-3)*f(1,2,-1)', 'f(-3,5,6)*f(-2,1,2)*f(-1,-2,-3)*f(3,4,-1)', 'f(-3,5,6)*f(-2,1,3)*f(-1,-2,-3)*f(2,4,-1)', 'f(-3,5,6)*f(-2,1,4)*f(-1,-2,-3)*f(2,3,-1)', 'f(-3,5,6)*f(-2,2,3)*f(-1,-2,-3)*f(1,4,-1)', 'f(-3,5,6)*f(-2,2,4)*f(-1,-2,-3)*f(1,3,-1)', 'f(-3,5,6)*f(-2,3,4)*f(-1,-2,-3)*f(1,2,-1)' ],
              lorentz = [ L.VVVVVV1, L.VVVVVV10, L.VVVVVV11, L.VVVVVV12, L.VVVVVV13, L.VVVVVV14, L.VVVVVV15, L.VVVVVV16, L.VVVVVV17, L.VVVVVV18, L.VVVVVV19, L.VVVVVV2, L.VVVVVV20, L.VVVVVV21, L.VVVVVV22, L.VVVVVV23, L.VVVVVV24, L.VVVVVV25, L.VVVVVV26, L.VVVVVV27, L.VVVVVV28, L.VVVVVV29, L.VVVVVV3, L.VVVVVV30, L.VVVVVV31, L.VVVVVV32, L.VVVVVV33, L.VVVVVV34, L.VVVVVV35, L.VVVVVV36, L.VVVVVV37, L.VVVVVV38, L.VVVVVV39, L.VVVVVV4, L.VVVVVV40, L.VVVVVV41, L.VVVVVV42, L.VVVVVV43, L.VVVVVV44, L.VVVVVV45, L.VVVVVV46, L.VVVVVV47, L.VVVVVV48, L.VVVVVV49, L.VVVVVV5, L.VVVVVV50, L.VVVVVV51, L.VVVVVV52, L.VVVVVV54, L.VVVVVV55, L.VVVVVV56, L.VVVVVV57, L.VVVVVV58, L.VVVVVV59, L.VVVVVV6, L.VVVVVV60, L.VVVVVV61, L.VVVVVV7, L.VVVVVV8, L.VVVVVV9 ],
              couplings = {(41,58):C.GC_160,(47,59):C.GC_159,(53,7):C.GC_160,(35,57):C.GC_159,(46,14):C.GC_160,(52,17):C.GC_159,(34,2):C.GC_160,(40,10):C.GC_159,(51,37):C.GC_160,(33,4):C.GC_159,(39,21):C.GC_160,(45,30):C.GC_159,(17,57):C.GC_160,(23,2):C.GC_159,(29,4):C.GC_160,(11,58):C.GC_159,(22,10):C.GC_160,(28,21):C.GC_159,(10,59):C.GC_160,(16,14):C.GC_159,(27,30):C.GC_160,(9,7):C.GC_159,(15,17):C.GC_160,(21,37):C.GC_159,(59,0):C.GC_160,(65,11):C.GC_159,(71,44):C.GC_160,(64,12):C.GC_160,(70,23):C.GC_159,(58,16):C.GC_159,(69,31):C.GC_160,(57,19):C.GC_160,(63,39):C.GC_159,(5,0):C.GC_159,(20,16):C.GC_160,(26,19):C.GC_159,(4,11):C.GC_160,(14,12):C.GC_159,(25,39):C.GC_160,(3,44):C.GC_159,(13,23):C.GC_160,(19,31):C.GC_159,(77,22):C.GC_159,(83,33):C.GC_160,(76,1):C.GC_160,(82,8):C.GC_159,(81,40):C.GC_160,(75,35):C.GC_159,(2,22):C.GC_160,(8,1):C.GC_159,(24,35):C.GC_160,(1,33):C.GC_159,(7,8):C.GC_160,(18,40):C.GC_159,(89,54):C.GC_160,(88,6):C.GC_159,(87,25):C.GC_160,(0,54):C.GC_159,(6,6):C.GC_160,(12,25):C.GC_159,(62,15):C.GC_160,(68,18):C.GC_159,(56,13):C.GC_159,(67,38):C.GC_160,(55,24):C.GC_160,(61,32):C.GC_159,(44,13):C.GC_160,(50,24):C.GC_159,(38,15):C.GC_159,(49,32):C.GC_160,(37,18):C.GC_160,(43,38):C.GC_159,(74,3):C.GC_160,(80,5):C.GC_159,(79,34):C.GC_160,(73,42):C.GC_159,(32,3):C.GC_159,(48,42):C.GC_160,(31,5):C.GC_160,(42,34):C.GC_159,(86,9):C.GC_159,(85,20):C.GC_160,(30,9):C.GC_160,(36,20):C.GC_159,(78,41):C.GC_160,(72,36):C.GC_159,(66,36):C.GC_160,(60,41):C.GC_159,(65,43):C.GC_157,(71,46):C.GC_158,(77,46):C.GC_157,(83,43):C.GC_158,(41,28):C.GC_157,(53,50):C.GC_157,(76,50):C.GC_158,(88,28):C.GC_158,(35,29):C.GC_157,(52,53):C.GC_157,(64,53):C.GC_158,(87,29):C.GC_158,(34,52):C.GC_158,(40,51):C.GC_158,(69,51):C.GC_157,(81,52):C.GC_157,(17,29):C.GC_158,(23,52):C.GC_157,(80,52):C.GC_158,(86,29):C.GC_157,(11,28):C.GC_158,(22,51):C.GC_157,(68,51):C.GC_158,(85,28):C.GC_157,(9,50):C.GC_158,(15,53):C.GC_158,(61,53):C.GC_157,(73,50):C.GC_157,(4,43):C.GC_158,(14,53):C.GC_157,(49,53):C.GC_158,(78,43):C.GC_157,(3,46):C.GC_157,(19,51):C.GC_158,(37,51):C.GC_157,(72,46):C.GC_158,(2,46):C.GC_158,(8,50):C.GC_157,(48,50):C.GC_158,(66,46):C.GC_157,(1,43):C.GC_157,(18,52):C.GC_158,(31,52):C.GC_157,(60,43):C.GC_158,(6,28):C.GC_157,(12,29):C.GC_157,(30,29):C.GC_158,(36,28):C.GC_158,(47,48):C.GC_157,(82,48):C.GC_158,(46,55):C.GC_157,(70,55):C.GC_158,(33,56):C.GC_158,(39,49):C.GC_158,(63,49):C.GC_157,(75,56):C.GC_157,(29,56):C.GC_157,(74,56):C.GC_158,(28,49):C.GC_157,(62,49):C.GC_158,(10,48):C.GC_158,(16,55):C.GC_158,(67,55):C.GC_157,(79,48):C.GC_157,(25,49):C.GC_158,(38,49):C.GC_157,(13,55):C.GC_157,(43,55):C.GC_158,(24,56):C.GC_158,(32,56):C.GC_157,(7,48):C.GC_157,(42,48):C.GC_158,(84,26):C.GC_160,(54,26):C.GC_159,(59,27):C.GC_157,(89,27):C.GC_158,(51,45):C.GC_157,(58,45):C.GC_158,(21,45):C.GC_158,(55,45):C.GC_157,(5,27):C.GC_158,(20,45):C.GC_157,(50,45):C.GC_158,(84,27):C.GC_157,(0,27):C.GC_157,(54,27):C.GC_158,(45,47):C.GC_158,(57,47):C.GC_157,(27,47):C.GC_157,(56,47):C.GC_158,(26,47):C.GC_158,(44,47):C.GC_157})

V_13 = Vertex(name = 'V_13',
              particles = [ P.a, P.a, P.W__minus__, P.W__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVVV15, L.VVVV5, L.VVVV8 ],
              couplings = {(0,1):C.GC_285,(0,0):C.GC_284,(0,2):C.GC_5})

V_14 = Vertex(name = 'V_14',
              particles = [ P.a, P.a, P.W__minus__, P.W__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVVV8 ],
              couplings = {(0,0):C.GC_476})

V_15 = Vertex(name = 'V_15',
              particles = [ P.a, P.a, P.W__minus__, P.W__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVVV8 ],
              couplings = {(0,0):C.GC_477})

V_16 = Vertex(name = 'V_16',
              particles = [ P.a, P.a, P.W__minus__, P.W__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVVV8 ],
              couplings = {(0,0):C.GC_497})

V_17 = Vertex(name = 'V_17',
              particles = [ P.a, P.a, P.W__minus__, P.W__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVVV8 ],
              couplings = {(0,0):C.GC_537})

V_18 = Vertex(name = 'V_18',
              particles = [ P.a, P.W__minus__, P.W__plus__, P.Z ],
              color = [ '1' ],
              lorentz = [ L.VVVV11, L.VVVV14, L.VVVV7 ],
              couplings = {(0,2):C.GC_140,(0,1):C.GC_139,(0,0):C.GC_206})

V_19 = Vertex(name = 'V_19',
              particles = [ P.a, P.W__minus__, P.W__plus__, P.Z ],
              color = [ '1' ],
              lorentz = [ L.VVVV11 ],
              couplings = {(0,0):C.GC_564})

V_20 = Vertex(name = 'V_20',
              particles = [ P.a, P.W__minus__, P.W__plus__, P.Z ],
              color = [ '1' ],
              lorentz = [ L.VVVV11 ],
              couplings = {(0,0):C.GC_565})

V_21 = Vertex(name = 'V_21',
              particles = [ P.a, P.W__minus__, P.W__plus__, P.Z ],
              color = [ '1' ],
              lorentz = [ L.VVVV11 ],
              couplings = {(0,0):C.GC_536})

V_22 = Vertex(name = 'V_22',
              particles = [ P.a, P.W__minus__, P.W__plus__, P.Z ],
              color = [ '1' ],
              lorentz = [ L.VVVV11 ],
              couplings = {(0,0):C.GC_539})

V_23 = Vertex(name = 'V_23',
              particles = [ P.a, P.a, P.a, P.W__minus__, P.W__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVVVV16, L.VVVVV8 ],
              couplings = {(0,1):C.GC_287,(0,0):C.GC_286})

V_24 = Vertex(name = 'V_24',
              particles = [ P.W__minus__, P.W__minus__, P.W__plus__, P.W__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVVV15, L.VVVV5, L.VVVV8 ],
              couplings = {(0,1):C.GC_243,(0,0):C.GC_241,(0,2):C.GC_175})

V_25 = Vertex(name = 'V_25',
              particles = [ P.W__minus__, P.W__minus__, P.W__plus__, P.W__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVVV8 ],
              couplings = {(0,0):C.GC_490})

V_26 = Vertex(name = 'V_26',
              particles = [ P.W__minus__, P.W__minus__, P.W__plus__, P.W__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVVV8 ],
              couplings = {(0,0):C.GC_492})

V_27 = Vertex(name = 'V_27',
              particles = [ P.a, P.a, P.W__minus__, P.W__plus__, P.Z ],
              color = [ '1' ],
              lorentz = [ L.VVVVV26, L.VVVVV55 ],
              couplings = {(0,0):C.GC_143,(0,1):C.GC_142})

V_28 = Vertex(name = 'V_28',
              particles = [ P.a, P.W__minus__, P.W__minus__, P.W__plus__, P.W__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVVVV32, L.VVVVV57 ],
              couplings = {(0,0):C.GC_248,(0,1):C.GC_246})

V_29 = Vertex(name = 'V_29',
              particles = [ P.a, P.a, P.W__minus__, P.W__minus__, P.W__plus__, P.W__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVVVVV53 ],
              couplings = {(0,0):C.GC_250})

V_30 = Vertex(name = 'V_30',
              particles = [ P.W__minus__, P.W__plus__, P.Z, P.Z ],
              color = [ '1' ],
              lorentz = [ L.VVVV15, L.VVVV5, L.VVVV8 ],
              couplings = {(0,1):C.GC_244,(0,0):C.GC_242,(0,2):C.GC_177})

V_31 = Vertex(name = 'V_31',
              particles = [ P.W__minus__, P.W__plus__, P.Z, P.Z ],
              color = [ '1' ],
              lorentz = [ L.VVVV8 ],
              couplings = {(0,0):C.GC_498})

V_32 = Vertex(name = 'V_32',
              particles = [ P.W__minus__, P.W__plus__, P.Z, P.Z ],
              color = [ '1' ],
              lorentz = [ L.VVVV8 ],
              couplings = {(0,0):C.GC_499})

V_33 = Vertex(name = 'V_33',
              particles = [ P.W__minus__, P.W__plus__, P.Z, P.Z ],
              color = [ '1' ],
              lorentz = [ L.VVVV8 ],
              couplings = {(0,0):C.GC_500})

V_34 = Vertex(name = 'V_34',
              particles = [ P.W__minus__, P.W__plus__, P.Z, P.Z ],
              color = [ '1' ],
              lorentz = [ L.VVVV8 ],
              couplings = {(0,0):C.GC_538})

V_35 = Vertex(name = 'V_35',
              particles = [ P.a, P.W__minus__, P.W__plus__, P.Z, P.Z ],
              color = [ '1' ],
              lorentz = [ L.VVVVV38, L.VVVVV45 ],
              couplings = {(0,0):C.GC_249,(0,1):C.GC_247})

V_36 = Vertex(name = 'V_36',
              particles = [ P.W__minus__, P.W__minus__, P.W__plus__, P.W__plus__, P.Z ],
              color = [ '1' ],
              lorentz = [ L.VVVVV27, L.VVVVV52 ],
              couplings = {(0,0):C.GC_183,(0,1):C.GC_181})

V_37 = Vertex(name = 'V_37',
              particles = [ P.a, P.W__minus__, P.W__minus__, P.W__plus__, P.W__plus__, P.Z ],
              color = [ '1' ],
              lorentz = [ L.VVVVVV62 ],
              couplings = {(0,0):C.GC_185})

V_38 = Vertex(name = 'V_38',
              particles = [ P.W__minus__, P.W__minus__, P.W__plus__, P.W__plus__, P.Z, P.Z ],
              color = [ '1' ],
              lorentz = [ L.VVVVVV53 ],
              couplings = {(0,0):C.GC_173})

V_39 = Vertex(name = 'V_39',
              particles = [ P.W__minus__, P.W__plus__, P.Z, P.Z, P.Z ],
              color = [ '1' ],
              lorentz = [ L.VVVVV39, L.VVVVV56 ],
              couplings = {(0,0):C.GC_184,(0,1):C.GC_182})

V_40 = Vertex(name = 'V_40',
              particles = [ P.a, P.a, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS2, L.VVSS4 ],
              couplings = {(0,0):C.GC_107,(0,1):C.GC_106})

V_41 = Vertex(name = 'V_41',
              particles = [ P.a, P.a, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS2, L.VVSS4 ],
              couplings = {(0,0):C.GC_278,(0,1):C.GC_275})

V_42 = Vertex(name = 'V_42',
              particles = [ P.a, P.a, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS2, L.VVSS4 ],
              couplings = {(0,0):C.GC_296,(0,1):C.GC_295})

V_43 = Vertex(name = 'V_43',
              particles = [ P.a, P.a, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVS2, L.VVS4 ],
              couplings = {(0,0):C.GC_360,(0,1):C.GC_299})

V_44 = Vertex(name = 'V_44',
              particles = [ P.a, P.a, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVS2, L.VVS4 ],
              couplings = {(0,0):C.GC_449,(0,1):C.GC_359})

V_45 = Vertex(name = 'V_45',
              particles = [ P.a, P.a, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVS2, L.VVS4 ],
              couplings = {(0,0):C.GC_456,(0,1):C.GC_446})

V_46 = Vertex(name = 'V_46',
              particles = [ P.a, P.a, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVS4 ],
              couplings = {(0,0):C.GC_455})

V_47 = Vertex(name = 'V_47',
              particles = [ P.g, P.g, P.H, P.H ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.VVSS2, L.VVSS4 ],
              couplings = {(0,0):C.GC_37,(0,1):C.GC_36})

V_48 = Vertex(name = 'V_48',
              particles = [ P.g, P.g, P.H ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.VVS2, L.VVS4, L.VVS5, L.VVS6, L.VVS7 ],
              couplings = {(0,0):C.GC_348,(0,1):C.GC_300,(0,3):C.GC_313,(0,2):C.GC_309,(0,4):C.GC_304})

V_49 = Vertex(name = 'V_49',
              particles = [ P.g, P.g, P.H ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.VVS4 ],
              couplings = {(0,0):C.GC_347})

V_50 = Vertex(name = 'V_50',
              particles = [ P.W__minus__, P.W__plus__, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS1, L.VVSS3, L.VVSS4 ],
              couplings = {(0,0):C.GC_39,(0,2):C.GC_38,(0,1):C.GC_174})

V_51 = Vertex(name = 'V_51',
              particles = [ P.W__minus__, P.W__plus__, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS3 ],
              couplings = {(0,0):C.GC_487})

V_52 = Vertex(name = 'V_52',
              particles = [ P.W__minus__, P.W__plus__, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS3 ],
              couplings = {(0,0):C.GC_488})

V_53 = Vertex(name = 'V_53',
              particles = [ P.W__minus__, P.W__plus__, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS3 ],
              couplings = {(0,0):C.GC_489})

V_54 = Vertex(name = 'V_54',
              particles = [ P.W__minus__, P.W__plus__, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS3 ],
              couplings = {(0,0):C.GC_491})

V_55 = Vertex(name = 'V_55',
              particles = [ P.W__minus__, P.W__plus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVS1, L.VVS3, L.VVS4 ],
              couplings = {(0,0):C.GC_350,(0,2):C.GC_349,(0,1):C.GC_391})

V_56 = Vertex(name = 'V_56',
              particles = [ P.W__minus__, P.W__plus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVS3 ],
              couplings = {(0,0):C.GC_554})

V_57 = Vertex(name = 'V_57',
              particles = [ P.W__minus__, P.W__plus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVS3 ],
              couplings = {(0,0):C.GC_555})

V_58 = Vertex(name = 'V_58',
              particles = [ P.W__minus__, P.W__plus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVS3 ],
              couplings = {(0,0):C.GC_556})

V_59 = Vertex(name = 'V_59',
              particles = [ P.W__minus__, P.W__plus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVS3 ],
              couplings = {(0,0):C.GC_557})

V_60 = Vertex(name = 'V_60',
              particles = [ P.a, P.Z, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS1, L.VVSS4 ],
              couplings = {(0,0):C.GC_298,(0,1):C.GC_297})

V_61 = Vertex(name = 'V_61',
              particles = [ P.a, P.Z, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS1, L.VVSS4 ],
              couplings = {(0,0):C.GC_273,(0,1):C.GC_272})

V_62 = Vertex(name = 'V_62',
              particles = [ P.a, P.Z, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS1, L.VVSS4 ],
              couplings = {(0,0):C.GC_279,(0,1):C.GC_274})

V_63 = Vertex(name = 'V_63',
              particles = [ P.a, P.Z, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVS1, L.VVS4 ],
              couplings = {(0,0):C.GC_563,(0,1):C.GC_303})

V_64 = Vertex(name = 'V_64',
              particles = [ P.a, P.Z, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVS1, L.VVS4 ],
              couplings = {(0,0):C.GC_444,(0,1):C.GC_562})

V_65 = Vertex(name = 'V_65',
              particles = [ P.a, P.Z, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVS1, L.VVS4 ],
              couplings = {(0,0):C.GC_450,(0,1):C.GC_443})

V_66 = Vertex(name = 'V_66',
              particles = [ P.a, P.Z, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVS4 ],
              couplings = {(0,0):C.GC_445})

V_67 = Vertex(name = 'V_67',
              particles = [ P.Z, P.Z, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS2, L.VVSS3, L.VVSS4 ],
              couplings = {(0,0):C.GC_109,(0,2):C.GC_108,(0,1):C.GC_176})

V_68 = Vertex(name = 'V_68',
              particles = [ P.Z, P.Z, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS2, L.VVSS3, L.VVSS4 ],
              couplings = {(0,0):C.GC_277,(0,2):C.GC_276,(0,1):C.GC_493})

V_69 = Vertex(name = 'V_69',
              particles = [ P.Z, P.Z, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS2, L.VVSS3, L.VVSS4 ],
              couplings = {(0,0):C.GC_294,(0,2):C.GC_293,(0,1):C.GC_494})

V_70 = Vertex(name = 'V_70',
              particles = [ P.Z, P.Z, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS3 ],
              couplings = {(0,0):C.GC_495})

V_71 = Vertex(name = 'V_71',
              particles = [ P.Z, P.Z, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS3 ],
              couplings = {(0,0):C.GC_496})

V_72 = Vertex(name = 'V_72',
              particles = [ P.Z, P.Z, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVS2, L.VVS3, L.VVS4 ],
              couplings = {(0,0):C.GC_362,(0,2):C.GC_361,(0,1):C.GC_392})

V_73 = Vertex(name = 'V_73',
              particles = [ P.Z, P.Z, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVS2, L.VVS3, L.VVS4 ],
              couplings = {(0,0):C.GC_448,(0,2):C.GC_447,(0,1):C.GC_558})

V_74 = Vertex(name = 'V_74',
              particles = [ P.Z, P.Z, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVS2, L.VVS3, L.VVS4 ],
              couplings = {(0,0):C.GC_454,(0,2):C.GC_453,(0,1):C.GC_559})

V_75 = Vertex(name = 'V_75',
              particles = [ P.Z, P.Z, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVS3 ],
              couplings = {(0,0):C.GC_560})

V_76 = Vertex(name = 'V_76',
              particles = [ P.Z, P.Z, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVS3 ],
              couplings = {(0,0):C.GC_561})

V_77 = Vertex(name = 'V_77',
              particles = [ P.g, P.g, P.g, P.H, P.H ],
              color = [ 'f(1,2,3)' ],
              lorentz = [ L.VVVSS3, L.VVVSS6 ],
              couplings = {(0,0):C.GC_149,(0,1):C.GC_148})

V_78 = Vertex(name = 'V_78',
              particles = [ P.g, P.g, P.g, P.H ],
              color = [ 'f(1,2,3)' ],
              lorentz = [ L.VVVS10, L.VVVS3, L.VVVS6, L.VVVS7, L.VVVS8, L.VVVS9 ],
              couplings = {(0,1):C.GC_384,(0,4):C.GC_305,(0,0):C.GC_314,(0,5):C.GC_310,(0,3):C.GC_307,(0,2):C.GC_301})

V_79 = Vertex(name = 'V_79',
              particles = [ P.g, P.g, P.g, P.H ],
              color = [ 'f(1,2,3)' ],
              lorentz = [ L.VVVS6 ],
              couplings = {(0,0):C.GC_383})

V_80 = Vertex(name = 'V_80',
              particles = [ P.a, P.W__minus__, P.W__plus__, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVVSS1, L.VVVSS3, L.VVVSS4, L.VVVSS6 ],
              couplings = {(0,1):C.GC_135,(0,0):C.GC_234,(0,3):C.GC_132,(0,2):C.GC_233})

V_81 = Vertex(name = 'V_81',
              particles = [ P.a, P.W__minus__, P.W__plus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVVS1, L.VVVS3, L.VVVS4, L.VVVS6 ],
              couplings = {(0,1):C.GC_376,(0,0):C.GC_429,(0,3):C.GC_373,(0,2):C.GC_428})

V_82 = Vertex(name = 'V_82',
              particles = [ P.W__minus__, P.W__plus__, P.Z, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVVSS2, L.VVVSS3, L.VVVSS5, L.VVVSS6 ],
              couplings = {(0,1):C.GC_235,(0,0):C.GC_134,(0,3):C.GC_232,(0,2):C.GC_133})

V_83 = Vertex(name = 'V_83',
              particles = [ P.W__minus__, P.W__plus__, P.Z, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVVS2, L.VVVS3, L.VVVS5, L.VVVS6 ],
              couplings = {(0,1):C.GC_430,(0,0):C.GC_375,(0,3):C.GC_427,(0,2):C.GC_374})

V_84 = Vertex(name = 'V_84',
              particles = [ P.H, P.H, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSSS1, L.SSSS2, L.SSSS3 ],
              couplings = {(0,0):C.GC_9,(0,2):C.GC_34,(0,1):C.GC_35})

V_85 = Vertex(name = 'V_85',
              particles = [ P.H, P.H, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSSS1 ],
              couplings = {(0,0):C.GC_459})

V_86 = Vertex(name = 'V_86',
              particles = [ P.H, P.H, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSSS1 ],
              couplings = {(0,0):C.GC_479})

V_87 = Vertex(name = 'V_87',
              particles = [ P.H, P.H, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSSS1 ],
              couplings = {(0,0):C.GC_480})

V_88 = Vertex(name = 'V_88',
              particles = [ P.H, P.H, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSSS1 ],
              couplings = {(0,0):C.GC_481})

V_89 = Vertex(name = 'V_89',
              particles = [ P.H, P.H, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSSS1 ],
              couplings = {(0,0):C.GC_482})

V_90 = Vertex(name = 'V_90',
              particles = [ P.H, P.H, P.H, P.H, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSSSSS1 ],
              couplings = {(0,0):C.GC_33})

V_91 = Vertex(name = 'V_91',
              particles = [ P.H, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSS1, L.SSS2, L.SSS3 ],
              couplings = {(0,0):C.GC_336,(0,2):C.GC_345,(0,1):C.GC_346})

V_92 = Vertex(name = 'V_92',
              particles = [ P.H, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSS1 ],
              couplings = {(0,0):C.GC_549})

V_93 = Vertex(name = 'V_93',
              particles = [ P.H, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSS1 ],
              couplings = {(0,0):C.GC_550})

V_94 = Vertex(name = 'V_94',
              particles = [ P.H, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSS1 ],
              couplings = {(0,0):C.GC_551})

V_95 = Vertex(name = 'V_95',
              particles = [ P.H, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSS1 ],
              couplings = {(0,0):C.GC_552})

V_96 = Vertex(name = 'V_96',
              particles = [ P.H, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSS1 ],
              couplings = {(0,0):C.GC_553})

V_97 = Vertex(name = 'V_97',
              particles = [ P.H, P.H, P.H, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSSSS1 ],
              couplings = {(0,0):C.GC_344})

V_98 = Vertex(name = 'V_98',
              particles = [ P.g, P.g, P.g, P.g, P.H, P.H ],
              color = [ 'f(-1,1,2)*f(3,4,-1)', 'f(-1,1,3)*f(2,4,-1)', 'f(-1,1,4)*f(2,3,-1)' ],
              lorentz = [ L.VVVVSS1, L.VVVVSS3, L.VVVVSS4 ],
              couplings = {(1,1):C.GC_156,(0,0):C.GC_156,(2,2):C.GC_156})

V_99 = Vertex(name = 'V_99',
              particles = [ P.g, P.g, P.g, P.g, P.H ],
              color = [ 'f(-1,1,2)*f(3,4,-1)', 'f(-1,1,3)*f(2,4,-1)', 'f(-1,1,4)*f(2,3,-1)' ],
              lorentz = [ L.VVVVS1, L.VVVVS10, L.VVVVS11, L.VVVVS12, L.VVVVS13, L.VVVVS14, L.VVVVS15, L.VVVVS16, L.VVVVS17, L.VVVVS19, L.VVVVS2, L.VVVVS3, L.VVVVS4, L.VVVVS7, L.VVVVS8 ],
              couplings = {(2,5):C.GC_306,(2,8):C.GC_315,(1,4):C.GC_306,(1,9):C.GC_315,(2,6):C.GC_312,(0,11):C.GC_308,(0,12):C.GC_316,(1,7):C.GC_312,(0,3):C.GC_311,(1,2):C.GC_308,(2,1):C.GC_308,(0,10):C.GC_306,(1,13):C.GC_302,(0,0):C.GC_302,(2,14):C.GC_302})

V_100 = Vertex(name = 'V_100',
               particles = [ P.g, P.g, P.g, P.g, P.H ],
               color = [ 'f(-1,1,2)*f(3,4,-1)', 'f(-1,1,3)*f(2,4,-1)', 'f(-1,1,4)*f(2,3,-1)' ],
               lorentz = [ L.VVVVS1, L.VVVVS7, L.VVVVS8 ],
               couplings = {(1,1):C.GC_387,(0,0):C.GC_387,(2,2):C.GC_387})

V_101 = Vertex(name = 'V_101',
               particles = [ P.a, P.a, P.W__minus__, P.W__plus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVVVSS2 ],
               couplings = {(0,0):C.GC_141})

V_102 = Vertex(name = 'V_102',
               particles = [ P.a, P.a, P.W__minus__, P.W__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVVVS6 ],
               couplings = {(0,0):C.GC_380})

V_103 = Vertex(name = 'V_103',
               particles = [ P.W__minus__, P.W__minus__, P.W__plus__, P.W__plus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVVVSS2 ],
               couplings = {(0,0):C.GC_178})

V_104 = Vertex(name = 'V_104',
               particles = [ P.W__minus__, P.W__minus__, P.W__plus__, P.W__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVVVS6 ],
               couplings = {(0,0):C.GC_393})

V_105 = Vertex(name = 'V_105',
               particles = [ P.a, P.W__minus__, P.W__plus__, P.Z, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVVVSS5 ],
               couplings = {(0,0):C.GC_245})

V_106 = Vertex(name = 'V_106',
               particles = [ P.a, P.W__minus__, P.W__plus__, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVVVS9 ],
               couplings = {(0,0):C.GC_436})

V_107 = Vertex(name = 'V_107',
               particles = [ P.Z, P.Z, P.H, P.H, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_179})

V_108 = Vertex(name = 'V_108',
               particles = [ P.Z, P.Z, P.H, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSS1 ],
               couplings = {(0,0):C.GC_394})

V_109 = Vertex(name = 'V_109',
               particles = [ P.W__minus__, P.W__plus__, P.Z, P.Z, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVVVSS2 ],
               couplings = {(0,0):C.GC_180})

V_110 = Vertex(name = 'V_110',
               particles = [ P.W__minus__, P.W__plus__, P.Z, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVVVS6 ],
               couplings = {(0,0):C.GC_395})

V_111 = Vertex(name = 'V_111',
               particles = [ P.H, P.H, P.H, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.SSSS1 ],
               couplings = {(0,0):C.GC_165})

V_112 = Vertex(name = 'V_112',
               particles = [ P.H, P.H, P.H1, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.SSSS1 ],
               couplings = {(0,0):C.GC_170})

V_113 = Vertex(name = 'V_113',
               particles = [ P.H, P.H1, P.H1, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.SSSS1 ],
               couplings = {(0,0):C.GC_171})

V_114 = Vertex(name = 'V_114',
               particles = [ P.H1, P.H1, P.H1, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.SSSS1 ],
               couplings = {(0,0):C.GC_172})

V_115 = Vertex(name = 'V_115',
               particles = [ P.H, P.H, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.SSS1 ],
               couplings = {(0,0):C.GC_388})

V_116 = Vertex(name = 'V_116',
               particles = [ P.H, P.H1, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.SSS1 ],
               couplings = {(0,0):C.GC_389})

V_117 = Vertex(name = 'V_117',
               particles = [ P.H1, P.H1, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.SSS1 ],
               couplings = {(0,0):C.GC_390})

V_118 = Vertex(name = 'V_118',
               particles = [ P.a, P.W__minus__, P.W1__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVV5 ],
               couplings = {(0,0):C.GC_162})

V_119 = Vertex(name = 'V_119',
               particles = [ P.a, P.W1__minus__, P.W1__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVV5 ],
               couplings = {(0,0):C.GC_167})

V_120 = Vertex(name = 'V_120',
               particles = [ P.a, P.W1__minus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVV5 ],
               couplings = {(0,0):C.GC_162})

V_121 = Vertex(name = 'V_121',
               particles = [ P.W__minus__, P.W1__plus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSS3 ],
               couplings = {(0,0):C.GC_186})

V_122 = Vertex(name = 'V_122',
               particles = [ P.W__minus__, P.W1__plus__, P.H, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.VVSS3 ],
               couplings = {(0,0):C.GC_190})

V_123 = Vertex(name = 'V_123',
               particles = [ P.W__minus__, P.W1__plus__, P.H1, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.VVSS3 ],
               couplings = {(0,0):C.GC_194})

V_124 = Vertex(name = 'V_124',
               particles = [ P.W__minus__, P.W1__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVS3 ],
               couplings = {(0,0):C.GC_396})

V_125 = Vertex(name = 'V_125',
               particles = [ P.W__minus__, P.W1__plus__, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.VVS3 ],
               couplings = {(0,0):C.GC_398})

V_126 = Vertex(name = 'V_126',
               particles = [ P.a, P.a, P.W__minus__, P.W1__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVV8 ],
               couplings = {(0,0):C.GC_163})

V_127 = Vertex(name = 'V_127',
               particles = [ P.W__minus__, P.W1__plus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.VVV5 ],
               couplings = {(0,0):C.GC_254})

V_128 = Vertex(name = 'V_128',
               particles = [ P.W__minus__, P.W1__plus__, P.Z1 ],
               color = [ '1' ],
               lorentz = [ L.VVV5 ],
               couplings = {(0,0):C.GC_258})

V_129 = Vertex(name = 'V_129',
               particles = [ P.W1__minus__, P.W1__plus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSS3 ],
               couplings = {(0,0):C.GC_190})

V_130 = Vertex(name = 'V_130',
               particles = [ P.W1__minus__, P.W1__plus__, P.H, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.VVSS3 ],
               couplings = {(0,0):C.GC_194})

V_131 = Vertex(name = 'V_131',
               particles = [ P.W1__minus__, P.W1__plus__, P.H1, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.VVSS3 ],
               couplings = {(0,0):C.GC_198})

V_132 = Vertex(name = 'V_132',
               particles = [ P.W1__minus__, P.W1__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVS3 ],
               couplings = {(0,0):C.GC_398})

V_133 = Vertex(name = 'V_133',
               particles = [ P.W1__minus__, P.W1__plus__, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.VVS3 ],
               couplings = {(0,0):C.GC_400})

V_134 = Vertex(name = 'V_134',
               particles = [ P.a, P.a, P.W1__minus__, P.W1__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVV8 ],
               couplings = {(0,0):C.GC_168})

V_135 = Vertex(name = 'V_135',
               particles = [ P.W1__minus__, P.W1__plus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.VVV5 ],
               couplings = {(0,0):C.GC_258})

V_136 = Vertex(name = 'V_136',
               particles = [ P.W1__minus__, P.W1__plus__, P.Z1 ],
               color = [ '1' ],
               lorentz = [ L.VVV5 ],
               couplings = {(0,0):C.GC_261})

V_137 = Vertex(name = 'V_137',
               particles = [ P.W__minus__, P.W__minus__, P.W1__plus__, P.W1__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVV8 ],
               couplings = {(0,0):C.GC_191})

V_138 = Vertex(name = 'V_138',
               particles = [ P.W__minus__, P.W1__minus__, P.W1__plus__, P.W1__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVV8 ],
               couplings = {(0,0):C.GC_195})

V_139 = Vertex(name = 'V_139',
               particles = [ P.W1__minus__, P.W1__minus__, P.W1__plus__, P.W1__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVV8 ],
               couplings = {(0,0):C.GC_199})

V_140 = Vertex(name = 'V_140',
               particles = [ P.W__minus__, P.W__plus__, P.H, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.VVSS3 ],
               couplings = {(0,0):C.GC_186})

V_141 = Vertex(name = 'V_141',
               particles = [ P.W__minus__, P.W__plus__, P.H1, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.VVSS3 ],
               couplings = {(0,0):C.GC_190})

V_142 = Vertex(name = 'V_142',
               particles = [ P.W__minus__, P.W__plus__, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.VVS3 ],
               couplings = {(0,0):C.GC_396})

V_143 = Vertex(name = 'V_143',
               particles = [ P.W__minus__, P.W__plus__, P.Z1 ],
               color = [ '1' ],
               lorentz = [ L.VVV5 ],
               couplings = {(0,0):C.GC_254})

V_144 = Vertex(name = 'V_144',
               particles = [ P.W1__minus__, P.W__plus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSS3 ],
               couplings = {(0,0):C.GC_186})

V_145 = Vertex(name = 'V_145',
               particles = [ P.W1__minus__, P.W__plus__, P.H, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.VVSS3 ],
               couplings = {(0,0):C.GC_190})

V_146 = Vertex(name = 'V_146',
               particles = [ P.W1__minus__, P.W__plus__, P.H1, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.VVSS3 ],
               couplings = {(0,0):C.GC_194})

V_147 = Vertex(name = 'V_147',
               particles = [ P.W1__minus__, P.W__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVS3 ],
               couplings = {(0,0):C.GC_396})

V_148 = Vertex(name = 'V_148',
               particles = [ P.W1__minus__, P.W__plus__, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.VVS3 ],
               couplings = {(0,0):C.GC_398})

V_149 = Vertex(name = 'V_149',
               particles = [ P.a, P.a, P.W1__minus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVV8 ],
               couplings = {(0,0):C.GC_163})

V_150 = Vertex(name = 'V_150',
               particles = [ P.W1__minus__, P.W__plus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.VVV5 ],
               couplings = {(0,0):C.GC_254})

V_151 = Vertex(name = 'V_151',
               particles = [ P.W1__minus__, P.W__plus__, P.Z1 ],
               color = [ '1' ],
               lorentz = [ L.VVV5 ],
               couplings = {(0,0):C.GC_258})

V_152 = Vertex(name = 'V_152',
               particles = [ P.W__minus__, P.W__minus__, P.W1__plus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVV8 ],
               couplings = {(0,0):C.GC_187})

V_153 = Vertex(name = 'V_153',
               particles = [ P.W__minus__, P.W1__minus__, P.W1__plus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVV8 ],
               couplings = {(0,0):C.GC_191})

V_154 = Vertex(name = 'V_154',
               particles = [ P.W1__minus__, P.W1__minus__, P.W1__plus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVV8 ],
               couplings = {(0,0):C.GC_195})

V_155 = Vertex(name = 'V_155',
               particles = [ P.W__minus__, P.W1__minus__, P.W__plus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVV8 ],
               couplings = {(0,0):C.GC_187})

V_156 = Vertex(name = 'V_156',
               particles = [ P.W1__minus__, P.W1__minus__, P.W__plus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVV8 ],
               couplings = {(0,0):C.GC_191})

V_157 = Vertex(name = 'V_157',
               particles = [ P.a, P.W__minus__, P.W1__plus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.VVVV11 ],
               couplings = {(0,0):C.GC_255})

V_158 = Vertex(name = 'V_158',
               particles = [ P.a, P.W1__minus__, P.W1__plus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.VVVV11 ],
               couplings = {(0,0):C.GC_259})

V_159 = Vertex(name = 'V_159',
               particles = [ P.a, P.W1__minus__, P.W__plus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.VVVV11 ],
               couplings = {(0,0):C.GC_255})

V_160 = Vertex(name = 'V_160',
               particles = [ P.Z, P.Z, P.H, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.VVSS3 ],
               couplings = {(0,0):C.GC_188})

V_161 = Vertex(name = 'V_161',
               particles = [ P.Z, P.Z, P.H1, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.VVSS3 ],
               couplings = {(0,0):C.GC_192})

V_162 = Vertex(name = 'V_162',
               particles = [ P.Z, P.Z, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.VVS3 ],
               couplings = {(0,0):C.GC_397})

V_163 = Vertex(name = 'V_163',
               particles = [ P.W__minus__, P.W1__plus__, P.Z, P.Z ],
               color = [ '1' ],
               lorentz = [ L.VVVV8 ],
               couplings = {(0,0):C.GC_189})

V_164 = Vertex(name = 'V_164',
               particles = [ P.W1__minus__, P.W1__plus__, P.Z, P.Z ],
               color = [ '1' ],
               lorentz = [ L.VVVV8 ],
               couplings = {(0,0):C.GC_193})

V_165 = Vertex(name = 'V_165',
               particles = [ P.W1__minus__, P.W__plus__, P.Z, P.Z ],
               color = [ '1' ],
               lorentz = [ L.VVVV8 ],
               couplings = {(0,0):C.GC_189})

V_166 = Vertex(name = 'V_166',
               particles = [ P.a, P.W__minus__, P.W1__plus__, P.Z1 ],
               color = [ '1' ],
               lorentz = [ L.VVVV11 ],
               couplings = {(0,0):C.GC_259})

V_167 = Vertex(name = 'V_167',
               particles = [ P.a, P.W1__minus__, P.W1__plus__, P.Z1 ],
               color = [ '1' ],
               lorentz = [ L.VVVV11 ],
               couplings = {(0,0):C.GC_262})

V_168 = Vertex(name = 'V_168',
               particles = [ P.a, P.W__minus__, P.W__plus__, P.Z1 ],
               color = [ '1' ],
               lorentz = [ L.VVVV11 ],
               couplings = {(0,0):C.GC_255})

V_169 = Vertex(name = 'V_169',
               particles = [ P.a, P.W1__minus__, P.W__plus__, P.Z1 ],
               color = [ '1' ],
               lorentz = [ L.VVVV11 ],
               couplings = {(0,0):C.GC_259})

V_170 = Vertex(name = 'V_170',
               particles = [ P.Z, P.Z1, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSS3 ],
               couplings = {(0,0):C.GC_188})

V_171 = Vertex(name = 'V_171',
               particles = [ P.Z, P.Z1, P.H, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.VVSS3 ],
               couplings = {(0,0):C.GC_192})

V_172 = Vertex(name = 'V_172',
               particles = [ P.Z, P.Z1, P.H1, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.VVSS3 ],
               couplings = {(0,0):C.GC_196})

V_173 = Vertex(name = 'V_173',
               particles = [ P.Z, P.Z1, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVS3 ],
               couplings = {(0,0):C.GC_397})

V_174 = Vertex(name = 'V_174',
               particles = [ P.Z, P.Z1, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.VVS3 ],
               couplings = {(0,0):C.GC_399})

V_175 = Vertex(name = 'V_175',
               particles = [ P.W__minus__, P.W1__plus__, P.Z, P.Z1 ],
               color = [ '1' ],
               lorentz = [ L.VVVV8 ],
               couplings = {(0,0):C.GC_193})

V_176 = Vertex(name = 'V_176',
               particles = [ P.W1__minus__, P.W1__plus__, P.Z, P.Z1 ],
               color = [ '1' ],
               lorentz = [ L.VVVV8 ],
               couplings = {(0,0):C.GC_197})

V_177 = Vertex(name = 'V_177',
               particles = [ P.W__minus__, P.W__plus__, P.Z, P.Z1 ],
               color = [ '1' ],
               lorentz = [ L.VVVV8 ],
               couplings = {(0,0):C.GC_189})

V_178 = Vertex(name = 'V_178',
               particles = [ P.W1__minus__, P.W__plus__, P.Z, P.Z1 ],
               color = [ '1' ],
               lorentz = [ L.VVVV8 ],
               couplings = {(0,0):C.GC_193})

V_179 = Vertex(name = 'V_179',
               particles = [ P.Z1, P.Z1, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSS3 ],
               couplings = {(0,0):C.GC_192})

V_180 = Vertex(name = 'V_180',
               particles = [ P.Z1, P.Z1, P.H, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.VVSS3 ],
               couplings = {(0,0):C.GC_196})

V_181 = Vertex(name = 'V_181',
               particles = [ P.Z1, P.Z1, P.H1, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.VVSS3 ],
               couplings = {(0,0):C.GC_200})

V_182 = Vertex(name = 'V_182',
               particles = [ P.Z1, P.Z1, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVS3 ],
               couplings = {(0,0):C.GC_399})

V_183 = Vertex(name = 'V_183',
               particles = [ P.Z1, P.Z1, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.VVS3 ],
               couplings = {(0,0):C.GC_401})

V_184 = Vertex(name = 'V_184',
               particles = [ P.W__minus__, P.W1__plus__, P.Z1, P.Z1 ],
               color = [ '1' ],
               lorentz = [ L.VVVV8 ],
               couplings = {(0,0):C.GC_197})

V_185 = Vertex(name = 'V_185',
               particles = [ P.W1__minus__, P.W1__plus__, P.Z1, P.Z1 ],
               color = [ '1' ],
               lorentz = [ L.VVVV8 ],
               couplings = {(0,0):C.GC_201})

V_186 = Vertex(name = 'V_186',
               particles = [ P.W__minus__, P.W__plus__, P.Z1, P.Z1 ],
               color = [ '1' ],
               lorentz = [ L.VVVV8 ],
               couplings = {(0,0):C.GC_193})

V_187 = Vertex(name = 'V_187',
               particles = [ P.W1__minus__, P.W__plus__, P.Z1, P.Z1 ],
               color = [ '1' ],
               lorentz = [ L.VVVV8 ],
               couplings = {(0,0):C.GC_197})

V_188 = Vertex(name = 'V_188',
               particles = [ P.b__tilde__, P.b, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV9 ],
               couplings = {(0,0):C.GC_1,(0,2):C.GC_354,(0,1):C.GC_353})

V_189 = Vertex(name = 'V_189',
               particles = [ P.b__tilde__, P.b, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV9 ],
               couplings = {(0,0):C.GC_462,(0,2):C.GC_440,(0,1):C.GC_439})

V_190 = Vertex(name = 'V_190',
               particles = [ P.b__tilde__, P.b, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_472})

V_191 = Vertex(name = 'V_191',
               particles = [ P.b__tilde__, P.b, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_483})

V_192 = Vertex(name = 'V_192',
               particles = [ P.b__tilde__, P.b, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_530})

V_193 = Vertex(name = 'V_193',
               particles = [ P.d__tilde__, P.d, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV9 ],
               couplings = {(0,0):C.GC_1,(0,2):C.GC_709,(0,1):C.GC_708})

V_194 = Vertex(name = 'V_194',
               particles = [ P.d__tilde__, P.d, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV9 ],
               couplings = {(0,0):C.GC_462,(0,2):C.GC_725,(0,1):C.GC_724})

V_195 = Vertex(name = 'V_195',
               particles = [ P.d__tilde__, P.d, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_472})

V_196 = Vertex(name = 'V_196',
               particles = [ P.d__tilde__, P.d, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_483})

V_197 = Vertex(name = 'V_197',
               particles = [ P.d__tilde__, P.d, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_530})

V_198 = Vertex(name = 'V_198',
               particles = [ P.s__tilde__, P.s, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV9 ],
               couplings = {(0,0):C.GC_1,(0,2):C.GC_965,(0,1):C.GC_964})

V_199 = Vertex(name = 'V_199',
               particles = [ P.s__tilde__, P.s, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV9 ],
               couplings = {(0,0):C.GC_462,(0,2):C.GC_981,(0,1):C.GC_980})

V_200 = Vertex(name = 'V_200',
               particles = [ P.s__tilde__, P.s, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_472})

V_201 = Vertex(name = 'V_201',
               particles = [ P.s__tilde__, P.s, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_483})

V_202 = Vertex(name = 'V_202',
               particles = [ P.s__tilde__, P.s, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_530})

V_203 = Vertex(name = 'V_203',
               particles = [ P.e__plus__, P.e__minus__, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV9 ],
               couplings = {(0,0):C.GC_4,(0,2):C.GC_792,(0,1):C.GC_791})

V_204 = Vertex(name = 'V_204',
               particles = [ P.e__plus__, P.e__minus__, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV9 ],
               couplings = {(0,0):C.GC_464,(0,2):C.GC_806,(0,1):C.GC_805})

V_205 = Vertex(name = 'V_205',
               particles = [ P.e__plus__, P.e__minus__, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_475})

V_206 = Vertex(name = 'V_206',
               particles = [ P.e__plus__, P.e__minus__, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_485})

V_207 = Vertex(name = 'V_207',
               particles = [ P.e__plus__, P.e__minus__, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_532})

V_208 = Vertex(name = 'V_208',
               particles = [ P.mu__plus__, P.mu__minus__, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV9 ],
               couplings = {(0,0):C.GC_4,(0,2):C.GC_872,(0,1):C.GC_871})

V_209 = Vertex(name = 'V_209',
               particles = [ P.mu__plus__, P.mu__minus__, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV9 ],
               couplings = {(0,0):C.GC_464,(0,2):C.GC_886,(0,1):C.GC_885})

V_210 = Vertex(name = 'V_210',
               particles = [ P.mu__plus__, P.mu__minus__, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_475})

V_211 = Vertex(name = 'V_211',
               particles = [ P.mu__plus__, P.mu__minus__, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_485})

V_212 = Vertex(name = 'V_212',
               particles = [ P.mu__plus__, P.mu__minus__, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_532})

V_213 = Vertex(name = 'V_213',
               particles = [ P.ta__plus__, P.ta__minus__, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV9 ],
               couplings = {(0,0):C.GC_4,(0,2):C.GC_1077,(0,1):C.GC_1076})

V_214 = Vertex(name = 'V_214',
               particles = [ P.ta__plus__, P.ta__minus__, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV9 ],
               couplings = {(0,0):C.GC_464,(0,2):C.GC_1091,(0,1):C.GC_1090})

V_215 = Vertex(name = 'V_215',
               particles = [ P.ta__plus__, P.ta__minus__, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_475})

V_216 = Vertex(name = 'V_216',
               particles = [ P.ta__plus__, P.ta__minus__, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_485})

V_217 = Vertex(name = 'V_217',
               particles = [ P.ta__plus__, P.ta__minus__, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_532})

V_218 = Vertex(name = 'V_218',
               particles = [ P.t__tilde__, P.t, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV9 ],
               couplings = {(0,0):C.GC_2,(0,2):C.GC_358,(0,1):C.GC_357})

V_219 = Vertex(name = 'V_219',
               particles = [ P.t__tilde__, P.t, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV9 ],
               couplings = {(0,0):C.GC_463,(0,2):C.GC_452,(0,1):C.GC_451})

V_220 = Vertex(name = 'V_220',
               particles = [ P.t__tilde__, P.t, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_473})

V_221 = Vertex(name = 'V_221',
               particles = [ P.t__tilde__, P.t, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_484})

V_222 = Vertex(name = 'V_222',
               particles = [ P.t__tilde__, P.t, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_531})

V_223 = Vertex(name = 'V_223',
               particles = [ P.u__tilde__, P.u, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV9 ],
               couplings = {(0,0):C.GC_2,(0,2):C.GC_1166,(0,1):C.GC_1165})

V_224 = Vertex(name = 'V_224',
               particles = [ P.u__tilde__, P.u, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV9 ],
               couplings = {(0,0):C.GC_463,(0,2):C.GC_1189,(0,1):C.GC_1188})

V_225 = Vertex(name = 'V_225',
               particles = [ P.u__tilde__, P.u, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_473})

V_226 = Vertex(name = 'V_226',
               particles = [ P.u__tilde__, P.u, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_484})

V_227 = Vertex(name = 'V_227',
               particles = [ P.u__tilde__, P.u, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_531})

V_228 = Vertex(name = 'V_228',
               particles = [ P.c__tilde__, P.c, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV9 ],
               couplings = {(0,0):C.GC_2,(0,2):C.GC_622,(0,1):C.GC_621})

V_229 = Vertex(name = 'V_229',
               particles = [ P.c__tilde__, P.c, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV9 ],
               couplings = {(0,0):C.GC_463,(0,2):C.GC_645,(0,1):C.GC_644})

V_230 = Vertex(name = 'V_230',
               particles = [ P.c__tilde__, P.c, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_473})

V_231 = Vertex(name = 'V_231',
               particles = [ P.c__tilde__, P.c, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_484})

V_232 = Vertex(name = 'V_232',
               particles = [ P.c__tilde__, P.c, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_531})

V_233 = Vertex(name = 'V_233',
               particles = [ P.b__tilde__, P.b, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV3, L.FFV5, L.FFV8, L.FFV9 ],
               couplings = {(0,0):C.GC_263,(0,2):C.GC_204,(0,4):C.GC_510,(0,3):C.GC_508,(0,5):C.GC_356,(0,1):C.GC_355})

V_234 = Vertex(name = 'V_234',
               particles = [ P.b__tilde__, P.b, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV3, L.FFV9 ],
               couplings = {(0,0):C.GC_466,(0,2):C.GC_519,(0,3):C.GC_438,(0,1):C.GC_437})

V_235 = Vertex(name = 'V_235',
               particles = [ P.b__tilde__, P.b, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_540,(0,1):C.GC_521})

V_236 = Vertex(name = 'V_236',
               particles = [ P.b__tilde__, P.b, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_543,(0,1):C.GC_523})

V_237 = Vertex(name = 'V_237',
               particles = [ P.b__tilde__, P.b, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_546,(0,1):C.GC_527})

V_238 = Vertex(name = 'V_238',
               particles = [ P.d__tilde__, P.d, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV3, L.FFV5, L.FFV8, L.FFV9 ],
               couplings = {(0,0):C.GC_263,(0,2):C.GC_204,(0,4):C.GC_510,(0,3):C.GC_509,(0,5):C.GC_711,(0,1):C.GC_710})

V_239 = Vertex(name = 'V_239',
               particles = [ P.d__tilde__, P.d, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV3, L.FFV9 ],
               couplings = {(0,0):C.GC_466,(0,2):C.GC_515,(0,3):C.GC_723,(0,1):C.GC_722})

V_240 = Vertex(name = 'V_240',
               particles = [ P.d__tilde__, P.d, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_540,(0,1):C.GC_517})

V_241 = Vertex(name = 'V_241',
               particles = [ P.d__tilde__, P.d, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_543,(0,1):C.GC_519})

V_242 = Vertex(name = 'V_242',
               particles = [ P.d__tilde__, P.d, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_546,(0,1):C.GC_527})

V_243 = Vertex(name = 'V_243',
               particles = [ P.s__tilde__, P.s, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV3, L.FFV5, L.FFV8, L.FFV9 ],
               couplings = {(0,0):C.GC_263,(0,2):C.GC_204,(0,4):C.GC_510,(0,3):C.GC_509,(0,5):C.GC_967,(0,1):C.GC_966})

V_244 = Vertex(name = 'V_244',
               particles = [ P.s__tilde__, P.s, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV3, L.FFV9 ],
               couplings = {(0,0):C.GC_466,(0,2):C.GC_515,(0,3):C.GC_979,(0,1):C.GC_978})

V_245 = Vertex(name = 'V_245',
               particles = [ P.s__tilde__, P.s, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_540,(0,1):C.GC_517})

V_246 = Vertex(name = 'V_246',
               particles = [ P.s__tilde__, P.s, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_543,(0,1):C.GC_519})

V_247 = Vertex(name = 'V_247',
               particles = [ P.s__tilde__, P.s, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_546,(0,1):C.GC_527})

V_248 = Vertex(name = 'V_248',
               particles = [ P.e__plus__, P.e__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV3, L.FFV5, L.FFV6, L.FFV9 ],
               couplings = {(0,0):C.GC_265,(0,2):C.GC_204,(0,4):C.GC_513,(0,3):C.GC_514,(0,5):C.GC_794,(0,1):C.GC_793})

V_249 = Vertex(name = 'V_249',
               particles = [ P.e__plus__, P.e__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV3, L.FFV9 ],
               couplings = {(0,0):C.GC_469,(0,2):C.GC_518,(0,3):C.GC_804,(0,1):C.GC_803})

V_250 = Vertex(name = 'V_250',
               particles = [ P.e__plus__, P.e__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_542,(0,1):C.GC_527})

V_251 = Vertex(name = 'V_251',
               particles = [ P.e__plus__, P.e__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_545})

V_252 = Vertex(name = 'V_252',
               particles = [ P.e__plus__, P.e__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_548})

V_253 = Vertex(name = 'V_253',
               particles = [ P.mu__plus__, P.mu__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV3, L.FFV5, L.FFV6, L.FFV9 ],
               couplings = {(0,0):C.GC_265,(0,2):C.GC_204,(0,4):C.GC_513,(0,3):C.GC_514,(0,5):C.GC_874,(0,1):C.GC_873})

V_254 = Vertex(name = 'V_254',
               particles = [ P.mu__plus__, P.mu__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV3, L.FFV9 ],
               couplings = {(0,0):C.GC_469,(0,2):C.GC_518,(0,3):C.GC_884,(0,1):C.GC_883})

V_255 = Vertex(name = 'V_255',
               particles = [ P.mu__plus__, P.mu__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_542,(0,1):C.GC_527})

V_256 = Vertex(name = 'V_256',
               particles = [ P.mu__plus__, P.mu__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_545})

V_257 = Vertex(name = 'V_257',
               particles = [ P.mu__plus__, P.mu__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_548})

V_258 = Vertex(name = 'V_258',
               particles = [ P.ta__plus__, P.ta__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV3, L.FFV5, L.FFV6, L.FFV9 ],
               couplings = {(0,0):C.GC_265,(0,2):C.GC_204,(0,4):C.GC_513,(0,3):C.GC_514,(0,5):C.GC_1079,(0,1):C.GC_1078})

V_259 = Vertex(name = 'V_259',
               particles = [ P.ta__plus__, P.ta__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV3, L.FFV9 ],
               couplings = {(0,0):C.GC_469,(0,2):C.GC_518,(0,3):C.GC_1089,(0,1):C.GC_1088})

V_260 = Vertex(name = 'V_260',
               particles = [ P.ta__plus__, P.ta__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_542,(0,1):C.GC_527})

V_261 = Vertex(name = 'V_261',
               particles = [ P.ta__plus__, P.ta__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_545})

V_262 = Vertex(name = 'V_262',
               particles = [ P.ta__plus__, P.ta__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_548})

V_263 = Vertex(name = 'V_263',
               particles = [ P.t__tilde__, P.t, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV3, L.FFV5, L.FFV7, L.FFV9 ],
               couplings = {(0,0):C.GC_264,(0,2):C.GC_203,(0,4):C.GC_512,(0,3):C.GC_524,(0,5):C.GC_369,(0,1):C.GC_367})

V_264 = Vertex(name = 'V_264',
               particles = [ P.t__tilde__, P.t, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV3, L.FFV9 ],
               couplings = {(0,0):C.GC_467,(0,2):C.GC_520,(0,3):C.GC_442,(0,1):C.GC_441})

V_265 = Vertex(name = 'V_265',
               particles = [ P.t__tilde__, P.t, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_541,(0,1):C.GC_521})

V_266 = Vertex(name = 'V_266',
               particles = [ P.t__tilde__, P.t, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_544,(0,1):C.GC_522})

V_267 = Vertex(name = 'V_267',
               particles = [ P.t__tilde__, P.t, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_547,(0,1):C.GC_526})

V_268 = Vertex(name = 'V_268',
               particles = [ P.u__tilde__, P.u, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV3, L.FFV5, L.FFV7, L.FFV9 ],
               couplings = {(0,0):C.GC_264,(0,2):C.GC_203,(0,4):C.GC_512,(0,3):C.GC_525,(0,5):C.GC_1175,(0,1):C.GC_1173})

V_269 = Vertex(name = 'V_269',
               particles = [ P.u__tilde__, P.u, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV3, L.FFV9 ],
               couplings = {(0,0):C.GC_467,(0,2):C.GC_515,(0,3):C.GC_1187,(0,1):C.GC_1186})

V_270 = Vertex(name = 'V_270',
               particles = [ P.u__tilde__, P.u, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_541,(0,1):C.GC_516})

V_271 = Vertex(name = 'V_271',
               particles = [ P.u__tilde__, P.u, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_544,(0,1):C.GC_520})

V_272 = Vertex(name = 'V_272',
               particles = [ P.u__tilde__, P.u, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_547,(0,1):C.GC_526})

V_273 = Vertex(name = 'V_273',
               particles = [ P.c__tilde__, P.c, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV3, L.FFV5, L.FFV7, L.FFV9 ],
               couplings = {(0,0):C.GC_264,(0,2):C.GC_203,(0,4):C.GC_512,(0,3):C.GC_525,(0,5):C.GC_631,(0,1):C.GC_629})

V_274 = Vertex(name = 'V_274',
               particles = [ P.c__tilde__, P.c, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV3, L.FFV9 ],
               couplings = {(0,0):C.GC_467,(0,2):C.GC_515,(0,3):C.GC_643,(0,1):C.GC_642})

V_275 = Vertex(name = 'V_275',
               particles = [ P.c__tilde__, P.c, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_541,(0,1):C.GC_516})

V_276 = Vertex(name = 'V_276',
               particles = [ P.c__tilde__, P.c, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_544,(0,1):C.GC_520})

V_277 = Vertex(name = 'V_277',
               particles = [ P.c__tilde__, P.c, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_547,(0,1):C.GC_526})

V_278 = Vertex(name = 'V_278',
               particles = [ P.b__tilde__, P.b, P.g ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV9 ],
               couplings = {(0,0):C.GC_6,(0,2):C.GC_338,(0,1):C.GC_337})

V_279 = Vertex(name = 'V_279',
               particles = [ P.d__tilde__, P.d, P.g ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV9 ],
               couplings = {(0,0):C.GC_6,(0,2):C.GC_702,(0,1):C.GC_701})

V_280 = Vertex(name = 'V_280',
               particles = [ P.s__tilde__, P.s, P.g ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV9 ],
               couplings = {(0,0):C.GC_6,(0,2):C.GC_958,(0,1):C.GC_957})

V_281 = Vertex(name = 'V_281',
               particles = [ P.t__tilde__, P.t, P.g ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV9 ],
               couplings = {(0,0):C.GC_6,(0,2):C.GC_352,(0,1):C.GC_351})

V_282 = Vertex(name = 'V_282',
               particles = [ P.u__tilde__, P.u, P.g ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV9 ],
               couplings = {(0,0):C.GC_6,(0,2):C.GC_1168,(0,1):C.GC_1167})

V_283 = Vertex(name = 'V_283',
               particles = [ P.c__tilde__, P.c, P.g ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFV1, L.FFV2, L.FFV9 ],
               couplings = {(0,0):C.GC_6,(0,2):C.GC_624,(0,1):C.GC_623})

V_284 = Vertex(name = 'V_284',
               particles = [ P.b__tilde__, P.t, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10, L.FFV3, L.FFV4, L.FFV5 ],
               couplings = {(0,2):C.GC_366,(0,0):C.GC_341,(0,1):C.GC_202,(0,3):C.GC_504})

V_285 = Vertex(name = 'V_285',
               particles = [ P.b__tilde__, P.t, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10, L.FFV3, L.FFV4, L.FFV5 ],
               couplings = {(0,2):C.GC_368,(0,0):C.GC_343,(0,1):C.GC_502,(0,3):C.GC_506})

V_286 = Vertex(name = 'V_286',
               particles = [ P.b__tilde__, P.t, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_503})

V_287 = Vertex(name = 'V_287',
               particles = [ P.b__tilde__, P.t, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_507})

V_288 = Vertex(name = 'V_288',
               particles = [ P.b__tilde__, P.t, P.W__minus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,1):C.GC_212,(0,0):C.GC_211})

V_289 = Vertex(name = 'V_289',
               particles = [ P.b__tilde__, P.t, P.W__minus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS2 ],
               couplings = {(0,0):C.GC_214})

V_290 = Vertex(name = 'V_290',
               particles = [ P.b__tilde__, P.t, P.W__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS3, L.FFVS4, L.FFVS5, L.FFVS7 ],
               couplings = {(0,1):C.GC_119,(0,3):C.GC_23,(0,2):C.GC_407,(0,0):C.GC_406})

V_291 = Vertex(name = 'V_291',
               particles = [ P.b__tilde__, P.t, P.W__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS4, L.FFVS5, L.FFVS7 ],
               couplings = {(0,0):C.GC_121,(0,2):C.GC_25,(0,1):C.GC_409})

V_292 = Vertex(name = 'V_292',
               particles = [ P.d__tilde__, P.u, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10, L.FFV3, L.FFV4, L.FFV5 ],
               couplings = {(0,2):C.GC_1172,(0,0):C.GC_705,(0,1):C.GC_202,(0,3):C.GC_1224})

V_293 = Vertex(name = 'V_293',
               particles = [ P.d__tilde__, P.u, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10, L.FFV3, L.FFV4, L.FFV5 ],
               couplings = {(0,2):C.GC_1174,(0,0):C.GC_707,(0,1):C.GC_501,(0,3):C.GC_1226})

V_294 = Vertex(name = 'V_294',
               particles = [ P.d__tilde__, P.u, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_502})

V_295 = Vertex(name = 'V_295',
               particles = [ P.d__tilde__, P.u, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_507})

V_296 = Vertex(name = 'V_296',
               particles = [ P.s__tilde__, P.c, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10, L.FFV3, L.FFV4, L.FFV5 ],
               couplings = {(0,2):C.GC_628,(0,0):C.GC_961,(0,1):C.GC_202,(0,3):C.GC_1016})

V_297 = Vertex(name = 'V_297',
               particles = [ P.s__tilde__, P.c, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10, L.FFV3, L.FFV4, L.FFV5 ],
               couplings = {(0,2):C.GC_630,(0,0):C.GC_963,(0,1):C.GC_501,(0,3):C.GC_1018})

V_298 = Vertex(name = 'V_298',
               particles = [ P.s__tilde__, P.c, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_502})

V_299 = Vertex(name = 'V_299',
               particles = [ P.s__tilde__, P.c, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_507})

V_300 = Vertex(name = 'V_300',
               particles = [ P.d__tilde__, P.u, P.W__minus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,1):C.GC_1218,(0,0):C.GC_209})

V_301 = Vertex(name = 'V_301',
               particles = [ P.d__tilde__, P.u, P.W__minus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS2 ],
               couplings = {(0,0):C.GC_1220})

V_302 = Vertex(name = 'V_302',
               particles = [ P.s__tilde__, P.c, P.W__minus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,1):C.GC_1010,(0,0):C.GC_209})

V_303 = Vertex(name = 'V_303',
               particles = [ P.s__tilde__, P.c, P.W__minus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS2 ],
               couplings = {(0,0):C.GC_1012})

V_304 = Vertex(name = 'V_304',
               particles = [ P.d__tilde__, P.u, P.W__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS3, L.FFVS4, L.FFVS5, L.FFVS7 ],
               couplings = {(0,1):C.GC_1146,(0,3):C.GC_657,(0,2):C.GC_1221,(0,0):C.GC_404})

V_305 = Vertex(name = 'V_305',
               particles = [ P.d__tilde__, P.u, P.W__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS4, L.FFVS5, L.FFVS7 ],
               couplings = {(0,0):C.GC_1148,(0,2):C.GC_659,(0,1):C.GC_1223})

V_306 = Vertex(name = 'V_306',
               particles = [ P.s__tilde__, P.c, P.W__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS3, L.FFVS4, L.FFVS5, L.FFVS7 ],
               couplings = {(0,1):C.GC_602,(0,3):C.GC_913,(0,2):C.GC_1013,(0,0):C.GC_404})

V_307 = Vertex(name = 'V_307',
               particles = [ P.s__tilde__, P.c, P.W__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS4, L.FFVS5, L.FFVS7 ],
               couplings = {(0,0):C.GC_604,(0,2):C.GC_915,(0,1):C.GC_1015})

V_308 = Vertex(name = 'V_308',
               particles = [ P.e__plus__, P.ve, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV10, L.FFV3 ],
               couplings = {(0,0):C.GC_788,(0,1):C.GC_202})

V_309 = Vertex(name = 'V_309',
               particles = [ P.e__plus__, P.ve, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV10, L.FFV3 ],
               couplings = {(0,0):C.GC_790,(0,1):C.GC_507})

V_310 = Vertex(name = 'V_310',
               particles = [ P.mu__plus__, P.vm, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV10, L.FFV3 ],
               couplings = {(0,0):C.GC_868,(0,1):C.GC_202})

V_311 = Vertex(name = 'V_311',
               particles = [ P.mu__plus__, P.vm, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV10, L.FFV3 ],
               couplings = {(0,0):C.GC_870,(0,1):C.GC_507})

V_312 = Vertex(name = 'V_312',
               particles = [ P.ta__plus__, P.vt, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV10, L.FFV3 ],
               couplings = {(0,0):C.GC_1073,(0,1):C.GC_202})

V_313 = Vertex(name = 'V_313',
               particles = [ P.ta__plus__, P.vt, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV10, L.FFV3 ],
               couplings = {(0,0):C.GC_1075,(0,1):C.GC_507})

V_314 = Vertex(name = 'V_314',
               particles = [ P.e__plus__, P.ve, P.W__minus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_210})

V_315 = Vertex(name = 'V_315',
               particles = [ P.mu__plus__, P.vm, P.W__minus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_210})

V_316 = Vertex(name = 'V_316',
               particles = [ P.ta__plus__, P.vt, P.W__minus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_210})

V_317 = Vertex(name = 'V_317',
               particles = [ P.e__plus__, P.ve, P.W__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS3, L.FFVS7 ],
               couplings = {(0,1):C.GC_751,(0,0):C.GC_405})

V_318 = Vertex(name = 'V_318',
               particles = [ P.e__plus__, P.ve, P.W__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS7 ],
               couplings = {(0,0):C.GC_753})

V_319 = Vertex(name = 'V_319',
               particles = [ P.mu__plus__, P.vm, P.W__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS3, L.FFVS7 ],
               couplings = {(0,1):C.GC_831,(0,0):C.GC_405})

V_320 = Vertex(name = 'V_320',
               particles = [ P.mu__plus__, P.vm, P.W__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS7 ],
               couplings = {(0,0):C.GC_833})

V_321 = Vertex(name = 'V_321',
               particles = [ P.ta__plus__, P.vt, P.W__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS3, L.FFVS7 ],
               couplings = {(0,1):C.GC_1036,(0,0):C.GC_405})

V_322 = Vertex(name = 'V_322',
               particles = [ P.ta__plus__, P.vt, P.W__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS7 ],
               couplings = {(0,0):C.GC_1038})

V_323 = Vertex(name = 'V_323',
               particles = [ P.t__tilde__, P.b, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10, L.FFV3, L.FFV4, L.FFV5 ],
               couplings = {(0,2):C.GC_342,(0,0):C.GC_365,(0,1):C.GC_202,(0,3):C.GC_505})

V_324 = Vertex(name = 'V_324',
               particles = [ P.t__tilde__, P.b, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10, L.FFV3, L.FFV4, L.FFV5 ],
               couplings = {(0,2):C.GC_343,(0,0):C.GC_368,(0,1):C.GC_502,(0,3):C.GC_506})

V_325 = Vertex(name = 'V_325',
               particles = [ P.t__tilde__, P.b, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_503})

V_326 = Vertex(name = 'V_326',
               particles = [ P.t__tilde__, P.b, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_507})

V_327 = Vertex(name = 'V_327',
               particles = [ P.t__tilde__, P.b, P.W__plus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,1):C.GC_213,(0,0):C.GC_211})

V_328 = Vertex(name = 'V_328',
               particles = [ P.t__tilde__, P.b, P.W__plus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS2 ],
               couplings = {(0,0):C.GC_214})

V_329 = Vertex(name = 'V_329',
               particles = [ P.t__tilde__, P.b, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS3, L.FFVS4, L.FFVS5, L.FFVS7 ],
               couplings = {(0,1):C.GC_24,(0,3):C.GC_118,(0,2):C.GC_408,(0,0):C.GC_406})

V_330 = Vertex(name = 'V_330',
               particles = [ P.t__tilde__, P.b, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS4, L.FFVS5, L.FFVS7 ],
               couplings = {(0,0):C.GC_25,(0,2):C.GC_121,(0,1):C.GC_409})

V_331 = Vertex(name = 'V_331',
               particles = [ P.u__tilde__, P.d, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10, L.FFV3, L.FFV4, L.FFV5 ],
               couplings = {(0,2):C.GC_706,(0,0):C.GC_1171,(0,1):C.GC_202,(0,3):C.GC_1225})

V_332 = Vertex(name = 'V_332',
               particles = [ P.u__tilde__, P.d, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10, L.FFV3, L.FFV4, L.FFV5 ],
               couplings = {(0,2):C.GC_707,(0,0):C.GC_1174,(0,1):C.GC_501,(0,3):C.GC_1226})

V_333 = Vertex(name = 'V_333',
               particles = [ P.u__tilde__, P.d, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_502})

V_334 = Vertex(name = 'V_334',
               particles = [ P.u__tilde__, P.d, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_507})

V_335 = Vertex(name = 'V_335',
               particles = [ P.c__tilde__, P.s, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10, L.FFV3, L.FFV4, L.FFV5 ],
               couplings = {(0,2):C.GC_962,(0,0):C.GC_627,(0,1):C.GC_202,(0,3):C.GC_1017})

V_336 = Vertex(name = 'V_336',
               particles = [ P.c__tilde__, P.s, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10, L.FFV3, L.FFV4, L.FFV5 ],
               couplings = {(0,2):C.GC_963,(0,0):C.GC_630,(0,1):C.GC_501,(0,3):C.GC_1018})

V_337 = Vertex(name = 'V_337',
               particles = [ P.c__tilde__, P.s, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_502})

V_338 = Vertex(name = 'V_338',
               particles = [ P.c__tilde__, P.s, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_507})

V_339 = Vertex(name = 'V_339',
               particles = [ P.u__tilde__, P.d, P.W__plus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,1):C.GC_1219,(0,0):C.GC_209})

V_340 = Vertex(name = 'V_340',
               particles = [ P.u__tilde__, P.d, P.W__plus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS2 ],
               couplings = {(0,0):C.GC_1220})

V_341 = Vertex(name = 'V_341',
               particles = [ P.c__tilde__, P.s, P.W__plus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,1):C.GC_1011,(0,0):C.GC_209})

V_342 = Vertex(name = 'V_342',
               particles = [ P.c__tilde__, P.s, P.W__plus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS2 ],
               couplings = {(0,0):C.GC_1012})

V_343 = Vertex(name = 'V_343',
               particles = [ P.u__tilde__, P.d, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS3, L.FFVS4, L.FFVS5, L.FFVS7 ],
               couplings = {(0,1):C.GC_658,(0,3):C.GC_1145,(0,2):C.GC_1222,(0,0):C.GC_404})

V_344 = Vertex(name = 'V_344',
               particles = [ P.u__tilde__, P.d, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS4, L.FFVS5, L.FFVS7 ],
               couplings = {(0,0):C.GC_659,(0,2):C.GC_1148,(0,1):C.GC_1223})

V_345 = Vertex(name = 'V_345',
               particles = [ P.c__tilde__, P.s, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS3, L.FFVS4, L.FFVS5, L.FFVS7 ],
               couplings = {(0,1):C.GC_914,(0,3):C.GC_601,(0,2):C.GC_1014,(0,0):C.GC_404})

V_346 = Vertex(name = 'V_346',
               particles = [ P.c__tilde__, P.s, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS4, L.FFVS5, L.FFVS7 ],
               couplings = {(0,0):C.GC_915,(0,2):C.GC_604,(0,1):C.GC_1015})

V_347 = Vertex(name = 'V_347',
               particles = [ P.ve__tilde__, P.e__minus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV3, L.FFV4 ],
               couplings = {(0,1):C.GC_789,(0,0):C.GC_202})

V_348 = Vertex(name = 'V_348',
               particles = [ P.ve__tilde__, P.e__minus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV3, L.FFV4 ],
               couplings = {(0,1):C.GC_790,(0,0):C.GC_507})

V_349 = Vertex(name = 'V_349',
               particles = [ P.vm__tilde__, P.mu__minus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV3, L.FFV4 ],
               couplings = {(0,1):C.GC_869,(0,0):C.GC_202})

V_350 = Vertex(name = 'V_350',
               particles = [ P.vm__tilde__, P.mu__minus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV3, L.FFV4 ],
               couplings = {(0,1):C.GC_870,(0,0):C.GC_507})

V_351 = Vertex(name = 'V_351',
               particles = [ P.vt__tilde__, P.ta__minus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV3, L.FFV4 ],
               couplings = {(0,1):C.GC_1074,(0,0):C.GC_202})

V_352 = Vertex(name = 'V_352',
               particles = [ P.vt__tilde__, P.ta__minus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV3, L.FFV4 ],
               couplings = {(0,1):C.GC_1075,(0,0):C.GC_507})

V_353 = Vertex(name = 'V_353',
               particles = [ P.ve__tilde__, P.e__minus__, P.W__plus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_210})

V_354 = Vertex(name = 'V_354',
               particles = [ P.vm__tilde__, P.mu__minus__, P.W__plus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_210})

V_355 = Vertex(name = 'V_355',
               particles = [ P.vt__tilde__, P.ta__minus__, P.W__plus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_210})

V_356 = Vertex(name = 'V_356',
               particles = [ P.ve__tilde__, P.e__minus__, P.W__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS3, L.FFVS4 ],
               couplings = {(0,1):C.GC_752,(0,0):C.GC_405})

V_357 = Vertex(name = 'V_357',
               particles = [ P.ve__tilde__, P.e__minus__, P.W__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS4 ],
               couplings = {(0,0):C.GC_753})

V_358 = Vertex(name = 'V_358',
               particles = [ P.vm__tilde__, P.mu__minus__, P.W__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS3, L.FFVS4 ],
               couplings = {(0,1):C.GC_832,(0,0):C.GC_405})

V_359 = Vertex(name = 'V_359',
               particles = [ P.vm__tilde__, P.mu__minus__, P.W__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS4 ],
               couplings = {(0,0):C.GC_833})

V_360 = Vertex(name = 'V_360',
               particles = [ P.vt__tilde__, P.ta__minus__, P.W__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS3, L.FFVS4 ],
               couplings = {(0,1):C.GC_1037,(0,0):C.GC_405})

V_361 = Vertex(name = 'V_361',
               particles = [ P.vt__tilde__, P.ta__minus__, P.W__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS4 ],
               couplings = {(0,0):C.GC_1038})

V_362 = Vertex(name = 'V_362',
               particles = [ P.b__tilde__, P.b, P.Z, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_224,(0,1):C.GC_215})

V_363 = Vertex(name = 'V_363',
               particles = [ P.b__tilde__, P.b, P.Z, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_226})

V_364 = Vertex(name = 'V_364',
               particles = [ P.b__tilde__, P.b, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2, L.FFVS3, L.FFVS5, L.FFVS6 ],
               couplings = {(0,1):C.GC_419,(0,2):C.GC_410,(0,3):C.GC_103,(0,0):C.GC_102})

V_365 = Vertex(name = 'V_365',
               particles = [ P.b__tilde__, P.b, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2, L.FFVS3, L.FFVS6 ],
               couplings = {(0,1):C.GC_421,(0,2):C.GC_267,(0,0):C.GC_266})

V_366 = Vertex(name = 'V_366',
               particles = [ P.d__tilde__, P.d, P.Z, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_218,(0,1):C.GC_216})

V_367 = Vertex(name = 'V_367',
               particles = [ P.d__tilde__, P.d, P.Z, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_220})

V_368 = Vertex(name = 'V_368',
               particles = [ P.s__tilde__, P.s, P.Z, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_218,(0,1):C.GC_216})

V_369 = Vertex(name = 'V_369',
               particles = [ P.s__tilde__, P.s, P.Z, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_220})

V_370 = Vertex(name = 'V_370',
               particles = [ P.d__tilde__, P.d, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2, L.FFVS3, L.FFVS5, L.FFVS6 ],
               couplings = {(0,1):C.GC_413,(0,2):C.GC_411,(0,3):C.GC_685,(0,0):C.GC_684})

V_371 = Vertex(name = 'V_371',
               particles = [ P.d__tilde__, P.d, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2, L.FFVS3, L.FFVS6 ],
               couplings = {(0,1):C.GC_415,(0,2):C.GC_698,(0,0):C.GC_697})

V_372 = Vertex(name = 'V_372',
               particles = [ P.s__tilde__, P.s, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2, L.FFVS3, L.FFVS5, L.FFVS6 ],
               couplings = {(0,1):C.GC_413,(0,2):C.GC_411,(0,3):C.GC_941,(0,0):C.GC_940})

V_373 = Vertex(name = 'V_373',
               particles = [ P.s__tilde__, P.s, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2, L.FFVS3, L.FFVS6 ],
               couplings = {(0,1):C.GC_415,(0,2):C.GC_954,(0,0):C.GC_953})

V_374 = Vertex(name = 'V_374',
               particles = [ P.e__plus__, P.e__minus__, P.Z, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_221,(0,1):C.GC_217})

V_375 = Vertex(name = 'V_375',
               particles = [ P.e__plus__, P.e__minus__, P.Z, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_223})

V_376 = Vertex(name = 'V_376',
               particles = [ P.mu__plus__, P.mu__minus__, P.Z, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_221,(0,1):C.GC_217})

V_377 = Vertex(name = 'V_377',
               particles = [ P.mu__plus__, P.mu__minus__, P.Z, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_223})

V_378 = Vertex(name = 'V_378',
               particles = [ P.ta__plus__, P.ta__minus__, P.Z, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_221,(0,1):C.GC_217})

V_379 = Vertex(name = 'V_379',
               particles = [ P.ta__plus__, P.ta__minus__, P.Z, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_223})

V_380 = Vertex(name = 'V_380',
               particles = [ P.e__plus__, P.e__minus__, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS2, L.FFVS3, L.FFVS5, L.FFVS6 ],
               couplings = {(0,1):C.GC_416,(0,2):C.GC_412,(0,3):C.GC_772,(0,0):C.GC_771})

V_381 = Vertex(name = 'V_381',
               particles = [ P.e__plus__, P.e__minus__, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS2, L.FFVS3, L.FFVS6 ],
               couplings = {(0,1):C.GC_418,(0,2):C.GC_783,(0,0):C.GC_782})

V_382 = Vertex(name = 'V_382',
               particles = [ P.mu__plus__, P.mu__minus__, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS2, L.FFVS3, L.FFVS5, L.FFVS6 ],
               couplings = {(0,1):C.GC_416,(0,2):C.GC_412,(0,3):C.GC_852,(0,0):C.GC_851})

V_383 = Vertex(name = 'V_383',
               particles = [ P.mu__plus__, P.mu__minus__, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS2, L.FFVS3, L.FFVS6 ],
               couplings = {(0,1):C.GC_418,(0,2):C.GC_863,(0,0):C.GC_862})

V_384 = Vertex(name = 'V_384',
               particles = [ P.ta__plus__, P.ta__minus__, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS2, L.FFVS3, L.FFVS5, L.FFVS6 ],
               couplings = {(0,1):C.GC_416,(0,2):C.GC_412,(0,3):C.GC_1057,(0,0):C.GC_1056})

V_385 = Vertex(name = 'V_385',
               particles = [ P.ta__plus__, P.ta__minus__, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS2, L.FFVS3, L.FFVS6 ],
               couplings = {(0,1):C.GC_418,(0,2):C.GC_1068,(0,0):C.GC_1067})

V_386 = Vertex(name = 'V_386',
               particles = [ P.t__tilde__, P.t, P.Z, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_224,(0,1):C.GC_227})

V_387 = Vertex(name = 'V_387',
               particles = [ P.t__tilde__, P.t, P.Z, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_225})

V_388 = Vertex(name = 'V_388',
               particles = [ P.t__tilde__, P.t, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2, L.FFVS3, L.FFVS5, L.FFVS6 ],
               couplings = {(0,1):C.GC_419,(0,2):C.GC_422,(0,3):C.GC_122,(0,0):C.GC_120})

V_389 = Vertex(name = 'V_389',
               particles = [ P.t__tilde__, P.t, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2, L.FFVS3, L.FFVS6 ],
               couplings = {(0,1):C.GC_420,(0,2):C.GC_271,(0,0):C.GC_270})

V_390 = Vertex(name = 'V_390',
               particles = [ P.u__tilde__, P.u, P.Z, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_218,(0,1):C.GC_228})

V_391 = Vertex(name = 'V_391',
               particles = [ P.u__tilde__, P.u, P.Z, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_219})

V_392 = Vertex(name = 'V_392',
               particles = [ P.c__tilde__, P.c, P.Z, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_218,(0,1):C.GC_228})

V_393 = Vertex(name = 'V_393',
               particles = [ P.c__tilde__, P.c, P.Z, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_219})

V_394 = Vertex(name = 'V_394',
               particles = [ P.u__tilde__, P.u, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2, L.FFVS3, L.FFVS5, L.FFVS6 ],
               couplings = {(0,1):C.GC_413,(0,2):C.GC_423,(0,3):C.GC_1149,(0,0):C.GC_1147})

V_395 = Vertex(name = 'V_395',
               particles = [ P.u__tilde__, P.u, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2, L.FFVS3, L.FFVS6 ],
               couplings = {(0,1):C.GC_414,(0,2):C.GC_1162,(0,0):C.GC_1161})

V_396 = Vertex(name = 'V_396',
               particles = [ P.c__tilde__, P.c, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2, L.FFVS3, L.FFVS5, L.FFVS6 ],
               couplings = {(0,1):C.GC_413,(0,2):C.GC_423,(0,3):C.GC_605,(0,0):C.GC_603})

V_397 = Vertex(name = 'V_397',
               particles = [ P.c__tilde__, P.c, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2, L.FFVS3, L.FFVS6 ],
               couplings = {(0,1):C.GC_414,(0,2):C.GC_618,(0,0):C.GC_617})

V_398 = Vertex(name = 'V_398',
               particles = [ P.ve__tilde__, P.ve, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_203})

V_399 = Vertex(name = 'V_399',
               particles = [ P.ve__tilde__, P.ve, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_511})

V_400 = Vertex(name = 'V_400',
               particles = [ P.ve__tilde__, P.ve, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_518})

V_401 = Vertex(name = 'V_401',
               particles = [ P.ve__tilde__, P.ve, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_526})

V_402 = Vertex(name = 'V_402',
               particles = [ P.vm__tilde__, P.vm, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_203})

V_403 = Vertex(name = 'V_403',
               particles = [ P.vm__tilde__, P.vm, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_511})

V_404 = Vertex(name = 'V_404',
               particles = [ P.vm__tilde__, P.vm, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_518})

V_405 = Vertex(name = 'V_405',
               particles = [ P.vm__tilde__, P.vm, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_526})

V_406 = Vertex(name = 'V_406',
               particles = [ P.vt__tilde__, P.vt, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_203})

V_407 = Vertex(name = 'V_407',
               particles = [ P.vt__tilde__, P.vt, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_511})

V_408 = Vertex(name = 'V_408',
               particles = [ P.vt__tilde__, P.vt, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_518})

V_409 = Vertex(name = 'V_409',
               particles = [ P.vt__tilde__, P.vt, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_526})

V_410 = Vertex(name = 'V_410',
               particles = [ P.ve__tilde__, P.ve, P.Z, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_221})

V_411 = Vertex(name = 'V_411',
               particles = [ P.ve__tilde__, P.ve, P.Z, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_222})

V_412 = Vertex(name = 'V_412',
               particles = [ P.vm__tilde__, P.vm, P.Z, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_221})

V_413 = Vertex(name = 'V_413',
               particles = [ P.vm__tilde__, P.vm, P.Z, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_222})

V_414 = Vertex(name = 'V_414',
               particles = [ P.vt__tilde__, P.vt, P.Z, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_221})

V_415 = Vertex(name = 'V_415',
               particles = [ P.vt__tilde__, P.vt, P.Z, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_222})

V_416 = Vertex(name = 'V_416',
               particles = [ P.ve__tilde__, P.ve, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_416})

V_417 = Vertex(name = 'V_417',
               particles = [ P.ve__tilde__, P.ve, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_417})

V_418 = Vertex(name = 'V_418',
               particles = [ P.vm__tilde__, P.vm, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_416})

V_419 = Vertex(name = 'V_419',
               particles = [ P.vm__tilde__, P.vm, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_417})

V_420 = Vertex(name = 'V_420',
               particles = [ P.vt__tilde__, P.vt, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_416})

V_421 = Vertex(name = 'V_421',
               particles = [ P.vt__tilde__, P.vt, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_417})

V_422 = Vertex(name = 'V_422',
               particles = [ P.t1__tilde__, P.t1, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_166})

V_423 = Vertex(name = 'V_423',
               particles = [ P.t__tilde__, P.t1, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_161})

V_424 = Vertex(name = 'V_424',
               particles = [ P.t1__tilde__, P.t, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_161})

V_425 = Vertex(name = 'V_425',
               particles = [ P.t1__tilde__, P.t1, P.g ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_169})

V_426 = Vertex(name = 'V_426',
               particles = [ P.t__tilde__, P.t1, P.g ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_164})

V_427 = Vertex(name = 'V_427',
               particles = [ P.t1__tilde__, P.t, P.g ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_164})

V_428 = Vertex(name = 'V_428',
               particles = [ P.b__tilde__, P.t1, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_251})

V_429 = Vertex(name = 'V_429',
               particles = [ P.b__tilde__, P.t1, P.W1__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_256})

V_430 = Vertex(name = 'V_430',
               particles = [ P.b__tilde__, P.t, P.W1__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_251})

V_431 = Vertex(name = 'V_431',
               particles = [ P.d__tilde__, P.u, P.W1__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_251})

V_432 = Vertex(name = 'V_432',
               particles = [ P.s__tilde__, P.c, P.W1__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_251})

V_433 = Vertex(name = 'V_433',
               particles = [ P.e__plus__, P.ve, P.W1__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_251})

V_434 = Vertex(name = 'V_434',
               particles = [ P.mu__plus__, P.vm, P.W1__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_251})

V_435 = Vertex(name = 'V_435',
               particles = [ P.ta__plus__, P.vt, P.W1__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_251})

V_436 = Vertex(name = 'V_436',
               particles = [ P.t1__tilde__, P.b, P.W1__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_256})

V_437 = Vertex(name = 'V_437',
               particles = [ P.t__tilde__, P.b, P.W1__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_251})

V_438 = Vertex(name = 'V_438',
               particles = [ P.u__tilde__, P.d, P.W1__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_251})

V_439 = Vertex(name = 'V_439',
               particles = [ P.c__tilde__, P.s, P.W1__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_251})

V_440 = Vertex(name = 'V_440',
               particles = [ P.ve__tilde__, P.e__minus__, P.W1__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_251})

V_441 = Vertex(name = 'V_441',
               particles = [ P.vm__tilde__, P.mu__minus__, P.W1__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_251})

V_442 = Vertex(name = 'V_442',
               particles = [ P.vt__tilde__, P.ta__minus__, P.W1__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_251})

V_443 = Vertex(name = 'V_443',
               particles = [ P.t1__tilde__, P.b, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_251})

V_444 = Vertex(name = 'V_444',
               particles = [ P.t1__tilde__, P.t1, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_291,(0,1):C.GC_257})

V_445 = Vertex(name = 'V_445',
               particles = [ P.t__tilde__, P.t1, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_289,(0,1):C.GC_252})

V_446 = Vertex(name = 'V_446',
               particles = [ P.t1__tilde__, P.t, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_289,(0,1):C.GC_252})

V_447 = Vertex(name = 'V_447',
               particles = [ P.b__tilde__, P.b, P.Z1 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_288,(0,1):C.GC_253})

V_448 = Vertex(name = 'V_448',
               particles = [ P.d__tilde__, P.d, P.Z1 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_288,(0,1):C.GC_253})

V_449 = Vertex(name = 'V_449',
               particles = [ P.s__tilde__, P.s, P.Z1 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_288,(0,1):C.GC_253})

V_450 = Vertex(name = 'V_450',
               particles = [ P.e__plus__, P.e__minus__, P.Z1 ],
               color = [ '1' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_290,(0,1):C.GC_253})

V_451 = Vertex(name = 'V_451',
               particles = [ P.mu__plus__, P.mu__minus__, P.Z1 ],
               color = [ '1' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_290,(0,1):C.GC_253})

V_452 = Vertex(name = 'V_452',
               particles = [ P.ta__plus__, P.ta__minus__, P.Z1 ],
               color = [ '1' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_290,(0,1):C.GC_253})

V_453 = Vertex(name = 'V_453',
               particles = [ P.t1__tilde__, P.t1, P.Z1 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_292,(0,1):C.GC_260})

V_454 = Vertex(name = 'V_454',
               particles = [ P.t__tilde__, P.t1, P.Z1 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_291,(0,1):C.GC_257})

V_455 = Vertex(name = 'V_455',
               particles = [ P.t1__tilde__, P.t, P.Z1 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_291,(0,1):C.GC_257})

V_456 = Vertex(name = 'V_456',
               particles = [ P.t__tilde__, P.t, P.Z1 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_289,(0,1):C.GC_252})

V_457 = Vertex(name = 'V_457',
               particles = [ P.u__tilde__, P.u, P.Z1 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_289,(0,1):C.GC_252})

V_458 = Vertex(name = 'V_458',
               particles = [ P.c__tilde__, P.c, P.Z1 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1, L.FFV3 ],
               couplings = {(0,0):C.GC_289,(0,1):C.GC_252})

V_459 = Vertex(name = 'V_459',
               particles = [ P.ve__tilde__, P.ve, P.Z1 ],
               color = [ '1' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_252})

V_460 = Vertex(name = 'V_460',
               particles = [ P.vm__tilde__, P.vm, P.Z1 ],
               color = [ '1' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_252})

V_461 = Vertex(name = 'V_461',
               particles = [ P.vt__tilde__, P.vt, P.Z1 ],
               color = [ '1' ],
               lorentz = [ L.FFV3 ],
               couplings = {(0,0):C.GC_252})

V_462 = Vertex(name = 'V_462',
               particles = [ P.b__tilde__, P.b, P.H, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSSS1, L.FFSSS2 ],
               couplings = {(0,0):C.GC_16,(0,1):C.GC_17})

V_463 = Vertex(name = 'V_463',
               particles = [ P.b__tilde__, P.b, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_339,(0,1):C.GC_340})

V_464 = Vertex(name = 'V_464',
               particles = [ P.b__tilde__, P.b, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS1, L.FFS2 ],
               couplings = {(0,0):C.GC_457,(0,1):C.GC_458})

V_465 = Vertex(name = 'V_465',
               particles = [ P.b__tilde__, P.b, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_566})

V_466 = Vertex(name = 'V_466',
               particles = [ P.b__tilde__, P.b, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_568})

V_467 = Vertex(name = 'V_467',
               particles = [ P.b__tilde__, P.b, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_569})

V_468 = Vertex(name = 'V_468',
               particles = [ P.b__tilde__, P.b, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_570})

V_469 = Vertex(name = 'V_469',
               particles = [ P.b__tilde__, P.b, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_571})

V_470 = Vertex(name = 'V_470',
               particles = [ P.t__tilde__, P.t, P.H, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSSS1, L.FFSSS2 ],
               couplings = {(0,0):C.GC_110,(0,1):C.GC_111})

V_471 = Vertex(name = 'V_471',
               particles = [ P.t__tilde__, P.t, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_363,(0,1):C.GC_364})

V_472 = Vertex(name = 'V_472',
               particles = [ P.t__tilde__, P.t, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS1, L.FFS2 ],
               couplings = {(0,0):C.GC_460,(0,1):C.GC_461})

V_473 = Vertex(name = 'V_473',
               particles = [ P.t__tilde__, P.t, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_1025})

V_474 = Vertex(name = 'V_474',
               particles = [ P.t__tilde__, P.t, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_1029})

V_475 = Vertex(name = 'V_475',
               particles = [ P.t__tilde__, P.t, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_1030})

V_476 = Vertex(name = 'V_476',
               particles = [ P.t__tilde__, P.t, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_1031})

V_477 = Vertex(name = 'V_477',
               particles = [ P.t__tilde__, P.t, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_1032})

V_478 = Vertex(name = 'V_478',
               particles = [ P.d__tilde__, P.d, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS1, L.FFS2 ],
               couplings = {(0,0):C.GC_726,(0,1):C.GC_652})

V_479 = Vertex(name = 'V_479',
               particles = [ P.d__tilde__, P.d, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_727})

V_480 = Vertex(name = 'V_480',
               particles = [ P.d__tilde__, P.d, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_728})

V_481 = Vertex(name = 'V_481',
               particles = [ P.d__tilde__, P.d, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_729})

V_482 = Vertex(name = 'V_482',
               particles = [ P.d__tilde__, P.d, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_730})

V_483 = Vertex(name = 'V_483',
               particles = [ P.d__tilde__, P.d, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_731})

V_484 = Vertex(name = 'V_484',
               particles = [ P.s__tilde__, P.s, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS1, L.FFS2 ],
               couplings = {(0,0):C.GC_982,(0,1):C.GC_908})

V_485 = Vertex(name = 'V_485',
               particles = [ P.s__tilde__, P.s, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_983})

V_486 = Vertex(name = 'V_486',
               particles = [ P.s__tilde__, P.s, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_984})

V_487 = Vertex(name = 'V_487',
               particles = [ P.s__tilde__, P.s, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_985})

V_488 = Vertex(name = 'V_488',
               particles = [ P.s__tilde__, P.s, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_986})

V_489 = Vertex(name = 'V_489',
               particles = [ P.s__tilde__, P.s, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_987})

V_490 = Vertex(name = 'V_490',
               particles = [ P.d__tilde__, P.d, P.H, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSSS1, L.FFSSS2 ],
               couplings = {(0,0):C.GC_655,(0,1):C.GC_656})

V_491 = Vertex(name = 'V_491',
               particles = [ P.s__tilde__, P.s, P.H, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSSS1, L.FFSSS2 ],
               couplings = {(0,0):C.GC_911,(0,1):C.GC_912})

V_492 = Vertex(name = 'V_492',
               particles = [ P.d__tilde__, P.d, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_703,(0,1):C.GC_704})

V_493 = Vertex(name = 'V_493',
               particles = [ P.s__tilde__, P.s, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_959,(0,1):C.GC_960})

V_494 = Vertex(name = 'V_494',
               particles = [ P.e__plus__, P.e__minus__, P.H, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFSSS1, L.FFSSS2 ],
               couplings = {(0,0):C.GC_749,(0,1):C.GC_750})

V_495 = Vertex(name = 'V_495',
               particles = [ P.mu__plus__, P.mu__minus__, P.H, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFSSS1, L.FFSSS2 ],
               couplings = {(0,0):C.GC_829,(0,1):C.GC_830})

V_496 = Vertex(name = 'V_496',
               particles = [ P.ta__plus__, P.ta__minus__, P.H, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFSSS1, L.FFSSS2 ],
               couplings = {(0,0):C.GC_1034,(0,1):C.GC_1035})

V_497 = Vertex(name = 'V_497',
               particles = [ P.e__plus__, P.e__minus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_786,(0,1):C.GC_787})

V_498 = Vertex(name = 'V_498',
               particles = [ P.mu__plus__, P.mu__minus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_866,(0,1):C.GC_867})

V_499 = Vertex(name = 'V_499',
               particles = [ P.ta__plus__, P.ta__minus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_1071,(0,1):C.GC_1072})

V_500 = Vertex(name = 'V_500',
               particles = [ P.e__plus__, P.e__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFS1, L.FFS2 ],
               couplings = {(0,0):C.GC_807,(0,1):C.GC_748})

V_501 = Vertex(name = 'V_501',
               particles = [ P.e__plus__, P.e__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_808})

V_502 = Vertex(name = 'V_502',
               particles = [ P.e__plus__, P.e__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_809})

V_503 = Vertex(name = 'V_503',
               particles = [ P.e__plus__, P.e__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_810})

V_504 = Vertex(name = 'V_504',
               particles = [ P.e__plus__, P.e__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_811})

V_505 = Vertex(name = 'V_505',
               particles = [ P.e__plus__, P.e__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_812})

V_506 = Vertex(name = 'V_506',
               particles = [ P.mu__plus__, P.mu__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFS1, L.FFS2 ],
               couplings = {(0,0):C.GC_887,(0,1):C.GC_828})

V_507 = Vertex(name = 'V_507',
               particles = [ P.mu__plus__, P.mu__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_888})

V_508 = Vertex(name = 'V_508',
               particles = [ P.mu__plus__, P.mu__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_889})

V_509 = Vertex(name = 'V_509',
               particles = [ P.mu__plus__, P.mu__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_890})

V_510 = Vertex(name = 'V_510',
               particles = [ P.mu__plus__, P.mu__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_891})

V_511 = Vertex(name = 'V_511',
               particles = [ P.mu__plus__, P.mu__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_892})

V_512 = Vertex(name = 'V_512',
               particles = [ P.ta__plus__, P.ta__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFS1, L.FFS2 ],
               couplings = {(0,0):C.GC_1092,(0,1):C.GC_1033})

V_513 = Vertex(name = 'V_513',
               particles = [ P.ta__plus__, P.ta__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_1093})

V_514 = Vertex(name = 'V_514',
               particles = [ P.ta__plus__, P.ta__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_1094})

V_515 = Vertex(name = 'V_515',
               particles = [ P.ta__plus__, P.ta__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_1095})

V_516 = Vertex(name = 'V_516',
               particles = [ P.ta__plus__, P.ta__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_1096})

V_517 = Vertex(name = 'V_517',
               particles = [ P.ta__plus__, P.ta__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_1097})

V_518 = Vertex(name = 'V_518',
               particles = [ P.u__tilde__, P.u, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS1, L.FFS2 ],
               couplings = {(0,0):C.GC_1194,(0,1):C.GC_1116})

V_519 = Vertex(name = 'V_519',
               particles = [ P.u__tilde__, P.u, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_1190})

V_520 = Vertex(name = 'V_520',
               particles = [ P.u__tilde__, P.u, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_1191})

V_521 = Vertex(name = 'V_521',
               particles = [ P.u__tilde__, P.u, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_1192})

V_522 = Vertex(name = 'V_522',
               particles = [ P.u__tilde__, P.u, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_1193})

V_523 = Vertex(name = 'V_523',
               particles = [ P.u__tilde__, P.u, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_1195})

V_524 = Vertex(name = 'V_524',
               particles = [ P.c__tilde__, P.c, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS1, L.FFS2 ],
               couplings = {(0,0):C.GC_650,(0,1):C.GC_572})

V_525 = Vertex(name = 'V_525',
               particles = [ P.c__tilde__, P.c, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_646})

V_526 = Vertex(name = 'V_526',
               particles = [ P.c__tilde__, P.c, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_647})

V_527 = Vertex(name = 'V_527',
               particles = [ P.c__tilde__, P.c, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_648})

V_528 = Vertex(name = 'V_528',
               particles = [ P.c__tilde__, P.c, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_649})

V_529 = Vertex(name = 'V_529',
               particles = [ P.c__tilde__, P.c, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_651})

V_530 = Vertex(name = 'V_530',
               particles = [ P.u__tilde__, P.u, P.H, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSSS1, L.FFSSS2 ],
               couplings = {(0,0):C.GC_1143,(0,1):C.GC_1144})

V_531 = Vertex(name = 'V_531',
               particles = [ P.c__tilde__, P.c, P.H, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSSS1, L.FFSSS2 ],
               couplings = {(0,0):C.GC_599,(0,1):C.GC_600})

V_532 = Vertex(name = 'V_532',
               particles = [ P.u__tilde__, P.u, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_1169,(0,1):C.GC_1170})

V_533 = Vertex(name = 'V_533',
               particles = [ P.c__tilde__, P.c, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_625,(0,1):C.GC_626})

V_534 = Vertex(name = 'V_534',
               particles = [ P.b__tilde__, P.b, P.H1 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_567})

V_535 = Vertex(name = 'V_535',
               particles = [ P.t1__tilde__, P.t1, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_1027})

V_536 = Vertex(name = 'V_536',
               particles = [ P.t1__tilde__, P.t1, P.H1 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_1028})

V_537 = Vertex(name = 'V_537',
               particles = [ P.t__tilde__, P.t1, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_1026})

V_538 = Vertex(name = 'V_538',
               particles = [ P.t__tilde__, P.t1, P.H1 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_1027})

V_539 = Vertex(name = 'V_539',
               particles = [ P.t1__tilde__, P.t, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_1026})

V_540 = Vertex(name = 'V_540',
               particles = [ P.t1__tilde__, P.t, P.H1 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_1027})

V_541 = Vertex(name = 'V_541',
               particles = [ P.t__tilde__, P.t, P.H1 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_1026})

V_542 = Vertex(name = 'V_542',
               particles = [ P.d__tilde__, P.d, P.H1 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_691})

V_543 = Vertex(name = 'V_543',
               particles = [ P.s__tilde__, P.s, P.H1 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_947})

V_544 = Vertex(name = 'V_544',
               particles = [ P.e__plus__, P.e__minus__, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_776})

V_545 = Vertex(name = 'V_545',
               particles = [ P.mu__plus__, P.mu__minus__, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_856})

V_546 = Vertex(name = 'V_546',
               particles = [ P.ta__plus__, P.ta__minus__, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_1061})

V_547 = Vertex(name = 'V_547',
               particles = [ P.u__tilde__, P.u, P.H1 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_1155})

V_548 = Vertex(name = 'V_548',
               particles = [ P.c__tilde__, P.c, P.H1 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               couplings = {(0,0):C.GC_611})

V_549 = Vertex(name = 'V_549',
               particles = [ P.b__tilde__, P.b, P.a, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2, L.FFVS6 ],
               couplings = {(0,1):C.GC_101,(0,0):C.GC_100})

V_550 = Vertex(name = 'V_550',
               particles = [ P.b__tilde__, P.b, P.a, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2, L.FFVS6 ],
               couplings = {(0,1):C.GC_269,(0,0):C.GC_268})

V_551 = Vertex(name = 'V_551',
               particles = [ P.e__plus__, P.e__minus__, P.a, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS2, L.FFVS6 ],
               couplings = {(0,1):C.GC_770,(0,0):C.GC_769})

V_552 = Vertex(name = 'V_552',
               particles = [ P.e__plus__, P.e__minus__, P.a, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS2, L.FFVS6 ],
               couplings = {(0,1):C.GC_785,(0,0):C.GC_784})

V_553 = Vertex(name = 'V_553',
               particles = [ P.mu__plus__, P.mu__minus__, P.a, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS2, L.FFVS6 ],
               couplings = {(0,1):C.GC_850,(0,0):C.GC_849})

V_554 = Vertex(name = 'V_554',
               particles = [ P.mu__plus__, P.mu__minus__, P.a, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS2, L.FFVS6 ],
               couplings = {(0,1):C.GC_865,(0,0):C.GC_864})

V_555 = Vertex(name = 'V_555',
               particles = [ P.ta__plus__, P.ta__minus__, P.a, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS2, L.FFVS6 ],
               couplings = {(0,1):C.GC_1055,(0,0):C.GC_1054})

V_556 = Vertex(name = 'V_556',
               particles = [ P.ta__plus__, P.ta__minus__, P.a, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS2, L.FFVS6 ],
               couplings = {(0,1):C.GC_1070,(0,0):C.GC_1069})

V_557 = Vertex(name = 'V_557',
               particles = [ P.t__tilde__, P.t, P.a, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2, L.FFVS6 ],
               couplings = {(0,2):C.GC_105,(0,1):C.GC_104,(0,0):C.GC_280})

V_558 = Vertex(name = 'V_558',
               particles = [ P.t__tilde__, P.t, P.a, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS6 ],
               couplings = {(0,0):C.GC_281})

V_559 = Vertex(name = 'V_559',
               particles = [ P.b__tilde__, P.b, P.g, P.H ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFVS2, L.FFVS6 ],
               couplings = {(0,1):C.GC_15,(0,0):C.GC_14})

V_560 = Vertex(name = 'V_560',
               particles = [ P.b__tilde__, P.b, P.g, P.g, P.H ],
               color = [ 'f(-1,3,4)*T(-1,2,1)' ],
               lorentz = [ L.FFVVS1, L.FFVVS3 ],
               couplings = {(0,1):C.GC_145,(0,0):C.GC_144})

V_561 = Vertex(name = 'V_561',
               particles = [ P.b__tilde__, P.b, P.g, P.g ],
               color = [ 'f(-1,3,4)*T(-1,2,1)' ],
               lorentz = [ L.FFVV1, L.FFVV3 ],
               couplings = {(0,1):C.GC_382,(0,0):C.GC_381})

V_562 = Vertex(name = 'V_562',
               particles = [ P.t__tilde__, P.t, P.g, P.H ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFVS2, L.FFVS6 ],
               couplings = {(0,1):C.GC_99,(0,0):C.GC_98})

V_563 = Vertex(name = 'V_563',
               particles = [ P.t__tilde__, P.t, P.g, P.g, P.H ],
               color = [ 'f(-1,3,4)*T(-1,2,1)' ],
               lorentz = [ L.FFVVS1, L.FFVVS3 ],
               couplings = {(0,1):C.GC_151,(0,0):C.GC_150})

V_564 = Vertex(name = 'V_564',
               particles = [ P.t__tilde__, P.t, P.g, P.g ],
               color = [ 'f(-1,3,4)*T(-1,2,1)' ],
               lorentz = [ L.FFVV1, L.FFVV3 ],
               couplings = {(0,1):C.GC_386,(0,0):C.GC_385})

V_565 = Vertex(name = 'V_565',
               particles = [ P.d__tilde__, P.d, P.a, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2, L.FFVS6 ],
               couplings = {(0,1):C.GC_683,(0,0):C.GC_682})

V_566 = Vertex(name = 'V_566',
               particles = [ P.d__tilde__, P.d, P.a, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2, L.FFVS6 ],
               couplings = {(0,1):C.GC_700,(0,0):C.GC_699})

V_567 = Vertex(name = 'V_567',
               particles = [ P.s__tilde__, P.s, P.a, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2, L.FFVS6 ],
               couplings = {(0,1):C.GC_939,(0,0):C.GC_938})

V_568 = Vertex(name = 'V_568',
               particles = [ P.s__tilde__, P.s, P.a, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2, L.FFVS6 ],
               couplings = {(0,1):C.GC_956,(0,0):C.GC_955})

V_569 = Vertex(name = 'V_569',
               particles = [ P.d__tilde__, P.d, P.g, P.H ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFVS2, L.FFVS6 ],
               couplings = {(0,1):C.GC_654,(0,0):C.GC_653})

V_570 = Vertex(name = 'V_570',
               particles = [ P.s__tilde__, P.s, P.g, P.H ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFVS2, L.FFVS6 ],
               couplings = {(0,1):C.GC_910,(0,0):C.GC_909})

V_571 = Vertex(name = 'V_571',
               particles = [ P.d__tilde__, P.d, P.g, P.g, P.H ],
               color = [ 'f(-1,3,4)*T(-1,2,1)' ],
               lorentz = [ L.FFVVS1, L.FFVVS3 ],
               couplings = {(0,1):C.GC_690,(0,0):C.GC_689})

V_572 = Vertex(name = 'V_572',
               particles = [ P.s__tilde__, P.s, P.g, P.g, P.H ],
               color = [ 'f(-1,3,4)*T(-1,2,1)' ],
               lorentz = [ L.FFVVS1, L.FFVVS3 ],
               couplings = {(0,1):C.GC_946,(0,0):C.GC_945})

V_573 = Vertex(name = 'V_573',
               particles = [ P.d__tilde__, P.d, P.g, P.g ],
               color = [ 'f(-1,3,4)*T(-1,2,1)' ],
               lorentz = [ L.FFVV1, L.FFVV3 ],
               couplings = {(0,1):C.GC_716,(0,0):C.GC_715})

V_574 = Vertex(name = 'V_574',
               particles = [ P.s__tilde__, P.s, P.g, P.g ],
               color = [ 'f(-1,3,4)*T(-1,2,1)' ],
               lorentz = [ L.FFVV1, L.FFVV3 ],
               couplings = {(0,1):C.GC_972,(0,0):C.GC_971})

V_575 = Vertex(name = 'V_575',
               particles = [ P.u__tilde__, P.u, P.a, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2, L.FFVS6 ],
               couplings = {(0,1):C.GC_1140,(0,0):C.GC_1139})

V_576 = Vertex(name = 'V_576',
               particles = [ P.u__tilde__, P.u, P.a, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2, L.FFVS6 ],
               couplings = {(0,1):C.GC_1164,(0,0):C.GC_1163})

V_577 = Vertex(name = 'V_577',
               particles = [ P.c__tilde__, P.c, P.a, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2, L.FFVS6 ],
               couplings = {(0,1):C.GC_596,(0,0):C.GC_595})

V_578 = Vertex(name = 'V_578',
               particles = [ P.c__tilde__, P.c, P.a, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2, L.FFVS6 ],
               couplings = {(0,1):C.GC_620,(0,0):C.GC_619})

V_579 = Vertex(name = 'V_579',
               particles = [ P.u__tilde__, P.u, P.g, P.H ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFVS2, L.FFVS6 ],
               couplings = {(0,1):C.GC_1142,(0,0):C.GC_1141})

V_580 = Vertex(name = 'V_580',
               particles = [ P.c__tilde__, P.c, P.g, P.H ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFVS2, L.FFVS6 ],
               couplings = {(0,1):C.GC_598,(0,0):C.GC_597})

V_581 = Vertex(name = 'V_581',
               particles = [ P.u__tilde__, P.u, P.g, P.g, P.H ],
               color = [ 'f(-1,3,4)*T(-1,2,1)' ],
               lorentz = [ L.FFVVS1, L.FFVVS3 ],
               couplings = {(0,1):C.GC_1154,(0,0):C.GC_1153})

V_582 = Vertex(name = 'V_582',
               particles = [ P.c__tilde__, P.c, P.g, P.g, P.H ],
               color = [ 'f(-1,3,4)*T(-1,2,1)' ],
               lorentz = [ L.FFVVS1, L.FFVVS3 ],
               couplings = {(0,1):C.GC_610,(0,0):C.GC_609})

V_583 = Vertex(name = 'V_583',
               particles = [ P.u__tilde__, P.u, P.g, P.g ],
               color = [ 'f(-1,3,4)*T(-1,2,1)' ],
               lorentz = [ L.FFVV1, L.FFVV3 ],
               couplings = {(0,1):C.GC_1180,(0,0):C.GC_1179})

V_584 = Vertex(name = 'V_584',
               particles = [ P.c__tilde__, P.c, P.g, P.g ],
               color = [ 'f(-1,3,4)*T(-1,2,1)' ],
               lorentz = [ L.FFVV1, L.FFVV3 ],
               couplings = {(0,1):C.GC_636,(0,0):C.GC_635})

V_585 = Vertex(name = 'V_585',
               particles = [ P.t__tilde__, P.b, P.a, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS2, L.FFVVS4 ],
               couplings = {(0,0):C.GC_129,(0,1):C.GC_136})

V_586 = Vertex(name = 'V_586',
               particles = [ P.t__tilde__, P.b, P.a, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS2, L.FFVVS4 ],
               couplings = {(0,0):C.GC_131,(0,1):C.GC_138})

V_587 = Vertex(name = 'V_587',
               particles = [ P.t__tilde__, P.b, P.a, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV2, L.FFVV4 ],
               couplings = {(0,0):C.GC_370,(0,1):C.GC_377})

V_588 = Vertex(name = 'V_588',
               particles = [ P.t__tilde__, P.b, P.a, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV2, L.FFVV4 ],
               couplings = {(0,0):C.GC_372,(0,1):C.GC_379})

V_589 = Vertex(name = 'V_589',
               particles = [ P.b__tilde__, P.b, P.W__minus__, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS1, L.FFVVS3 ],
               couplings = {(0,1):C.GC_208,(0,0):C.GC_207})

V_590 = Vertex(name = 'V_590',
               particles = [ P.b__tilde__, P.b, P.W__minus__, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV1, L.FFVV3 ],
               couplings = {(0,1):C.GC_403,(0,0):C.GC_402})

V_591 = Vertex(name = 'V_591',
               particles = [ P.t__tilde__, P.b, P.W__plus__, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS2, L.FFVVS4 ],
               couplings = {(0,0):C.GC_229,(0,1):C.GC_237})

V_592 = Vertex(name = 'V_592',
               particles = [ P.t__tilde__, P.b, P.W__plus__, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS2, L.FFVVS4 ],
               couplings = {(0,0):C.GC_230,(0,1):C.GC_239})

V_593 = Vertex(name = 'V_593',
               particles = [ P.t__tilde__, P.b, P.W__plus__, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV2, L.FFVV4 ],
               couplings = {(0,0):C.GC_424,(0,1):C.GC_432})

V_594 = Vertex(name = 'V_594',
               particles = [ P.t__tilde__, P.b, P.W__plus__, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV2, L.FFVV4 ],
               couplings = {(0,0):C.GC_425,(0,1):C.GC_434})

V_595 = Vertex(name = 'V_595',
               particles = [ P.ve__tilde__, P.e__minus__, P.a, P.W__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVVS2 ],
               couplings = {(0,0):C.GC_773})

V_596 = Vertex(name = 'V_596',
               particles = [ P.ve__tilde__, P.e__minus__, P.a, P.W__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVVS2 ],
               couplings = {(0,0):C.GC_775})

V_597 = Vertex(name = 'V_597',
               particles = [ P.vm__tilde__, P.mu__minus__, P.a, P.W__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVVS2 ],
               couplings = {(0,0):C.GC_853})

V_598 = Vertex(name = 'V_598',
               particles = [ P.vm__tilde__, P.mu__minus__, P.a, P.W__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVVS2 ],
               couplings = {(0,0):C.GC_855})

V_599 = Vertex(name = 'V_599',
               particles = [ P.vt__tilde__, P.ta__minus__, P.a, P.W__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVVS2 ],
               couplings = {(0,0):C.GC_1058})

V_600 = Vertex(name = 'V_600',
               particles = [ P.vt__tilde__, P.ta__minus__, P.a, P.W__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVVS2 ],
               couplings = {(0,0):C.GC_1060})

V_601 = Vertex(name = 'V_601',
               particles = [ P.ve__tilde__, P.e__minus__, P.a, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVV2 ],
               couplings = {(0,0):C.GC_795})

V_602 = Vertex(name = 'V_602',
               particles = [ P.ve__tilde__, P.e__minus__, P.a, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVV2 ],
               couplings = {(0,0):C.GC_797})

V_603 = Vertex(name = 'V_603',
               particles = [ P.vm__tilde__, P.mu__minus__, P.a, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVV2 ],
               couplings = {(0,0):C.GC_875})

V_604 = Vertex(name = 'V_604',
               particles = [ P.vm__tilde__, P.mu__minus__, P.a, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVV2 ],
               couplings = {(0,0):C.GC_877})

V_605 = Vertex(name = 'V_605',
               particles = [ P.vt__tilde__, P.ta__minus__, P.a, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVV2 ],
               couplings = {(0,0):C.GC_1080})

V_606 = Vertex(name = 'V_606',
               particles = [ P.vt__tilde__, P.ta__minus__, P.a, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVV2 ],
               couplings = {(0,0):C.GC_1082})

V_607 = Vertex(name = 'V_607',
               particles = [ P.e__plus__, P.e__minus__, P.W__minus__, P.W__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVVS1, L.FFVVS3 ],
               couplings = {(0,1):C.GC_778,(0,0):C.GC_777})

V_608 = Vertex(name = 'V_608',
               particles = [ P.mu__plus__, P.mu__minus__, P.W__minus__, P.W__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVVS1, L.FFVVS3 ],
               couplings = {(0,1):C.GC_858,(0,0):C.GC_857})

V_609 = Vertex(name = 'V_609',
               particles = [ P.ta__plus__, P.ta__minus__, P.W__minus__, P.W__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVVS1, L.FFVVS3 ],
               couplings = {(0,1):C.GC_1063,(0,0):C.GC_1062})

V_610 = Vertex(name = 'V_610',
               particles = [ P.e__plus__, P.e__minus__, P.W__minus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVV1, L.FFVV3 ],
               couplings = {(0,1):C.GC_799,(0,0):C.GC_798})

V_611 = Vertex(name = 'V_611',
               particles = [ P.mu__plus__, P.mu__minus__, P.W__minus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVV1, L.FFVV3 ],
               couplings = {(0,1):C.GC_879,(0,0):C.GC_878})

V_612 = Vertex(name = 'V_612',
               particles = [ P.ta__plus__, P.ta__minus__, P.W__minus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVV1, L.FFVV3 ],
               couplings = {(0,1):C.GC_1084,(0,0):C.GC_1083})

V_613 = Vertex(name = 'V_613',
               particles = [ P.ve__tilde__, P.e__minus__, P.W__plus__, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVVS2 ],
               couplings = {(0,0):C.GC_779})

V_614 = Vertex(name = 'V_614',
               particles = [ P.ve__tilde__, P.e__minus__, P.W__plus__, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVVS2 ],
               couplings = {(0,0):C.GC_780})

V_615 = Vertex(name = 'V_615',
               particles = [ P.vm__tilde__, P.mu__minus__, P.W__plus__, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVVS2 ],
               couplings = {(0,0):C.GC_859})

V_616 = Vertex(name = 'V_616',
               particles = [ P.vm__tilde__, P.mu__minus__, P.W__plus__, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVVS2 ],
               couplings = {(0,0):C.GC_860})

V_617 = Vertex(name = 'V_617',
               particles = [ P.vt__tilde__, P.ta__minus__, P.W__plus__, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVVS2 ],
               couplings = {(0,0):C.GC_1064})

V_618 = Vertex(name = 'V_618',
               particles = [ P.vt__tilde__, P.ta__minus__, P.W__plus__, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVVS2 ],
               couplings = {(0,0):C.GC_1065})

V_619 = Vertex(name = 'V_619',
               particles = [ P.ve__tilde__, P.e__minus__, P.W__plus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFVV2 ],
               couplings = {(0,0):C.GC_800})

V_620 = Vertex(name = 'V_620',
               particles = [ P.ve__tilde__, P.e__minus__, P.W__plus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFVV2 ],
               couplings = {(0,0):C.GC_801})

V_621 = Vertex(name = 'V_621',
               particles = [ P.vm__tilde__, P.mu__minus__, P.W__plus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFVV2 ],
               couplings = {(0,0):C.GC_880})

V_622 = Vertex(name = 'V_622',
               particles = [ P.vm__tilde__, P.mu__minus__, P.W__plus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFVV2 ],
               couplings = {(0,0):C.GC_881})

V_623 = Vertex(name = 'V_623',
               particles = [ P.vt__tilde__, P.ta__minus__, P.W__plus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFVV2 ],
               couplings = {(0,0):C.GC_1085})

V_624 = Vertex(name = 'V_624',
               particles = [ P.vt__tilde__, P.ta__minus__, P.W__plus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFVV2 ],
               couplings = {(0,0):C.GC_1086})

V_625 = Vertex(name = 'V_625',
               particles = [ P.b__tilde__, P.t, P.a, P.W__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS2, L.FFVVS4 ],
               couplings = {(0,0):C.GC_136,(0,1):C.GC_129})

V_626 = Vertex(name = 'V_626',
               particles = [ P.b__tilde__, P.t, P.a, P.W__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS2, L.FFVVS4 ],
               couplings = {(0,0):C.GC_137,(0,1):C.GC_130})

V_627 = Vertex(name = 'V_627',
               particles = [ P.b__tilde__, P.t, P.a, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV2, L.FFVV4 ],
               couplings = {(0,0):C.GC_377,(0,1):C.GC_370})

V_628 = Vertex(name = 'V_628',
               particles = [ P.b__tilde__, P.t, P.a, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV2, L.FFVV4 ],
               couplings = {(0,0):C.GC_378,(0,1):C.GC_371})

V_629 = Vertex(name = 'V_629',
               particles = [ P.t__tilde__, P.t, P.W__minus__, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS1, L.FFVVS3 ],
               couplings = {(0,1):C.GC_238,(0,0):C.GC_236})

V_630 = Vertex(name = 'V_630',
               particles = [ P.t__tilde__, P.t, P.W__minus__, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV1, L.FFVV3 ],
               couplings = {(0,1):C.GC_433,(0,0):C.GC_431})

V_631 = Vertex(name = 'V_631',
               particles = [ P.b__tilde__, P.t, P.W__minus__, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS2, L.FFVVS4 ],
               couplings = {(0,0):C.GC_237,(0,1):C.GC_229})

V_632 = Vertex(name = 'V_632',
               particles = [ P.b__tilde__, P.t, P.W__minus__, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS2, L.FFVVS4 ],
               couplings = {(0,0):C.GC_240,(0,1):C.GC_231})

V_633 = Vertex(name = 'V_633',
               particles = [ P.b__tilde__, P.t, P.W__minus__, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV2, L.FFVV4 ],
               couplings = {(0,0):C.GC_432,(0,1):C.GC_424})

V_634 = Vertex(name = 'V_634',
               particles = [ P.b__tilde__, P.t, P.W__minus__, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV2, L.FFVV4 ],
               couplings = {(0,0):C.GC_435,(0,1):C.GC_426})

V_635 = Vertex(name = 'V_635',
               particles = [ P.u__tilde__, P.d, P.a, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS2, L.FFVVS4 ],
               couplings = {(0,0):C.GC_686,(0,1):C.GC_1150})

V_636 = Vertex(name = 'V_636',
               particles = [ P.u__tilde__, P.d, P.a, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS2, L.FFVVS4 ],
               couplings = {(0,0):C.GC_688,(0,1):C.GC_1152})

V_637 = Vertex(name = 'V_637',
               particles = [ P.c__tilde__, P.s, P.a, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS2, L.FFVVS4 ],
               couplings = {(0,0):C.GC_942,(0,1):C.GC_606})

V_638 = Vertex(name = 'V_638',
               particles = [ P.c__tilde__, P.s, P.a, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS2, L.FFVVS4 ],
               couplings = {(0,0):C.GC_944,(0,1):C.GC_608})

V_639 = Vertex(name = 'V_639',
               particles = [ P.u__tilde__, P.d, P.a, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV2, L.FFVV4 ],
               couplings = {(0,0):C.GC_712,(0,1):C.GC_1176})

V_640 = Vertex(name = 'V_640',
               particles = [ P.u__tilde__, P.d, P.a, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV2, L.FFVV4 ],
               couplings = {(0,0):C.GC_714,(0,1):C.GC_1178})

V_641 = Vertex(name = 'V_641',
               particles = [ P.c__tilde__, P.s, P.a, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV2, L.FFVV4 ],
               couplings = {(0,0):C.GC_968,(0,1):C.GC_632})

V_642 = Vertex(name = 'V_642',
               particles = [ P.c__tilde__, P.s, P.a, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV2, L.FFVV4 ],
               couplings = {(0,0):C.GC_970,(0,1):C.GC_634})

V_643 = Vertex(name = 'V_643',
               particles = [ P.d__tilde__, P.d, P.W__minus__, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS1, L.FFVVS3 ],
               couplings = {(0,1):C.GC_693,(0,0):C.GC_692})

V_644 = Vertex(name = 'V_644',
               particles = [ P.s__tilde__, P.s, P.W__minus__, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS1, L.FFVVS3 ],
               couplings = {(0,1):C.GC_949,(0,0):C.GC_948})

V_645 = Vertex(name = 'V_645',
               particles = [ P.d__tilde__, P.d, P.W__minus__, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV1, L.FFVV3 ],
               couplings = {(0,1):C.GC_718,(0,0):C.GC_717})

V_646 = Vertex(name = 'V_646',
               particles = [ P.s__tilde__, P.s, P.W__minus__, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV1, L.FFVV3 ],
               couplings = {(0,1):C.GC_974,(0,0):C.GC_973})

V_647 = Vertex(name = 'V_647',
               particles = [ P.u__tilde__, P.d, P.W__plus__, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS2, L.FFVVS4 ],
               couplings = {(0,0):C.GC_694,(0,1):C.GC_1157})

V_648 = Vertex(name = 'V_648',
               particles = [ P.u__tilde__, P.d, P.W__plus__, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS2, L.FFVVS4 ],
               couplings = {(0,0):C.GC_695,(0,1):C.GC_1159})

V_649 = Vertex(name = 'V_649',
               particles = [ P.c__tilde__, P.s, P.W__plus__, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS2, L.FFVVS4 ],
               couplings = {(0,0):C.GC_950,(0,1):C.GC_613})

V_650 = Vertex(name = 'V_650',
               particles = [ P.c__tilde__, P.s, P.W__plus__, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS2, L.FFVVS4 ],
               couplings = {(0,0):C.GC_951,(0,1):C.GC_615})

V_651 = Vertex(name = 'V_651',
               particles = [ P.u__tilde__, P.d, P.W__plus__, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV2, L.FFVV4 ],
               couplings = {(0,0):C.GC_719,(0,1):C.GC_1182})

V_652 = Vertex(name = 'V_652',
               particles = [ P.u__tilde__, P.d, P.W__plus__, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV2, L.FFVV4 ],
               couplings = {(0,0):C.GC_720,(0,1):C.GC_1184})

V_653 = Vertex(name = 'V_653',
               particles = [ P.c__tilde__, P.s, P.W__plus__, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV2, L.FFVV4 ],
               couplings = {(0,0):C.GC_975,(0,1):C.GC_638})

V_654 = Vertex(name = 'V_654',
               particles = [ P.c__tilde__, P.s, P.W__plus__, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV2, L.FFVV4 ],
               couplings = {(0,0):C.GC_976,(0,1):C.GC_640})

V_655 = Vertex(name = 'V_655',
               particles = [ P.d__tilde__, P.u, P.a, P.W__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS2, L.FFVVS4 ],
               couplings = {(0,0):C.GC_1150,(0,1):C.GC_686})

V_656 = Vertex(name = 'V_656',
               particles = [ P.d__tilde__, P.u, P.a, P.W__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS2, L.FFVVS4 ],
               couplings = {(0,0):C.GC_1151,(0,1):C.GC_687})

V_657 = Vertex(name = 'V_657',
               particles = [ P.s__tilde__, P.c, P.a, P.W__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS2, L.FFVVS4 ],
               couplings = {(0,0):C.GC_606,(0,1):C.GC_942})

V_658 = Vertex(name = 'V_658',
               particles = [ P.s__tilde__, P.c, P.a, P.W__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS2, L.FFVVS4 ],
               couplings = {(0,0):C.GC_607,(0,1):C.GC_943})

V_659 = Vertex(name = 'V_659',
               particles = [ P.d__tilde__, P.u, P.a, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV2, L.FFVV4 ],
               couplings = {(0,0):C.GC_1176,(0,1):C.GC_712})

V_660 = Vertex(name = 'V_660',
               particles = [ P.d__tilde__, P.u, P.a, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV2, L.FFVV4 ],
               couplings = {(0,0):C.GC_1177,(0,1):C.GC_713})

V_661 = Vertex(name = 'V_661',
               particles = [ P.s__tilde__, P.c, P.a, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV2, L.FFVV4 ],
               couplings = {(0,0):C.GC_632,(0,1):C.GC_968})

V_662 = Vertex(name = 'V_662',
               particles = [ P.s__tilde__, P.c, P.a, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV2, L.FFVV4 ],
               couplings = {(0,0):C.GC_633,(0,1):C.GC_969})

V_663 = Vertex(name = 'V_663',
               particles = [ P.u__tilde__, P.u, P.W__minus__, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS1, L.FFVVS3 ],
               couplings = {(0,1):C.GC_1158,(0,0):C.GC_1156})

V_664 = Vertex(name = 'V_664',
               particles = [ P.c__tilde__, P.c, P.W__minus__, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS1, L.FFVVS3 ],
               couplings = {(0,1):C.GC_614,(0,0):C.GC_612})

V_665 = Vertex(name = 'V_665',
               particles = [ P.u__tilde__, P.u, P.W__minus__, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV1, L.FFVV3 ],
               couplings = {(0,1):C.GC_1183,(0,0):C.GC_1181})

V_666 = Vertex(name = 'V_666',
               particles = [ P.c__tilde__, P.c, P.W__minus__, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV1, L.FFVV3 ],
               couplings = {(0,1):C.GC_639,(0,0):C.GC_637})

V_667 = Vertex(name = 'V_667',
               particles = [ P.d__tilde__, P.u, P.W__minus__, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS2, L.FFVVS4 ],
               couplings = {(0,0):C.GC_1157,(0,1):C.GC_694})

V_668 = Vertex(name = 'V_668',
               particles = [ P.d__tilde__, P.u, P.W__minus__, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS2, L.FFVVS4 ],
               couplings = {(0,0):C.GC_1160,(0,1):C.GC_696})

V_669 = Vertex(name = 'V_669',
               particles = [ P.s__tilde__, P.c, P.W__minus__, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS2, L.FFVVS4 ],
               couplings = {(0,0):C.GC_613,(0,1):C.GC_950})

V_670 = Vertex(name = 'V_670',
               particles = [ P.s__tilde__, P.c, P.W__minus__, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS2, L.FFVVS4 ],
               couplings = {(0,0):C.GC_616,(0,1):C.GC_952})

V_671 = Vertex(name = 'V_671',
               particles = [ P.d__tilde__, P.u, P.W__minus__, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV2, L.FFVV4 ],
               couplings = {(0,0):C.GC_1182,(0,1):C.GC_719})

V_672 = Vertex(name = 'V_672',
               particles = [ P.d__tilde__, P.u, P.W__minus__, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV2, L.FFVV4 ],
               couplings = {(0,0):C.GC_1185,(0,1):C.GC_721})

V_673 = Vertex(name = 'V_673',
               particles = [ P.s__tilde__, P.c, P.W__minus__, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV2, L.FFVV4 ],
               couplings = {(0,0):C.GC_638,(0,1):C.GC_975})

V_674 = Vertex(name = 'V_674',
               particles = [ P.s__tilde__, P.c, P.W__minus__, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV2, L.FFVV4 ],
               couplings = {(0,0):C.GC_641,(0,1):C.GC_977})

V_675 = Vertex(name = 'V_675',
               particles = [ P.e__plus__, P.ve, P.a, P.W__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVVS4 ],
               couplings = {(0,0):C.GC_773})

V_676 = Vertex(name = 'V_676',
               particles = [ P.e__plus__, P.ve, P.a, P.W__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVVS4 ],
               couplings = {(0,0):C.GC_774})

V_677 = Vertex(name = 'V_677',
               particles = [ P.mu__plus__, P.vm, P.a, P.W__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVVS4 ],
               couplings = {(0,0):C.GC_853})

V_678 = Vertex(name = 'V_678',
               particles = [ P.mu__plus__, P.vm, P.a, P.W__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVVS4 ],
               couplings = {(0,0):C.GC_854})

V_679 = Vertex(name = 'V_679',
               particles = [ P.ta__plus__, P.vt, P.a, P.W__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVVS4, L.FFVVS5 ],
               couplings = {(0,0):C.GC_1058,(0,1):C.GC_1059})

V_680 = Vertex(name = 'V_680',
               particles = [ P.e__plus__, P.ve, P.a, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVV4 ],
               couplings = {(0,0):C.GC_795})

V_681 = Vertex(name = 'V_681',
               particles = [ P.e__plus__, P.ve, P.a, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVV4 ],
               couplings = {(0,0):C.GC_796})

V_682 = Vertex(name = 'V_682',
               particles = [ P.mu__plus__, P.vm, P.a, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVV4 ],
               couplings = {(0,0):C.GC_875})

V_683 = Vertex(name = 'V_683',
               particles = [ P.mu__plus__, P.vm, P.a, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVV4 ],
               couplings = {(0,0):C.GC_876})

V_684 = Vertex(name = 'V_684',
               particles = [ P.ta__plus__, P.vt, P.a, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVV4 ],
               couplings = {(0,0):C.GC_1080})

V_685 = Vertex(name = 'V_685',
               particles = [ P.ta__plus__, P.vt, P.a, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVV4 ],
               couplings = {(0,0):C.GC_1081})

V_686 = Vertex(name = 'V_686',
               particles = [ P.e__plus__, P.ve, P.W__minus__, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVVS4 ],
               couplings = {(0,0):C.GC_779})

V_687 = Vertex(name = 'V_687',
               particles = [ P.e__plus__, P.ve, P.W__minus__, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVVS4 ],
               couplings = {(0,0):C.GC_781})

V_688 = Vertex(name = 'V_688',
               particles = [ P.mu__plus__, P.vm, P.W__minus__, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVVS4 ],
               couplings = {(0,0):C.GC_859})

V_689 = Vertex(name = 'V_689',
               particles = [ P.mu__plus__, P.vm, P.W__minus__, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVVS4 ],
               couplings = {(0,0):C.GC_861})

V_690 = Vertex(name = 'V_690',
               particles = [ P.ta__plus__, P.vt, P.W__minus__, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVVS4 ],
               couplings = {(0,0):C.GC_1064})

V_691 = Vertex(name = 'V_691',
               particles = [ P.ta__plus__, P.vt, P.W__minus__, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVVS4 ],
               couplings = {(0,0):C.GC_1066})

V_692 = Vertex(name = 'V_692',
               particles = [ P.e__plus__, P.ve, P.W__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFVV4 ],
               couplings = {(0,0):C.GC_800})

V_693 = Vertex(name = 'V_693',
               particles = [ P.e__plus__, P.ve, P.W__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFVV4 ],
               couplings = {(0,0):C.GC_802})

V_694 = Vertex(name = 'V_694',
               particles = [ P.mu__plus__, P.vm, P.W__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFVV4 ],
               couplings = {(0,0):C.GC_880})

V_695 = Vertex(name = 'V_695',
               particles = [ P.mu__plus__, P.vm, P.W__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFVV4 ],
               couplings = {(0,0):C.GC_882})

V_696 = Vertex(name = 'V_696',
               particles = [ P.ta__plus__, P.vt, P.W__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFVV4 ],
               couplings = {(0,0):C.GC_1085})

V_697 = Vertex(name = 'V_697',
               particles = [ P.ta__plus__, P.vt, P.W__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFVV4 ],
               couplings = {(0,0):C.GC_1087})

V_698 = Vertex(name = 'V_698',
               particles = [ P.b__tilde__, P.b, P.b__tilde__, P.b ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF13, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF17, L.FFFF18, L.FFFF3, L.FFFF4 ],
               couplings = {(1,6):C.GC_79,(3,6):C.GC_80,(0,7):C.GC_79,(2,7):C.GC_80,(3,1):C.GC_63,(1,2):C.GC_10,(0,4):C.GC_62,(2,4):C.GC_63,(0,5):C.GC_10,(0,0):C.GC_62,(2,0):C.GC_63,(1,3):C.GC_62,(3,3):C.GC_63,(1,1):C.GC_62})

V_699 = Vertex(name = 'V_699',
               particles = [ P.b__tilde__, P.b, P.d__tilde__, P.d ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF13, L.FFFF14, L.FFFF16, L.FFFF17, L.FFFF18, L.FFFF4 ],
               couplings = {(0,5):C.GC_67,(2,5):C.GC_68,(0,0):C.GC_18,(2,0):C.GC_19,(1,2):C.GC_660,(3,2):C.GC_663,(1,1):C.GC_661,(3,1):C.GC_664,(0,3):C.GC_64,(2,3):C.GC_65,(0,4):C.GC_11,(2,4):C.GC_12})

V_700 = Vertex(name = 'V_700',
               particles = [ P.b__tilde__, P.b, P.d__tilde__, P.d ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF14, L.FFFF16, L.FFFF4 ],
               couplings = {(0,2):C.GC_70,(2,2):C.GC_73,(1,1):C.GC_662,(3,1):C.GC_665,(1,0):C.GC_662,(3,0):C.GC_665})

V_701 = Vertex(name = 'V_701',
               particles = [ P.b__tilde__, P.b, P.s__tilde__, P.s ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF13, L.FFFF14, L.FFFF16, L.FFFF17, L.FFFF18, L.FFFF4 ],
               couplings = {(0,5):C.GC_67,(2,5):C.GC_68,(0,0):C.GC_18,(2,0):C.GC_19,(1,2):C.GC_916,(3,2):C.GC_919,(1,1):C.GC_917,(3,1):C.GC_920,(0,3):C.GC_64,(2,3):C.GC_65,(0,4):C.GC_11,(2,4):C.GC_12})

V_702 = Vertex(name = 'V_702',
               particles = [ P.b__tilde__, P.b, P.s__tilde__, P.s ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF14, L.FFFF16, L.FFFF4 ],
               couplings = {(0,2):C.GC_70,(2,2):C.GC_73,(1,1):C.GC_918,(3,1):C.GC_921,(1,0):C.GC_918,(3,0):C.GC_921})

V_703 = Vertex(name = 'V_703',
               particles = [ P.d__tilde__, P.d, P.d__tilde__, P.d ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF13, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF17, L.FFFF18, L.FFFF3, L.FFFF4 ],
               couplings = {(1,6):C.GC_43,(3,6):C.GC_44,(0,7):C.GC_43,(2,7):C.GC_44,(0,0):C.GC_40,(2,0):C.GC_41,(1,3):C.GC_40,(3,3):C.GC_41,(1,1):C.GC_40,(3,1):C.GC_41,(1,2):C.GC_26,(3,2):C.GC_27,(0,4):C.GC_40,(2,4):C.GC_41,(0,5):C.GC_26,(2,5):C.GC_27})

V_704 = Vertex(name = 'V_704',
               particles = [ P.d__tilde__, P.d, P.d__tilde__, P.d ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF3, L.FFFF4 ],
               couplings = {(1,0):C.GC_46,(3,0):C.GC_49,(0,1):C.GC_46,(2,1):C.GC_49})

V_705 = Vertex(name = 'V_705',
               particles = [ P.d__tilde__, P.d, P.s__tilde__, P.s ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'T(-1,2,1)*T(-1,4,3)' ],
               lorentz = [ L.FFFF13, L.FFFF17, L.FFFF18, L.FFFF4 ],
               couplings = {(0,3):C.GC_43,(1,3):C.GC_44,(0,0):C.GC_40,(1,0):C.GC_41,(0,1):C.GC_40,(1,1):C.GC_41,(0,2):C.GC_26,(1,2):C.GC_27})

V_706 = Vertex(name = 'V_706',
               particles = [ P.d__tilde__, P.d, P.s__tilde__, P.s ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'T(-1,2,1)*T(-1,4,3)' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_46,(1,0):C.GC_49})

V_707 = Vertex(name = 'V_707',
               particles = [ P.s__tilde__, P.s, P.s__tilde__, P.s ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF13, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF17, L.FFFF18, L.FFFF3, L.FFFF4 ],
               couplings = {(1,6):C.GC_43,(3,6):C.GC_44,(0,7):C.GC_43,(2,7):C.GC_44,(0,0):C.GC_40,(2,0):C.GC_41,(1,3):C.GC_40,(3,3):C.GC_41,(1,1):C.GC_40,(3,1):C.GC_41,(1,2):C.GC_26,(3,2):C.GC_27,(0,4):C.GC_40,(2,4):C.GC_41,(0,5):C.GC_26,(2,5):C.GC_27})

V_708 = Vertex(name = 'V_708',
               particles = [ P.s__tilde__, P.s, P.s__tilde__, P.s ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF3, L.FFFF4 ],
               couplings = {(1,0):C.GC_46,(3,0):C.GC_49,(0,1):C.GC_46,(2,1):C.GC_49})

V_709 = Vertex(name = 'V_709',
               particles = [ P.e__plus__, P.e__minus__, P.b__tilde__, P.b ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF13, L.FFFF17, L.FFFF18, L.FFFF4, L.FFFF9 ],
               couplings = {(0,4):C.GC_75,(0,1):C.GC_66,(0,2):C.GC_20,(0,3):C.GC_13,(0,5):C.GC_754,(0,0):C.GC_755})

V_710 = Vertex(name = 'V_710',
               particles = [ P.e__plus__, P.e__minus__, P.b__tilde__, P.b ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF4, L.FFFF9 ],
               couplings = {(0,1):C.GC_77,(0,2):C.GC_756,(0,0):C.GC_756})

V_711 = Vertex(name = 'V_711',
               particles = [ P.mu__plus__, P.mu__minus__, P.b__tilde__, P.b ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF13, L.FFFF17, L.FFFF18, L.FFFF4, L.FFFF9 ],
               couplings = {(0,4):C.GC_75,(0,1):C.GC_66,(0,2):C.GC_20,(0,3):C.GC_13,(0,5):C.GC_834,(0,0):C.GC_835})

V_712 = Vertex(name = 'V_712',
               particles = [ P.mu__plus__, P.mu__minus__, P.b__tilde__, P.b ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF4, L.FFFF9 ],
               couplings = {(0,1):C.GC_77,(0,2):C.GC_836,(0,0):C.GC_836})

V_713 = Vertex(name = 'V_713',
               particles = [ P.ta__plus__, P.ta__minus__, P.b__tilde__, P.b ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF13, L.FFFF17, L.FFFF18, L.FFFF4, L.FFFF9 ],
               couplings = {(0,4):C.GC_75,(0,1):C.GC_66,(0,2):C.GC_20,(0,3):C.GC_13,(0,5):C.GC_1039,(0,0):C.GC_1040})

V_714 = Vertex(name = 'V_714',
               particles = [ P.ta__plus__, P.ta__minus__, P.b__tilde__, P.b ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF4, L.FFFF9 ],
               couplings = {(0,1):C.GC_77,(0,2):C.GC_1041,(0,0):C.GC_1041})

V_715 = Vertex(name = 'V_715',
               particles = [ P.e__plus__, P.e__minus__, P.d__tilde__, P.d ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF13, L.FFFF17, L.FFFF18, L.FFFF4, L.FFFF9 ],
               couplings = {(0,4):C.GC_55,(0,1):C.GC_42,(0,2):C.GC_53,(0,3):C.GC_28,(0,5):C.GC_825,(0,0):C.GC_826})

V_716 = Vertex(name = 'V_716',
               particles = [ P.e__plus__, P.e__minus__, P.d__tilde__, P.d ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF4, L.FFFF9 ],
               couplings = {(0,1):C.GC_57,(0,2):C.GC_827,(0,0):C.GC_827})

V_717 = Vertex(name = 'V_717',
               particles = [ P.mu__plus__, P.mu__minus__, P.d__tilde__, P.d ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF13, L.FFFF17, L.FFFF18, L.FFFF4, L.FFFF9 ],
               couplings = {(0,4):C.GC_55,(0,1):C.GC_42,(0,2):C.GC_53,(0,3):C.GC_28,(0,5):C.GC_905,(0,0):C.GC_906})

V_718 = Vertex(name = 'V_718',
               particles = [ P.mu__plus__, P.mu__minus__, P.d__tilde__, P.d ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF4, L.FFFF9 ],
               couplings = {(0,1):C.GC_57,(0,2):C.GC_907,(0,0):C.GC_907})

V_719 = Vertex(name = 'V_719',
               particles = [ P.ta__plus__, P.ta__minus__, P.d__tilde__, P.d ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF13, L.FFFF17, L.FFFF18, L.FFFF4, L.FFFF9 ],
               couplings = {(0,4):C.GC_55,(0,1):C.GC_42,(0,2):C.GC_53,(0,3):C.GC_28,(0,5):C.GC_1110,(0,0):C.GC_1111})

V_720 = Vertex(name = 'V_720',
               particles = [ P.ta__plus__, P.ta__minus__, P.d__tilde__, P.d ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF4, L.FFFF9 ],
               couplings = {(0,1):C.GC_57,(0,2):C.GC_1112,(0,0):C.GC_1112})

V_721 = Vertex(name = 'V_721',
               particles = [ P.e__plus__, P.e__minus__, P.s__tilde__, P.s ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF13, L.FFFF17, L.FFFF18, L.FFFF4, L.FFFF9 ],
               couplings = {(0,4):C.GC_55,(0,1):C.GC_42,(0,2):C.GC_53,(0,3):C.GC_28,(0,5):C.GC_1019,(0,0):C.GC_1020})

V_722 = Vertex(name = 'V_722',
               particles = [ P.e__plus__, P.e__minus__, P.s__tilde__, P.s ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF4, L.FFFF9 ],
               couplings = {(0,1):C.GC_57,(0,2):C.GC_1021,(0,0):C.GC_1021})

V_723 = Vertex(name = 'V_723',
               particles = [ P.mu__plus__, P.mu__minus__, P.s__tilde__, P.s ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF13, L.FFFF17, L.FFFF18, L.FFFF4, L.FFFF9 ],
               couplings = {(0,4):C.GC_55,(0,1):C.GC_42,(0,2):C.GC_53,(0,3):C.GC_28,(0,5):C.GC_1022,(0,0):C.GC_1023})

V_724 = Vertex(name = 'V_724',
               particles = [ P.mu__plus__, P.mu__minus__, P.s__tilde__, P.s ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF4, L.FFFF9 ],
               couplings = {(0,1):C.GC_57,(0,2):C.GC_1024,(0,0):C.GC_1024})

V_725 = Vertex(name = 'V_725',
               particles = [ P.ta__plus__, P.ta__minus__, P.s__tilde__, P.s ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF13, L.FFFF17, L.FFFF18, L.FFFF4, L.FFFF9 ],
               couplings = {(0,4):C.GC_55,(0,1):C.GC_42,(0,2):C.GC_53,(0,3):C.GC_28,(0,5):C.GC_1113,(0,0):C.GC_1114})

V_726 = Vertex(name = 'V_726',
               particles = [ P.ta__plus__, P.ta__minus__, P.s__tilde__, P.s ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF4, L.FFFF9 ],
               couplings = {(0,1):C.GC_57,(0,2):C.GC_1115,(0,0):C.GC_1115})

V_727 = Vertex(name = 'V_727',
               particles = [ P.e__plus__, P.e__minus__, P.e__plus__, P.e__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFFF13, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF17, L.FFFF18, L.FFFF3, L.FFFF4 ],
               couplings = {(0,6):C.GC_59,(0,7):C.GC_59,(0,0):C.GC_54,(0,3):C.GC_54,(0,1):C.GC_54,(0,2):C.GC_29,(0,4):C.GC_54,(0,5):C.GC_29})

V_728 = Vertex(name = 'V_728',
               particles = [ P.e__plus__, P.e__minus__, P.e__plus__, P.e__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFFF3, L.FFFF4 ],
               couplings = {(0,0):C.GC_60,(0,1):C.GC_60})

V_729 = Vertex(name = 'V_729',
               particles = [ P.mu__plus__, P.e__minus__, P.e__plus__, P.mu__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF3, L.FFFF4 ],
               couplings = {(0,3):C.GC_59,(0,4):C.GC_60,(0,2):C.GC_54,(0,0):C.GC_54,(0,1):C.GC_29})

V_730 = Vertex(name = 'V_730',
               particles = [ P.ta__plus__, P.e__minus__, P.e__plus__, P.ta__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF3, L.FFFF4 ],
               couplings = {(0,3):C.GC_59,(0,4):C.GC_60,(0,2):C.GC_54,(0,0):C.GC_54,(0,1):C.GC_29})

V_731 = Vertex(name = 'V_731',
               particles = [ P.mu__plus__, P.mu__minus__, P.mu__plus__, P.mu__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFFF13, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF17, L.FFFF18, L.FFFF3, L.FFFF4 ],
               couplings = {(0,6):C.GC_59,(0,7):C.GC_59,(0,0):C.GC_54,(0,3):C.GC_54,(0,1):C.GC_54,(0,2):C.GC_29,(0,4):C.GC_54,(0,5):C.GC_29})

V_732 = Vertex(name = 'V_732',
               particles = [ P.mu__plus__, P.mu__minus__, P.mu__plus__, P.mu__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFFF3, L.FFFF4 ],
               couplings = {(0,0):C.GC_60,(0,1):C.GC_60})

V_733 = Vertex(name = 'V_733',
               particles = [ P.ta__plus__, P.mu__minus__, P.mu__plus__, P.ta__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF3, L.FFFF4 ],
               couplings = {(0,3):C.GC_59,(0,4):C.GC_60,(0,2):C.GC_54,(0,0):C.GC_54,(0,1):C.GC_29})

V_734 = Vertex(name = 'V_734',
               particles = [ P.ta__plus__, P.ta__minus__, P.ta__plus__, P.ta__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFFF13, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF17, L.FFFF18, L.FFFF3, L.FFFF4 ],
               couplings = {(0,6):C.GC_59,(0,7):C.GC_59,(0,0):C.GC_54,(0,3):C.GC_54,(0,1):C.GC_54,(0,2):C.GC_29,(0,4):C.GC_54,(0,5):C.GC_29})

V_735 = Vertex(name = 'V_735',
               particles = [ P.ta__plus__, P.ta__minus__, P.ta__plus__, P.ta__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFFF3, L.FFFF4 ],
               couplings = {(0,0):C.GC_60,(0,1):C.GC_60})

V_736 = Vertex(name = 'V_736',
               particles = [ P.ve__tilde__, P.e__minus__, P.b__tilde__, P.t ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF2, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,2):C.GC_78,(0,5):C.GC_761,(0,3):C.GC_763,(0,4):C.GC_763,(0,1):C.GC_758,(0,0):C.GC_755})

V_737 = Vertex(name = 'V_737',
               particles = [ P.ve__tilde__, P.e__minus__, P.b__tilde__, P.t ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF2, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,4):C.GC_767,(0,2):C.GC_766,(0,3):C.GC_766,(0,1):C.GC_760,(0,0):C.GC_756})

V_738 = Vertex(name = 'V_738',
               particles = [ P.vm__tilde__, P.mu__minus__, P.b__tilde__, P.t ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF2, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,2):C.GC_78,(0,5):C.GC_841,(0,3):C.GC_843,(0,4):C.GC_843,(0,1):C.GC_838,(0,0):C.GC_835})

V_739 = Vertex(name = 'V_739',
               particles = [ P.vm__tilde__, P.mu__minus__, P.b__tilde__, P.t ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF2, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,4):C.GC_847,(0,2):C.GC_846,(0,3):C.GC_846,(0,1):C.GC_840,(0,0):C.GC_836})

V_740 = Vertex(name = 'V_740',
               particles = [ P.vt__tilde__, P.ta__minus__, P.b__tilde__, P.t ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF2, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,2):C.GC_78,(0,5):C.GC_1046,(0,3):C.GC_1048,(0,4):C.GC_1048,(0,1):C.GC_1043,(0,0):C.GC_1040})

V_741 = Vertex(name = 'V_741',
               particles = [ P.vt__tilde__, P.ta__minus__, P.b__tilde__, P.t ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF2, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,4):C.GC_1052,(0,2):C.GC_1051,(0,3):C.GC_1051,(0,1):C.GC_1045,(0,0):C.GC_1041})

V_742 = Vertex(name = 'V_742',
               particles = [ P.e__plus__, P.e__minus__, P.t__tilde__, P.t ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF12, L.FFFF13, L.FFFF17, L.FFFF18, L.FFFF19, L.FFFF2, L.FFFF20, L.FFFF21, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,8):C.GC_75,(0,11):C.GC_764,(0,9):C.GC_762,(0,10):C.GC_762,(0,5):C.GC_757,(0,1):C.GC_66,(0,2):C.GC_114,(0,3):C.GC_97,(0,7):C.GC_761,(0,4):C.GC_763,(0,6):C.GC_763,(0,0):C.GC_758})

V_743 = Vertex(name = 'V_743',
               particles = [ P.e__plus__, P.e__minus__, P.t__tilde__, P.t ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF12, L.FFFF19, L.FFFF2, L.FFFF20, L.FFFF21, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,5):C.GC_76,(0,8):C.GC_768,(0,6):C.GC_765,(0,7):C.GC_765,(0,2):C.GC_759,(0,4):C.GC_768,(0,1):C.GC_765,(0,3):C.GC_765,(0,0):C.GC_759})

V_744 = Vertex(name = 'V_744',
               particles = [ P.mu__plus__, P.mu__minus__, P.t__tilde__, P.t ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF12, L.FFFF13, L.FFFF17, L.FFFF18, L.FFFF19, L.FFFF2, L.FFFF20, L.FFFF21, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,8):C.GC_75,(0,11):C.GC_844,(0,9):C.GC_842,(0,10):C.GC_842,(0,5):C.GC_837,(0,1):C.GC_66,(0,2):C.GC_114,(0,3):C.GC_97,(0,7):C.GC_841,(0,4):C.GC_843,(0,6):C.GC_843,(0,0):C.GC_838})

V_745 = Vertex(name = 'V_745',
               particles = [ P.mu__plus__, P.mu__minus__, P.t__tilde__, P.t ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF12, L.FFFF19, L.FFFF2, L.FFFF20, L.FFFF21, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,5):C.GC_76,(0,8):C.GC_848,(0,6):C.GC_845,(0,7):C.GC_845,(0,2):C.GC_839,(0,4):C.GC_848,(0,1):C.GC_845,(0,3):C.GC_845,(0,0):C.GC_839})

V_746 = Vertex(name = 'V_746',
               particles = [ P.ta__plus__, P.ta__minus__, P.t__tilde__, P.t ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF12, L.FFFF13, L.FFFF17, L.FFFF18, L.FFFF19, L.FFFF2, L.FFFF20, L.FFFF21, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,8):C.GC_75,(0,11):C.GC_1049,(0,9):C.GC_1047,(0,10):C.GC_1047,(0,5):C.GC_1042,(0,1):C.GC_66,(0,2):C.GC_114,(0,3):C.GC_97,(0,7):C.GC_1046,(0,4):C.GC_1048,(0,6):C.GC_1048,(0,0):C.GC_1043})

V_747 = Vertex(name = 'V_747',
               particles = [ P.ta__plus__, P.ta__minus__, P.t__tilde__, P.t ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF12, L.FFFF19, L.FFFF2, L.FFFF20, L.FFFF21, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,5):C.GC_76,(0,8):C.GC_1053,(0,6):C.GC_1050,(0,7):C.GC_1050,(0,2):C.GC_1044,(0,4):C.GC_1053,(0,1):C.GC_1050,(0,3):C.GC_1050,(0,0):C.GC_1044})

V_748 = Vertex(name = 'V_748',
               particles = [ P.u__tilde__, P.d, P.b__tilde__, P.t ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF2, L.FFFF4 ],
               couplings = {(0,7):C.GC_71,(2,7):C.GC_74,(1,0):C.GC_667,(3,0):C.GC_671,(0,6):C.GC_674,(2,6):C.GC_678,(1,5):C.GC_661,(3,5):C.GC_664,(1,3):C.GC_1117,(3,3):C.GC_1120,(1,4):C.GC_1213,(3,4):C.GC_1216,(1,1):C.GC_1131,(3,1):C.GC_1135,(0,2):C.GC_1124,(2,2):C.GC_1128})

V_749 = Vertex(name = 'V_749',
               particles = [ P.u__tilde__, P.d, P.b__tilde__, P.t ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF2 ],
               couplings = {(1,0):C.GC_669,(3,0):C.GC_673,(0,6):C.GC_676,(2,6):C.GC_680,(1,5):C.GC_662,(3,5):C.GC_665,(1,3):C.GC_1119,(3,3):C.GC_1122,(1,4):C.GC_1214,(3,4):C.GC_1217,(1,1):C.GC_1134,(3,1):C.GC_1138,(0,2):C.GC_1125,(2,2):C.GC_1129})

V_750 = Vertex(name = 'V_750',
               particles = [ P.c__tilde__, P.s, P.b__tilde__, P.t ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF2, L.FFFF4 ],
               couplings = {(0,7):C.GC_71,(2,7):C.GC_74,(1,0):C.GC_923,(3,0):C.GC_927,(0,6):C.GC_930,(2,6):C.GC_934,(1,5):C.GC_917,(3,5):C.GC_920,(1,3):C.GC_573,(3,3):C.GC_576,(1,4):C.GC_1005,(3,4):C.GC_1008,(1,1):C.GC_587,(3,1):C.GC_591,(0,2):C.GC_580,(2,2):C.GC_584})

V_751 = Vertex(name = 'V_751',
               particles = [ P.c__tilde__, P.s, P.b__tilde__, P.t ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF2 ],
               couplings = {(1,0):C.GC_925,(3,0):C.GC_929,(0,6):C.GC_932,(2,6):C.GC_936,(1,5):C.GC_918,(3,5):C.GC_921,(1,3):C.GC_575,(3,3):C.GC_578,(1,4):C.GC_1006,(3,4):C.GC_1009,(1,1):C.GC_590,(3,1):C.GC_594,(0,2):C.GC_581,(2,2):C.GC_585})

V_752 = Vertex(name = 'V_752',
               particles = [ P.t__tilde__, P.b, P.b__tilde__, P.t ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF2, L.FFFF3 ],
               couplings = {(1,7):C.GC_79,(3,7):C.GC_80,(1,0):C.GC_84,(3,0):C.GC_88,(0,6):C.GC_83,(2,6):C.GC_87,(1,5):C.GC_62,(3,5):C.GC_63,(1,3):C.GC_81,(3,3):C.GC_82,(1,4):C.GC_93,(3,4):C.GC_94,(1,1):C.GC_83,(3,1):C.GC_87,(0,2):C.GC_84,(2,2):C.GC_88})

V_753 = Vertex(name = 'V_753',
               particles = [ P.t__tilde__, P.b, P.b__tilde__, P.t ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF2 ],
               couplings = {(1,0):C.GC_86,(3,0):C.GC_90,(0,3):C.GC_85,(2,3):C.GC_89,(1,1):C.GC_86,(3,1):C.GC_90,(0,2):C.GC_85,(2,2):C.GC_89})

V_754 = Vertex(name = 'V_754',
               particles = [ P.t__tilde__, P.d, P.d__tilde__, P.t ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF2, L.FFFF3 ],
               couplings = {(1,7):C.GC_67,(3,7):C.GC_68,(1,0):C.GC_675,(3,0):C.GC_679,(0,6):C.GC_666,(2,6):C.GC_670,(1,5):C.GC_64,(3,5):C.GC_65,(1,3):C.GC_112,(3,3):C.GC_113,(1,4):C.GC_95,(3,4):C.GC_96,(1,1):C.GC_674,(3,1):C.GC_678,(0,2):C.GC_667,(2,2):C.GC_671})

V_755 = Vertex(name = 'V_755',
               particles = [ P.t__tilde__, P.d, P.d__tilde__, P.t ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF2, L.FFFF3 ],
               couplings = {(1,4):C.GC_69,(3,4):C.GC_72,(1,0):C.GC_677,(3,0):C.GC_681,(0,3):C.GC_668,(2,3):C.GC_672,(1,1):C.GC_677,(3,1):C.GC_681,(0,2):C.GC_668,(2,2):C.GC_672})

V_756 = Vertex(name = 'V_756',
               particles = [ P.t__tilde__, P.s, P.s__tilde__, P.t ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF2, L.FFFF3 ],
               couplings = {(1,7):C.GC_67,(3,7):C.GC_68,(1,0):C.GC_931,(3,0):C.GC_935,(0,6):C.GC_922,(2,6):C.GC_926,(1,5):C.GC_64,(3,5):C.GC_65,(1,3):C.GC_112,(3,3):C.GC_113,(1,4):C.GC_95,(3,4):C.GC_96,(1,1):C.GC_930,(3,1):C.GC_934,(0,2):C.GC_923,(2,2):C.GC_927})

V_757 = Vertex(name = 'V_757',
               particles = [ P.t__tilde__, P.s, P.s__tilde__, P.t ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF2, L.FFFF3 ],
               couplings = {(1,4):C.GC_69,(3,4):C.GC_72,(1,0):C.GC_933,(3,0):C.GC_937,(0,3):C.GC_924,(2,3):C.GC_928,(1,1):C.GC_933,(3,1):C.GC_937,(0,2):C.GC_924,(2,2):C.GC_928})

V_758 = Vertex(name = 'V_758',
               particles = [ P.t__tilde__, P.t, P.t__tilde__, P.t ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF13, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF17, L.FFFF18, L.FFFF3, L.FFFF4 ],
               couplings = {(1,6):C.GC_79,(3,6):C.GC_80,(0,7):C.GC_79,(2,7):C.GC_80,(0,0):C.GC_81,(2,0):C.GC_82,(1,3):C.GC_81,(3,3):C.GC_82,(1,1):C.GC_81,(3,1):C.GC_82,(1,2):C.GC_115,(0,4):C.GC_81,(2,4):C.GC_82,(0,5):C.GC_115})

V_759 = Vertex(name = 'V_759',
               particles = [ P.ve__tilde__, P.e__minus__, P.d__tilde__, P.u ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF2, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,2):C.GC_58,(0,5):C.GC_1231,(0,3):C.GC_1233,(0,4):C.GC_1233,(0,1):C.GC_1228,(0,0):C.GC_826})

V_760 = Vertex(name = 'V_760',
               particles = [ P.ve__tilde__, P.e__minus__, P.d__tilde__, P.u ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF2, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,4):C.GC_1237,(0,2):C.GC_1236,(0,3):C.GC_1236,(0,1):C.GC_1230,(0,0):C.GC_827})

V_761 = Vertex(name = 'V_761',
               particles = [ P.ve__tilde__, P.e__minus__, P.s__tilde__, P.c ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF2, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,2):C.GC_58,(0,5):C.GC_817,(0,3):C.GC_819,(0,4):C.GC_819,(0,1):C.GC_814,(0,0):C.GC_1020})

V_762 = Vertex(name = 'V_762',
               particles = [ P.ve__tilde__, P.e__minus__, P.s__tilde__, P.c ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF2, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,4):C.GC_823,(0,2):C.GC_822,(0,3):C.GC_822,(0,1):C.GC_816,(0,0):C.GC_1021})

V_763 = Vertex(name = 'V_763',
               particles = [ P.vm__tilde__, P.mu__minus__, P.d__tilde__, P.u ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF2, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,2):C.GC_58,(0,5):C.GC_1243,(0,3):C.GC_1245,(0,4):C.GC_1245,(0,1):C.GC_1240,(0,0):C.GC_906})

V_764 = Vertex(name = 'V_764',
               particles = [ P.vm__tilde__, P.mu__minus__, P.d__tilde__, P.u ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF2, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,4):C.GC_1249,(0,2):C.GC_1248,(0,3):C.GC_1248,(0,1):C.GC_1242,(0,0):C.GC_907})

V_765 = Vertex(name = 'V_765',
               particles = [ P.vm__tilde__, P.mu__minus__, P.s__tilde__, P.c ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF2, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,2):C.GC_58,(0,5):C.GC_897,(0,3):C.GC_899,(0,4):C.GC_899,(0,1):C.GC_894,(0,0):C.GC_1023})

V_766 = Vertex(name = 'V_766',
               particles = [ P.vm__tilde__, P.mu__minus__, P.s__tilde__, P.c ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF2, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,4):C.GC_903,(0,2):C.GC_902,(0,3):C.GC_902,(0,1):C.GC_896,(0,0):C.GC_1024})

V_767 = Vertex(name = 'V_767',
               particles = [ P.vt__tilde__, P.ta__minus__, P.d__tilde__, P.u ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF2, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,2):C.GC_58,(0,5):C.GC_1271,(0,3):C.GC_1273,(0,4):C.GC_1273,(0,1):C.GC_1268,(0,0):C.GC_1111})

V_768 = Vertex(name = 'V_768',
               particles = [ P.vt__tilde__, P.ta__minus__, P.d__tilde__, P.u ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF2, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,4):C.GC_1277,(0,2):C.GC_1276,(0,3):C.GC_1276,(0,1):C.GC_1270,(0,0):C.GC_1112})

V_769 = Vertex(name = 'V_769',
               particles = [ P.vt__tilde__, P.ta__minus__, P.s__tilde__, P.c ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF2, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF7, L.FFFF8 ],
               couplings = {(0,4):C.GC_1104,(0,2):C.GC_58,(0,6):C.GC_1102,(0,3):C.GC_1104,(0,5):C.GC_1107,(0,1):C.GC_1099,(0,0):C.GC_1114})

V_770 = Vertex(name = 'V_770',
               particles = [ P.vt__tilde__, P.ta__minus__, P.s__tilde__, P.c ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF11, L.FFFF2, L.FFFF6, L.FFFF8 ],
               couplings = {(0,2):C.GC_1107,(0,3):C.GC_1108,(0,1):C.GC_1101,(0,0):C.GC_1115})

V_771 = Vertex(name = 'V_771',
               particles = [ P.e__plus__, P.e__minus__, P.u__tilde__, P.u ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF12, L.FFFF13, L.FFFF17, L.FFFF18, L.FFFF19, L.FFFF2, L.FFFF20, L.FFFF21, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,8):C.GC_55,(0,11):C.GC_1234,(0,9):C.GC_1232,(0,10):C.GC_1232,(0,5):C.GC_1227,(0,1):C.GC_42,(0,2):C.GC_61,(0,3):C.GC_30,(0,7):C.GC_1231,(0,4):C.GC_1233,(0,6):C.GC_1233,(0,0):C.GC_1228})

V_772 = Vertex(name = 'V_772',
               particles = [ P.e__plus__, P.e__minus__, P.u__tilde__, P.u ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF12, L.FFFF19, L.FFFF2, L.FFFF20, L.FFFF21, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,5):C.GC_56,(0,8):C.GC_1238,(0,6):C.GC_1235,(0,7):C.GC_1235,(0,2):C.GC_1229,(0,4):C.GC_1238,(0,1):C.GC_1235,(0,3):C.GC_1235,(0,0):C.GC_1229})

V_773 = Vertex(name = 'V_773',
               particles = [ P.e__plus__, P.e__minus__, P.c__tilde__, P.c ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF12, L.FFFF13, L.FFFF17, L.FFFF18, L.FFFF19, L.FFFF2, L.FFFF20, L.FFFF21, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,8):C.GC_55,(0,11):C.GC_820,(0,9):C.GC_818,(0,10):C.GC_818,(0,5):C.GC_813,(0,1):C.GC_42,(0,2):C.GC_61,(0,3):C.GC_30,(0,7):C.GC_817,(0,4):C.GC_819,(0,6):C.GC_819,(0,0):C.GC_814})

V_774 = Vertex(name = 'V_774',
               particles = [ P.e__plus__, P.e__minus__, P.c__tilde__, P.c ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF12, L.FFFF19, L.FFFF2, L.FFFF20, L.FFFF21, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,5):C.GC_56,(0,8):C.GC_824,(0,6):C.GC_821,(0,7):C.GC_821,(0,2):C.GC_815,(0,4):C.GC_824,(0,1):C.GC_821,(0,3):C.GC_821,(0,0):C.GC_815})

V_775 = Vertex(name = 'V_775',
               particles = [ P.mu__plus__, P.mu__minus__, P.u__tilde__, P.u ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF12, L.FFFF13, L.FFFF17, L.FFFF18, L.FFFF19, L.FFFF2, L.FFFF20, L.FFFF21, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,8):C.GC_55,(0,11):C.GC_1246,(0,9):C.GC_1244,(0,10):C.GC_1244,(0,5):C.GC_1239,(0,1):C.GC_42,(0,2):C.GC_61,(0,3):C.GC_30,(0,7):C.GC_1243,(0,4):C.GC_1245,(0,6):C.GC_1245,(0,0):C.GC_1240})

V_776 = Vertex(name = 'V_776',
               particles = [ P.mu__plus__, P.mu__minus__, P.u__tilde__, P.u ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF12, L.FFFF19, L.FFFF2, L.FFFF20, L.FFFF21, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,5):C.GC_56,(0,8):C.GC_1250,(0,6):C.GC_1247,(0,7):C.GC_1247,(0,2):C.GC_1241,(0,4):C.GC_1250,(0,1):C.GC_1247,(0,3):C.GC_1247,(0,0):C.GC_1241})

V_777 = Vertex(name = 'V_777',
               particles = [ P.mu__plus__, P.mu__minus__, P.c__tilde__, P.c ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF12, L.FFFF13, L.FFFF17, L.FFFF18, L.FFFF19, L.FFFF2, L.FFFF20, L.FFFF21, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,8):C.GC_55,(0,11):C.GC_900,(0,9):C.GC_898,(0,10):C.GC_898,(0,5):C.GC_893,(0,1):C.GC_42,(0,2):C.GC_61,(0,3):C.GC_30,(0,7):C.GC_897,(0,4):C.GC_899,(0,6):C.GC_899,(0,0):C.GC_894})

V_778 = Vertex(name = 'V_778',
               particles = [ P.mu__plus__, P.mu__minus__, P.c__tilde__, P.c ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF12, L.FFFF19, L.FFFF2, L.FFFF20, L.FFFF21, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,5):C.GC_56,(0,8):C.GC_904,(0,6):C.GC_901,(0,7):C.GC_901,(0,2):C.GC_895,(0,4):C.GC_904,(0,1):C.GC_901,(0,3):C.GC_901,(0,0):C.GC_895})

V_779 = Vertex(name = 'V_779',
               particles = [ P.ta__plus__, P.ta__minus__, P.u__tilde__, P.u ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF12, L.FFFF13, L.FFFF17, L.FFFF18, L.FFFF19, L.FFFF2, L.FFFF20, L.FFFF21, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,8):C.GC_55,(0,11):C.GC_1274,(0,9):C.GC_1272,(0,10):C.GC_1272,(0,5):C.GC_1267,(0,1):C.GC_42,(0,2):C.GC_61,(0,3):C.GC_30,(0,7):C.GC_1271,(0,4):C.GC_1273,(0,6):C.GC_1273,(0,0):C.GC_1268})

V_780 = Vertex(name = 'V_780',
               particles = [ P.ta__plus__, P.ta__minus__, P.u__tilde__, P.u ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF12, L.FFFF19, L.FFFF2, L.FFFF20, L.FFFF21, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,5):C.GC_56,(0,8):C.GC_1278,(0,6):C.GC_1275,(0,7):C.GC_1275,(0,2):C.GC_1269,(0,4):C.GC_1278,(0,1):C.GC_1275,(0,3):C.GC_1275,(0,0):C.GC_1269})

V_781 = Vertex(name = 'V_781',
               particles = [ P.ta__plus__, P.ta__minus__, P.c__tilde__, P.c ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF12, L.FFFF13, L.FFFF17, L.FFFF18, L.FFFF19, L.FFFF2, L.FFFF20, L.FFFF21, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,8):C.GC_55,(0,11):C.GC_1105,(0,9):C.GC_1103,(0,10):C.GC_1103,(0,5):C.GC_1098,(0,1):C.GC_42,(0,2):C.GC_61,(0,3):C.GC_30,(0,7):C.GC_1102,(0,4):C.GC_1104,(0,6):C.GC_1104,(0,0):C.GC_1099})

V_782 = Vertex(name = 'V_782',
               particles = [ P.ta__plus__, P.ta__minus__, P.c__tilde__, P.c ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF12, L.FFFF19, L.FFFF2, L.FFFF20, L.FFFF21, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,5):C.GC_56,(0,8):C.GC_1109,(0,6):C.GC_1106,(0,7):C.GC_1106,(0,2):C.GC_1100,(0,4):C.GC_1109,(0,1):C.GC_1106,(0,3):C.GC_1106,(0,0):C.GC_1100})

V_783 = Vertex(name = 'V_783',
               particles = [ P.t__tilde__, P.b, P.d__tilde__, P.u ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF2, L.FFFF4 ],
               couplings = {(0,7):C.GC_71,(2,7):C.GC_74,(1,0):C.GC_1132,(3,0):C.GC_1136,(0,6):C.GC_1123,(2,6):C.GC_1127,(1,5):C.GC_660,(3,5):C.GC_663,(1,3):C.GC_1118,(3,3):C.GC_1121,(1,4):C.GC_1212,(3,4):C.GC_1215,(1,1):C.GC_666,(3,1):C.GC_670,(0,2):C.GC_675,(2,2):C.GC_679})

V_784 = Vertex(name = 'V_784',
               particles = [ P.t__tilde__, P.b, P.d__tilde__, P.u ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF2 ],
               couplings = {(1,0):C.GC_1134,(3,0):C.GC_1138,(0,6):C.GC_1125,(2,6):C.GC_1129,(1,5):C.GC_662,(3,5):C.GC_665,(1,3):C.GC_1119,(3,3):C.GC_1122,(1,4):C.GC_1214,(3,4):C.GC_1217,(1,1):C.GC_669,(3,1):C.GC_673,(0,2):C.GC_676,(2,2):C.GC_680})

V_785 = Vertex(name = 'V_785',
               particles = [ P.t__tilde__, P.b, P.s__tilde__, P.c ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF2, L.FFFF4 ],
               couplings = {(0,7):C.GC_71,(2,7):C.GC_74,(1,0):C.GC_588,(3,0):C.GC_592,(0,6):C.GC_579,(2,6):C.GC_583,(1,5):C.GC_916,(3,5):C.GC_919,(1,3):C.GC_574,(3,3):C.GC_577,(1,4):C.GC_1004,(3,4):C.GC_1007,(1,1):C.GC_922,(3,1):C.GC_926,(0,2):C.GC_931,(2,2):C.GC_935})

V_786 = Vertex(name = 'V_786',
               particles = [ P.t__tilde__, P.b, P.s__tilde__, P.c ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF2 ],
               couplings = {(1,0):C.GC_590,(3,0):C.GC_594,(0,6):C.GC_581,(2,6):C.GC_585,(1,5):C.GC_918,(3,5):C.GC_921,(1,3):C.GC_575,(3,3):C.GC_578,(1,4):C.GC_1006,(3,4):C.GC_1009,(1,1):C.GC_925,(3,1):C.GC_929,(0,2):C.GC_932,(2,2):C.GC_936})

V_787 = Vertex(name = 'V_787',
               particles = [ P.u__tilde__, P.b, P.b__tilde__, P.u ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF2, L.FFFF3 ],
               couplings = {(1,7):C.GC_67,(3,7):C.GC_68,(1,0):C.GC_1124,(3,0):C.GC_1128,(0,6):C.GC_1131,(2,6):C.GC_1135,(1,5):C.GC_18,(3,5):C.GC_19,(1,3):C.GC_91,(3,3):C.GC_92,(1,4):C.GC_21,(3,4):C.GC_22,(1,1):C.GC_1123,(3,1):C.GC_1127,(0,2):C.GC_1132,(2,2):C.GC_1136})

V_788 = Vertex(name = 'V_788',
               particles = [ P.u__tilde__, P.b, P.b__tilde__, P.u ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF2, L.FFFF3 ],
               couplings = {(1,4):C.GC_69,(3,4):C.GC_72,(1,0):C.GC_1126,(3,0):C.GC_1130,(0,3):C.GC_1133,(2,3):C.GC_1137,(1,1):C.GC_1126,(3,1):C.GC_1130,(0,2):C.GC_1133,(2,2):C.GC_1137})

V_789 = Vertex(name = 'V_789',
               particles = [ P.c__tilde__, P.b, P.b__tilde__, P.c ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF2, L.FFFF3 ],
               couplings = {(1,7):C.GC_67,(3,7):C.GC_68,(1,0):C.GC_580,(3,0):C.GC_584,(0,6):C.GC_587,(2,6):C.GC_591,(1,5):C.GC_18,(3,5):C.GC_19,(1,3):C.GC_91,(3,3):C.GC_92,(1,4):C.GC_21,(3,4):C.GC_22,(1,1):C.GC_579,(3,1):C.GC_583,(0,2):C.GC_588,(2,2):C.GC_592})

V_790 = Vertex(name = 'V_790',
               particles = [ P.c__tilde__, P.b, P.b__tilde__, P.c ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF2, L.FFFF3 ],
               couplings = {(1,4):C.GC_69,(3,4):C.GC_72,(1,0):C.GC_582,(3,0):C.GC_586,(0,3):C.GC_589,(2,3):C.GC_593,(1,1):C.GC_582,(3,1):C.GC_586,(0,2):C.GC_589,(2,2):C.GC_593})

V_791 = Vertex(name = 'V_791',
               particles = [ P.t__tilde__, P.t, P.u__tilde__, P.u ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF13, L.FFFF14, L.FFFF16, L.FFFF17, L.FFFF18, L.FFFF4 ],
               couplings = {(0,5):C.GC_67,(2,5):C.GC_68,(0,0):C.GC_112,(2,0):C.GC_113,(1,2):C.GC_1117,(3,2):C.GC_1120,(1,1):C.GC_1118,(3,1):C.GC_1121,(0,3):C.GC_91,(2,3):C.GC_92,(0,4):C.GC_116,(2,4):C.GC_117})

V_792 = Vertex(name = 'V_792',
               particles = [ P.t__tilde__, P.t, P.u__tilde__, P.u ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF14, L.FFFF16, L.FFFF4 ],
               couplings = {(0,2):C.GC_70,(2,2):C.GC_73,(1,1):C.GC_1119,(3,1):C.GC_1122,(1,0):C.GC_1119,(3,0):C.GC_1122})

V_793 = Vertex(name = 'V_793',
               particles = [ P.t__tilde__, P.t, P.c__tilde__, P.c ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF13, L.FFFF14, L.FFFF16, L.FFFF17, L.FFFF18, L.FFFF4 ],
               couplings = {(0,5):C.GC_67,(2,5):C.GC_68,(0,0):C.GC_112,(2,0):C.GC_113,(1,2):C.GC_573,(3,2):C.GC_576,(1,1):C.GC_574,(3,1):C.GC_577,(0,3):C.GC_91,(2,3):C.GC_92,(0,4):C.GC_116,(2,4):C.GC_117})

V_794 = Vertex(name = 'V_794',
               particles = [ P.t__tilde__, P.t, P.c__tilde__, P.c ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF14, L.FFFF16, L.FFFF4 ],
               couplings = {(0,2):C.GC_70,(2,2):C.GC_73,(1,1):C.GC_575,(3,1):C.GC_578,(1,0):C.GC_575,(3,0):C.GC_578})

V_795 = Vertex(name = 'V_795',
               particles = [ P.u__tilde__, P.d, P.d__tilde__, P.u ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF2, L.FFFF3, L.FFFF4 ],
               couplings = {(1,7):C.GC_43,(3,7):C.GC_44,(0,8):C.GC_47,(2,8):C.GC_50,(1,0):C.GC_1197,(3,0):C.GC_1205,(0,6):C.GC_1196,(2,6):C.GC_1204,(1,5):C.GC_40,(3,5):C.GC_41,(1,3):C.GC_51,(3,3):C.GC_52,(1,4):C.GC_123,(3,4):C.GC_124,(1,1):C.GC_1196,(3,1):C.GC_1204,(0,2):C.GC_1197,(2,2):C.GC_1205})

V_796 = Vertex(name = 'V_796',
               particles = [ P.u__tilde__, P.d, P.d__tilde__, P.u ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF2, L.FFFF3 ],
               couplings = {(1,4):C.GC_45,(3,4):C.GC_48,(1,0):C.GC_1199,(3,0):C.GC_1207,(0,3):C.GC_1198,(2,3):C.GC_1206,(1,1):C.GC_1199,(3,1):C.GC_1207,(0,2):C.GC_1198,(2,2):C.GC_1206})

V_797 = Vertex(name = 'V_797',
               particles = [ P.u__tilde__, P.d, P.d__tilde__, P.u ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF2 ],
               couplings = {(1,0):C.GC_1201,(3,0):C.GC_1209,(0,3):C.GC_1200,(2,3):C.GC_1208,(1,1):C.GC_1200,(3,1):C.GC_1208,(0,2):C.GC_1201,(2,2):C.GC_1209})

V_798 = Vertex(name = 'V_798',
               particles = [ P.u__tilde__, P.d, P.d__tilde__, P.u ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF2 ],
               couplings = {(1,0):C.GC_1203,(3,0):C.GC_1211,(0,3):C.GC_1202,(2,3):C.GC_1210,(1,1):C.GC_1203,(3,1):C.GC_1211,(0,2):C.GC_1202,(2,2):C.GC_1210})

V_799 = Vertex(name = 'V_799',
               particles = [ P.c__tilde__, P.d, P.d__tilde__, P.c ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF2, L.FFFF3 ],
               couplings = {(1,7):C.GC_43,(3,7):C.GC_44,(1,0):C.GC_737,(3,0):C.GC_745,(0,6):C.GC_732,(2,6):C.GC_740,(1,5):C.GC_40,(3,5):C.GC_41,(1,3):C.GC_51,(3,3):C.GC_52,(1,4):C.GC_123,(3,4):C.GC_124,(1,1):C.GC_736,(3,1):C.GC_744,(0,2):C.GC_733,(2,2):C.GC_741})

V_800 = Vertex(name = 'V_800',
               particles = [ P.c__tilde__, P.d, P.d__tilde__, P.c ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF2, L.FFFF3 ],
               couplings = {(1,4):C.GC_45,(3,4):C.GC_48,(1,0):C.GC_739,(3,0):C.GC_747,(0,3):C.GC_734,(2,3):C.GC_742,(1,1):C.GC_739,(3,1):C.GC_747,(0,2):C.GC_734,(2,2):C.GC_742})

V_801 = Vertex(name = 'V_801',
               particles = [ P.u__tilde__, P.d, P.s__tilde__, P.c ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF2, L.FFFF4 ],
               couplings = {(0,4):C.GC_47,(2,4):C.GC_50,(1,0):C.GC_733,(3,0):C.GC_741,(0,3):C.GC_736,(2,3):C.GC_744,(1,1):C.GC_1251,(3,1):C.GC_1259,(0,2):C.GC_1256,(2,2):C.GC_1264})

V_802 = Vertex(name = 'V_802',
               particles = [ P.u__tilde__, P.d, P.s__tilde__, P.c ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF2 ],
               couplings = {(1,0):C.GC_735,(3,0):C.GC_743,(0,3):C.GC_738,(2,3):C.GC_746,(1,1):C.GC_1254,(3,1):C.GC_1262,(0,2):C.GC_1257,(2,2):C.GC_1265})

V_803 = Vertex(name = 'V_803',
               particles = [ P.c__tilde__, P.s, P.d__tilde__, P.u ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF2, L.FFFF4 ],
               couplings = {(0,4):C.GC_47,(2,4):C.GC_50,(1,0):C.GC_1252,(3,0):C.GC_1260,(0,3):C.GC_1255,(2,3):C.GC_1263,(1,1):C.GC_732,(3,1):C.GC_740,(0,2):C.GC_737,(2,2):C.GC_745})

V_804 = Vertex(name = 'V_804',
               particles = [ P.c__tilde__, P.s, P.d__tilde__, P.u ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF2 ],
               couplings = {(1,0):C.GC_1254,(3,0):C.GC_1262,(0,3):C.GC_1257,(2,3):C.GC_1265,(1,1):C.GC_735,(3,1):C.GC_743,(0,2):C.GC_738,(2,2):C.GC_746})

V_805 = Vertex(name = 'V_805',
               particles = [ P.u__tilde__, P.s, P.s__tilde__, P.u ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF2, L.FFFF3 ],
               couplings = {(1,7):C.GC_43,(3,7):C.GC_44,(1,0):C.GC_1256,(3,0):C.GC_1264,(0,6):C.GC_1251,(2,6):C.GC_1259,(1,5):C.GC_40,(3,5):C.GC_41,(1,3):C.GC_51,(3,3):C.GC_52,(1,4):C.GC_123,(3,4):C.GC_124,(1,1):C.GC_1255,(3,1):C.GC_1263,(0,2):C.GC_1252,(2,2):C.GC_1260})

V_806 = Vertex(name = 'V_806',
               particles = [ P.u__tilde__, P.s, P.s__tilde__, P.u ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF2, L.FFFF3 ],
               couplings = {(1,4):C.GC_45,(3,4):C.GC_48,(1,0):C.GC_1258,(3,0):C.GC_1266,(0,3):C.GC_1253,(2,3):C.GC_1261,(1,1):C.GC_1258,(3,1):C.GC_1266,(0,2):C.GC_1253,(2,2):C.GC_1261})

V_807 = Vertex(name = 'V_807',
               particles = [ P.c__tilde__, P.s, P.s__tilde__, P.c ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF2, L.FFFF3, L.FFFF4 ],
               couplings = {(1,7):C.GC_43,(3,7):C.GC_44,(0,8):C.GC_47,(2,8):C.GC_50,(1,0):C.GC_989,(3,0):C.GC_997,(0,6):C.GC_988,(2,6):C.GC_996,(1,5):C.GC_40,(3,5):C.GC_41,(1,3):C.GC_51,(3,3):C.GC_52,(1,4):C.GC_123,(3,4):C.GC_124,(1,1):C.GC_988,(3,1):C.GC_996,(0,2):C.GC_989,(2,2):C.GC_997})

V_808 = Vertex(name = 'V_808',
               particles = [ P.c__tilde__, P.s, P.s__tilde__, P.c ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF2, L.FFFF3 ],
               couplings = {(1,4):C.GC_45,(3,4):C.GC_48,(1,0):C.GC_991,(3,0):C.GC_999,(0,3):C.GC_990,(2,3):C.GC_998,(1,1):C.GC_991,(3,1):C.GC_999,(0,2):C.GC_990,(2,2):C.GC_998})

V_809 = Vertex(name = 'V_809',
               particles = [ P.c__tilde__, P.s, P.s__tilde__, P.c ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF2 ],
               couplings = {(1,0):C.GC_993,(3,0):C.GC_1001,(0,3):C.GC_992,(2,3):C.GC_1000,(1,1):C.GC_992,(3,1):C.GC_1000,(0,2):C.GC_993,(2,2):C.GC_1001})

V_810 = Vertex(name = 'V_810',
               particles = [ P.c__tilde__, P.s, P.s__tilde__, P.c ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF12, L.FFFF2 ],
               couplings = {(1,0):C.GC_995,(3,0):C.GC_1003,(0,3):C.GC_994,(2,3):C.GC_1002,(1,1):C.GC_995,(3,1):C.GC_1003,(0,2):C.GC_994,(2,2):C.GC_1002})

V_811 = Vertex(name = 'V_811',
               particles = [ P.u__tilde__, P.u, P.u__tilde__, P.u ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF13, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF17, L.FFFF18, L.FFFF3, L.FFFF4 ],
               couplings = {(1,6):C.GC_43,(3,6):C.GC_44,(0,7):C.GC_43,(2,7):C.GC_44,(0,0):C.GC_51,(2,0):C.GC_52,(1,3):C.GC_51,(3,3):C.GC_52,(1,1):C.GC_51,(3,1):C.GC_52,(1,2):C.GC_125,(3,2):C.GC_126,(0,4):C.GC_51,(2,4):C.GC_52,(0,5):C.GC_125,(2,5):C.GC_126})

V_812 = Vertex(name = 'V_812',
               particles = [ P.u__tilde__, P.u, P.u__tilde__, P.u ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF3, L.FFFF4 ],
               couplings = {(1,0):C.GC_46,(3,0):C.GC_49,(0,1):C.GC_46,(2,1):C.GC_49})

V_813 = Vertex(name = 'V_813',
               particles = [ P.u__tilde__, P.u, P.c__tilde__, P.c ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'T(-1,2,1)*T(-1,4,3)' ],
               lorentz = [ L.FFFF13, L.FFFF17, L.FFFF18, L.FFFF4 ],
               couplings = {(0,3):C.GC_43,(1,3):C.GC_44,(0,0):C.GC_51,(1,0):C.GC_52,(0,1):C.GC_51,(1,1):C.GC_52,(0,2):C.GC_125,(1,2):C.GC_126})

V_814 = Vertex(name = 'V_814',
               particles = [ P.u__tilde__, P.u, P.c__tilde__, P.c ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'T(-1,2,1)*T(-1,4,3)' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_46,(1,0):C.GC_49})

V_815 = Vertex(name = 'V_815',
               particles = [ P.c__tilde__, P.c, P.c__tilde__, P.c ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF13, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF17, L.FFFF18, L.FFFF3, L.FFFF4 ],
               couplings = {(1,6):C.GC_43,(3,6):C.GC_44,(0,7):C.GC_43,(2,7):C.GC_44,(0,0):C.GC_51,(2,0):C.GC_52,(1,3):C.GC_51,(3,3):C.GC_52,(1,1):C.GC_51,(3,1):C.GC_52,(1,2):C.GC_125,(3,2):C.GC_126,(0,4):C.GC_51,(2,4):C.GC_52,(0,5):C.GC_125,(2,5):C.GC_126})

V_816 = Vertex(name = 'V_816',
               particles = [ P.c__tilde__, P.c, P.c__tilde__, P.c ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF3, L.FFFF4 ],
               couplings = {(1,0):C.GC_46,(3,0):C.GC_49,(0,1):C.GC_46,(2,1):C.GC_49})

V_817 = Vertex(name = 'V_817',
               particles = [ P.t__tilde__, P.b, P.e__plus__, P.ve ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF11, L.FFFF12, L.FFFF19, L.FFFF20, L.FFFF21, L.FFFF4 ],
               couplings = {(0,5):C.GC_78,(0,4):C.GC_764,(0,2):C.GC_762,(0,3):C.GC_762,(0,0):C.GC_754,(0,1):C.GC_757})

V_818 = Vertex(name = 'V_818',
               particles = [ P.t__tilde__, P.b, P.e__plus__, P.ve ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF11, L.FFFF12, L.FFFF19, L.FFFF20, L.FFFF21 ],
               couplings = {(0,4):C.GC_767,(0,2):C.GC_766,(0,3):C.GC_766,(0,0):C.GC_756,(0,1):C.GC_760})

V_819 = Vertex(name = 'V_819',
               particles = [ P.t__tilde__, P.b, P.mu__plus__, P.vm ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF11, L.FFFF12, L.FFFF19, L.FFFF20, L.FFFF21, L.FFFF4 ],
               couplings = {(0,5):C.GC_78,(0,4):C.GC_844,(0,2):C.GC_842,(0,3):C.GC_842,(0,0):C.GC_834,(0,1):C.GC_837})

V_820 = Vertex(name = 'V_820',
               particles = [ P.t__tilde__, P.b, P.mu__plus__, P.vm ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF11, L.FFFF12, L.FFFF19, L.FFFF20, L.FFFF21 ],
               couplings = {(0,4):C.GC_847,(0,2):C.GC_846,(0,3):C.GC_846,(0,0):C.GC_836,(0,1):C.GC_840})

V_821 = Vertex(name = 'V_821',
               particles = [ P.t__tilde__, P.b, P.ta__plus__, P.vt ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF11, L.FFFF12, L.FFFF19, L.FFFF20, L.FFFF21, L.FFFF4 ],
               couplings = {(0,5):C.GC_78,(0,4):C.GC_1049,(0,2):C.GC_1047,(0,3):C.GC_1047,(0,0):C.GC_1039,(0,1):C.GC_1042})

V_822 = Vertex(name = 'V_822',
               particles = [ P.t__tilde__, P.b, P.ta__plus__, P.vt ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF11, L.FFFF12, L.FFFF19, L.FFFF20, L.FFFF21 ],
               couplings = {(0,4):C.GC_1052,(0,2):C.GC_1051,(0,3):C.GC_1051,(0,0):C.GC_1041,(0,1):C.GC_1045})

V_823 = Vertex(name = 'V_823',
               particles = [ P.u__tilde__, P.d, P.e__plus__, P.ve ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF11, L.FFFF12, L.FFFF19, L.FFFF20, L.FFFF21, L.FFFF4 ],
               couplings = {(0,5):C.GC_58,(0,4):C.GC_1234,(0,2):C.GC_1232,(0,3):C.GC_1232,(0,0):C.GC_825,(0,1):C.GC_1227})

V_824 = Vertex(name = 'V_824',
               particles = [ P.u__tilde__, P.d, P.e__plus__, P.ve ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF11, L.FFFF12, L.FFFF19, L.FFFF20, L.FFFF21 ],
               couplings = {(0,4):C.GC_1237,(0,2):C.GC_1236,(0,3):C.GC_1236,(0,0):C.GC_827,(0,1):C.GC_1230})

V_825 = Vertex(name = 'V_825',
               particles = [ P.u__tilde__, P.d, P.mu__plus__, P.vm ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF11, L.FFFF12, L.FFFF19, L.FFFF20, L.FFFF21, L.FFFF4 ],
               couplings = {(0,5):C.GC_58,(0,4):C.GC_1246,(0,2):C.GC_1244,(0,3):C.GC_1244,(0,0):C.GC_905,(0,1):C.GC_1239})

V_826 = Vertex(name = 'V_826',
               particles = [ P.u__tilde__, P.d, P.mu__plus__, P.vm ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF11, L.FFFF12, L.FFFF19, L.FFFF20, L.FFFF21 ],
               couplings = {(0,4):C.GC_1249,(0,2):C.GC_1248,(0,3):C.GC_1248,(0,0):C.GC_907,(0,1):C.GC_1242})

V_827 = Vertex(name = 'V_827',
               particles = [ P.u__tilde__, P.d, P.ta__plus__, P.vt ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF11, L.FFFF12, L.FFFF19, L.FFFF20, L.FFFF21, L.FFFF4 ],
               couplings = {(0,5):C.GC_58,(0,4):C.GC_1274,(0,2):C.GC_1272,(0,3):C.GC_1272,(0,0):C.GC_1110,(0,1):C.GC_1267})

V_828 = Vertex(name = 'V_828',
               particles = [ P.u__tilde__, P.d, P.ta__plus__, P.vt ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF11, L.FFFF12, L.FFFF19, L.FFFF20, L.FFFF21 ],
               couplings = {(0,4):C.GC_1277,(0,2):C.GC_1276,(0,3):C.GC_1276,(0,0):C.GC_1112,(0,1):C.GC_1270})

V_829 = Vertex(name = 'V_829',
               particles = [ P.c__tilde__, P.s, P.e__plus__, P.ve ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF11, L.FFFF12, L.FFFF19, L.FFFF20, L.FFFF21, L.FFFF4 ],
               couplings = {(0,5):C.GC_58,(0,4):C.GC_820,(0,2):C.GC_818,(0,3):C.GC_818,(0,0):C.GC_1019,(0,1):C.GC_813})

V_830 = Vertex(name = 'V_830',
               particles = [ P.c__tilde__, P.s, P.e__plus__, P.ve ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF11, L.FFFF12, L.FFFF19, L.FFFF20, L.FFFF21 ],
               couplings = {(0,4):C.GC_823,(0,2):C.GC_822,(0,3):C.GC_822,(0,0):C.GC_1021,(0,1):C.GC_816})

V_831 = Vertex(name = 'V_831',
               particles = [ P.c__tilde__, P.s, P.mu__plus__, P.vm ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF11, L.FFFF12, L.FFFF19, L.FFFF20, L.FFFF21, L.FFFF4 ],
               couplings = {(0,5):C.GC_58,(0,4):C.GC_900,(0,2):C.GC_898,(0,3):C.GC_898,(0,0):C.GC_1022,(0,1):C.GC_893})

V_832 = Vertex(name = 'V_832',
               particles = [ P.c__tilde__, P.s, P.mu__plus__, P.vm ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF11, L.FFFF12, L.FFFF19, L.FFFF20, L.FFFF21 ],
               couplings = {(0,4):C.GC_903,(0,2):C.GC_902,(0,3):C.GC_902,(0,0):C.GC_1024,(0,1):C.GC_896})

V_833 = Vertex(name = 'V_833',
               particles = [ P.c__tilde__, P.s, P.ta__plus__, P.vt ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF11, L.FFFF12, L.FFFF19, L.FFFF20, L.FFFF21, L.FFFF4 ],
               couplings = {(0,5):C.GC_58,(0,4):C.GC_1105,(0,2):C.GC_1103,(0,3):C.GC_1103,(0,0):C.GC_1113,(0,1):C.GC_1098})

V_834 = Vertex(name = 'V_834',
               particles = [ P.c__tilde__, P.s, P.ta__plus__, P.vt ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF11, L.FFFF12, L.FFFF19, L.FFFF20, L.FFFF21 ],
               couplings = {(0,4):C.GC_1108,(0,2):C.GC_1107,(0,3):C.GC_1107,(0,0):C.GC_1115,(0,1):C.GC_1101})

V_835 = Vertex(name = 'V_835',
               particles = [ P.b__tilde__, P.b, P.ve__tilde__, P.ve ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_75,(0,0):C.GC_20})

V_836 = Vertex(name = 'V_836',
               particles = [ P.b__tilde__, P.b, P.ve__tilde__, P.ve ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_76})

V_837 = Vertex(name = 'V_837',
               particles = [ P.b__tilde__, P.b, P.vm__tilde__, P.vm ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_75,(0,0):C.GC_20})

V_838 = Vertex(name = 'V_838',
               particles = [ P.b__tilde__, P.b, P.vm__tilde__, P.vm ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_76})

V_839 = Vertex(name = 'V_839',
               particles = [ P.b__tilde__, P.b, P.vt__tilde__, P.vt ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_75,(0,0):C.GC_20})

V_840 = Vertex(name = 'V_840',
               particles = [ P.b__tilde__, P.b, P.vt__tilde__, P.vt ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_76})

V_841 = Vertex(name = 'V_841',
               particles = [ P.d__tilde__, P.d, P.ve__tilde__, P.ve ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_55,(0,0):C.GC_53})

V_842 = Vertex(name = 'V_842',
               particles = [ P.d__tilde__, P.d, P.ve__tilde__, P.ve ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_56})

V_843 = Vertex(name = 'V_843',
               particles = [ P.d__tilde__, P.d, P.vm__tilde__, P.vm ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_55,(0,0):C.GC_53})

V_844 = Vertex(name = 'V_844',
               particles = [ P.d__tilde__, P.d, P.vm__tilde__, P.vm ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_56})

V_845 = Vertex(name = 'V_845',
               particles = [ P.d__tilde__, P.d, P.vt__tilde__, P.vt ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_55,(0,0):C.GC_53})

V_846 = Vertex(name = 'V_846',
               particles = [ P.d__tilde__, P.d, P.vt__tilde__, P.vt ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_56})

V_847 = Vertex(name = 'V_847',
               particles = [ P.s__tilde__, P.s, P.ve__tilde__, P.ve ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_55,(0,0):C.GC_53})

V_848 = Vertex(name = 'V_848',
               particles = [ P.s__tilde__, P.s, P.ve__tilde__, P.ve ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_56})

V_849 = Vertex(name = 'V_849',
               particles = [ P.s__tilde__, P.s, P.vm__tilde__, P.vm ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_55,(0,0):C.GC_53})

V_850 = Vertex(name = 'V_850',
               particles = [ P.s__tilde__, P.s, P.vm__tilde__, P.vm ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_56})

V_851 = Vertex(name = 'V_851',
               particles = [ P.s__tilde__, P.s, P.vt__tilde__, P.vt ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_55,(0,0):C.GC_53})

V_852 = Vertex(name = 'V_852',
               particles = [ P.s__tilde__, P.s, P.vt__tilde__, P.vt ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_56})

V_853 = Vertex(name = 'V_853',
               particles = [ P.t__tilde__, P.t, P.ve__tilde__, P.ve ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_75,(0,0):C.GC_114})

V_854 = Vertex(name = 'V_854',
               particles = [ P.t__tilde__, P.t, P.ve__tilde__, P.ve ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_77})

V_855 = Vertex(name = 'V_855',
               particles = [ P.t__tilde__, P.t, P.vm__tilde__, P.vm ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_75,(0,0):C.GC_114})

V_856 = Vertex(name = 'V_856',
               particles = [ P.t__tilde__, P.t, P.vm__tilde__, P.vm ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_77})

V_857 = Vertex(name = 'V_857',
               particles = [ P.t__tilde__, P.t, P.vt__tilde__, P.vt ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_75,(0,0):C.GC_114})

V_858 = Vertex(name = 'V_858',
               particles = [ P.t__tilde__, P.t, P.vt__tilde__, P.vt ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_77})

V_859 = Vertex(name = 'V_859',
               particles = [ P.u__tilde__, P.u, P.ve__tilde__, P.ve ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_55,(0,0):C.GC_61})

V_860 = Vertex(name = 'V_860',
               particles = [ P.u__tilde__, P.u, P.ve__tilde__, P.ve ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_57})

V_861 = Vertex(name = 'V_861',
               particles = [ P.u__tilde__, P.u, P.vm__tilde__, P.vm ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_55,(0,0):C.GC_61})

V_862 = Vertex(name = 'V_862',
               particles = [ P.u__tilde__, P.u, P.vm__tilde__, P.vm ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_57})

V_863 = Vertex(name = 'V_863',
               particles = [ P.u__tilde__, P.u, P.vt__tilde__, P.vt ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_55,(0,0):C.GC_61})

V_864 = Vertex(name = 'V_864',
               particles = [ P.u__tilde__, P.u, P.vt__tilde__, P.vt ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_57})

V_865 = Vertex(name = 'V_865',
               particles = [ P.c__tilde__, P.c, P.ve__tilde__, P.ve ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_55,(0,0):C.GC_61})

V_866 = Vertex(name = 'V_866',
               particles = [ P.c__tilde__, P.c, P.ve__tilde__, P.ve ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_57})

V_867 = Vertex(name = 'V_867',
               particles = [ P.c__tilde__, P.c, P.vm__tilde__, P.vm ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_55,(0,0):C.GC_61})

V_868 = Vertex(name = 'V_868',
               particles = [ P.c__tilde__, P.c, P.vm__tilde__, P.vm ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_57})

V_869 = Vertex(name = 'V_869',
               particles = [ P.c__tilde__, P.c, P.vt__tilde__, P.vt ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_55,(0,0):C.GC_61})

V_870 = Vertex(name = 'V_870',
               particles = [ P.c__tilde__, P.c, P.vt__tilde__, P.vt ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_57})

V_871 = Vertex(name = 'V_871',
               particles = [ P.e__plus__, P.e__minus__, P.ve__tilde__, P.ve ],
               color = [ '1' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_59,(0,0):C.GC_54})

V_872 = Vertex(name = 'V_872',
               particles = [ P.e__plus__, P.e__minus__, P.ve__tilde__, P.ve ],
               color = [ '1' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_60})

V_873 = Vertex(name = 'V_873',
               particles = [ P.e__plus__, P.e__minus__, P.vm__tilde__, P.vm ],
               color = [ '1' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_59,(0,0):C.GC_54})

V_874 = Vertex(name = 'V_874',
               particles = [ P.e__plus__, P.e__minus__, P.vt__tilde__, P.vt ],
               color = [ '1' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_59,(0,0):C.GC_54})

V_875 = Vertex(name = 'V_875',
               particles = [ P.mu__plus__, P.e__minus__, P.ve__tilde__, P.vm ],
               color = [ '1' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_60})

V_876 = Vertex(name = 'V_876',
               particles = [ P.ta__plus__, P.e__minus__, P.ve__tilde__, P.vt ],
               color = [ '1' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_60})

V_877 = Vertex(name = 'V_877',
               particles = [ P.e__plus__, P.mu__minus__, P.vm__tilde__, P.ve ],
               color = [ '1' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_60})

V_878 = Vertex(name = 'V_878',
               particles = [ P.mu__plus__, P.mu__minus__, P.ve__tilde__, P.ve ],
               color = [ '1' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_59,(0,0):C.GC_54})

V_879 = Vertex(name = 'V_879',
               particles = [ P.mu__plus__, P.mu__minus__, P.vm__tilde__, P.vm ],
               color = [ '1' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_59,(0,0):C.GC_54})

V_880 = Vertex(name = 'V_880',
               particles = [ P.mu__plus__, P.mu__minus__, P.vm__tilde__, P.vm ],
               color = [ '1' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_60})

V_881 = Vertex(name = 'V_881',
               particles = [ P.mu__plus__, P.mu__minus__, P.vt__tilde__, P.vt ],
               color = [ '1' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_59,(0,0):C.GC_54})

V_882 = Vertex(name = 'V_882',
               particles = [ P.ta__plus__, P.mu__minus__, P.vm__tilde__, P.vt ],
               color = [ '1' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_60})

V_883 = Vertex(name = 'V_883',
               particles = [ P.e__plus__, P.ta__minus__, P.vt__tilde__, P.ve ],
               color = [ '1' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_60})

V_884 = Vertex(name = 'V_884',
               particles = [ P.mu__plus__, P.ta__minus__, P.vt__tilde__, P.vm ],
               color = [ '1' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_60})

V_885 = Vertex(name = 'V_885',
               particles = [ P.ta__plus__, P.ta__minus__, P.ve__tilde__, P.ve ],
               color = [ '1' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_59,(0,0):C.GC_54})

V_886 = Vertex(name = 'V_886',
               particles = [ P.ta__plus__, P.ta__minus__, P.vm__tilde__, P.vm ],
               color = [ '1' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_59,(0,0):C.GC_54})

V_887 = Vertex(name = 'V_887',
               particles = [ P.ta__plus__, P.ta__minus__, P.vt__tilde__, P.vt ],
               color = [ '1' ],
               lorentz = [ L.FFFF13, L.FFFF4 ],
               couplings = {(0,1):C.GC_59,(0,0):C.GC_54})

V_888 = Vertex(name = 'V_888',
               particles = [ P.ta__plus__, P.ta__minus__, P.vt__tilde__, P.vt ],
               color = [ '1' ],
               lorentz = [ L.FFFF4 ],
               couplings = {(0,0):C.GC_60})

V_889 = Vertex(name = 'V_889',
               particles = [ P.ve__tilde__, P.ve, P.ve__tilde__, P.ve ],
               color = [ '1' ],
               lorentz = [ L.FFFF3, L.FFFF4 ],
               couplings = {(0,0):C.GC_59,(0,1):C.GC_59})

V_890 = Vertex(name = 'V_890',
               particles = [ P.ve__tilde__, P.ve, P.ve__tilde__, P.ve ],
               color = [ '1' ],
               lorentz = [ L.FFFF3, L.FFFF4 ],
               couplings = {(0,0):C.GC_60,(0,1):C.GC_60})

V_891 = Vertex(name = 'V_891',
               particles = [ P.vm__tilde__, P.ve, P.ve__tilde__, P.vm ],
               color = [ '1' ],
               lorentz = [ L.FFFF3, L.FFFF4 ],
               couplings = {(0,0):C.GC_59,(0,1):C.GC_60})

V_892 = Vertex(name = 'V_892',
               particles = [ P.vt__tilde__, P.ve, P.ve__tilde__, P.vt ],
               color = [ '1' ],
               lorentz = [ L.FFFF3, L.FFFF4 ],
               couplings = {(0,0):C.GC_59,(0,1):C.GC_60})

V_893 = Vertex(name = 'V_893',
               particles = [ P.vm__tilde__, P.vm, P.vm__tilde__, P.vm ],
               color = [ '1' ],
               lorentz = [ L.FFFF3, L.FFFF4 ],
               couplings = {(0,0):C.GC_59,(0,1):C.GC_59})

V_894 = Vertex(name = 'V_894',
               particles = [ P.vm__tilde__, P.vm, P.vm__tilde__, P.vm ],
               color = [ '1' ],
               lorentz = [ L.FFFF3, L.FFFF4 ],
               couplings = {(0,0):C.GC_60,(0,1):C.GC_60})

V_895 = Vertex(name = 'V_895',
               particles = [ P.vt__tilde__, P.vm, P.vm__tilde__, P.vt ],
               color = [ '1' ],
               lorentz = [ L.FFFF3, L.FFFF4 ],
               couplings = {(0,0):C.GC_59,(0,1):C.GC_60})

V_896 = Vertex(name = 'V_896',
               particles = [ P.vt__tilde__, P.vt, P.vt__tilde__, P.vt ],
               color = [ '1' ],
               lorentz = [ L.FFFF3, L.FFFF4 ],
               couplings = {(0,0):C.GC_59,(0,1):C.GC_59})

V_897 = Vertex(name = 'V_897',
               particles = [ P.vt__tilde__, P.vt, P.vt__tilde__, P.vt ],
               color = [ '1' ],
               lorentz = [ L.FFFF3, L.FFFF4 ],
               couplings = {(0,0):C.GC_60,(0,1):C.GC_60})

V_898 = Vertex(name = 'V_898',
               particles = [ P.a, P.a, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.VVS4 ],
               couplings = {(0,0):C.GC_317})

V_899 = Vertex(name = 'V_899',
               particles = [ P.g, P.g, P.H1 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.VVS4, L.VVS5, L.VVS6, L.VVS7 ],
               couplings = {(0,0):C.GC_318,(0,2):C.GC_331,(0,1):C.GC_327,(0,3):C.GC_322})

V_900 = Vertex(name = 'V_900',
               particles = [ P.a, P.Z, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.VVS4 ],
               couplings = {(0,0):C.GC_321})

V_901 = Vertex(name = 'V_901',
               particles = [ P.a, P.Z1, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVS4 ],
               couplings = {(0,0):C.GC_321})

V_902 = Vertex(name = 'V_902',
               particles = [ P.a, P.Z1, P.H1 ],
               color = [ '1' ],
               lorentz = [ L.VVS4 ],
               couplings = {(0,0):C.GC_335})

V_903 = Vertex(name = 'V_903',
               particles = [ P.g, P.g, P.g, P.H1 ],
               color = [ 'f(1,2,3)' ],
               lorentz = [ L.VVVS10, L.VVVS6, L.VVVS7, L.VVVS8, L.VVVS9 ],
               couplings = {(0,3):C.GC_323,(0,0):C.GC_332,(0,4):C.GC_328,(0,2):C.GC_325,(0,1):C.GC_319})

V_904 = Vertex(name = 'V_904',
               particles = [ P.g, P.g, P.g, P.g, P.H1 ],
               color = [ 'f(-1,1,2)*f(-1,3,4)', 'f(-1,1,3)*f(-1,2,4)', 'f(-1,1,4)*f(-1,2,3)' ],
               lorentz = [ L.VVVVS1, L.VVVVS10, L.VVVVS11, L.VVVVS12, L.VVVVS13, L.VVVVS14, L.VVVVS15, L.VVVVS16, L.VVVVS18, L.VVVVS2, L.VVVVS20, L.VVVVS3, L.VVVVS5, L.VVVVS7, L.VVVVS8 ],
               couplings = {(2,5):C.GC_324,(2,8):C.GC_333,(1,4):C.GC_324,(1,10):C.GC_333,(2,6):C.GC_330,(0,11):C.GC_326,(0,12):C.GC_334,(1,7):C.GC_330,(0,3):C.GC_329,(1,2):C.GC_326,(2,1):C.GC_326,(0,9):C.GC_324,(1,13):C.GC_320,(0,0):C.GC_320,(2,14):C.GC_320})

