import pyves
from pprint import pprint

p1 = pyves.Pointf(1,2,3)
p2 = pyves.Pointf(0,1,2)
print("p1+p2", p1+p2)

print()
pos1 = pyves.Positioni(1,2,3)
pos2 = pyves.Positioni(1,1,1)
print("pos1", pos1)
print("pos2", pos2)
print("pos1+pos2", pos1+pos2)
print("pos1", pos1)
pos1.translation(pos2)
print("pos1.translation(pos2)", pos1)
pos1 += pos2
print("pos1 += pos2", pos1)
print("pos2", pos2)

print()
pos3 = pyves.Positioni(pyves.Pointi(9,9,9))
print("pos3", pos3)
pos4 = pyves.Positioni(pyves.Positioni(9,9,9))
print("pos3", pos3)