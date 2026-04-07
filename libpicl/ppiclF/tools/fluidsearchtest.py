import math

# these are fortran (1 based) indexes
# Indicies 1-3: Centroid x,y,z position
# Indicies 4-6: Max cell dx,dy,dz based on any vertex combination
# Index      7: Cell Volume 
class FluidCell:
    def __init__(self, homeRank: int, index: int, centroidPos: tuple[float,float,float], maxDX: tuple[float,float,float], volume: float):
        self.homeRank: int = homeRank
        self.index: int = index # fortran based indexing
        self.centroidPos: tuple[float,float,float] = centroidPos
        self.maxDX: tuple[float,float,float] = maxDX
        self.volume: float = volume
    def __str__(self):
        # return "(" + str(self.homeRank) + "," + str(self.index) + ")["+str(self.centroidPos) + "," + str(self.maxDX) + "," + str(self.volume) + "]"
        return "(" + str(self.homeRank) + "," + str(self.index) + ")"
    def __repr__(self):
        return str(self)
class Particle:
    def __init__(self, pos: tuple[float,float,float]):
        self.pos: tuple[float,float,float] = pos

def LoadFluidCells(rank, fileName) -> list[FluidCell]:
    with open(fileName, mode="r") as f:
        lines = f.readlines()
    # print(len(lines))
    # print(lines[0])
    lineFloats = [[float(v) for v in line.split()] for line in lines]
    # print(lineFloats[0])
    cells = [FluidCell(rank, i + 1, tuple(line[0:3]), tuple(line[3:6]), line[6]) for i, line in enumerate(lineFloats)]
    return cells


AllCells: list[FluidCell] = []

for i in range(16):
    name = "rank" + str(i).zfill(3) + "_fluid_grid.txt"
    AllCells += LoadFluidCells(i, name)

print("Total: " + str(len(AllCells)))

testParticle = Particle((0.13736237586016914, 0.0069842089924479689, 0.0071062286821715502))
nearest: list[tuple[int, FluidCell, float, bool, tuple[float,float,float]]] = [] # index(in python array), Cell, distance^2, inside cell 
ppiclf_interp_dchk = (0.0030120495000000302, 0.0025000000050000009,0.0025000000050000009)
dSQ_check_components = [ppiclf_interp_dchk[l] ** 2 for l in range(3)]
# print(dSQ_check_components)
for ie in range(len(AllCells)):
    # print(AllCells[ie].centroidPos[1])
    # print(testParticle.pos[1])
    dSQ_components = [(AllCells[ie].centroidPos[l] - testParticle.pos[l])**2 for l in range(3)]
    dSQ = sum(dSQ_components)
    if any([c > dSQ_check_components[l] for l, c in enumerate(dSQ_components)]):
        continue
    inCell = all([abs(AllCells[ie].centroidPos[l] - testParticle.pos[l]) < ((ppiclf_interp_dchk[l]/1.5) *0.5) for l in range(3)])

    nearest.append((ie, AllCells[ie], dSQ, inCell))
nearest.sort(key=lambda x: x[2])
print("NNearest: " + str(len(nearest)))
for i in range(len(nearest)):
    print(nearest[i])
    print("")