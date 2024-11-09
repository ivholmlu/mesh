import meshio

msh = meshio.read("data/simple.msh")
print(msh)

points = msh.points
cells = msh.cells

print("Checking triangle 221")
print("Type: ", cells[1].type)
print("Data: ", cells[1].data[221])


