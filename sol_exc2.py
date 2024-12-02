from main import Mesh

msh = Mesh("data/simple.msh")
msh.find_neighbours()

print(msh._cells[4])
print(msh._cells[189])
print(msh._cells[222]) 
