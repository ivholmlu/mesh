import meshio
from abc import ABC, abstractmethod
from typing import Union

class CellFactory:
    """Facotry class for creating cells"""
    def __init__(self):
        self.cell_types = {}

    def __register__(self, cell_type, cell_class):
        self.cell_types[cell_type] = cell_class
    
    def __call__(self, cell_type, id, points):
        return self.cell_types[cell_type](id, points)

class Point:
    """Object representing a point in 2D space.
    """
    def __init__(self, id, x, y):
        self.id = id
        self._coords = (x, y)
    
    def __str__(self):
        return f"Point {self.id} at x: {self._coords[0]:.2f} y: {self._coords[1]:.2f}" 
    
    def x(self):
        return self._coords[0]
    
    def y(self):
        return self._coords[1]
    
class Cell(ABC):
    """Abstract class representing a cell in 2D space.

    Args:
        ABC (_type_): 
    """
    def __init__(self, id, points = []):
        self._id = id
        self._points = points
        self.neighbours = []
    
    def __str__(self):
        return f"{self.__class__} {self._id} with points {self._points} and neighbours {self.neighbours}"
    
    def points(self):
        return set(self._points)

    @abstractmethod
    def get_neighbours(self):
        pass

    @abstractmethod
    def store_neighbours(self, cells):
        pass


class Triangle(Cell):
    """
    Object representing a triangle in 2D space.
    """
    def __init__(self, id, points):
        super().__init__(id, points)
    
    def get_neighbours(self, cells):
        pass

    def store_neighbours(self, cells):
        self_points = self.points()

        for idx, cell in enumerate(cells):
            if len(self.neighbours) == 3:
                break
            matching_points = self_points & cell.points()
            if len(matching_points) == 2:
                self.neighbours.append(idx)

class Line(Cell):
    """
    Object representing a line in 2D space.
    """
    def __init__(self, id, points = []):
        super().__init__(id, points)
    
    def store_neighbours(self, cells):
        self_points = self.points()

        for id, cell in enumerate(cells):
            matches = self_points & cell.points() # determine how many points cells share
            if (len(matches) == 1 and len(cell._points) == 2) or (len(matches) == 2 and id != self._id):
                if(len(self.neighbours) == 3):
                    break
                self.neighbours.append(id)
    
    def get_neighbours(self, cells):
        pass

class Mesh:
    """
    Object representing a mesh in 2D space.
    """
    def __init__(self, mesh_path):
        self.mesh = meshio.read(mesh_path)
        points = self.mesh.points
        self.cells = self.mesh.cells
        self._points = []
        for idx, point in enumerate(points):
            self._points.append(Point(idx, point[0], point[1]))
        
        factory = CellFactory()
        factory.__register__("triangle", Triangle)
        factory.__register__("line", Line)
        self._cells = []
        for cell in self.cells:
            cell_type = cell.type
            cell_points = cell.data
            for points in cell_points:
                self._cells.append(factory(cell_type, len(self._cells), points))
    
    def __str__(self):
        return f"Mesh with {len(self.points)} points and {self._cells} cells."
    
    def find_neighbours(self):
        for cell in self._cells:
            cell.store_neighbours(self._cells)
        

msh = Mesh("data/simple.msh")
msh.find_neighbours()

ngh_id = msh._cells[4].neighbours[2]

print(msh._cells[4])
print(msh._cells[189])
print(msh._cells[222])
print(msh._cells[226])
print(msh._cells[224])