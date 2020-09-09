import unittest
import pyves
import numpy as np



class MainTest(unittest.TestCase):
    def test_cell_methods(self):
        box = pyves.BoxPBC(50,50,50)
        cell = pyves.Cell([0,1,0.5], [3,3,3], box)
        bb = cell.bounding_box
        self.assertAlmostEqual(bb.min[0], 0-pyves.CellBoundOffset()[0])
        self.assertAlmostEqual(bb.min[1], 1-pyves.CellBoundOffset()[1])
        self.assertAlmostEqual(bb.min[2], 0.5-pyves.CellBoundOffset()[2])
        bb.min = [0.1,1,0.5]
        self.assertAlmostEqual(bb.min[0], 0.1)
        self.assertAlmostEqual(bb.min[1], 1)
        self.assertAlmostEqual(bb.min[2], 0.5)
        cell.min = [0.1,1,0.5]
        self.assertAlmostEqual(bb.min[0], 0.1-pyves.CellBoundOffset()[0])
        self.assertAlmostEqual(bb.min[1], 1-pyves.CellBoundOffset()[1])
        self.assertAlmostEqual(bb.min[2], 0.5-pyves.CellBoundOffset()[2])
        
        box = pyves.BoxPBC(12,12,12)
        cell1 = pyves.Cell(min=[0,0,0], max=[3,3,3], box=box)
        cell2 = pyves.Cell(max=[6,6,6], min=[3,3,3], box=box)
        cell3 = pyves.Cell(min=[9,9,9], max=[12,12,12], box=box)
        cell4 = pyves.Cell(min=[3,0,0], max=[6,3,3], box=box)
        cell5 = pyves.Cell(min=[6,6,6], max=[9,9,9], box=box)
        self.assertFalse(cell1.isNeighbourOf(cell1))
        self.assertFalse(cell2.isNeighbourOf(cell2))
        self.assertFalse(cell3.isNeighbourOf(cell3))
        self.assertFalse(cell4.isNeighbourOf(cell4))
        self.assertTrue(cell1.isNeighbourOf(cell2))
        self.assertTrue(cell1.isNeighbourOf(cell3))
        self.assertTrue(cell1.isNeighbourOf(cell4))
        self.assertTrue(cell1.isNeighbourOf(cell4))
        self.assertFalse(cell1.isNeighbourOf(cell5))
        box = pyves.BoxPBC(13,13,13)
        cell1 = pyves.Cell(min=[0,0,0], max=[3,3,3], box=box)
        cell3 = pyves.Cell(min=[9,9,9], max=[12,12,12], box=box)
        self.assertFalse(cell1.isNeighbourOf(cell3))
        



if __name__ == '__main__':
    unittest.main()