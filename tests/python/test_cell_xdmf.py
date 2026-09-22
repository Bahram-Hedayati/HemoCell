"""Regression coverage for scalar particle metadata that crashed ParaView.

Run: python3 -m unittest discover -s tests/python -p 'test_*.py'
"""
import importlib.util
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET

import h5py

spec = importlib.util.spec_from_file_location(
    'cell_xdmf', Path(__file__).resolve().parents[2] / 'scripts/CellHDF5toXMF.py')
converter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(converter)


class ParticleAttributeMetadata(unittest.TestCase):
    def test_scalar_vector_and_matrix_shapes(self):
        with tempfile.TemporaryDirectory() as directory:
            with h5py.File(Path(directory) / 'cells.h5', 'w') as handle:
                handle.create_dataset('Position', shape=(4, 3), dtype='f4')
                handle.create_dataset('Cell Id', shape=(4, 1), dtype='i4')
                handle.create_dataset('Flat scalar', shape=(4,), dtype='f4')
                handle.create_dataset('Matrix', shape=(4, 2, 2), dtype='f4')
                data = converter.updateDictForXDMFStringsCell(handle, {})
                data['pathToHDF5'] = 'cells.h5'
                expected = {'Position': ('Vector', '4 3'),
                            'Cell Id': ('Scalar', '4 1'),
                            'Flat scalar': ('Scalar', '4'),
                            'Matrix': ('Matrix', '4 2 2')}
                for item in converter.iteratePossibleDataSetsDict(data):
                    xml = converter.createH5AttibuteCell(converter.XMLIndentation()) % item
                    element = ET.fromstring(xml)
                    kind, dimensions = expected[element.attrib['Name']]
                    self.assertEqual(element.attrib['AttributeType'], kind)
                    self.assertEqual(element.attrib['Center'], 'Node')
                    self.assertEqual(element.find('DataItem').attrib['Dimensions'].strip(), dimensions)


if __name__ == '__main__':
    unittest.main()
