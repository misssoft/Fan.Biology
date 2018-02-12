'''
#Test tree_2_matrix.py
#python -m unittest tests.test_tree-2_matrix (under Fan.Biology)
'''
import unittest
from src import tree_2_matrix

class TestTree2Matrix(unittest.TestCase):
    def test_read_path_to_dictionary(self):
        input_file = 'E:\\Share\\Code\\Fan.Biology\\tests\\test_tree.txt'
        item, dictionary, weighed_dictionary = tree_2_matrix.read_path_to_dictionary(input_file)
        self.assertEqual(4, item)
        self.assertEqual(6, len(dictionary))
        self.assertEqual(6, len(weighed_dictionary))

    def test_find_path(self):
        input_file = 'E:\\Share\\Code\\Fan.Biology\\tests\\test_tree.txt'
        item, dictionary, weighed_dictionary = tree_2_matrix.read_path_to_dictionary(input_file)
        path = tree_2_matrix.find_path(dictionary,0,1)
        self.assertEqual(3, len(path))

    def test_get_weight_of_path(self):
        input_file = 'E:\\Share\\Code\\Fan.Biology\\tests\\test_tree.txt'
        item, dictionary, weighed_dictionary = tree_2_matrix.read_path_to_dictionary(input_file)
        path = tree_2_matrix.find_path(dictionary, 0, 1)
        weight = tree_2_matrix.get_weight_of_path(weighed_dictionary,path)
        self.assertEqual(13,weight)