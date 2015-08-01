import gatools
import numpy
import numpy.testing as nptest
import unittest
import tempfile

class gatoolsTest(unittest.TestCase):
    def test_LoadLabels(self):
        labels = gatools.LoadLabels('TestResources/labels.txt')
         
        self.assertEqual(labels, ['AB', 'C', 'DEF'])

    def test_SaveLabels(self):
        labels = ['AB', 'C', 'DEF']
        
        with tempfile.NamedTemporaryFile() as tempFile:
            gatools.SaveLabels(tempFile.name, labels)
            
            content = tempFile.readlines()
            self.assertEqual(content, ['AB\n', 'C\n', 'DEF'])
            
    def test_Save_and_LoadLabels(self):
        labels = ['AB', 'C', 'DEF']
        with tempfile.NamedTemporaryFile() as tempFile:
            gatools.SaveLabels(tempFile.name, labels)
            loaded_labels = gatools.LoadLabels(tempFile.name)
            
            self.assertEqual(loaded_labels, labels)
    
    def test_LoadPartition(self):
        partition = gatools.LoadPartition('TestResources/matrix.txt')
        
        expected = [[0, 1, 2], [2, 0, 3], [1, 2, 0]]
        self.assertEqual(partition, expected)
    
    def test_SavePartition(self):
        partition = [[0, 1] , [2, 3]]
        with tempfile.NamedTemporaryFile() as tempFile:
            gatools.SavePartition(tempFile.name, partition)
            
            content = tempFile.readlines()
            self.assertEqual(content, ['0 1\n', '2 3'])
    
    def test_LoadFromPajek(self):
        pass
    
    def test_Save2Pajek(self):
        pass
    
    def test_ExtractSubmatrix(self):
        testmat = numpy.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], int)
        
        # TODO: fix it- for me this runs into an error and I cannot fix it since I don't grasp the purpose of this function:
        # IndexError: index 3 is out of bounds for axis 0 with size 3
        result = gatools.ExtractSubmatrix(testmat, [0, 3])
        expected = numpy.array([[1, 4], [13, 16]])
        
        nptest.assert_array_equal(result, expected) 
    
    def test_SymmetriseMatrix(self):
        matrix = numpy.array([[1, 2], [0, 1]])
        nptest.assert_equal(gatools.SymmetriseMatrix(matrix), numpy.array([[1, 1], [1, 1]]))
    
    def test_LaplacianMatrix(self):
        pass
    
    def test_CleanPaths(self):
        paths = [[1, 2, 3], [3, 2, 1], [2, 4, 7, 8], [8, 7, 4, 2]]
        gatools.CleanPaths(paths)
        
        self.assertEqual([[1, 2, 3], [2, 4, 7, 8]], paths)
    
    def test_HammingDistance(self):
        array1 = numpy.array([[1, 2], [3, 4]])
        array2 = numpy.array([[1, 0], [3, 4]])
        array3 = numpy.array([[5, 6], [7, 8]])
        array4 = numpy.array([1, 2, 3, 4])
        
        # test identity
        self.assertEqual(0, gatools.HammingDistance(array1, array1))
        
        # test little difference
        self.assertEqual(0.25, gatools.HammingDistance(array1, array2))
        
        # test total difference
        self.assertEqual(1, gatools.HammingDistance(array1, array3))
        
        # test normalization
        self.assertEqual((1, 0), gatools.HammingDistance(array1, array3, normed=True))
        
        try:
            gatools.HammingDistance(array1, array4)
            self.fail()
        except:
            pass
        
    def test_NonZeroMin(self):
        array1 = numpy.array([-1, 0, 1, 2])
        self.assertEqual(-1, gatools.NonZeroMin(array1))
        
        array2 = numpy.array([0, 1, 2])
        self.assertEqual(1, gatools.NonZeroMin(array2))
    
    def test_CumulativeDistribution(self):
        pass
    
    def test_Factorial(self):
        self.assertEqual(1, gatools.Factorial(0))
        self.assertEqual(1, gatools.Factorial(1))
        self.assertEqual(2, gatools.Factorial(2))
        self.assertEqual(24, gatools.Factorial(4))
        
        try:
            gatools.Factorial(-1)
            self.fail()
        except:
            pass
    
    def test_BinomialCoefficient(self):
        # test trivial combinations returning 1
        self.assertEqual(1, gatools.BinomialCoefficient(1, 1))
        self.assertEqual(1, gatools.BinomialCoefficient(2, 2))
        self.assertEqual(1, gatools.BinomialCoefficient(17, 0))
        
        # test efficiency for m = 1, n != 1
        self.assertEqual(12345, gatools.BinomialCoefficient(12345, 1))
        
        # test some non-trivial combination
        self.assertEqual(10, gatools.BinomialCoefficient(5, 3))
        
        # test if algorithm is efficient on something big
        self.assertEqual(500, gatools.BinomialCoefficient(500, 499))
    
    def test_StdDeviation(self):
        data1 = numpy.array([1, 1, 1])
        data2 = numpy.array([-2, -2, 2, 2])
        
        self.assertEqual((1, 0), gatools.StdDeviation(data1))
        self.assertEqual((0, 2), gatools.StdDeviation(data2))
    
    def test_Quartiles(self):
        # trivial example
        data = numpy.array([1, 1, 1])
        self.assertEqual((1, 1, 1), gatools.Quartiles(data))
        
        # easy example
        data = numpy.array([4, 3, 2, 1, 0])
        self.assertEqual((1, 2, 3), gatools.Quartiles(data))
    
    def test_AllPermutations(self):
        permutations = gatools.AllPermutations([1, 2, 3])
        self.assertEqual([(1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1)], permutations)
        
        # test on boundary conditions
        self.assertEqual([], gatools.AllPermutations([]))
        self.assertEqual([1], gatools.AllPermutations([1]))
        
    
    def test_AllCombinations(self):
        combinations = gatools.AllCombinations((1, 2, 3, 4), 2)
        self.assertEqual([(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)], combinations)
        
        combinations = gatools.AllCombinations((1, 2, 3, 4), 3)
        self.assertEqual([(1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4)], combinations)
        
    
    def test_AllBipartitions(self):
        expected = [
            ((1,), (2, 3)),
            ((2,), (1, 3)),
            ((3,), (1, 2))
        ]
        self.assertEqual(expected, gatools.AllBipartitions((1, 2, 3)))
        
        expected = [
            ((1,), (2, 3, 4)),
            ((2,), (1, 3, 4)),
            ((3,), (1, 2, 4)),
            ((4,), (1, 2, 3)),
            ((1, 2), (3, 4)),
            ((1, 3), (2, 4)),
            ((1, 4), (2, 3))
        ]
        
        self.assertEqual(expected, gatools.AllBipartitions((1, 2, 3, 4)))
        
    def test_MeanCorrelation(self):
        pass

if __name__ == '__main__':
    unittest.main()