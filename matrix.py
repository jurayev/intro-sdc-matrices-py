import math
from math import sqrt
import numbers

def zeroes(height, width):
        """Creates a matrix of zeroes.
        
        Args:
            height: matrix rows number.
            width: matrix columns number.
         
        Returns:
            A matrix object. 
            
            For example:
            
            > z = zeroes(3, 3)
            > z
              Matrix([[0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0]])
        """
        g = [[0.0 for _ in range(width)] for __ in range(height)]
        return Matrix(g)

def identity(n):
        """Creates a n x n identity matrix.
        
        Args:
            n: - matrix dimension.
        
        Returns:
            A matrix object. 
            
            For example:
            
            > I = identity(3)
            > I
              Matrix([[1, 0.0, 0.0],
                      [0.0, 1, 0.0],
                      [0.0, 0.0, 1]])
        """
        I = zeroes(n, n)
        for i in range(n):
            I.g[i][i] = 1.0
        return I

def vectors_product(vector_a, vector_b):
    product = [vector_a[i] * vector_b[i] for i in range(len(vector_a))]
    return sum(product)
            
class Matrix(object):

    # Constructor
    def __init__(self, grid):
        self.g = grid
        self.h = len(grid)
        self.w = len(grid[0])

    #
    # Primary matrix math methods
    #############################
 
    def determinant(self):
        """Calculates the determinant of a 1x1 or 2x2 matrix.
                     
        Returns:
            A determinant of the given matrix object. For example:

            > m = Matrix([[2, 2], 
                          [3, 1]])
            > print(m.inverse)
              -4

        Raises:
            ValueError: Cannot calculate determinant of non-square matrix.
            NotImplementedError: Calculating determinant not implemented for matrices largerer than 2x2.
        """    
        if not self.is_square():
            raise ValueError("Cannot calculate determinant of non-square matrix.")
        if self.h > 2:
            raise NotImplementedError("Calculating determinant not implemented for matrices largerer than 2x2.")
        det_g = []    
        if self.h == 1:
            det_g = self.g[0][0]
        else:
            det_g = self.g[0][0] * self.g[1][1] - self.g[0][1] * self.g[1][0]
        return det_g 

    def trace(self):
        """Calculates the trace of a matrix.
        
        Returns:
            A sum of diagonal entries for the given matrix object. 
            
            For example:
            
            > m = Matrix([[100, 5, 3],
                           [6, 100, 3],
                           [6, 1, 100]])
            > print(m.trace())
              300  
            
        Raises:
            ValueError: Cannot calculate the trace of a non-square matrix.
        """
        if not self.is_square():
            raise ValueError("Cannot calculate the trace of a non-square matrix.")
        return sum([self.g[i][i] for i in range(self.h)])

    def inverse(self):
        """Calculates the inverse of a 1x1 or 2x2 Matrix.
        
        Returns:
            An inverse matrix of the given matrix object. 
            
            For example:
            
            > m = Matrix([[2, 2], 
                          [3, 1]])
            > print(m.inverse())
              -0.25  0.5 
              0.75  -0.5 
            
        Raises:
            ValueError: Non-square Matrix does not have an inverse.
            NotImplementedError: Inversion not implemented for matrices larger than 2x2.
            ValueError: Determinant of the given matrix is 0. Non-inverse matrix.
        """
        if not self.is_square():
            raise ValueError("Non-square Matrix does not have an inverse.")
        if self.h > 2:
            raise NotImplementedError("Inversion not implemented for matrices larger than 2x2.")
        det_g = self.determinant()
        if det_g == 0:
                raise ValueError("Determinant of the given matrix is 0. Non-inverse matrix.")
        if self.h == 1:
            return Matrix([[1 / det_g]])
        return 1.0 / det_g * (self.trace() * identity(self.w) - self)   
        
    def T(self):
        """Transposes the given matrix.
        
        Returns: 
            A transposed copy of the given matrix. 
            
            For example:
            
            > m = Matrix([[4, 5, 3],
                          [6, 7, 1]])
            > print(m.T()) 
              4  6
              5  7
              3  1        
        """
        trans_g = zeroes(self.w, self.h)
        for i in range(self.w):
            for j in range(self.h):
                trans_g[i][j] = self.g[j][i]   
        return trans_g    

    def is_square(self):
        return self.h == self.w
    
    #
    # Begin Operator Overloading
    ############################
    def __getitem__(self, idx):
        """Defines the behavior of using square brackets [] on instances
        of this class.
        
        Returns:
            Value of the given position of the matrix grid.
            
            Example:

            > my_matrix = Matrix([ [1, 2], [3, 4] ])
            > my_matrix[0]
              [1, 2]

            > my_matrix[0][0]
              1
        """
        return self.g[idx]

    def __repr__(self):
        """Defines the behavior of calling print on an instance of this class.
        
        Returns:
            Formatted string representation.
        """
        s = ""
        for row in self.g:
            s += " ".join(["{} ".format(x) for x in row])
            s += "\n"
        return s

    def __add__(self, other):
        """Defines the behavior of the "+" addition operator.
        
        Returns:
            Sum of two matrix objects with the same dimentions. 
            
            For example:
            
            > matrix_a = Matrix([[4, 5],[6, 7]])
            > matrix_b = Matrix([[6, 9],[1, 3]])
            > matrix_a + matrix_b
              [[12, 14], [7, 10]]
        
        Raises:
            ValueError: Matrices can only be added if the dimensions are the same.
        """
        if self.h != other.h or self.w != other.w:
            raise ValueError("Matrices can only be added if the dimensions are the same.")
        add_product = zeroes(self.h, self.w)    
        for i in range(self.h):
            for j in range(self.w):
                add_product[i][j] = self.g[i][j] + other.g[i][j] 
        return add_product    

    def __neg__(self):
        """
        Defines the behavior of "-" negative operator (NOT subtraction).
        
        Returns:
            Matrix object with negative operator applied.
            
        For example:

        > my_matrix = Matrix([ [1, 2], [3, 4] ])
        > negative  = -my_matrix
        > print(negative)
          -1.0  -2.0
          -3.0  -4.0
        """
        neg_product = [[-self.g[i][j] for j in range(self.w)] for i in range(self.h)]
        return Matrix(neg_product)

    def __sub__(self, other):
        """
        Defines the behavior of "-" subtraction operator.
        
        Returns:
            Substraction product of two matrix objects with the same dimentions. 
            
            For example:
            
            > matrix_a = Matrix([[4, 5],[6, 7]]) 
            > matrix_b = Matrix([[6, 9],[1, -3]])
            > matrix_a - matrix_b
              [[-2, -4],[5, 10]]
        
        Raises:
            ValueError: Matrices can only be substracted if the dimensions are the same.
        """
        if self.h != other.h or self.w != other.w:
            raise ValueError("Matrices can only be substracted if the dimensions are the same.")
        sub_product = zeroes(self.h, self.w)    
        for i in range(self.h):
            for j in range(self.w):
                sub_product[i][j] = self.g[i][j] - other.g[i][j] 
        return sub_product

    def __mul__(self, other):
        """Defines the behavior of "*" operator (matrix multiplication)
        
        Returns:
            Product of matrix multiplication with the product size of self.h x other.w.
            
            For example:
            
            > matrix_a = Matrix([[4, 5, 1],[6, -7, 1]])
            > matrix_b = Matrix([[5, 2],[-4, 3],[-4, 1]])
            > print(matrix_a * matrix_b)
              -4  24 
              54  -8
            
        Raises:
            ValueError: Matrices can only be substracted if first matrix width equals second matrix height.
        """
        if self.w != other.h:
            raise ValueError("Matrices can only be substracted if first matrix width equals second matrix height.")
        t_other_g = other.T()
        mul_product = zeroes(self.h, other.w)
        for i in range(self.h):
            for j in range(other.w):
                mul_product[i][j] = vectors_product(self.g[i], t_other_g[j])
        return mul_product
    
    def __rmul__(self, other):
        """
        Called when the thing on the left of the * is not a matrix.
        
        Returns:
            Scalar multiplication of each matrix element with passed value.
            
            For example:

            > identity = Matrix([ [1,0], [0,1] ])
            > doubled  = 2 * identity
            > print(doubled)
              2.0  0.0
              0.0  2.0
              
        Raises:
            ValueError: Scalar multiplication supports instance of Number class only.
        """
        if isinstance(other, numbers.Number):
            rmul_product = zeroes(self.h, self.w)
            for i in range(self.h):
                for j in range(self.w):
                    rmul_product[i][j] = other * self.g[i][j]
            return rmul_product        
        raise ValueError("Scalar multiplication supports instance of Number class only.")
            
