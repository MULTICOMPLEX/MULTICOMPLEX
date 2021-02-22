import numpy as np

def faddeev_leverrier(A):
    '''
    Given an n x n matrix A, we return the coefficients of it's
    characteristic polynomial
    P(x) ^= det(xI - A) = a_0 * x^n + a_1 * x^(n - 1) + ... + a_n
    It is a property of P(x) that a_0 = 1, a_n = det(A)
    We return the list a = [a_0, a_1, ..., a_n]
    '''
    A = np.array(A) #Ensure we have a numpy array
    n = A.shape[0]
    assert A.shape[1] == n, 'Array must be square!'

    a = np.array([1.])
    Ak = np.array(A)
    for k in range(1, n + 1):
        ak = -Ak.trace() / k
        a = np.append(a, ak)
        Ak += np.diag(np.repeat(ak, n))
        Ak = np.dot(A, Ak)
    return a

if __name__ == '__main__':
  A = np.array([[1., 2.], [3., 4.]])
  print (faddeev_leverrier(A))
 
  A = [[2.,1.,4.,5.,1.,-3.],
       [-1.,0.,2.,9.,4.,5.],
       [3.,4.,-1.,-1.,-1.,-1.],
       [1.,2.,3.,4.,5.,6.],
       [-6.,-5.,-4.,-3.,-2.,-1.],
       [1.,-1.,10.,-10.,100.,-100.]]
  #[1, 97, -132, 778, 59654, 190191, 189700] from wolframalpha
  
  print (faddeev_leverrier(A))