```python
from cvxopt import matrix, solvers
Q = 2*matrix([ [2, .5], [.5, 1] ])
p = matrix([1.0, 1.0])
G = matrix([[-1.0,0.0],[0.0,-1.0]])
h = matrix([0.0,0.0])
A = matrix([1.0, 1.0], (1,2))
b = matrix(1.0)
sol=solvers.qp(Q, p, G, h, A, b)
```

         pcost       dcost       gap    pres   dres
     0:  1.8889e+00  7.7778e-01  1e+00  3e-16  2e+00
     1:  1.8769e+00  1.8320e+00  4e-02  2e-16  6e-02
     2:  1.8750e+00  1.8739e+00  1e-03  2e-16  5e-04
     3:  1.8750e+00  1.8750e+00  1e-05  1e-16  5e-06
     4:  1.8750e+00  1.8750e+00  1e-07  1e-16  5e-08
    Optimal solution found.



```python
print(sol['x'])
```

    [ 2.50e-01]
    [ 7.50e-01]
    



```python
import numpy as np
from cvxopt import matrix, solvers
```


```python
dmatrix = np.array([[1956, 1425,  966],[1534,916,709],[1333 ,1274 ,1058],[2551 , 576, 1489],[673, 2174 ,1233 ],[890, 1194 ,1460]])
```


```python
dmatrix = np.array([[1956, 1425,  966],[1534,916,709],[1333 ,1274 ,1058]])
```


```python
dmatrix = dmatrix / np.sum(dmatrix,axis=1)[:, np.newaxis]
```


```python
print(dmatrix)
```

    [[0.44996549 0.32781228 0.22222222]
     [0.48559671 0.28996518 0.22443811]
     [0.36371078 0.34761255 0.28867667]]



```python
Qmatrix = np.zeros([9,9])
Qmatrix[0,0] = np.sum(dmatrix[:2,0] * dmatrix[:2,0])
Qmatrix[3,3] = Qmatrix[0,0]
Qmatrix[6,6] = Qmatrix[0,0]
Qmatrix[1,1] = np.sum(dmatrix[:2,1] * dmatrix[:2,1])
Qmatrix[4,4] = Qmatrix[1,1]
Qmatrix[7,7] = Qmatrix[1,1]
Qmatrix[2,2] = np.sum(dmatrix[:2,2] * dmatrix[:2,2])
Qmatrix[5,5] = Qmatrix[2,2]
Qmatrix[8,8] = Qmatrix[2,2]

Qmatrix[0,1] = np.sum(dmatrix[:2,0] * dmatrix[:2,1])
Qmatrix[0,2] = np.sum(dmatrix[:2,0] * dmatrix[:2,2])
Qmatrix[1,2] = np.sum(dmatrix[:2,1] * dmatrix[:2,2])
Qmatrix[3,4] = Qmatrix[0,1]
Qmatrix[3,5] = Qmatrix[0,2]
Qmatrix[4,5] = Qmatrix[1,2]
Qmatrix[6,7] = Qmatrix[0,1]
Qmatrix[6,8] = Qmatrix[0,2]
Qmatrix[7,8] = Qmatrix[1,2]

for i in range(1,9):
    for j in range(i):
        Qmatrix[i,j] = Qmatrix[j,i]
        
print(np.round(Qmatrix,2))
```

    [[0.44 0.29 0.21 0.   0.   0.   0.   0.   0.  ]
     [0.29 0.19 0.14 0.   0.   0.   0.   0.   0.  ]
     [0.21 0.14 0.1  0.   0.   0.   0.   0.   0.  ]
     [0.   0.   0.   0.44 0.29 0.21 0.   0.   0.  ]
     [0.   0.   0.   0.29 0.19 0.14 0.   0.   0.  ]
     [0.   0.   0.   0.21 0.14 0.1  0.   0.   0.  ]
     [0.   0.   0.   0.   0.   0.   0.44 0.29 0.21]
     [0.   0.   0.   0.   0.   0.   0.29 0.19 0.14]
     [0.   0.   0.   0.   0.   0.   0.21 0.14 0.1 ]]



```python
pvector = np.zeros(9)
pvector[0] = (-2)*np.sum(dmatrix[:2,0] * dmatrix[1:3,0])
pvector[1] = (-2)*np.sum(dmatrix[:2,1] * dmatrix[1:3,0])
pvector[2] = (-2)*np.sum(dmatrix[:2,2] * dmatrix[1:3,0])
pvector[3] = (-2)*np.sum(dmatrix[:2,0] * dmatrix[1:3,1])
pvector[4] = (-2)*np.sum(dmatrix[:2,1] * dmatrix[1:3,1])
pvector[5] = (-2)*np.sum(dmatrix[:2,2] * dmatrix[1:3,1])
pvector[6] = (-2)*np.sum(dmatrix[:2,0] * dmatrix[1:3,2])
pvector[7] = (-2)*np.sum(dmatrix[:2,1] * dmatrix[1:3,2])
pvector[8] = (-2)*np.sum(dmatrix[:2,2] * dmatrix[1:3,2])
print(np.round(pvector,2))
```

    [-0.79 -0.53 -0.38 -0.6  -0.39 -0.28 -0.48 -0.31 -0.23]



```python
Amatrix = np.zeros([3,9])
Amatrix[0,:] = np.array([1,0,0]*3)
Amatrix[1,:] = np.array([0,1,0]*3)
Amatrix[2,:] = np.array([0,0,1]*3)
print(Amatrix)
```

    [[1. 0. 0. 1. 0. 0. 1. 0. 0.]
     [0. 1. 0. 0. 1. 0. 0. 1. 0.]
     [0. 0. 1. 0. 0. 1. 0. 0. 1.]]



```python
Q = 2*matrix(Qmatrix)
p = matrix(pvector)
A = matrix(Amatrix)
b = matrix(np.array([1.0,1.0,1.0]),(3,1))
G = matrix(-np.eye(9))
h = matrix(np.zeros(9))
```


```python
sol=solvers.qp(P = Q, q = p,G=G,h=h , A=A, b=b)
```

         pcost       dcost       gap    pres   dres
     0: -6.9062e-01 -3.8739e+00  3e+00  1e-16  3e+00
     1: -6.9429e-01 -7.9617e-01  1e-01  3e-16  9e-02
     2: -6.9672e-01 -7.0763e-01  1e-02  1e-16  9e-03
     3: -7.0046e-01 -7.0240e-01  2e-03  1e-16  9e-04
     4: -7.0096e-01 -7.0148e-01  5e-04  2e-16  2e-04
     5: -7.0125e-01 -7.0128e-01  3e-05  2e-16  1e-06
     6: -7.0126e-01 -7.0126e-01  1e-06  2e-16  2e-08
     7: -7.0126e-01 -7.0126e-01  1e-08  3e-16  2e-10
    Optimal solution found.



```python
print(sol['x'])
```

    [ 1.56e-07]
    [ 1.00e+00]
    [ 5.17e-01]
    [ 4.51e-01]
    [ 8.24e-08]
    [ 4.83e-01]
    [ 5.49e-01]
    [ 1.56e-07]
    [ 1.32e-04]
    



```python
np.round(np.transpose(np.reshape(np.array(sol['x']),(3,3))),3)
```




    array([[0.   , 0.451, 0.549],
           [1.   , 0.   , 0.   ],
           [0.517, 0.483, 0.   ]])




```python

```
