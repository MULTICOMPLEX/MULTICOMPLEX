import numpy as geek 
  
#Working on 1D 
arr = geek.arange(5) 
print("arr : \n", arr) 
  
repetitions = 2
a = geek.repeat(arr, repetitions) 
print("\nRepeating arr 2 times : \n", a) 
print("Shape : ", a.shape) 
  
repetitions = 3
a = geek.repeat(arr, repetitions) 
print("\nRepeating arr 3 times : \n", a) 
# [0 0 0 ..., 4 4 4] means [0 0 0 1 1 1 2 2 2 3 3 3 4 4 4] 
# since it was long output, so it uses [ ... ] 
print("Shape : ", a.shape) 