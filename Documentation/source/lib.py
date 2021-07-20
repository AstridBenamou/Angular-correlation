import numpy as np



"""Usefull functions """

def make_big_box(data):
  """
  Returns data extended to 0.5 on every side and repeats points inside 

  :param data: array of shape = [n,3] for n points in box 
  :type data: numpy array

  :return: bigger box 
  :rtype: numpy array of size [m,3] (m > n)

  """
  new_data= data
  m = np.array([[0,0,0], [0,0,1], [0,1,0], [1,0,0], [1,1,0], [1,0,1], [0,1,1], [1,1,1], [0,0,-1], [0,-1,0], [-1,0,0], [-1,-1,0], [-1,0,1], [0,-1,-1], [-1,-1,-1], [-1,1,0], [-1,0,1], [1,-1,0], [1,0,-1], [0,1,-1], [0,-1,1], [-1,-1,1], [-1,1,-1], [1,-1,-1], [1,1,-1], [1,-1,1], [-1,1,1]])
  for i in range(len(m)):
      new = (data) + m[i]
      new = np.delete(new, np.where((new > 1.5) | (new < -0.5))[0], axis=0)
      new_data = np.concatenate((new_data, new), axis=0) 
  return (new_data)


