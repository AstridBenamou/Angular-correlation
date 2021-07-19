import numpy as np
import math
import scipy as sp
from scipy import stats
import matplotlib.pyplot as plt
import scipy.stats as stats
import plotly.express as px
import plotly.graph_objects as go
from scipy.spatial import cKDTree


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

def make_bins(data, r_max, dr):
  """
  Returns array of bins containing for each the number of points seperated by a distance between r and dr. 

  :param data: array of shape = [n,3] containing coordinates of points in box
  :type data: numpy array 
  :param r_max: distance max to search for points separation 
  :type r_max: float 
  :param dr: size of bins (thickness of layers in which to search for points in each bins)
  :type dr: float 

  :returns: Bins chart with for each i the number of points separated by a distance between i*dr and (i+1)*dr
  :rtype: numpy array 

  """
  new_data = make_big_box(data) 
  tree = cKDTree(new_data)
  b = 2*np.ones(len(data))
  bins = np.empty(math.floor(r_max/dr))
  for i in range(math.floor(r_max/dr)):
      a = tree.query_ball_point(data, (i+1)*dr, return_length = True)
      bins[i] = sum(a) - sum(b)
      b = a
  return(bins)

def normalize(bins, r_max, dr, n):
  """
  Returns 
  """
  volumes = (np.arange(math.floor(r_max/dr))*dr)**2*4*np.pi*dr
  volumes[0] = 4*np.pi*(dr**3)/3
  return (np.divide(bins, volumes*n))


def pair_correlation(data:np.array, r_max:np.float, dr:np.float)->list:
  """
  Returns the pair correlation of a sample of points in a box.

  :param data: array of shape = [n,3] containing coordinates of points in box
  :type data: numpy array 
  :param r_max: distance max to search for points separation 
  :type r_max: float 
  :param dr: size of bins (thickness of layers in which to search for points in each bins)
  :type dr: float 

  Returns:
      list: [numpy array of the different radius considered, numpy array with the normalized bins]
  """
  bins = make_bins(data, r_max, dr)
  bins = normalize(bins, r_max, dr, len(data))
  r = np.arange(0, r_max, dr)
  return([r, bins])

def trace_pair_correlation(abscisse, correlation, distribution_name, background, chose_color, sample_size):
  """
  Makes figure from the result of pair_correlation function. 

  :param abscisse: first return of pair_correlation function, used as x-axis for our figure.
  :type abscisse: numpy.array
  :param correlation: second return of pair_correlation function, used as y-axis for our figure.
  :type correlation: numpy.array
  :param distribution_name: Distribution that will be indicated in the title of the figure.
  :type distribution_name: string
  :param background: color of the background of the figure
  :type background: string - css color
  :param chose_color: color of the dots and axis
  :type chose_color: string - css color
  :param sample_size: size of the initial data
  :type sample_size: int

  :returns: figure 
  :rtype: plotly.graph_objects.Figure
  """
  layout={
  "title":"Pair correlation function of a " + distribution_name + " distribution",
  "xaxis_title":"radius",
  "yaxis_title":"correlation for a sample of "+ str(sample_size) + " variables",
  "font":dict(
      family="Courier New, monospace",
      size=14,
      color= chose_color
  ),
  "plot_bgcolor":background
  }
  
  fig = go.Figure(
      data=[go.Scatter(
          x = abscisse,
          y = correlation, 
          mode='markers',
          marker_color = chose_color
      )],
      layout=layout
      
  )    
  
  return(fig)

def arc_length(point_1, point_2, r):
  """Returns arc length between two points on a sphere of radius r.
  
  :param point_1: two angular coordinates of the point on the sphere 
  :type point_1: numpy array 
  :param point_2: two angular coordinates of the point on the sphere
  :type point_2: numpy array
  :param r: radius of the sphere 
  :type r: float 

  :returns: arc length between the two points
  :rtype: float
  """
  delta_sig = np.arccos(np.sin(point_1[1])*np.sin(point_2[1])+np.cos(point_1[1])*np.cos(point_2[1])*np.cos(abs(point_1[2]-point_2[2])))
  return(delta_sig*r)


def angular_pair_correlation(data, r, d):
  """
  Returns number of pairs of points separated by an arc length of r for each r 

  :param data: angular coordinates of each point on sphere 
  :type: numpy array of shape (n,2)
  :param r: radius of the sphere on which the points are 
  :type: float 
  :param d: interval between two values of bins we want to try
  :type: float

  :returns: array of bins in which we put the number of pairs of point separated by a certain arc length. 
  :rtype: numpy array 
  """
  bins = np.zeros(math.floor((np.pi*r)/d)+1)
  for i in range(len(data)):
      for j in range(i+1, len(data)):
          l = arc_length(data[i],data[j], r)
          bins[math.floor(l/d)] += 1
  #we need to normalize  divide each bins by the perimeter of the circle of radius r 
  return(bins)

