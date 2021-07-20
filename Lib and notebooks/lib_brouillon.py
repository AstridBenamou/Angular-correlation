# ----------------------Calculate without peridoic distance----------------------
def norme_R3(x,y):
  return (np.sqrt((x[0]-y[:, 0])**2 + (x[1]-y[:, 1])**2 + (x[2]-y[:, 2])**2))

def near_pts(index, df, r):
  """
  Returns all points in a radius r of the center particule [index] 

  :param index: index of the center particule
  :type index: integer
  :param df: dataframe with all coordinates of the particules 
  :type df: pd.DataFrame
  :param r: radius of sphere in with we want to get the particules 
  :type r: float

  :returns: particules in a radius r of the center particule 
  :rtype: pd.DataFrame

  """
  x,y,z = df.iloc[index,]
  df = df.drop(index, axis = 0)
  neighbours = df[(abs(df['x']- x) <= r) & (abs(df['y']- y) <= r) & (abs(df['z']- z) <= r)]
  return (neighbours)

def near_pts0(df, r):
  """
  Returns all points in a radius r of the center of the box

  :param df: dataframe with all coordinates of the particules 
  :type df: pd.DataFrame
  :param r: radius of sphere in with we want to get the particules 
  :type r: float

  :returns: particules in a radius r of the center of the box 
  :rtype: pd.DataFrame

  """
  neighbours = df[(abs(df['x'] - 0.5) <= r) & (abs(df['y'] - 0.5) <= r) & (abs(df['z'] - 0.5) <= r)]
  return (neighbours)

def pts_at_r(indice, df, r, dr):
  """
  Returns number of points at a distance between r and r+dr from the center particule [indice] 

  :param index: index of the center particule
  :type index: integer
  :param df: dataframe with all coordinates of the particules 
  :type df: pd.DataFrame
  :param r: radius of sphere in with we want to get the particules 
  :type r: float
  :param dr: precision of the distance 
  :type dr: float 

  :returns: particules in a radius r of the center particule 
  :rtype: pd.DataFrame

  """
  x,y,z = df.iloc[indice,]
  df = df.drop(indice, axis = 0)
  n1 = df[(abs(df['x']- x-r) <= dr) & (abs(df['y']- y - r) <= dr) & (abs(df['z']- z - r) <= dr)]
  n2 = df[(abs(df['x']- x + r) <= dr) & (abs(df['y']- y + r) <= dr) & (abs(df['z']- z + r) <= dr)]
  return (pd.concate(n1,n2))

# ----------------------Calculate with peridoc distance----------------------
def dist_period(x1,x2):

  d1 = (x1-x2)**2
  d2 = (1 - x2 + x1)**2
  d3 = (1 + x2 - x1)**2
  return (min(d1,d2,d3))

def norm_period(x,y):

  m = np.array([[0,0,0], [0,0,1], [0,1,0], [1,0,0], [1,1,0], [1,0,1], [0,1,1], [1,1,1], [0,0,-1], [0,-1,0], [-1,0,0], [-1,-1,0], [-1,0,1], [0,-1,-1], [-1,-1,-1], [-1,1,0], [-1,0,1], [1,-1,0], [1,0,-1], [0,1,-1], [0,-1,1], [-1,-1,1], [-1,1,-1], [1,-1,-1], [1,1,-1], [1,-1,1], [-1,1,1]])
  d = norme_R3(x, y + m)
  res = np.min(d)
  return (res)

def repartition_pts(data,r_max,delta):
  bins = np.zeros(math.floor(r_max/delta)*2)
  for i in range(len(data)):
    for j in range(i+1,len(data)): 
      d = norm_period(data[j],data[i])
      bins[math.floor(d/delta)] += 1
  return(bins)

def near_pts_period(index, df, r):
  """Returns all points in a radius r of the center particule [index] considering the periodic distance.

  :param index: index of the center particule
  :type index: integer
  :param df: matrix nx3 with all coordinates of the particules 
  :type df: numpy.array
  :param r: radius of sphere in with we want to get the particules 
  :type r: float

  :returns: particules in a radius r of the center particule considering the box is repeated on each side
  :rtype: pd.DataFrame

  """  
  c = data[index]
  data = np.delete(data, index, 0)
  df["r3"]=df[['x','y','z']].apply(lambda x : norme_R3_period(c,x),axis=1)
  df=df.loc[df["r3"]<=r,["x","y","z"]]  
  return (df)

def pts_at_r_period(indice, df, r, dr):
  """Returns number of points at a distance between r and r+dr from the center particule [indice] periodicaly (box repeats itself on each side)

  :param index: index of the center particule
  :type index: integer
  :param df: dataframe with all coordinates of the particules 
  :type df: pd.DataFrame
  :param r: radius of sphere in with we want to get the particules 
  :type r: float
  :param dr: precision of the distance 
  :type dr: float 

  :returns: particules in a radius r of the center particule 
  :rtype: pd.DataFrame

  """
  c = df.iloc[indice,]
  df = df.drop(indice, axis = 0)
  df["r3"]=df[['x','y','z']].apply(lambda x : norme_R3_period(c,x),axis=1)
  df=df.loc[(df["r3"]>=r) & ((df["r3"]<=r+dr)),["x","y","z"]] 
  return (df)



def correlation(df, dr, range_max):
  bins = np.zeros(math.floor(range_max/dr)+1) 
  for i in range(len(df)):
    for j in range(len(bins)):
      bins[j] += len(pts_at_r_period(i,df,j*dr,dr))  
  bins = np.divide(bins,sum(bins)*np.ones(len(bins)))
  shell = (np.arange(1,len(bins)+1))**2*dr*np.pi
  #bins = np.divide(bins, shell)

  return bins

def all_distances(df):
    list_dist = []
    for i in range(len(df)):
        for j in range(i, len(df)):
            list_dist.append(norme_R3_period(df.iloc[i],df.iloc[j]))
    return (list_dist)

# ---------------------- Temporary Functions ----------------------
#list of all the distances between pairs of points in df 
def distances(df):
  dist = []
  for i in range(len(df)):
    for j in range(len(df)): 
      d = norme_R3(df.iloc[j],df.iloc[i])
      dist.append(d)
  return(dist)



def Corr(rad,center,df):
  neighbours = near_pts(center, df, rad) #dataframe with points in sphere
  distance_list = distances(neighbours) #list of distances between all pairs of particules in sphere 
  distance_list = np.around(distance_list,2)
  uniq, counts = np.unique(distance_list, return_counts=True)
  vol_sphere = np.ones(len(counts))*(np.pi*(rad**3))
  return ([uniq, counts, np.divide(counts,vol_sphere)])

def Corr_0(rad, df):
  neighbours = CountPtsNear0(df, rad) #dataframe with points in sphere
  distance_list = distances(neighbours) #list of distances between all pairs of particules in sphere 
  distance_list = np.around(distance_list,3)
  uniq, counts = np.unique(distance_list, return_counts=True)
  vol_sphere = np.ones(len(counts))*(4*np.pi*(rad**3))/4
  return ([uniq, counts, np.divide(counts,vol_sphere)])

#list of all the distances between pairs of points in df 
def distances_period(df):
  dist = []
  for i in range(len(df)):
    for j in range(len(df)): 
      d = np.sqrt(dist_period(df.iloc[j,0],df.iloc[i,0])+dist_period(df.iloc[j,1],df.iloc[i,1])+dist_period(df.iloc[j,2],df.iloc[i,2]))
      dist.append(d)
  return(dist)

def Corr_0(rad, df, N):
  neighbours = CountPtsNear0(df, rad) #dataframe with points in sphere
  distance_list = distances_period(neighbours) #list of distances between all pairs of particules in sphere 
  distance_list = np.around(distance_list,2) #dr = 0.01
  uniq, counts = np.unique(distance_list, return_counts=True)
  counts = np.divide(counts, N*np.ones(len(counts)) )
  vol_sphere = np.ones(len(counts))*(4*np.pi*(rad**3))/3
  return ([uniq, counts, np.divide(counts,vol_sphere)])


#Correlation of gaussian distribution using method from overleaf