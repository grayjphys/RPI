import pickle
import pandas as pd

#
#(103, 103)
#(103, 103)
#(109, 103)
#(225, 225, 225)
#(225, 225, 225)
#(225, 225, 225)
#
with open("Dens_Mat.pkl",'rb') as f:
     data = pickle.load(f)
     if data.shape != (103,103):
          print("DM has dimensions: ",data.shape)

with open("KS_Mat.pkl",'rb') as f:
     data = pickle.load(f)
     if data.shape != (103,103):
          print("KSM has dimensions: ",data.shape)

with open("States.pkl",'rb') as f:
     data = pickle.load(f)
     if data.shape != (109,103):
          print("States Matrix has dimensions: ",data.shape)

with open("EX.pkl",'rb') as f:
     data = pickle.load(f)
     if data.shape != (225,225,225):
          print("EX has dimensions: ",data.shape)

with open("EY.pkl",'rb') as f:
     data = pickle.load(f)
     if data.shape != (225,225,225):
          print("EY has dimensions: ",data.shape)

with open("EZ.pkl",'rb') as f:
     data = pickle.load(f)
     if data.shape != (225,225,225):
          print("EZ has dimensions: ",data.shape)
