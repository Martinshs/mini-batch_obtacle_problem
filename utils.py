from fenics import *
import imageio     # To make the gif
import os          # To obtain the path of the Python document
import random 


# Needed to have a nested conditional for the symbolic C++ rep.
parameters["form_compiler"]["representation"] = "uflacs"

def smoothmax(r, eps=1e-4):
    """ Smooth approximation of x |--> max(x, 0). 
        Got it from: M. Hintermuller and I. Kopacka. 
        "A smooth penalty approach and a
         nonlinear multigrid algorithm for 
         elliptic MPECs.""
         Computational Optimization and 
         Applications, 50(1):111-145, 2011.
"""
    return conditional(gt(r, eps), r-eps/2, conditional(lt(r, 0), 0, r**2/(2*eps)))


############################
####  GIF GENERATOR   ######
############################

def make_gif(names_figure,gif_name='Simulation_1'):
    
    cwd = os.getcwd()
    newpath_gif = os.path.join(cwd,"gifs") 
    if not os.path.exists(newpath_gif): #We create the a folder called gif
        os.makedirs(newpath_gif)
        
    new_image_names = []
    Frames = 2                      # Increase this parameter to make the animation slower 
    for itera in range(len(names_figure)):
        for frame in range(Frames):
            new_image_names.append(names_figure[itera])    
        
    # We generate the gif with the library imageio
    gif = os.path.abspath(newpath_gif+'/'+gif_name+'.gif')
    with imageio.get_writer(gif, mode='I') as writer:
        for filename in new_image_names:
            image = imageio.imread(filename)
            writer.append_data(image)
            
            
############################
##### UNITY PARTITION ######
############################

            
class unity_partition_1d(UserExpression):
    def __init__(self, n, **kwargs):
        super().__init__(**kwargs) # This part is new!
        self.n = n

    def eval(self, values, x):
        if self.n%2==0 :
            if -1<=x[0]<=-0.1:
                values[0] = 1
            elif -0.1<=x[0]<=0.1:
                values[0] = -(x[0]+0.1)/0.2 +1
            else:
                values[0]=0
        else:
            if 0.1<=x[0]<=1:
                values[0] = 1
            elif -0.1<=x[0]<=0.1:
                values[0] = (x[0]+0.1)/0.2
            else:
                values[0]=0
        #return values[0]
    def value_shape(self):
        return (1,)
    
    
class unity_partition_2d(UserExpression):
    def __init__(self, n, **kwargs):
        super().__init__(**kwargs) # This part is new!
        self.n = n
    def li(self,r):
        if r<= -0.1:
            return 0
        elif -0.1<=r<=0.1:
            return (r+0.1)/0.2
        else:
            return 1     
    def eval(self, values, x):
        if self.n==0:
            values[0] = (1-self.li(x[0]))*(1-self.li(x[1]))
        elif self.n==1:
            values[0] = self.li(x[0])*(1-self.li(x[1]))
        elif self.n==2:
            values[0] = (1-self.li(x[0]))*self.li(x[1])
        else:
            values[0] = self.li(x[0])*self.li(x[1])
    def value_shape(self):
        return ()
    
##################
#### OBSTALES ####
##################
    
class two_mountains(UserExpression):
    def eval(self, value, x):
        cx = -0.5
        cy = 0
        cx2 = 0.5
        cy2 = 0
        deep = -1
        if -((x[0]-cx)**2)*4 - ((x[1]-cy)**2)*4-deep>0:
            value[0] = -((x[0]-cx)**2)*4 - ((x[1]-cy)**2)*4
        elif -((x[0]-cx2)**2)*4 - ((x[1]-cy2)**2)*4-deep>0:
            value[0] = -((x[0]-cx2)**2)*4 - ((x[1]-cy2)**2)*4
        else:
            value[0] = deep
    def value_shape(self):
            return ()


class obs_1d(UserExpression):
    def eval(self, value, x):
        if -1 <= x[0] <= -1/2:
            value[0] = sin(-pi*x[0]-pi/2)-1
        elif 1/2 <= x[0] <= 1:
            value[0] = sin(pi*x[0]-pi/2)-1
        else:
            value[0] = -1
    def value_shape(self):
        return () 
    
##############################
##### BATCHES GENERATOR ######
##############################
    
def batches_gen(len_batch, itera,n_seed,example):
    Batches=[]
    for i in range(itera):
        random.seed(n_seed+i)
        if example=='example_1':
            loc_bat = random.choices([0,1], weights=[1/2,1/2], k=len_batch)
        elif example=='example_2':
            loc_bat = random.choices([0,1,2,3], weights=[1/4,1/4,1/4,1/4], k=len_batch)
        Batches.append(loc_bat)
    return Batches