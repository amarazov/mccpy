import vtk
from scipy import *
from numpy.random import multivariate_normal
from vtk.util.colors import tomato, banana, gold, blue,red

from scipy import integrate, linalg

    
def Lorenz(w, t):
    
    S = 10.0; R = 28.0; B = 8.0/3.0

    x, y, z = w
    return array([S*(y-x), R*x-y-x*z, x*y-B*z])

        
def intLorenz(w,time):
    
    S = 10.0; R = 28.0; B = 8.0/3.0

    res=[]
    for w0 in w:
        trajectory = integrate.odeint(Lorenz, w0, time, args=() )
        res.append(trajectory[time.size-1,:])
    
    return res

def intSLorenz(w,time):
    
    trajectory = integrate.odeint(Lorenz, w, time, args=() )
    return trajectory[time.size-1,:]

def getParam(X):
    
    N=X.shape[0]    # num of points;
    mu=sum(X,0)/N   # midpoint of the cloud;
    X=X-mu
    Sigma=dot(X.T,X)/(N-1)  # unbiased estimate;
    
    return mu, Sigma

def getRadius(Sigma):
    
    la,U = linalg.eig(Sigma)
    
    return max(real(la)**0.5)

def DrawTrajectory():
    
    nbrPoints = 1000
    w0 = array([0.0, 1.0, 0.0])
    time = linspace(0.0, 40.0, nbrPoints) 
    trajectory = integrate.odeint(Lorenz, w0, time,args=())
     
    
     
    # Total number of points.
    numberOfInputPoints = nbrPoints
     
    # One spline for each direction.
    aSplineX = vtk.vtkCardinalSpline()
    aSplineY = vtk.vtkCardinalSpline()
    aSplineZ = vtk.vtkCardinalSpline()
     
    inputPoints = vtk.vtkPoints()
    for i in range(0, numberOfInputPoints):
        x = trajectory[i,0]
        y = trajectory[i,1]
        z = trajectory[i,2]
        aSplineX.AddPoint(i, trajectory[i,0])
        aSplineY.AddPoint(i, trajectory[i,1])
        aSplineZ.AddPoint(i, trajectory[i,2])
        inputPoints.InsertPoint(i, x, y, z)
     
    # The following section will create glyphs for the pivot points
    # in order to make the effect of the spline more clear.
     
    # Create a polydata to be glyphed.
    inputData = vtk.vtkPolyData()
    inputData.SetPoints(inputPoints)
     
    # Use sphere as glyph source.
    balls = vtk.vtkSphereSource()
    balls.SetRadius(.01)
    balls.SetPhiResolution(10)
    balls.SetThetaResolution(10)
     
    glyphPoints = vtk.vtkGlyph3D()
    glyphPoints.SetInputData(inputData)
    glyphPoints.SetSourceConnection(balls.GetOutputPort())
     
    glyphMapper = vtk.vtkPolyDataMapper()
    glyphMapper.SetInputConnection(glyphPoints.GetOutputPort())
     
    glyph = vtk.vtkActor()
    glyph.SetMapper(glyphMapper)
    glyph.GetProperty().SetDiffuseColor(tomato)
    glyph.GetProperty().SetSpecular(.3)
    glyph.GetProperty().SetSpecularPower(30)
     
    # Generate the polyline for the spline.
    points = vtk.vtkPoints()
    profileData = vtk.vtkPolyData()
     
    # Number of points on the spline
    numberOfOutputPoints = nbrPoints*10
     
    # Interpolate x, y and z by using the three spline filters and
    # create new points
    for i in range(0, numberOfOutputPoints):
        t = (numberOfInputPoints-1.0)/(numberOfOutputPoints-1.0)*i
        points.InsertPoint(i, aSplineX.Evaluate(t), aSplineY.Evaluate(t),
                           aSplineZ.Evaluate(t))
     
    # Create the polyline.
    lines = vtk.vtkCellArray()
    lines.InsertNextCell(numberOfOutputPoints)
    for i in range(0, numberOfOutputPoints):
        lines.InsertCellPoint(i)
     
    profileData.SetPoints(points)
    profileData.SetLines(lines)
     
    # Add thickness to the resulting line.
    profileTubes = vtk.vtkTubeFilter()
    profileTubes.SetNumberOfSides(8)
    profileTubes.SetInputData(profileData)
    profileTubes.SetRadius(.05)
     
    profileMapper = vtk.vtkPolyDataMapper()
    profileMapper.SetInputConnection(profileTubes.GetOutputPort())
     
    profile = vtk.vtkActor()
    profile.SetMapper(profileMapper)
    profile.GetProperty().SetDiffuseColor(gold)
    profile.GetProperty().SetSpecular(.3)
    profile.GetProperty().SetSpecularPower(30)
    
    return glyph, profile

class vtkTimerCallback():

    def __init__(self,duration,actors,points, ellipsoid, dt, sampling=0):
            self.timer_count = 0
            self.duration = duration
            self.actors = actors
            self.ell=ellipsoid
            self.dt = dt
            self.time=[0]
            self.radius=[]
            mu,S=getParam(points)
            self.radius.append(getRadius(S))
            self.points=points
            self.sampling=sampling
    def ResetActors(self):
        j=0
        for act in self.actors:
            act.SetPosition(self.points[j])
            j+=1;

    def execute(self,obj,event):
        iren = obj
            
        if self.timer_count < self.duration:
            j=0
            for act in self.actors:
                self.points[j,:] = intSLorenz(self.points[j,:],array([0.0, self.dt]))
                act.SetPosition(self.points[j])
                
                j+=1;
            
            mu,Sigma=getParam(self.points)
            if self.sampling:
                self.points = multivariate_normal(mu, Sigma,j)
            
            self.radius.append(getRadius(Sigma))
            self.time.append(self.timer_count*self.dt)
            self.ell.SetPosition(mu)
            iren.CreateTimer(1)
            iren.GetRenderWindow().Render()
            self.timer_count += 1
        else:
            iren.DestroyTimer()
            

