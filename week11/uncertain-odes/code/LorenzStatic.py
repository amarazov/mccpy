#!/usr/bin/env python
 
import vtk
from scipy import *
from scipy import integrate
from pylab import *
from vtk.util.colors import tomato, purple, blue
 
nbrPoints = 5000
 
def Lorenz(w, t, S, R, B):
    x, y, z = w
    return array([S*(y-x), R*x-y-x*z, x*y-B*z])
 
w0 = array([0.0, 1.0, 0.0])
time = linspace(0.0, 50.0, nbrPoints)
S = 10.0; R = 28.0; B = 8.0/3.0
 
# Conditions initiales
y0=array([-7.5, -3.6, 30.0])
 
# On resous le tout avec odeint de scipy qui est
# lui-meme un wrap de ODEPACK en fortran
trajectory = integrate.odeint(Lorenz, w0, time, args=(S, R, B))
 
# This will be used later to get random numbers.
math = vtk.vtkMath()
 
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
glyphPoints.SetInput(inputData)
glyphPoints.SetSource(balls.GetOutput())
 
glyphMapper = vtk.vtkPolyDataMapper()
glyphMapper.SetInputConnection(glyphPoints.GetOutputPort())
 
glyph = vtk.vtkActor()
glyph.SetMapper(glyphMapper)
glyph.GetProperty().SetDiffuseColor(purple)
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
profileTubes.SetInput(profileData)
profileTubes.SetRadius(.05)
 
profileMapper = vtk.vtkPolyDataMapper()
profileMapper.SetInputConnection(profileTubes.GetOutputPort())
 
profile = vtk.vtkActor()
profile.SetMapper(profileMapper)
profile.GetProperty().SetDiffuseColor(blue)
profile.GetProperty().SetSpecular(.3)
profile.GetProperty().SetSpecularPower(30)
 
# Now create the RenderWindow, Renderer and Interactor
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
 
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
 
# Add the actors
ren.AddActor(glyph)
ren.AddActor(profile)
ren.SetBackground(1,1,1)

renWin.SetSize(500, 500)
 
iren.Initialize()
renWin.Render()
iren.Start()

_win2imgFilter = vtk.vtkWindowToImageFilter()
_win2imgFilter.SetInput(renWin)
  
# save the image to file
_outWriter = vtk.vtkPNGWriter()
_outWriter.SetInput(_win2imgFilter.GetOutput())
_outWriter.SetFileName("Lorenz.png")
_outWriter.Write()