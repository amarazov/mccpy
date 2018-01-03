import numpy

from LorenzModule import *


dt=0.1
sigma=0.1
x0=array([0,1,0])
sampling=1


withEllipsoid=0

N=100



w=(sigma**2)*randn(N,3)+x0
# create a rendering window and renderer
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)

# create a renderwindowinteractor
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

#newPoints(w)
sourceEllipsoid = vtk.vtkSphereSource()
sourceEllipsoid.SetRadius(10*sigma)
sourceEllipsoid.SetThetaResolution(32)
sourceEllipsoid.SetPhiResolution(32)

mapperEllipsoid = vtk.vtkPolyDataMapper()
mapperEllipsoid.SetInputConnection(sourceEllipsoid.GetOutputPort())
actorEllipsoid = vtk.vtkActor()
actorEllipsoid.SetMapper(mapperEllipsoid)

actorEllipsoid.GetProperty().SetColor(red)
if withEllipsoid:
    actorEllipsoid.GetProperty().SetOpacity(0.5)
else:
    actorEllipsoid.GetProperty().SetOpacity(0.0)

ren.AddActor(actorEllipsoid)


source=[]
mapper=[]
actor=[]
# create source
for j in range(N):
    source.append(vtk.vtkSphereSource())
    source[j].SetCenter(w[j,:])
    source[j].SetRadius(.4)
  
    # mapper
    mapper.append(vtk.vtkPolyDataMapper())
    mapper[j].SetInputConnection(source[j].GetOutputPort())
  
    # actor
    actor.append(vtk.vtkActor())
    actor[j].SetMapper(mapper[j])
  
    # assign actor to the renderer
    ren.AddActor(actor[j])

# Add the actors

glyph, profile = DrawTrajectory()

ren.AddActor(glyph)
ren.AddActor(profile)
 
renWin.SetSize(1366, 600)
renWin.SetFullScreen(1)

# enable user interface interactor
iren.Initialize()
renWin.Render()

cb = vtkTimerCallback(duration=10000, actors=actor, points=w, ellipsoid = actorEllipsoid, 
                      dt=dt,sampling=sampling) 

iren.AddObserver('TimerEvent', cb.execute)

iren.CreateTimer(0)

iren.Start()

resFile=file('../RadiusResults.txt','a')

numpy.savetxt(resFile,(cb.time,cb.radius))

resFile.close()
