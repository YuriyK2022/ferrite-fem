# script-version: 2.0
# Catalyst state generated using paraview version 5.13.0
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 13

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1458, 550]
renderView1.AxesGrid = 'Grid Axes 3D Actor'
renderView1.CenterOfRotation = [0.5, 0.4891558394429695, 0.5003691871547652]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [-2.9553435395696757, 2.2083543744578202, -0.9850191975737316]
renderView1.CameraFocalPoint = [-0.2692311917812466, 0.5549176561502385, 1.408907585410504]
renderView1.CameraViewUp = [0.3452565551438977, 0.9074679957936758, 0.239373657156292]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 1.0248691033587511
renderView1.LegendGrid = 'Legend Grid Actor'
renderView1.PolarGrid = 'Polar Grid Actor'
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1458, 550)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'XML Unstructured Grid Reader'
hyperelasticityvtu = XMLUnstructuredGridReader(registrationName='hyperelasticity.vtu', FileName=['C:\\Users\\yuriy\\myJulia_proj\\ferrite\\1.2 Hyperelasticity\\hyperelasticity.vtu'])
hyperelasticityvtu.PointArrayStatus = ['u']
hyperelasticityvtu.TimeArray = 'None'

# create a new 'Warp By Vector'
warpByVector1 = WarpByVector(registrationName='WarpByVector1', Input=hyperelasticityvtu)
warpByVector1.Vectors = ['POINTS', 'u']
warpByVector1.ScaleFactor = 0.4522245412204498

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from warpByVector1
warpByVector1Display = Show(warpByVector1, renderView1, 'UnstructuredGridRepresentation')

# get 2D transfer function for 'u'
uTF2D = GetTransferFunction2D('u')

# get color transfer function/color map for 'u'
uLUT = GetColorTransferFunction('u')
uLUT.TransferFunction2D = uTF2D
uLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 0.17677669529663687, 0.865003, 0.865003, 0.865003, 0.35355339059327373, 0.705882, 0.0156863, 0.14902]
uLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'u'
uPWF = GetOpacityTransferFunction('u')
uPWF.Points = [0.0, 0.0, 0.5, 0.0, 0.35355339059327373, 1.0, 0.5, 0.0]
uPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
warpByVector1Display.Representation = 'Surface'
warpByVector1Display.ColorArrayName = ['POINTS', 'u']
warpByVector1Display.LookupTable = uLUT
warpByVector1Display.SelectNormalArray = 'None'
warpByVector1Display.SelectTangentArray = 'None'
warpByVector1Display.SelectTCoordArray = 'None'
warpByVector1Display.TextureTransform = 'Transform2'
warpByVector1Display.OSPRayScaleArray = 'u'
warpByVector1Display.OSPRayScaleFunction = 'Piecewise Function'
warpByVector1Display.Assembly = ''
warpByVector1Display.SelectedBlockSelectors = ['']
warpByVector1Display.SelectOrientationVectors = 'None'
warpByVector1Display.ScaleFactor = 0.12660756815739588
warpByVector1Display.SelectScaleArray = 'None'
warpByVector1Display.GlyphType = 'Arrow'
warpByVector1Display.GlyphTableIndexArray = 'None'
warpByVector1Display.GaussianRadius = 0.006330378407869794
warpByVector1Display.SetScaleArray = ['POINTS', 'u']
warpByVector1Display.ScaleTransferFunction = 'Piecewise Function'
warpByVector1Display.OpacityArray = ['POINTS', 'u']
warpByVector1Display.OpacityTransferFunction = 'Piecewise Function'
warpByVector1Display.DataAxesGrid = 'Grid Axes Representation'
warpByVector1Display.PolarAxes = 'Polar Axes Representation'
warpByVector1Display.ScalarOpacityFunction = uPWF
warpByVector1Display.ScalarOpacityUnitDistance = 0.11280144063101548
warpByVector1Display.OpacityArrayName = ['POINTS', 'u']
warpByVector1Display.SelectInputVectors = ['POINTS', 'u']
warpByVector1Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
warpByVector1Display.ScaleTransferFunction.Points = [-0.030505934463539334, 0.0, 0.5, 0.0, 0.01894796859039277, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
warpByVector1Display.OpacityTransferFunction.Points = [-0.030505934463539334, 0.0, 0.5, 0.0, 0.01894796859039277, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for uLUT in view renderView1
uLUTColorBar = GetScalarBar(uLUT, renderView1)
uLUTColorBar.WindowLocation = 'Any Location'
uLUTColorBar.Position = [0.020576131687242746, 0.6327272727272727]
uLUTColorBar.Title = 'u'
uLUTColorBar.ComponentTitle = 'Magnitude'
uLUTColorBar.ScalarBarLength = 0.33000000000000007

# set color bar visibility
uLUTColorBar.Visibility = 1

# show color legend
warpByVector1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity maps used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup animation scene, tracks and keyframes
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get animation track
warpByVector1ScaleFactorTrack = GetAnimationTrack('ScaleFactor', index=0, proxy=warpByVector1)

# create a new key frame
keyFrame10332 = CompositeKeyFrame()

# create a new key frame
keyFrame10333 = CompositeKeyFrame()
keyFrame10333.KeyTime = 1.0
keyFrame10333.KeyValues = [0.4522245412204498]

# initialize the animation scene
warpByVector1ScaleFactorTrack.KeyFrames = [keyFrame10332, keyFrame10333]

# initialize the animation scene
keyFrame10333.KeyTime = 1.0
keyFrame10333.KeyValues = [0.4522245412204498]

# get time animation track
timeAnimationCue1 = GetTimeTrack()

# initialize the animation scene

# initialize the animation scene

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# initialize the timekeeper

# initialize the animation track

# initialize the animation track
warpByVector1ScaleFactorTrack.KeyFrames = [keyFrame10332, keyFrame10333]

# get animation scene
animationScene1 = GetAnimationScene()

# initialize the animation scene
animationScene1.ViewModules = renderView1
animationScene1.Cues = [timeAnimationCue1, warpByVector1ScaleFactorTrack]
animationScene1.AnimationTime = 1.0
animationScene1.NumberOfFrames = 150

# ----------------------------------------------------------------
# restore active source
SetActiveSource(warpByVector1)
# ----------------------------------------------------------------

# ------------------------------------------------------------------------------
# Catalyst options
from paraview import catalyst
options = catalyst.Options()
options.GlobalTrigger = 'Time Step'
options.CatalystLiveTrigger = 'Time Step'

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    # Code for non in-situ environments; if executing in post-processing
    # i.e. non-Catalyst mode, let's generate extracts using Catalyst options
    SaveExtractsUsingCatalystOptions(options)
