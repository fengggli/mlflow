#### import the simple module from the paraview
# this generate images from hdf5 sources
# paraview /usr/lib/paraview-5.1/site-packages/paraview/simple.py
#import sys
#sys.path.append('/usr/lib/paraview-5.1/site-packages/paraview/simple.py')
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

for timestep in range(1,9):
#timestep=1

#input and output path
#vtipath='data/isotropic_201_201_1_t_'+ str(timestep)+'.vtk'
    vtipath='data/isotropic_201_201_1_t_' + str(timestep) +'.vti'
    imagepath='data/result_images/isotropic_201_201_1_t_'+ str(timestep) + '.png'
    # create a new 'XML Image Data Reader'
    isotropic_201_201_1_v_pvti = XMLImageDataReader(FileName=vtipath)
    isotropic_201_201_1_v_pvti.PointArrayStatus = ['Pressure', 'Velocity']

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    #renderView1.ViewSize = [784, 784]

    # reset view to fit data
    renderView1.ResetCamera()

    # reset view to fit data
    renderView1.ResetCamera()

    # get color transfer function/color map for 'Pressure'
    pressureLUT = GetColorTransferFunction('Pressure')

    # show data in view
    isotropic_201_201_1_v_pvtiDisplay = Show(isotropic_201_201_1_v_pvti, renderView1)
    # trace defaults for the display properties.
    isotropic_201_201_1_v_pvtiDisplay.Representation = 'Slice'
    isotropic_201_201_1_v_pvtiDisplay.ColorArrayName = ['POINTS', 'Pressure']
    isotropic_201_201_1_v_pvtiDisplay.LookupTable = pressureLUT
    isotropic_201_201_1_v_pvtiDisplay.OSPRayScaleArray = 'Pressure'
    isotropic_201_201_1_v_pvtiDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    isotropic_201_201_1_v_pvtiDisplay.GlyphType = 'Arrow'
    isotropic_201_201_1_v_pvtiDisplay.ScalarOpacityUnitDistance = 0.05074636015285223
    isotropic_201_201_1_v_pvtiDisplay.Slice = 100

    # reset view to fit data
    renderView1.ResetCamera()

    #changing interaction mode based on data extents
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [10000.0, 0.6135923, 0.6135923]
    renderView1.CameraFocalPoint = [0.0, 0.6135923, 0.6135923]
    renderView1.CameraViewUp = [0.0, 0.0, 1.0]

    # show color bar/color legend
    #isotropic_201_201_1_v_pvtiDisplay.SetScalarBarVisibility(renderView1, True)

    # get opacity transfer function/opacity map for 'Pressure'
    #pressurePWF = GetOpacityTransferFunction('Pressure')

    # reset view to fit data
    renderView1.ResetCamera()

    # create a new 'Glyph'
    glyph1 = Glyph(Input=isotropic_201_201_1_v_pvti,
        GlyphType='Arrow')
    glyph1.Scalars = ['POINTS', 'None']
    glyph1.Vectors = ['POINTS', 'Velocity']
    glyph1.ScaleFactor = 0.12271846
    glyph1.GlyphTransform = 'Transform2'

    # Properties modified on glyph1
    glyph1.ScaleMode = 'vector'
    glyph1.ScaleFactor = 0.06510780153626579

    # show data in view
    glyph1Display = Show(glyph1, renderView1)
    # trace defaults for the display properties.
    glyph1Display.ColorArrayName = ['POINTS', 'Pressure']
    glyph1Display.LookupTable = pressureLUT
    glyph1Display.OSPRayScaleArray = 'GlyphVector'
    glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    glyph1Display.GlyphType = 'Arrow'
    glyph1Display.SetScaleArray = ['POINTS', 'Pressure']
    glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
    glyph1Display.OpacityArray = ['POINTS', 'Pressure']
    glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'

    # show color bar/color legend
    #glyph1Display.SetScalarBarVisibility(renderView1, True)

    # set scalar coloring
    #ColorBy(glyph1Display, ('POINTS', 'GlyphVector'))

    # rescale color and/or opacity maps used to include current data range
    #glyph1Display.RescaleTransferFunctionToDataRange(True)

    # show color bar/color legend
    #glyph1Display.SetScalarBarVisibility(renderView1, True)

    # get color transfer function/color map for 'GlyphVector'
    #glyphVectorLUT = GetColorTransferFunction('GlyphVector')

    # get opacity transfer function/opacity map for 'GlyphVector'
    #glyphVectorPWF = GetOpacityTransferFunction('GlyphVector')

    # turn off scalar coloring
    ColorBy(glyph1Display, None)

    # Properties modified on glyph1
    glyph1.ScaleFactor = 0.0527689378

    # current camera placement for renderView1
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [-3.3527306774660897, 0.6135923, 0.6135923]
    renderView1.CameraFocalPoint = [0.0, 0.6135923, 0.6135923]
    renderView1.CameraViewUp = [0.0, 0.0, -1.0]
    renderView1.CameraParallelScale = 0.8677505524277008

    # save screenshot
    SaveScreenshot(imagepath, magnification=1, quality=100, view=renderView1)
    print('image is saved in', imagepath

