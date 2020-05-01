This is a python module for post-processing images received from NOAA satellites. In particular, this
module will attempt to georeference the images.

Operations that currently seem to be working:

* Histogram stretch (global and adaptive)
* write (can write whole image, left channel only, right channel only, or both)
* Border (can draw border around the area of detected actual sat image)
* Denoise
* Show image (useful for playing with parameters)

Author: @tomaszmrugalski

# Usage

To process an image, you should call process_img and pass two parameters. First is the name
of the PNG file and the second parameter is a structure that explains what exactly you want
to be done with the image. Here's a minimal code that does image extraction and enhacements:

```python

from noaatools.imageproc import process_img
params = {
    "histogram": True,
    "histogram-adaptive": True,
    "border": True,
    "show": False,
    "write": True,
    "write-left": True,
    "write-right": False,
    "denoise": False,
    "georef": False
}

process_img("tests/1276.png", params)
```

# Georeferencing

Currently is being worked on. There's code in noaa_postproc/georef.py, but it doesn't work yet.


Example Cesium code:

```

    var p = Cesium.Cartesian3.fromElements(944166.043, 851659.627, 7099575.419)
    var ssp = Cesium.Cartesian3.fromDegrees(79.904783, 40.435743, 855.126588)

    console.log(p)
    console.log(ssp)

    viewer.entities.add({
    position : p,
    billboard : {
        image : '../images/whiteShapes.png',
        imageSubRegion : new Cesium.BoundingRectangle(49, 43, 18, 18),
        color : Cesium.Color.WHITE
    }
    });

    viewer.entities.add({
    position : ssp,
    billboard : {
        image : '../images/whiteShapes.png',
        imageSubRegion : new Cesium.BoundingRectangle(49, 43, 18, 18),
        color : Cesium.Color.YELLOW
    }
    });

```
