This is a python module for post-processing images received from NOAA satellites. In particular, this
module will attempt to georeference the images.

Operations that currently seem to be working:

* Histogram stretch (global and adaptive)
* Write to a file (can write whole image, left channel only, right channel only, or both)
* Border (can draw border around the area of detected actual sat image)
* Denoise
* Show image (useful for playing with parameters)

Operation under development:

* Georeferencing

Author: [Tomek Mrugalski](https://github.com/tomaszmrugalski/)

# Installation

Do the standard python thing:

```
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

# Running unit-tests

The unit-tests coverage is miserable so far. Nevertheless, you may do the following:

```
pip install pytest
python -m pytest -s -v
```

During the first run, an image for satnogs observation 1276 (about 3MB) will be downloaded. It's a
reasonably recent, good quality observation. I've decided to keep the repo clean and not pollute it
with unnecessary images.

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

Example code:

```
from noaatools.georef import georef

tle1 = '1 28654U 05018A   20098.54037539  .00000075  00000-0  65128-4 0  9992'
tle2 = '2 28654  99.0522 154.2797 0015184  73.2195 287.0641 14.12501077766909'

aos = '2020-04-12 09:01:03.063476'
los = '2020-04-12 09:17:06.466954'

filename = '1276.png'

georef(filename, tle1, tle2, aos, los)

# The output Cesium code will be generated to 1276.js
```

# Useful links

* GDAL + python: https://pcjericks.github.io/py-gdalogr-cookbook/
* Georeferencing in qgis: https://www.qgistutorials.com/en/docs/georeferencing_basics.html
* pymap3d has tons of useful coordinate transforms: https://geospace-code.github.io/pymap3d/index.html
* wxtoimg man page: https://wxtoimgrestored.xyz/support/wxtoimg.pdf
* Cesium support for displaying a single image is confusing. https://www.cesium.com/docs/cesiumjs-ref-doc/SingleTileImageryProvider.html
