This is a python module for post-processing images received from NOAA satellites. In particular, this
module will attempt to georeference the images.

Operations that currently seem to be working:

* Histogram stretch (global and adaptive)
* Write to a file (can write whole image, left channel only, right channel only, or both)
* Border (can draw border around the area of detected actual sat image)
* Denoise
* Show image (useful for playing with parameters)
* Georeferencing (experimental)

Author: [Tomek Mrugalski](https://github.com/tomaszmrugalski/)

## Installation

Do the standard python thing:

```
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Usage (image processing)

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

## Usage (georeferencing)

Currently is being worked on. There's code in noaa_postproc/georef.py, but it doesn't work yet.

Example code:

```
from noaatools import georef

# Specify the parameters
tle1 = '1 28654U 05018A   20098.54037539  .00000075  00000-0  65128-4 0  9992'
tle2 = '2 28654  99.0522 154.2797 0015184  73.2195 287.0641 14.12501077766909'

aos = '2020-04-12 09:01:03.063476'
los = '2020-04-12 09:17:06.466954'

filename = '1276.png'
method = georef.Method.PYMAP3D

# Call the routine to calculate the image corners.
d1, d2, aos_lla, los_lla, corner_ul, corner_ur, corner_ll, corner_lr = georef.georef(method, tle1, tle2, aos, los)

# You can then export the data:
export_js.export2cesium(outfile_js, satname, d1, d2, aos_lla, los_lla, corner_ul, corner_ur, corner_ll, corner_lr, tle1, tle2, method.name)
    export_csv.export2csv(outfile_csv, satname, d1, d2, aos_lla, los_lla, corner_ul, corner_ur, corner_ll, corner_lr, tle1, tle2, method.name)
```

## Using from the command line:

This currently supports georeferencing only.

```bash
python -m noaatools.noaatools --help

usage: noaatools [-h] [--aos AOS] [--los LOS] [--tle1 TLE1] [--tle2 TLE2]
                 [--method METHOD] [--satname SATNAME] [--file FILE]

optional arguments:
  -h, --help         show this help message and exit
  --aos AOS          Specify Acquision of Signal as timestamp (default:
                     2020-04-12Z09:01:03.063476)
  --los LOS          Specify Loss of signal as timestamp (default:
                     2020-04-12Z09:17:06.466954)
  --tle1 TLE1        First line of TLE (default: 1 28654U 05018A
                     20098.54037539 .00000075 00000-0 65128-4 0 9992)
  --tle2 TLE2        Second line of TLE (default: 2 28654 99.0522 154.2797
                     0015184 73.2195 287.0641 14.12501077766909)
  --method METHOD    Specify orbital position calculation method (SPHERICAL,
                     OBLATE, PYMAP3D, default: PYMAP3D)
  --satname SATNAME  Specifies satellite name (default: NOAA 18)
  --file FILE        Specifies output file pattern (default: data/1276.png)
```

The tool will export data to files similar to the specified. When you specify 1276.png, the actual output files
will be 1276.js and 1276.txt. By default this tool uses [satnogs observation 1276](https://satnogs.klub.com.pl/obs/1276)
as sample data. That way you can easily run it without any parameters to quickly see how it works:

```bash
python -m noaatools.noaatools
```

Currently there are three calculation method implemented:

- spherical (the simplest and least accurate method. It assumes Earth is completely round)
- oblate (this accounts for Earth oblateness, i.e. the bulge around equator)
- pymap3d (this method is using [PYMAP3D](https://pypi.org/project/pymap3d/)'s eci2ecef and ecef2geodetic functions)

Currently there are two export methods implemented:

- JavaScript - this generates a JavaScript code that is using Cesium. You can see the results easily by pasting the generated
  file on https://sandcastle.cesium.com
- CSV - Comma Separated File. This is a very flexible text format. You can import this data into OpenOffice, Microsoft Excel
  or many other tools and environments

## Running unit-tests

The unit-tests coverage is miserable so far. Nevertheless, you may do the following:

```
pip install pytest
python -m pytest -s -v
```

During the first run, an image for satnogs observation 1276 (about 3MB) will be downloaded. It's a
reasonably recent, good quality observation. I've decided to keep the repo clean and not pollute it
with unnecessary images.

## Useful links

* GDAL + python: https://pcjericks.github.io/py-gdalogr-cookbook/
* Georeferencing in qgis: https://www.qgistutorials.com/en/docs/georeferencing_basics.html
* pymap3d has tons of useful coordinate transforms: https://geospace-code.github.io/pymap3d/index.html
* wxtoimg man page: https://wxtoimgrestored.xyz/support/wxtoimg.pdf
* Cesium support for displaying a single image is confusing. https://www.cesium.com/docs/cesiumjs-ref-doc/SingleTileImageryProvider.html
* PROJ4 (recommended by @fivitti): https://github.com/pyproj4/pyproj
* https://github.com/Alexander-Barth/APTDecoder.jl
