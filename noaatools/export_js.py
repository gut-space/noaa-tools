# This file contains routines that export data to JavaScript using Cesium

from tletools import TLE
# This is needed to export orbit to CZML format (Cesium)
from poliastro.czml.extract_czml import CZMLExtractor
from astropy import time
from datetime import timezone, timedelta

def cesium_preamble():
    code = """
    // Code generated with noaa-tools.
    // This can be executed using Cesium. The easiest way to test it is to go to
    // https://sandcastle.cesium.com and paste this code there.

    var viewer = new Cesium.Viewer('cesiumContainer');
    var pinBuilder = new Cesium.PinBuilder();
    """

    return code

def export2cesium_point(lla, name, color = "RED"):
    """
    Exports specified LLA coordinates to filename, using name as a label.
    """

    code = """
    var questionPin = viewer.entities.add({
        name : '%s lat=%f,lon=%f',
        position : Cesium.Cartesian3.fromDegrees(%f, %f, %f),
        billboard : {
            image : pinBuilder.fromText('%s', Cesium.Color.%s, 48).toDataURL(),
            verticalOrigin : Cesium.VerticalOrigin.BOTTOM
        }
    });
    """ % (name, lla[0], lla[1], lla[1], lla[0], lla[2], name[0], color)

    return code

def export2cesium_tle(tle1, tle2, satname, aos, los):
    """
    Generates JS Cesium code that shows satellite trajectory between AOS and LOS.

    Parameters
    ==========
    tle1, tle2 - first and second line of TLE
    satname - to be displayed
    aos - Aquisition of Singal in Timedate format
    los - Loss of Signal in Timedate format
    """
    tle = TLE.from_lines(satname, tle1, tle2)
    orb = tle.to_orbit()

    sample_points = 10
    aos.replace(tzinfo=timezone.utc)
    aos += timedelta(hours=0)
    los += timedelta(hours=0)

    aos_astropy = time.Time(aos, scale="utc")
    los_astropy = time.Time(los, scale="utc")

    extractor = CZMLExtractor(aos_astropy, los_astropy, sample_points)
    extractor.add_orbit(orb, path_show=True, path_width=3, path_color=[125, 80, 120, 255], label_text=satname)

    txt =" var czml = [\n"
    for i in extractor.packets:
        txt += repr(i)
        txt += ",\n"
    txt += "];\n"

    # I don't understand what's exactly going on, but I suspect this is a poliastro 0.13.1 bug. When it exports data
    # to CZML, the timezone is messed up. I've fixed very similar problem in Poliastro. The patch for this been accepted
    # and will be released in upcoming 0.14.0.
    # This is how it looks like:  "availability": "2020-04-12T09:01:03Z/2020-04-12T09:17:06Z",
    # This is how it SHOULD look: "availability": "2020-04-12T09:01:03/2020-04-12T09:17:06",

    #txt = txt.replace("Z/", "/") # replace Z/ with /
    #txt = txt.replace('Z"', '"') # replace Z" with "

    txt += "var dataSourcePromise = viewer.dataSources.add(Cesium.CzmlDataSource.load(czml));"

    return txt

def expand3d(pt):
    """ Expand specified point to be full LLA (latitude, longitude, altitude), even if only lon, lat was specified """
    if (len(pt) == 3):
        return pt
    return pt[0], pt[1], 0.0

def export2cesium(outfile, imgfile, aos_ts, los_ts, aos_lla, los_lla,
                  corner_ul, corner_ur, corner_ll, corner_lr, tle1, tle2, text):
    """
    Exports all data to JavaScript that's usable with Cesium.

    Parameters
    ==========
    outfile - output filename to write content to
    imgfile - name of the input PNG file (currently not used yet)
    aos - aquisition of signal, in timedate format
    los - loss of signal, in timedate format
    lla_aos - an array of three coords that specify subsat point in LLA notation (longitude, lattitude, altitude)
    lla_aos - an array of three coords that specify subsat point in LLA notation (longitude, lattitude, altitude)
    tle1 - first line of TLE data
    tle2 - second line of TLE data
    """

    txt =  cesium_preamble()
    txt += export2cesium_tle(tle1, tle2, "satname", aos_ts, los_ts)
    txt += "\n\n"

    txt += export2cesium_point(aos_lla, "AOS(%s)" % text, "BLUE")
    txt += export2cesium_point(los_lla, "LOS(%s)" % text, "YELLOW")

    # Export corners
    if corner_ul:
        corner_ul = expand3d(corner_ul)
        txt += export2cesium_point(corner_ul, "Upper Left", "RED")
    if corner_ur:
        corner_ur = expand3d(corner_ur)
        txt += export2cesium_point(corner_ur, "Upper Right", "GREEN")
    if corner_ll:
        corner_ll = expand3d(corner_ll)
        txt += export2cesium_point(corner_ll, "Lower Left", "RED")
    if corner_lr:
        corner_lr = expand3d(corner_lr)
        txt += export2cesium_point(corner_lr, "Lower Right", "GREEN")

    f = open(outfile, "w")
    f.write(txt)
    f.close()

    print("Georeference data exported to %s" % outfile)
