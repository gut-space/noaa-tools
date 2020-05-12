from tletools import TLE
# This is needed to export orbit to CZML format (Cesium)
from poliastro.czml.extract_czml import CZMLExtractor
from astropy import time

def cesium_preamble():
    code = """
    // Code generated with noaa-tools.
    // This can be executed using Cesium. The easiest way to test it is to go to
    // https://sandcastle.cesium.com and paste this code there.

    var viewer = new Cesium.Viewer('cesiumContainer');
    var pinBuilder = new Cesium.PinBuilder();
    """

    return code

def export2cesium_point(lla, name):
    """
    Exports specified LLA coordinates to filename, using name as a label.
    """

    code = """
    var questionPin = viewer.entities.add({
        name : '%s',
        position : Cesium.Cartesian3.fromDegrees(%f, %f, %f),
        billboard : {
            image : pinBuilder.fromText('%s', Cesium.Color.RED, 48).toDataURL(),
            verticalOrigin : Cesium.VerticalOrigin.BOTTOM
        }
    });
    """ % (name, lla[0], lla[1], lla[2], name[0])

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

    txt = txt.replace("Z/", "/") # replace Z/ with /
    txt = txt.replace('Z"', '"') # replace Z" with "

    txt += "var dataSourcePromise = viewer.dataSources.add(Cesium.CzmlDataSource.load(czml));"

    return txt

def export2cesium(outfile, imgfile, aos, los, aos_list, los_list, methods, tle1, tle2):
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

    if len(methods) != len(aos_list) or len(methods) != len(los_list):
        raise Exception("Incorrect parameter size: len(methods) != len(aos_list) != len(los_list)")

    txt =  cesium_preamble()
    txt += export2cesium_tle(tle1, tle2, "satname", aos, los)
    txt += "\n\n"

    for i in range(0, len(methods)):
        txt += " // points %d out of %d" % (i, len(methods))
        txt += export2cesium_point(aos_list[i], "AOS:" + str(aos) + ", method " + methods[i])
        txt += export2cesium_point(los_list[i], "LOS:" + str(los) + ", method " + methods[i])

    f = open(outfile, "w")
    f.write(txt)
    f.close()

    print("Georeference data exported to %s" % outfile)