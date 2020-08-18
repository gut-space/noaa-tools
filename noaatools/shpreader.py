import shapefile

def read_shp(fname: str, filter_country: str):
    """Returns a list of lists of points"""

    # Read the shp file.
    shp = shapefile.Reader(fname)

    # Get both the shapes and records
    shrecs = shp.shapeRecords()

    countries = []

    # For each pair of shapes and records
    for c in shrecs:
        # check if it's the country We're looking for.
        # c.record has a ton of different records. This is .shp specific. In this particular case (we imported the data from
        # https://www.naturalearthdata.com), the administrative name is in column 8.
        if c.record[8] == filter_country or not filter_country:

            points = []
            # Now iterate over all points.
            for coord in c.shape.points:
                points.append([coord[0], coord[1]])

            countries.append(points)

    return countries

# Here's an example code that uses the function above:
#
# fname = "data/shp/countries3.shp"
# filter = "Poland"
# lines = read_shp(fname, filter)
#
# print("Got %d lines from %s file for %s: First pt is %.3f,%.3f" % (len(lines), fname, filter, lines[0][0], lines[0][1]) )
