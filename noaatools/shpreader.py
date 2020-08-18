import shapefile
import pprint

def read_shp(fname: str, filter_country: str):
    """Returns a list of lines. Each line is defined by [lat1, lon1, lat2, lon2] pair of floats."""

    # Read the shp file.
    shp = shapefile.Reader(fname)

    # Get both the shapes and records
    shrecs = shp.shapeRecords()

    lines = []

    # For each pair of shapes and records
    i = 0
    for c in shrecs:
        # check if it's the country We're looking for.
        # c.record has a ton of different records. This is .shp specific. In this particular case (we imported the data from
        # https://www.naturalearthdata.com), the administrative name is in column 8.
        if c.record[8] == filter_country or not filter_country:

            #print("Processing country %d: %s" % (i, c.record[8]))
            i += 1
            parts = c.shape.parts
            # For example, for Finland we get this: [0, 400, 406, 411, 422, 427, 433, 438, 446, 450, 458, 463]
            # It should be interpreted as 11 polygond. The first contains points from 0 to 399. The second from
            # 400 to 405, etc.
            if not len(c.shape.points):
                continue # Monaco is so small it doesn't have any points.
            prv = c.shape.points[0]
            index = 0
            for coord in c.shape.points:
                if (index not in parts):
                    lines.append([coord[0], coord[1], prv[0], prv[1]])
                prv = coord
                index += 1

    print("Returned %d lines" % len(lines))
    return lines

# Here's an example code that uses the function above:
#
# fname = "data/shp/countries3.shp"
# filter = "Poland"
# lines = read_shp(fname, filter)
#
# print("Got %d lines from %s file for %s: First pt is %.3f,%.3f" % (len(lines), fname, filter, lines[0][0], lines[0][1]) )
