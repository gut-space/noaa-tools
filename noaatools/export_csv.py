

def expand3d(pt):
    """ Expand specified point to be full LLA (latitude, longitude, altitude), even if only lon, lat was specified """
    if (len(pt) == 3):
        return pt
    return pt[0], pt[1], 0.0

def to_text(name, lla, comment):
    lla = expand3d(lla)
    return "%s, %11.7f, %11.7f, %11.7f, %s" % (name, lla[0], lla[1], lla[2], comment)

def export2csv(outfile, satname, aos_ts, los_ts, aos_lla, los_lla,
                  corner_ul, corner_ur, corner_ll, corner_lr, tle1, tle2, text):

    f = open(outfile, "w")

    txt =  "# Data generated using noaatools.\n"
    txt += "# %s\n" % satname
    txt += "# TLE:\n"
    txt += "# %s\n" % tle1
    txt += "# %s\n" % tle2
    txt += "object, longitude [deg], latitude [deg], altitude [km], comment\n"
    txt += to_text("AOS", aos_lla, aos_ts) + "\n"
    txt += to_text("LOS", los_lla, los_ts) + "\n"
    txt += to_text("Upper Left Corner", corner_ul, aos_ts) + "\n"
    txt += to_text("Upper Right Corner", corner_ur, aos_ts) + "\n"
    txt += to_text("Lower Left Corner", corner_ll, los_ts) + "\n"
    txt += to_text("lower Right Corner", corner_lr, los_ts) + "\n"

    f.write(txt)
    f.close()

    print("Georeference data exported in CSV format to %s" % outfile)
