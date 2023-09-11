""" Export georeference data to CSV format """

def expand3d(pt):
    """ Expand specified point to be full LLA (latitude, longitude, altitude), even if only lon, lat was specified """
    if len(pt) == 3:
        return pt
    return pt[0], pt[1], 0.0


def to_text(name, lla, comment):
    """ Convert specified latitude-longitude-altitude point to text """
    lla = expand3d(lla)
    return f"{name}, {lla[0]:11.7f}, {lla[1]:11.7f}, {lla[2]:11.7f}, {comment}"


def export2csv(outfile, satname, aos_ts, los_ts, aos_lla, los_lla,
               corner_ul, corner_ur, corner_ll, corner_lr, tle1, tle2, _):
    """ Exports georeference data to CSV format"""

    with open(outfile, "w", encoding="utf-8") as f:

        txt = "# Data generated using noaatools.\n"
        txt += f"# {satname}\n"
        txt += "# TLE:\n"
        txt += f"# {tle1}\n"
        txt += f"# {tle2}\n"
        txt += "object, latitude [deg], longitude [deg], altitude [km], comment\n"
        txt += to_text("AOS", aos_lla, aos_ts) + "\n"
        txt += to_text("LOS", los_lla, los_ts) + "\n"
        txt += to_text("Upper Left Corner", corner_ul, aos_ts) + "\n"
        txt += to_text("Upper Right Corner", corner_ur, aos_ts) + "\n"
        txt += to_text("Lower Left Corner", corner_ll, los_ts) + "\n"
        txt += to_text("lower Right Corner", corner_lr, los_ts) + "\n"

        f.write(txt)
        f.close()

    print(f"Georeference data exported in CSV format to {outfile}")
