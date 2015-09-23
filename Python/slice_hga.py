import sys
import os

def clearSlices(HGA_file):
    d = os.path.dirname(HGA_file)

    for dirpath, dnames, fnames in os.walk(d):
        for f in fnames:
            if f.startswith("%s." % os.path.basename(HGA_file)):
                print "REMOVING %s/%s" % (d,f)
                os.remove("%s/%s" % (d,f))


def sliceHGA(HGA_file):
    """ Reads in a HGA file, and then splits the file in slices.
        Each slice contains all lattice sites with a given x coordinate.
        The name of each slice file is <HGA_file>.<X_coordinate>
    """

    last_x = ""

    hga_source = open(HGA_file,"r")
    clearSlices(HGA_file)
    for line in hga_source:
        chunks = (line.strip()).split()
        if len(chunks)>2:
            x_coordinate = chunks[1]
            if x_coordinate != last_x:
                slice_file = open("%s.%s" % (HGA_file, x_coordinate), "a")
                last_x = x_coordinate
            slice_file.write(line)
          
            



if __name__ == '__main__':
    HGA_file = sys.argv[1]
    sliceHGA(HGA_file)
