#  Title:  coverage_summarize.py
#  Author: Tony Di Sera, Gabor Marth Lab, iobio project
#  Written in July 2015
#
#  DESCRIPTION
#  -----------
#  A simple python that takes the output from samtools mpileup and
#  produces the following output records:
#   specific_points
#   1.  Filter the records, keeping those that are the specific positions passed
#       in as an argument
#   reduced_points.
#   2.  Reduce the records down to the maximum points, passed in an argument.
#       To reduce the pont data, the app calculate the mean depth over
#       n positions.
#
#   ARGUMENTS
#   ---------
#   This app take 3 arguments, position dependent.  (TODO: make these args into
#   tag/value form and provide usage.)
#   1. maxpoints
#   2. the region (chromosome:start:end)
#   3. the specific positions (chromosome:start:end,chromosome:start:end, etc.)
#
#   Example:
#   >java -classpath ./ GetCoverage 1000 13:130000:150000 13:130044:1130045,13:140042:140043
#
#    maxpoints = 1000
#    region = chromosome 13 for region spanning 130000 through 150000
#    positions = 2 positions on chromosome 13; one at position 130044, the other at position 140042
#
#    IMPORTANT - The start position follows the samtools convention of zero-based half closed.
#                Note that the output follows the vcf convention, 1-based.
#
#   OUTPUT
#   ------
#   The output from this app a set of records.  The records are delimited by
#   newline (\n) and the fields are delimited by tabs. The first field is the
#   start position and next field is the depth. The data is returned in two
#   sections.  Each section is a header record follwed by detail records.
#
#   Example:
#
#   #specific_points
#   130044 40
#   140042 54
#   #reduced_points
#   130000 33
#   130100 44
#   130200 42

import sys



maxPoints = 0
refName = ""
regionStart = 0
regionEnd = 0
specificPoints = []

points = [];
pointsReserved = []
pointsRemaining = []

def splitPoints(points, specificPoints,  pointsReserved,  pointsRemaining):
    pointIdx = 0
    keepIdx = 0;
    while pointIdx < len(points):
        keepPoint = None
        if keepIdx < len(specificPoints):
            keepPoint = specificPoints[keepIdx]

        pointRow = list()
        pointRow = points[pointIdx]
        thePoint = pointRow[0]

        if keepPoint is not None and thePoint is not None and thePoint == keepPoint:
            pointsReserved.append(pointRow)
            pointIdx += 1
            keepIdx  += 1
        else:
            pointsRemaining.append(pointRow)
            pointIdx += 1

def reducePoints(regionStart, regionEnd, points, maxPoints):
    #
    # Zero fill the points that were not provided by samtools mpileup
    #
    regionSize = regionEnd - regionStart
    zeroFilledPoints = list();
    for i in range(regionSize+1):
        row = list()
        row.append(regionStart + i)
        row.append(0)
        zeroFilledPoints.append(row)

    for i in range(len(points)):
        startPos = points[i][0];

        # only keep mpileup entry if it is inside the expected region
        if startPos >= regionStart and startPos <= regionEnd:
            idx = startPos - regionStart
            zeroFilledPoints[idx][1] = points[i][1]

    # We are done if we are not reducing the number of points
    if maxPoints <= 1:
        return points

    factor = int(len(zeroFilledPoints) / maxPoints)
    modulo = len(zeroFilledPoints) % maxPoints


    # Don't bother reducing if we aren't reducing at least by a factor
    # of 2
    if factor <= 1:
        return points

    results = list()
    sum = 0
    avgWindow = None
    windowSize = factor
    remainingModulo = modulo

    finished = False
    i = 0
    while finished is False:
        if i >= len(zeroFilledPoints):
            finished = True
            continue

        anchorRow = zeroFilledPoints[i];

        # We need to spread the remaining (modulo) over
        # the chunks we are averaging.  Add 1 to each chunk
        # until we have consumed all of the remainder. Then
        # revert back to the average window size (the factor)
        if (remainingModulo > 0):
            windowSize = factor + 1
            remainingModulo -= 1
        else:
            windowSize = factor


        # Grab n rows.  This is the window that we will average
        # the points across
        avgWindow = list()
        for x in range(windowSize):
            if i+x >= len(zeroFilledPoints):
                finished = True
                continue
            avgWindow.append(zeroFilledPoints[i+x])
        windowSize = len(avgWindow)


        max = 9999999;
        min = 0;

        # Sum the depths in the window
        for j in range(windowSize):
            y = avgWindow[j][1];
            sum += y

            if y > max:
                max = y
            if y < min:
                min = y

        # Now calculate the median depth and associate it with
        # a data point that represents the span of the window
        average = int(sum / windowSize)
        resultRow = list()
        resultRow.append(anchorRow[0])
        resultRow.append(average)
        results.append(resultRow)

        sum = 0
        i += windowSize

    return results


def printPoints(pointsReserved, pointsReduced):
    print("#specific_points")
    for point in pointsReserved:
        print(str(point[0]) + "\t" + str(point[1]))
    print("#reduced_points")
    for point in pointsReduced:
        print(str(point[0]) + "\t" + str(point[1]))


if len(sys.argv) > 1:
    maxPoints = int(sys.argv[1])

if len(sys.argv) > 2:
    rangeString  = sys.argv[2]
    rangeTokens  = rangeString.split(":")
    refName      = rangeTokens[0]
    regionStart  = int(rangeTokens[1])
    regionEnd    = int(rangeTokens[2])

# Region arg looks like this 13:10000:20000
# Grab the start position and increment by 1.
# This will be the start position to match
# the the pileup records.
if len(sys.argv) > 3:
    specificPointsString  = sys.argv[3]
    for specificPoint in specificPointsString.split(","):
        tokens = specificPoint.split(":");
        if len(tokens) == 3:
            start = int(tokens[1])
            start = start + 1
            specificPoints.append(start)



# Read all of the output from mpileup and put in
# points array
points = []
for line in sys.stdin.readlines():
    if line.startswith("[mpileup]") or line.startswith("<mpileup>"):
        continue
    tokens = line.split("\t")
    row = list()
    row.append(int(tokens[1]))
    row.append(int(tokens[3]))

    points.append(row)


splitPoints(points, specificPoints, pointsReserved, pointsRemaining)
pointsReduced = reducePoints(regionStart, regionEnd, pointsRemaining, maxPoints)
printPoints(pointsReserved, pointsReduced)









