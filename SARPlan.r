####################################################
# SnowSAR Flight Planner
# JKing 16/08/2016
#
# Inputs:
# Ouputs:
####################################################
rm(list = ls()) #Clear workspace
graphics.off()

library(sp)  # vector data
library(raster)  # raster data
library(rgdal)  # input/output, projections
library(rgeos)  # geometry ops
library(spdep)  # spatial dependence
library(ggplot2)



##################################################################
getMinBBox <- function(xy) {
    stopifnot(is.matrix(xy), is.numeric(xy), nrow(xy) >= 2, ncol(xy) == 2)

    ## rotating calipers algorithm using the convex hull
    H    <- chull(xy)                    # hull indices, vertices ordered clockwise
    n    <- length(H)                    # number of hull vertices
    hull <- xy[H, ]                      # hull vertices

    ## unit basis vectors for all subspaces spanned by the hull edges
    hDir  <- diff(rbind(hull, hull[1,])) # account for circular hull vertices
    hLens <- sqrt(rowSums(hDir^2))       # length of basis vectors
    huDir <- diag(1/hLens) %*% hDir      # scaled to unit length

    ## unit basis vectors for the orthogonal subspaces
    ## rotation by 90 deg -> y' = x, x' = -y
    ouDir <- cbind(-huDir[ , 2], huDir[ , 1])

    ## project hull vertices on the subspaces spanned by the hull edges, and on
    ## the subspaces spanned by their orthogonal complements - in subspace coords
    projMat <- rbind(huDir, ouDir) %*% t(hull)

    ## range of projections and corresponding width/height of bounding rectangle
    rangeH  <- matrix(numeric(n*2), ncol=2)   # hull edge
    rangeO  <- matrix(numeric(n*2), ncol=2)   # orth subspace
    widths  <- numeric(n)
    heights <- numeric(n)
    for(i in seq(along=H)) {
        rangeH[i, ] <- range(projMat[  i, ])
        rangeO[i, ] <- range(projMat[n+i, ])  # orth subspace is in 2nd half
        widths[i]   <- abs(diff(rangeH[i, ]))
        heights[i]  <- abs(diff(rangeO[i, ]))
    }

    ## extreme projections for min-area rect in subspace coordinates
    eMin  <- which.min(widths*heights)   # hull edge leading to minimum-area
    hProj <- rbind(   rangeH[eMin, ], 0)
    oProj <- rbind(0, rangeO[eMin, ])

    ## move projections to rectangle corners
    hPts <- sweep(hProj, 1, oProj[ , 1], "+")
    oPts <- sweep(hProj, 1, oProj[ , 2], "+")

    ## corners in standard coordinates, rows = x,y, columns = corners
    ## in combined (4x2)-matrix: reverse point order to be usable in polygon()
    basis <- cbind(huDir[eMin, ], ouDir[eMin, ])  # basis formed by hull edge and orth
    hCorn <- basis %*% hPts
    oCorn <- basis %*% oPts
    pts   <- t(cbind(hCorn, oCorn[ , c(2, 1)]))

    return(list(pts=pts, width=widths[eMin], height=heights[eMin]))
}
#########################################################################################


workingDir = 'C:/Users/kingj/Documents/Projects/2016-2017/080816_SnowEx/Tools/SARPlan/geo_data'#getwd()
GCS = "+proj=longlat +ellps=WGS84" #CRS geographic coordniate system WGS84
PCS = "+proj=utm +zone=12 +south +units=m +ellps=WGS84" #CRS projected coordniate system UTM16N/WGS84
overlap = 25 #overlap in %
sarSwath = 400 #swath width in m

#Aircraft const
airSpeed = 248 #in NM/h; 248 = Wallops P3
airAltitude = 500 #Flight altitude in m; Not currently used
timeTurm = 0.1666 #Conservative estimate of the time per turn and align 0.1666 hours = 10 min

polyData <- readOGR(dsn = workingDir, layer = "GM_ROI_Small")
polyData <- spTransform(polyData, CRS(PCS))
polyDF = fortify(polyData)

xy = data.matrix(polyDF[1:2])
mbb <- getMinBBox(xy)       ## minimum bounding box
H   <- chull(xy)            ## convex hull

# plot original points, convex hull, and minimum bounding box
plot(xy, xlab="Easting (m)", ylab="Northing (m)", asp=1, type="n",
         xlim=range(c(xy[ , 1], mbb$pts[ , 1])),
         ylim=range(c(xy[ , 2], mbb$pts[ , 2])))
polygon(xy[H, ], col=NA)    ## show convex hull
polygon(mbb$pts, border="blue", lwd=2)
points(xy, pch=16, cex=1.5)

seqDist = rep(0,nrow(mbb$pts))
lineSeg = matrix(, nrow = 4, ncol = 3)
lineId = 0
segID = -1
for (pt in 1:nrow(mbb$pts)){
	if(pt == nrow(mbb$pts)){
	seqDist[pt] = sqrt((mbb$pts[1,1] - mbb$pts[pt,1])^2 + (mbb$pts[1,2] - mbb$pts[pt,2])^2)
	}
	else{
	seqDist[pt] = sqrt((mbb$pts[pt,1] - mbb$pts[(pt+1),1])^2 + (mbb$pts[pt,2] - mbb$pts[(pt+1),2])^2)
	}
	if (round(seqDist[pt])  == round(mbb$width)){
		segID = segID +  2
		lineId = lineId + 1
		lineSeg[segID:(segID+1),1] = lineId
		lineSeg[segID,2:3] = mbb$pts[pt,]
		if(pt == nrow(mbb$pts)){
				lineSeg[segID+1,2:3] = mbb$pts[1,]
			}
		else{
				lineSeg[segID+1,2:3] = mbb$pts[pt+1,]
			}	
	}
}

lineDF = data.frame(lineSeg)
colnames(lineDF)<-c("ID", "N", "E")

coordinates(lineDF) <- ~ N+E
proj4string(lineDF) <- CRS(PCS)

L1= Line(lineDF[1:2,])
L2 = Line(lineDF[3:4,])
 
Ls1 = Lines(list(L1), ID = "a")
Ls2 = Lines(list(L2), ID = "b")
SL = SpatialLines(list(Ls1, Ls2))
flightPoints = mbb$width/(sarSwath-(sarSwath*(overlap/100)))

ptsLs1 <- spsample(Ls1, flightPoints, type="regular")
ptsLs2 <- spsample(Ls2, flightPoints, type="regular")
proj4string(ptsLs1) <- CRS(PCS)
proj4string(ptsLs2) <- CRS(PCS)
points(ptsLs1)
points(ptsLs2)


numLines = length(ptsLs1)
 
flightList <- vector("list", numLines)
segLengthKm <- rep(NA, numLines)
segLengthNm <- rep(NA, numLines)
airTime <- rep(NA, numLines)

for(lns in 0:(numLines-1)) {
	flightList[[lns+1]] = Lines(Line(rbind(coordinates(ptsLs1[lns+1]), coordinates(ptsLs2[length(ptsLs2)-lns]))), ID=as.character(lns+1))
	#segLengthKm is approximate (ecludian distance)!
	segLengthKm[lns+1] = (sqrt((coordinates(ptsLs1[lns+1])[1] - coordinates(ptsLs2[length(ptsLs2)-lns])[1])^2+(coordinates(ptsLs1[lns+1])[2] - coordinates(ptsLs2[length(ptsLs2)-lns])[2])^2))/1000
	#segLengthNM[lns+1] = segLengthKm[lns+1]/1.852
	airTime[lns+1] = segLengthKm[lns+1]/(airSpeed*1.852*1000/3600)*1000/3600 #in hours
	}
	
totalAirTime = sum(airTime)+(numLines*timeTurm)
flightlines = SpatialLines(flightList)

footprintPoly <- gBuffer(flightlines, width=sarSwath/2, byid=TRUE, capStyle="FLAT", quadsegs=10)
footprintDf <- data.frame(ID=as.character(c(1:length(ptsLs1))))
footprintOut <- SpatialPolygonsDataFrame(footprintPoly,footprintDf)
writeOGR(obj=footprintOut, dsn=workingDir, layer="SARF12", driver="ESRI Shapefile")



plot(footprintOut, add=TRUE)

