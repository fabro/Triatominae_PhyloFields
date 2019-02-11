ConvexHull <- function(x, alpha=250000, increment=360, rnd=2)   {
    if (!require (sp)) stop("sp PACKAGE MISSING")
      if (!require (alphahull)) stop("sp PACKAGE MISSING")
        if (!require(maptools)) stop("maptools PACKAGE MISSING")
    if (!inherits(x, "SpatialPointsDataFrame") |  !inherits(x, "SpatialPoints") ) 
      stop(deparse(substitute(x)), " MUST BE A sp Points OBJECT")
        x.coords <- coordinates(x)    
    ahull2sp <- function(x, increment=360, rnd=10, proj4string=CRS(as.character(NA))){
      if (class(x) != "ahull")
        stop("x needs to be an ahull class object")
      xdf <- as.data.frame(x$arcs)
      xdf <- subset(xdf,xdf$r > 0)
      res <- NULL
      if (nrow(xdf) > 0){
        linesj <- list()
        prevx<-NULL
        prevy<-NULL
        j<-1
        for(i in 1:nrow(xdf)){
          rowi <- xdf[i,]
          v <- c(rowi$v.x, rowi$v.y)
          theta <- rowi$theta
          r <- rowi$r
          cc <- c(rowi$c1, rowi$c2)
          ipoints <- 2 + round(increment * (rowi$theta / 2),0)
          angles <- anglesArc(v, theta)
          seqang <- seq(angles[1], angles[2], length = ipoints)
          x <- round(cc[1] + r * cos(seqang),rnd)
          y <- round(cc[2] + r * sin(seqang),rnd)
          if (is.null(prevx)){
            prevx<-x
            prevy<-y
          } else if (x[1] == round(prevx[length(prevx)],rnd) && y[1] == round(prevy[length(prevy)],rnd)){
              if (i == nrow(xdf)){
                prevx<-append(prevx,x[2:ipoints])
                prevy<-append(prevy,y[2:ipoints])
                prevx[length(prevx)]<-prevx[1]
                prevy[length(prevy)]<-prevy[1]
                coordsj<-cbind(prevx,prevy)
                colnames(coordsj)<-NULL
                linej <- Line(coordsj)
                linesj[[j]] <- Lines(linej, ID = as.character(j))
              } else {
                prevx<-append(prevx,x[2:ipoints])
                prevy<-append(prevy,y[2:ipoints])
              }
            } else {
              prevx[length(prevx)]<-prevx[1]
              prevy[length(prevy)]<-prevy[1]
              coordsj<-cbind(prevx,prevy)
              colnames(coordsj)<-NULL
              linej <- Line(coordsj)
              linesj[[j]] <- Lines(linej, ID = as.character(j))
              j<-j+1
              prevx<-NULL
              prevy<-NULL
            }
        }
        lspl <- SpatialLines(linesj)
        lns <- slot(lspl, "lines")
        polys <- sapply(lns, function(x) { 
          crds <- slot(slot(x, "Lines")[[1]], "coords")
          identical(crds[1, ], crds[nrow(crds), ])
        }) 
        polyssl <- lspl[polys]
        list_of_Lines <- slot(polyssl, "lines")
        sppolys <- SpatialPolygons(list(Polygons(lapply(list_of_Lines, function(x) {
    Polygon(slot(slot(x, "Lines")[[1]], "coords")) }), ID = "1")), proj4string=proj4string)
        hid <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "ID"))
          areas <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "area"))
          df <- data.frame(hid,areas)
            names(df) <- c("HID","Area")
            rownames(df) <- df$HID
        res <- SpatialPolygonsDataFrame(sppolys, data=df)
       res <- res[which(res@data$Area > 0),]
      }  
      return(res)
  }
   a <- ahull(x.coords, alpha=alpha)
    return(ahull2sp(a, rnd=rnd)) 
}
