###function to interpolate the input coordinates to quiavalent distance
interpolateZ <- function(x,y,xout,space = 8){
    stopifnot(length(x)==length(y))
    if(any(duplicated(x))){
        ind = which(duplicated(x))
        x = x[-ind]
        y = y[-ind]
    }
    xInter.left = pmax(x - space,c(0,x[-length(x)])) ##shift "space" step to the left, if such shift is beyond the adjacent left probe, then we take the position of the left probe 
    xInter.right = pmin(x + space,c(xInter.left[-1],max(x)+space)) ##shift "space" step to the right, if such shift is beyond the ajacent already left shift probes, then we take the position of the left shift probes 
    xAdd = setdiff(xInter.right,xInter.left) ##the missing position that is in the right shift, but not in the left shift should be added and considered as the border that could take a value
    xInsert = c(xInter.left,xAdd)
    yInsert = c(y,rep(-1000,length(xAdd)))
    yInsert[is.na(yInsert)] = -1000 ##this step is very important to make sure that there is no NAs in y,otherwise the interpolate will have problems
    ord = order(xInsert)
    rv = approx(xInsert[ord], yInsert[ord], xout, method="constant", f=0)$y
    rv[rv== -1000] = NA
    return(rv)
}
raster.image = function(x, y, z, uniq=uniq,
        colRamp = colorRamp(brewer.pal(9, "YlGnBu")),
        width = 8,space = 8,...)   ## the typical spacing between neighbouring probes 
{
    
    rg = range(z, na.rm=TRUE)
    if( (rg[1]<0) || (rg[2]>1) )
        warning("rasterImage: 'z' contained values outside [0,1].\n")
    
    nx = nrow(z)
    ny = ncol(z)
    stopifnot(length(x)==nx, length(y)==ny, length(uniq)==nx)
    
    ## equi-distant steps in x and y direction
#    x.equi = seq(x[1], x[nx], length=nx) ##equi distance with the same number of the data points
    x.equi = seq(x[1], x[nx], by=width) ##equi distance with a fixed space witdh
    y.equi = seq(y[1], y[ny], length=ny)
    
    # check that values along y are equi-distant
    stopifnot(all(abs(y-y.equi)<1e-6))
    
    # interpolate z
#    z.equi = apply(z, 2, function(v) approx(x, y=v, xout=x.equi, method="constant", f=0.5)$y)
    z.equi = apply(z, 2, function(v) interpolateZ(x, y=v, xout=x.equi,space = space)) 
    
    mcol = colRamp(as.vector(z.equi)) / 256
    mcol[ uniq!=0, ] = 1
    indna <- is.na(mcol[,1])
    mcol[indna,] = 1
    rast = rgb(mcol[,1], mcol[,2], mcol[,3])
    dim(rast) = c(nrow(z.equi), ncol(z.equi))
    
    ## transpose and revert row order
    trsf = function(x) t(x)[seq(from=ncol(x), to=1, by=-1),, drop=FALSE]
    
    grid.raster(image = trsf(rast),
            y = y.equi[1], height = y.equi[length(y.equi)] - y.equi[1], hjust=0, 
            x = x.equi[1], width  = x.equi[length(x.equi)] - x.equi[1], vjust=0,
            default.units = "native", interpolate=FALSE)
    
}


