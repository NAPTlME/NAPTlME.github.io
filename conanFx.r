# Supporting functions for use in the Conan Map project for Solitude's Server
# Nathan Pratt
# 8/17/18
requiredPkgs = c("abind", "reshape2", "grid", "stringr", "raster")
missingPkgs = requiredPkgs[!(requiredPkgs %in% installed.packages()[,"Package"])]
if(length(missingPkgs)) install.packages(missingPkgs)
library(stringr)
library(abind)
library(reshape2)
library(grid)
library(raster)
# Set variable for array file
claimMatrixFile = "MetaData/mapMatrixBaseInfo.rds"
# Creates empty matrix to dimensions of map
CreateEmptyMapMatrix = function(rowNames,colNames){
  m = matrix(rep(0, times = length(rowNames) * length(colNames)), nrow = length(rowNames), ncol = length(colNames))
  rownames(m) = rowNames
  colnames(m) = colNames
  return(m)
}
CreateEmptyMapArray = function(rowNames, colNames, dNames = c("POI", "Out.Of.Play")){
  a = array(0, dim = c(length(rowNames), length(colNames), length(dNames)), dimnames = list(rowNames, colNames, dNames))
  return(a)
}
AddClanToMapArray = function(mapArray, clanName){
  if (clanName %in% dimnames(mapArray)[[3]]){
    warning(paste0(clanName, " already exists."))
  } else {
    rowNames = dimnames(mapArray)[[1]]
    colNames = dimnames(mapArray)[[2]]
    clanArray = CreateEmptyMapArray(rowNames, colNames, clanName)
    mapArray = abind(mapArray, clanArray, along = 3)
  }
  return(mapArray)
}
ConvertArrayValuesToT = function(m, rowcols, owner, claimableCoords, overwrite = F){
  rowcols = toupper(rowcols) #prevents case sensitivity
  lRowCols = lapply(rowcols, function(x) c(str_match(x, "[A-Z]+")[1,1], str_match(x, "[\\-\\d]+")))
  # check for valid cells
  invalidRows = sapply(lRowCols, function(x) !(x[1] %in% rownames(m)))
  invalidCols = sapply(lRowCols, function(x) !(x[2] %in% colnames(m)))
  invalidCells = invalidRows | invalidCols
  if (sum(invalidCells) > 0){
    warning(paste0("Invalid Cells detected (will ignore): ", paste0(rowcols[invalidCells], collapse = ", "), "."))
    lRowCols = lRowCols[!invalidCells] # use only valid cells
  }
  for (i in 1:length(lRowCols)){
    #warning(paste("Rowcols: ", rowcols[i]))
    #warning(str(claimableCoords))
    if (!overwrite){
      # Check if already owned, if not then assign
      if (sum(m[lRowCols[[i]][1], lRowCols[[i]][2],2:length(dimnames(m)[[3]])]) > 0){
        warning(paste0(rowcols[i], " cannot be assigned as it is already owned."))
      } else if (!(rowcols[i] %in% names(claimableCoords))){
        warning(paste0(rowcols[i], " cannot be assigned as it is not claimable."))
      } else {
        m[lRowCols[[i]][1], lRowCols[[i]][2], 1] = 0
        # check for case sensitivity in 
        m[lRowCols[[i]][1], lRowCols[[i]][2], owner] = 1
      }
    } else {
      # convert all other values to 0
      if (!(rowcols[i] %in% names(claimableCoords))){
        warning(paste0(rowcols[i], " cannot be assigned as it is not claimable."))
      } else {
        m[lRowCols[[i]][1], lRowCols[[i]][2], ] = 0
        m[lRowCols[[i]][1], lRowCols[[i]][2], owner] = 1 # overwrite
      }
    }
  }
  return(m)
}
GetDfFromArray = function(mapArray){
  # Get first owner found in matrix ### should only ever be one, but there's no reason to trust
  tmpArray = array(apply(mapArray, c(1, 2), function(x) 
    if(sum(as.numeric(x)) > 0){ 
      names(x)[which(as.numeric(x) == 1)[1]] 
    } else { 
      "None" 
    }), 
    dim = dim(mapArray), dimnames = list(dimnames(mapArray)[[1]], dimnames(mapArray)[[2]]))
  return(melt(tmpArray))
}
CreateBasePolygon = function(coords){
  startX = numeric(0)
  startY = numeric(0)
  endX = numeric(0)
  endY = numeric(0)
  for(coord in coords){
    coordY = which(LETTERS == str_match(toupper(coord), "[A-Z]")[1,1])
    coordX = as.numeric(str_match(coord, "[\\d\\-]+")[1,1])
    startX = c(startX, coordX - 0.5)
    endX = c(endX, coordX + 0.5)
    startY = c(startY, coordY - 0.5)
    endY = c(endY, coordY + 0.5)
  }
  returnVal = as.matrix(expand.grid(c(min(startY), max(endY)), c(min(startX), max(endX))))
  dimnames(returnVal)[[2]] = c("Y", "X")
  returnVal = returnVal[c(1,2,4,3), 2:1] # makes coords start in bottom left and go clockwise
  return(returnVal)
}
CreateCoordList = function(cells){
  cells = toupper(cells)
  returnVal = list()
  for (i in 1:length(cells)){
    mVal = CreateBasePolygon(cells[i])
    returnVal[[i]] = list(Polygons = list(mVal), Holes = list())
    #print(returnVal[[i]])
    #print("Names:")
    #print(names(returnVal[[i]]))
  }
  setNames(returnVal, cells)
}
GetExtent = function(coords){
  startX = numeric(0)
  startY = numeric(0)
  endX = numeric(0)
  endY = numeric(0)
  for(coord in coords){
    coordY = which(LETTERS == str_match(toupper(coord), "[A-Z]")[1,1])
    coordX = as.numeric(str_match(coord, "[\\d\\-]+")[1,1])
    startX = c(startX, coordX - 0.5)
    endX = c(endX, coordX + 0.5)
    startY = c(startY, coordY - 0.5)
    endY = c(endY, coordY + 0.5)
  }
  return(extent(x = min(startX), xmax = max(endX), ymin = min(startY), ymax = max(endY)))
}
plotZoomedImage = function(imag, coords){
  ex = GetExtent(coords)
  tmpImage = crop(imag, ex)
  xAxis = seq(from = ex@xmin + 0.5, to = ex@xmax - 0.5, by = 1)
  yAxis = seq(from = ex@ymin + 0.5, to = ex@ymax - 0.5, by = 1)
  plot(tmpImage, 1, xaxt = "n", yaxt = "n")
  axis(1, at = xAxis, labels = xAxis)
  axis(2, at = yAxis, labels = LETTERS[yAxis])
  # add vlines and hlines
  abline(v = seq(from = ex@xmin, to = ex@xmax, by = 1),  col= "white")
  abline(h = seq(from = ex@ymin, to = ex@ymax, by = 1), col = "white")
}
plotImg = function(imag, coords = NULL){ # expects a rasterBrick
  if (!is.null(coords)){
    imag = crop(imag, GetExtent(coords))
  }
  ex = extent(imag)
  plot(imag, 1, xaxt = "n", yaxt = "n")
  axis(1, at = 0:24, labels = 0:24)
  axis(2, at = 1:20, labels = LETTERS[1:20])
  abline(v = seq(from = ex@xmin, to = ex@xmax, by = 1), col = "white")
  abline(h = seq(from = ex@ymin, to = ex@ymax, by = 1), col = "white")
}
GetRoundedClickPoints = function(n = 80){
  return(round(click(n = n), 2))
}
setCoords = function(availCells, claimableCoords, OutOfPlayCoords, poiCoords, img){
  if (file.exists("ConanMap/cmMapCoordsTmp.rds")){
    file.remove("ConanMap/cmMapCoordsTmp.rds")
  }
  on.exit(saveRDS(list(Claimable = claimableCoords, OutofPlay = OutOfPlayCoords, POI = poiCoords), "ConanMap/cmMapCoordsTmp.rds"))
  for (cell in availCells){
    isClaimable = cell %in% names(claimableCoords)
    isOOP = cell %in% names(OutOfPlayCoords)
    isPOI = cell %in% names(poiCoords)
    print(paste0("Cell: ", cell))
    plotImg(img, cell)
    tmp = paste(paste0(c("Claimable: ", "OutofPlay: ", "POI: "), c(isClaimable, isOOP, isPOI), collapse = "\n"),
                 "Would you like to make any changes to the types of coords found in this cell? (y/n/stop/zOut/zIn): "
                 , sep = "\n")
    input = "zIn"
    while(input == "zIn" | input == "zOut"){
      input = readline(tmp)
      if (input == "zIn"){
        plotImg(img, cell)
      } else if (input == "zOut"){
        plotImg(img)
      }
    }
    if(tolower(input) == "stop"){
      return(list(claimableCoords, OutOfPlayCoords, poiCoords))
    }
    if(tolower(input) == "y"){
      tmp = modifyCells(cell, claimableCoords, OutOfPlayCoords, poiCoords)
      claimableCoords = tmp[[1]]
      OutOfPlayCoords = tmp[[2]]
      poiCoords = tmp[[3]]
      isClaimable = cell %in% names(claimableCoords)
      isOOP = cell %in% names(OutOfPlayCoords)
      isPOI = cell %in% names(poiCoords)
    }
    finished = F
    if (sum(c(isClaimable, isOOP, isPOI)) > 1){
      while(!finished){
        plotImg(img, cell)
        claimCorners = 0
        oopCorners = 0
        poiCorners = 0
        if(isClaimable){
          claimCorners = as.numeric(str_match_all(readline("Enter base corners for Claim polygon: "), "\\d+")[[1]][,1])
          if (length(claimCorners) == 0){
            claimCorners = 0
          }
        }
        if(isOOP){
          oopCorners = as.numeric(str_match_all(readline("Enter base corners for OOP polygon: "), "\\d+")[[1]][,1])
          if (length(oopCorners) == 0){
            oopCorners = 0
          }
        }
        if(isPOI){
          poiCorners = as.numeric(str_match_all(readline("Enter base corners for POI polygon: "), "\\d+")[[1]][,1])
          if (length(poiCorners) == 0){
            poiCorners = 0
          }
        }
        numLinesToDraw = 0
        while(numLinesToDraw == 0){
          input = as.numeric(readline("Enter an integer for the number of lines to draw (or 'stop' to return results up to now): "))
          if (tolower(input) == "stop"){
            return(list(claimableCoords, OutOfPlayCoords, poiCoords))
          }
          if (input > 0){
            numLinesToDraw = input
          }
        }
        eLines = lapply(1:numLinesToDraw, function(x) click(n = 80))
        # Claimable Lines
        claimableLinesToUse = numeric(0)
        if (isClaimable){
          while(length(claimableLinesToUse) == 0){
            input = readline(paste0("Enter lines to use for Claim or type 'show' with line indices to plot on map (1:", numLinesToDraw, "): "))
            plotImg(img, cell)
            if (str_detect(tolower(input), "show.*\\d")){
              sapply(as.numeric(str_match_all(input, "\\d+")[[1]][,1]), function(x) lines(eLines[[x]]))
            } else if (str_detect(input, "\\d+")){
              claimableLinesToUse = as.numeric(str_match_all(input, "\\d+")[[1]][,1])
            }
          }
        }
        # OOP Lines
        oopLinesToUse = numeric(0)
        if (isOOP){
          while(length(oopLinesToUse) == 0){
            input = readline(paste0("Enter lines to use for OOP or type 'show' with line indices to plot on map (1:", numLinesToDraw, "): "))
            plotImg(img, cell)
            if (str_detect(tolower(input), "show.*\\d")){
              sapply(as.numeric(str_match_all(input, "\\d+")[[1]][,1]), function(x) lines(eLines[[x]]))
            } else if (str_detect(input, "\\d+")){
              oopLinesToUse = as.numeric(str_match_all(input, "\\d+")[[1]][,1])
            }
          }
        }
        # POI Lines
        poiLinesToUse = numeric(0)
        if (isPOI){
          while(length(poiLinesToUse) == 0){
            input = readline(paste0("Enter lines to use for POI or type 'show' with line indices to plot on map (1:", numLinesToDraw, "): "))
            plotImg(img, cell)
            if (str_detect(tolower(input), "show.*\\d")){
              sapply(as.numeric(str_match_all(input, "\\d+")[[1]][,1]), function(x) lines(eLines[[x]]))
            } else if (str_detect(input, "\\d+")){
              poiLinesToUse = as.numeric(str_match_all(input, "\\d+")[[1]][,1])
            }
          }
        }
        # Distplay tmp polygons
        plotImg(img, cell)
        if (isClaimable){
          #print(paste("Claim Corners: ", paste0(claimCorners, collapse = ", ")))
          tmp = claimableCoords[[cell]]$Polygons[[1]][claimCorners,]
          #print(tmp)
          for (i in claimableLinesToUse){
            tmp = rbind(tmp, eLines[[i]])
          }
          polygon(tmp, col = adjustcolor("blue", alpha.f = 0.3), border = NA)
        }
        if (isOOP){
          #print(OutOfPlayCoords[[cell]]$Polygons[[1]])
          tmp = OutOfPlayCoords[[cell]]$Polygons[[1]][oopCorners,]
          #print(tmp)
          for (i in oopLinesToUse){
            tmp = rbind(tmp, eLines[[i]])
          }
          polygon(tmp, col = adjustcolor("red", alpha.f = 0.3), border = NA)
        }
        if (isPOI){
          tmp = poiCoords[[cell]]$Polygons[[1]][poiCorners,]
          for (i in poiLinesToUse){
            tmp = rbind(tmp, eLines[[i]])
          }
          polygon(tmp, col = adjustcolor("white", alpha.f = 0.3), border = NA)
        }
        #prompt user if finished
        input = readline("Finalize this cell? (y/n): ")
        finished = input == "y"
        if (finished){
          if (isClaimable){
            tmp = claimableCoords[[cell]]$Polygons[[1]][claimCorners,]
            for (i in claimableLinesToUse){
              tmp = rbind(tmp, eLines[[i]])
            }
            claimableCoords[[cell]]$Polygons[[1]] = tmp
          }
          if (isOOP){
            tmp = OutOfPlayCoords[[cell]]$Polygons[[1]][oopCorners,]
            for (i in oopLinesToUse){
              tmp = rbind(tmp, eLines[[i]])
            }
            OutOfPlayCoords[[cell]]$Polygons[[1]] = tmp
          }
          if(isPOI){
            tmp = poiCoords[[cell]]$Polygons[[1]][poiCorners,]
            for (i in poiLinesToUse){
              tmp = rbind(tmp, eLines[[i]])
            }
            poiCoords[[cell]]$Polygons[[1]] = tmp
          }
        }
      }
    }
  }
  return(list(claimableCoords, OutOfPlayCoords, poiCoords))
}
modifyCells = function(coord, claimableCoords, OutOfPlayCoords, poiCoords, index = 1){
  #prompt to remove first then add
  alreadyInCell = c(coord %in% names(claimableCoords), coord %in% names(OutOfPlayCoords), coord %in% names(poiCoords))
  allLists = c(1:3)
  input = readline(paste(paste0(allLists[alreadyInCell], ": ", c("Claimable", "Out of Play", "POI")[alreadyInCell],collapse = "\n"),
                          "Enter indices desired for removal separated by a ',': ", sep = "\n"))
  rmNums = as.numeric(str_match_all(input, "\\d+")[[1]][,1])
  if(any(!(rmNums %in% allLists[alreadyInCell]))){
    #warn and remove
    warning("Inputs should be between 1 and 3 and currently have coords in this cell")
    rmNums = rmNums[rmNums %in% allLists[alreadyInCell]]
  }
  if (1 %in% rmNums){
    claimableCoords[[coord]] = NULL
  }
  if (2 %in% rmNums){
    OutOfPlayCoords[[coord]] = NULL
  }
  if (3 %in% rmNums){
    poiCoords[[coord]] = NULL
  }
  notInCell = c(allLists[!alreadyInCell], rmNums)
  input = readline(paste(paste0(notInCell, ": ", c("Claimable", "Out of Play", "POI")[notInCell]),collapse = "\n",
                         "Enter indices desired for addition separated by a ',': ", sep = "\n"))
  adNums = as.numeric(str_match_all(input, "\\d+")[[1]][,1])
  #print(adNums)
  if (any(!(adNums %in% notInCell))){
    # warn and remove
    warning("Inputs should be between 1 and 3 and not currently have coordinates in this cell")
    adNums = adNums[adNums %in% notInCell]
  }
  if (1 %in% adNums){
    claimableCoords[[coord]] = list(Polygons = list(CreateBasePolygon(coord)),
                                    Holes = list())
  }
  if (2 %in% adNums){
    #print("2 in adnums")
    #print(names(OutOfPlayCoords))
    OutOfPlayCoords[[coord]] = list(Polygons = list(CreateBasePolygon(coord)),
                                   Holes = list())
    #print(names(OutOfPlayCoords))
  }
  if (3 %in% adNums){
    poiCoords[[coord]] = list(Polygons = list(CreateBasePolygon(coord)),
                              Holes = list())
  }
  return(list(claimableCoords, OutOfPlayCoords, poiCoords))
}