####################################################
# server.R File: Routines for PSHACalculatoR       #
#                                                  #
# For details and program execution see: runIT.R   #
# Michael Haas                                     #
# 15.01.2014                                       #
# mhaas@gfz-potsdam.de                             #      
####################################################

library(shiny)
library(sgeostat)

shinyServer(function(input, output) {

  ###########################################################################
  ###########################################################################
  ##                                                                       ##
  ##                  A)     Calculation Routines                          ##
  ##                                                                       ##
  ###########################################################################
  ###########################################################################
  
  
  
  ###########################################################
  # Get coordinates from UI                                 #
  ###########################################################
  getCrds <- reactive({
    # point source
    if(input$sourceType=="pt"){
      coords <- matrix(0,1,2)
      coords[1,1]=input$x1
      coords[1,2]=input$y1
      
    }else if(input$sourceType=="ft"){ #fault source
      
      # calculate segments with length sqrt(2) in fault
      # 1) get user input
      x <- c(input$x1,input$x2)
      y <- c(input$y1,input$y2)
      # 2) arange data and construct segments 
      dl <- 2 #2 points per segment of length sqrt(2)
      mat <- cbind(x[1]-x[2],y[1]-y[2])
      L1 <- norm(as.matrix(mat),'F')
      nseg <- floor(L1/dl) #number of segments 1 per length dl=2
      if (nseg < 1) nseg=1 #prevent error in plotfct
      strikeAngle <- atan2(abs(y[1]-y[2]),abs(x[1]-x[2])) #strike angle - 90?
      
      # determine which point north and east of the other
      if(y[1] < y[2]) {sgn1 <- 1 #point2 north of point1
      }else sgn1 <- -1 #point2 south of point1
      
      if(x[1] < x[2]) {sgn2 <- 1 #point2 east of point1
      }else sgn2 <- -1 #point2 west of point1
      
      # construct potential source points on fault
      coords <- matrix(0,nseg+1,2)
      coords[1,1]=x[1]
      coords[1,2]=y[1]
      for (i in 1:nseg){
        coords[i+1,1]=coords[i,1]+sgn2*dl*cos(strikeAngle)
        coords[i+1,2]=coords[i,2]+sgn1*dl*sin(strikeAngle)
      } # end of for
    }else{ # end of else if fault source
    
      #create a regular grid with 1km spacing over around defined area 
      #and store values which are inside the defined area
      
      x <- c(input$x1,input$x2,input$x3,input$x4)
      y <- c(input$y1,input$y2,input$y3,input$y4)
      
      #determine starting mesh size
      xVal <- min(x)
      yVal <- min(y)
      xMax <- max(x)
      yMax <- max(y)
      xDir <- xMax-xVal+1
      yDir <- yMax-yVal+1
      size <- xDir*yDir
            
      #create mesh  
      #start in lower left corner
      for(i in 1:(size)){
        if (in.chull(xVal, yVal,x,y)){
          if (exists("coords")){
          coords <- rbind(c(xVal,yVal),coords)
          }else{coords <- c(xVal,yVal)}#first entry
        }
        if (i%%xDir==0){
          xVal <- xVal-(xDir-1)
          yVal <- yVal+1 #increase y if end of line
        }else xVal <- xVal+1
      } #end for    
    } #end else for area
    return(coords)
  }) # end of reactive
  
  ###########################################################
  # Calculate source-site distances from Coordinates        #
  ###########################################################
  getDist <- reactive({
    #call coordinates function crds
    crds <- getCrds()
    # number of points to consider (depends on user input)
    npts2 <- length(crds[,1])
    #euclidian distances; Assumption: equal rates for all points (fault/area)    
    dist <- matrix(0,1,npts2)
    for (i in 1:npts2) {
      dist[i] <- norm(as.matrix(c(crds[i,1],crds[i,2])),'F')
    }
    dist <- as.vector(dist)
    return(dist)
  })
  
  ###########################################################
  # Create Histogram from Distances                         #                                                  #
  ###########################################################
  getHist <- reactive({
             #get Distances from coordinates
             dist <- getDist()
             # make histogram with distance prob
             nbreaks <- 11 #number of breaks between classes --> 12 classes in histogram
             #making histogram (absolut count)
             h <- hist(dist,breaks=nbreaks)
             return(h)
  })
  
  ###########################################################
  # Calculate Magnitude range from user input               #
  ###########################################################
  getMw <- reactive({
           dm <- 0.2
           Mw <- seq(input$MminGR,input$MmaxGR,dm)
  })
  
  ###########################################################
  # Calculate Gutenberg-Richter from user input             #
  ###########################################################
  getRates <- reactive({
      Mw <- getMw() #considered magnitude 
      MwRange <- Mw[length(Mw)]-Mw[1]
      
      #parameters
      alpha <- 2.303*input$aGR 
      beta <- 2.303*input$bGR
      nu <- exp(alpha-beta*Mw[1])    

      #rate for each magnitude
      rateM <- matrix(999,length(Mw),1)
      for (i in 1:length(Mw)){
      rateM[i] <- nu * (exp(-beta*(Mw[i]-Mw[1]))-
                     exp(-beta*(MwRange)))/(1-exp(-beta*(MwRange)))
      }
      
      #accumulate rates
      for (i in (length(Mw)-1):1){
        rateM[i] <- rateM[i+1] + rateM[i]
      }
      rateM <- as.vector(rateM)
      return(rateM)
  })
  
  ###################################################################
  # Calculation of PHA [cm/s?] (GMPE: Cornell 1979 in Kramer)       #
  ###################################################################
  # a) fixed Distance vector/Only one Magnitude
  getPHA <- reactive({
    R <- seq(1,300,1)
    lnPha <- 6.74+0.859*input$gmpeMw-1.80*log(R+25);
    return(lnPha)
  })
  # b) distances: midpoints of histogram distances/all magnitudes from magnitude range Mw
  getPHA2 <- reactive({
    #get user defined distances and Mw range
    R <- getHist()$mids
    Mw <- getMw()
    #initialize
    lnPha <- matrix(0,length(Mw),length(R))
    #calculate lnPha for each (Mw,R) pair
    for (i in 1:length(Mw)){
      for (j in 1:length(R)){
          lnPha[i,j] <- 6.74+0.859*Mw[i]-1.80*log(R[j]+25);          
      }
    }
    return(lnPha)
  })
  
  #########################################################################
  # Calculation of Sigma of PHA at 100km dist (Normally distributed)      #   
  #########################################################################
  getSig <- reactive({
    lnPha <- getPHA()
    x <- seq((lnPha[100]-5*input$gmpeSig),(lnPha[100]+5*input$gmpeSig),0.01)
    gauss <- 100*dnorm(x,lnPha[100],input$gmpeSig)
    for (i in 1:length(x)){
      if(exists("px")){
      px <- rbind(px,(100-gauss[i]))
      py <- rbind(py,exp(x[i]))
      }else{px <- 100-gauss[1]
            py <- exp(x[1])}
    }
    pxpy <- cbind(px,py)
    return(pxpy)
  })
  
  ###############################################################
  # Calculate Probability,that lev exceeded                     #
  ###############################################################
  
  getP <- reactive({
          lnPha <- getPHA()
          lev <- 1.2*exp(lnPha[100])
          P <- (1- pnorm(log(lev),lnPha[100],input$gmpeSig))
          return(P)
          })
  
  ##################################################
  # Calculate mean annual rates                    #
  ##################################################
  getRa <- reactive({
             #get rates from magnitude distribution
             Mw <- getMw() #user defined Mw range (vector)
             rateM <- getRates() #rates for this Mw range (vector)
             #get distances from coordinates and midpoints of corresponding histogram classes (vector)              
             dist <- getDist()
             n <- getHist()$counts
             n <- n/sum(n)
             #get lnPha values for the (Mw,Dist) pairs (matrix)
             lnPha <- getPHA2()
             
             #considered g-values as vector
             g <- 981*t(seq(0.01,1,0.02))
             
             #GR rate for Mmin
             nu=10.^(input$aGR-input$bGR*input$MminGR)
             
             #setup hazard matrix
             Sum <- matrix(0,1,length(g))
             #calculate S
             for (i in 1:length(Mw)){
                  for (j in 1:length(n)){
                       # Probability for each g value that g exceeded assuming normal distribution around lnPha with sigma 
                       P <- 1- pnorm(log(g),lnPha[i,j],input$gmpeSig)
                       # sum of ((rate of Mmin)*(Exceedance Prob of lnPha)*(Prob Mag exceeded)*(Prob Distance class)) 
                       # nu is the same for all sources here since we assume Mmin the same for all of them
                       Sum=Sum+nu*P*rateM[i]*n[j]
                       rm(P)
                  }
             }
             #return g and Sum in matrix
             return(rbind(g,Sum))
    })
  
  
  #########################################################################
  #########################################################################
  ##                                                                     ##  
  ##                    B)    Plotting Routines                          ##        
  ##                                                                     ##
  #########################################################################
  #########################################################################
  
  ###########################################################
  # plot source (depends on user input: Sources type/coords)#
  ###########################################################
  output$sourcePlot <- renderPlot({
    #call coordinates function crds
    crds <- getCrds()
    #determine limits +/- 1
    all <- rbind(c(0,0),crds)
    xMax <- max(all)+1
    yMax <- max(all)+1
    xMin <- min(all)-1
    yMin <- min(all)-1
    #plot point source
    if(input$sourceType=="pt"){ 
      plot(0,0,main="Source location",
           xlab="x",ylab="y",
           xlim=c(xMin,xMax),ylim=c(yMin,yMax),pch=1)#plot site location
      points(crds[1,1],crds[1,2],col="red",pch=2)
      legend("topright", c("Site","Source"), pch=c(1,2)) 
    }
    
    #plot fault source
    else if(input$sourceType=="ft"){ #fault source
           # plot the fault segments and the site location
           #site location
           plot(0,0,main="Source location",
                xlab="x",ylab="y",
                xlim=c(xMin,xMax),ylim=c(yMin,yMax),pch=1)#plot site location
           #points of segments
           points(crds[,1],crds[,2],col="red",pch=2)
           #segments
           for (i in 1:(length(crds[,1])-1)){
            segments(crds[i,1], crds[i,2], x1 = crds[i+1,1], y1 = crds[i+1,2],col="red")
           }
           #legend
           legend("topright", c("Site","Source"), pch=c(1,2))
    }
    
    # plot area source
    else{
    
        plot(0,0,main="Source location",
             xlab="x [km]",ylab="y [km]",
             xlim=c(xMin,xMax),ylim=c(yMin,yMax),pch=1)#plot site location
        polygon(c(input$x1,input$x2,input$x3,input$x4),
                c(input$y1,input$y2,input$y3,input$y4),
                col=gray(0.8))    
        points(crds[,1], crds[,2],col="red",
             #type="p",main="Source location",xlab="x",ylab="y",
             #xlim=c(xMin,xMax),ylim=c(yMin,yMax)),
             pch=2)
        legend("topright", c("Site","Source"), pch=c(1,2))
        }
  }) #end renderPlot

  ##########################################################################
  # Plot Distance Distribution (depends on user input: Sources type/coords)#
  ##########################################################################
  output$distPlot <- renderPlot({
      h <- getHist()
      #changing density to relative/probability(assumed source homogeneity)
      h$density <- h$counts/sum(h$counts)
      #ploting histogram
      #plot(h,freq=F,breaks=h$breaks,col="blue",
      plot(h,freq=F,col="blue",
           main="Distance Probabilities",
           xlab="Distance [km]",ylab="Probability",axes=F)
      axis(2)
      axis(1, at=h$mids)
  }) #end renderPlot
  
  #######################################################
  # Plot GR Magnitude Distribution
  #######################################################
  output$grPlot <- renderPlot({
    #Calculate rates
    rateM <- getRates()
    #plot the GR Magnitude Distribution according to user input
    dm <- 0.2 #Magnitude increment
    Mw <- paste(seq(input$MminGR,input$MmaxGR,dm))
    plt <- barplot(rateM,main="Magnitude Distribution",
                   xlab="Magnitude", ylab="Mean annual rate",
                   names.arg=Mw)
  }) #end renderPlot
  
  #######################################################
  # Plot PHA (GMPE: Cornell1979 in Kramer)              #
  #######################################################
  output$gmpePlot <- renderPlot({  
    # get values
    lnPha <- getPHA()
    sig <- getSig()
    
    #plot GMPE
    R <- seq(1,300,1)
    plot(R,exp(lnPha + input$gmpeSig),log="xy",pch="+",col="red",
         xlab="Distance [km]",ylab="PHA [cm/s?]",
         main="PHA (GMPE: Cornell 1979)")
    points(R,exp(lnPha),pch="+",)
    points(R,exp(lnPha - input$gmpeSig),col="red",pch="+",)
    #points(R,exp(lnPha),log="xy",pch="+",)
    #points(R,exp(lnPha - input$gmpeSig),log="xy",col="red",pch="+",)
    legend("topright", c("mean PHA","PHA +/- 1 sigma","120 % mean-PHA-Level"),pch=c("+","+","-"),col=c('black','red','blue'))
    
    #plot pdf
    points(sig[,1],sig[,2],col="green",pch="-")
    # horizontal to read PHA value of pdf peak
    lines(c(1,1000),exp(c(lnPha[100],lnPha[100])),col="black") 
    lines(c(100,100),c(1,10000),col="black")
    
    #plot line at 120% of mean PHA
    lev <- 1.2*exp(lnPha[100])
    lines(c(1,1000),c(lev,lev),col="blue")
    
    #mark part of pdf above lev
    #find function: values larger equal lev
    rfind <- function(x) seq(along=x)[x >= lev]
    i=rfind(sig[,2])
    xsel=sig[i,1]
    ysel=sig[i,2]
    #attach border values
    xsel=rbind(xsel,100) 
    ysel=rbind(ysel,(sig[i[1],2])) #first pha value above lev
    #add to plot
    polygon(xsel,ysel,col='green',border="green")
    
  }) #end renderPlot
  
  ######################################################
  # Show probability, that lev is exceeded             #
  ######################################################
  output$excP <- renderText({
    #show probability of 120% value exceeded
    paste("The probability, that the mean PHA value at 100 km distance
           is exceeded by more than 20% is: ",sprintf("%.2f",100*getP())," %")
  })
  
  ######################################################
  # Plot mean annual rate of exc.& prob. of exc in OP  #
  ######################################################
  output$annRatePlot <- renderPlot({
        #get mean annual rates
        meanAnRates <- getRa()
        #semilog plot of mean annual rate
        plot((meanAnRates[1,]/981),meanAnRates[2,],log="y",
             pch="-", xlab="PHA [g]", ylab="Mean annual rate of exceedance",
             main="Seismic Hazard at Site")
  })
  
  output$excProbPlot <- renderPlot({
          #Get mean annual rates of exceedance  
          meanAnRates <- getRa()
          #Poisson prob
          probPoiss=1-exp((-1)*meanAnRates[2,]*input$OP);
          
          #semilog plot of prob of exc in input$OP years
          
          plot((meanAnRates[1,]/981),probPoiss,pch="-",
                    log="y",xlab="PHA [g]", 
                    ylab=sprintf('Probability of exceedence in %.0f years',input$OP),
                    main="Poisson Model",axes=F)
          axis(1)
          axis(2, at=c(10^-5,10^-4,10^-3,10^-2,10^-1,10^0))
  })
    
}) #end server.R

