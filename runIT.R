#######################################################
#                                                     #
# Webapplication to calculate Seismic Hazard          # 
#                                                     #
# Original Matlab program by Dr. Dino Bindi           #
# bindi@gfz-potsdam.de                                #  
#                                                     #
# Translated to R by Michael Haas                     #
# 14.01.2014                                          #
# mhaas@gfz-potsdam.de                                #
#                                                     #
# Dependencies besides standard R:                    #
# - "shiny"-package: install.packages("shiny")        #
# - "raster"-package: install.packages("raster")      #
# - "maptools"-package: install.packages("maptools")  #
# - "sp"-package: install.packages("sp")              #
# - "sgeostat"-package: install.packages("sgeostat")  #
# - "pracma"-package: install.packages("pracma")      #
#                                                     #
# also required: Webbrowser                           #
#                                                     #
# Run the program by executing the code lines below   #
#                                                     #
# Note: Specify the path to the server.R and ui.R     # 
#       containing directory inthe variable           #
#       "path" below!                                 # 
#                                                     #
# Note2: Program is blocking.End execution by hitting # 
#        "esc"-key in  R-console                      #  
#                                                     #
#######################################################

rm(list=ls(all=TRUE))

library(shiny)
Sys.setlocale('LC_ALL','C') 
#specify path
#path <- "C:\\Users\\michael\\Desktop\\PhD\\Risk suite\\yurt" 
path <- "/home/mhaas/PhD/Routines/PSHACalculatoR/"

#uncomment line below to trace execution and errors:
#options(shiny.trace = TRUE)
runApp(appDir = path, port = 1234)
