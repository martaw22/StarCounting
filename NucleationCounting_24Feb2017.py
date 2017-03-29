# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 12:10:05 2016

@author: Marta
"""

#from multiprocessing import Pool
import os
from PIL import Image
import PIL.ImageOps
import numpy as np
import matplotlib.pylab as plt
from astropy.stats import mad_std
from photutils import datasets, daofind, aperture_photometry, CircularAperture
from scipy import stats
import xlwt
from xlutils.copy import copy 
from xlrd import open_workbook
from astropy.stats import sigma_clipped_stats


path1 = "/Users/Marta/Documents/Python/WarnerCell/22.6_27_16Mar2017/22.6_27_pics/"
path2 = "/Users/Marta/Documents/Python/WarnerCell/22.6_27_16Mar2017/22.6_27_grayscale/"
path3 = "/Users/Marta/Documents/Python/WarnerCell/22.6_27_16Mar2017/22.6_27_inverted/"

listing = os.listdir(path1) 

def process_fpath():
    # run the following code once per file inside 'listing'
    # with the full path of the file represented by item
    for item in listing:
        
        print("Checking file grayscale " + item)
        filename1, extension1 = os.path.splitext(path1+item) # separate the filename from the extension
        if extension1 == ".tif":
            print("Processing file grayscale " + item)
            if os.path.isfile(path1+item):
                im = Image.open(path1 + item)
                # "/Users/...../image.png"
                # f = image, e = png
                grayscale_image = PIL.ImageOps.grayscale(im) #convert to greyscale
                # print("Trying to save file to " + path2 + os.path.basename(filename) + '.tif')
                grayscale_image.save(path2 + os.path.basename(filename1) + '.tif') #save output somewhere else
        else:
            print("Skipping file grayscale: " + item)

    listing2 = os.listdir(path2)

    for item in listing2:
        
        print("Checking file inverted " + item)
        filename2, extension2 = os.path.splitext(path2+item) # separate the filename from the extension
        if extension2 == ".tif":
            print("Processing file inverted " + item)
            if os.path.isfile(path2+item):
                im2 = Image.open(path2 + item)
                # "/Users/...../image.png"
                # f = image, e = png
                inverted_image = PIL.ImageOps.invert(im2) #invert grayscale
                # print("Trying to save file to " + path2 + os.path.basename(filename) + '.tif')
                inverted_image.save(path3 + os.path.basename(filename2) + '.tif') #save output somewhere else
        else:
            print("Skipping file inverted: " + item)
    
    numberofnuclei_grayscale = []
    numberofnuclei_inverted = []
    backgroundmedian_grayscale = []
    backgroundmedian_inverted = []
    
    for item in listing2:
        print("Checking file for star counting " + item)
        filename2, extension2 = os.path.splitext(path2+item) # separate the filename from the extension
        if extension2 == ".tif":
            print("Processing file for star counting " + item)
            if os.path.isfile(path2 + item):
                full_image = plt.imread(path2 + item)
                # TODO: make sure it's greyscale first 
                # look up how to remove alpha channel from image?
                image = full_image[:,:].astype(float)

                #subtract median from image
                #median = np.median(image)
                #backgroundmedian_grayscale.append(median)
                #image -= np.median(image)
                
                mean, median, std = sigma_clipped_stats(image, sigma=3.0, iters=5)
                backgroundmedian_grayscale.append(mean)

                #get median absolute deviation
                bkg_sigma = mad_std(image)

                #use function to dectect stars. output (sources) is a table like object
                sources = daofind(image, fwhm=8, threshold=6*bkg_sigma)

                #get x-y coordinates of each "star" or light soruce
                positions = (sources['xcentroid'], sources['ycentroid'])

                #print number of detected sources
                num_stars_grayscale = len(sources['xcentroid'])
                numberofnuclei_grayscale.append(num_stars_grayscale)
                #print("Number of detected sources: %s" %num_stars)

                #define a circule around each source
                #apertures = CircularAperture(positions, r=4.)

                #find brightest source
                #phot_table = aperture_photometry(image, apertures)
                #brightest_source_id = phot_table['aperture_sum'].argmax()

                #plot image with blue circles around each detected "star"
                #plt.figure()
                #plt.imshow(image, cmap='gray_r', origin='lower')
                #apertures.plot(color='blue', lw=1.5, alpha=0.5)
        else:
            print("Skipping file for star counting: " + item)

    print("Number of detected sources from grayscale: %s" %numberofnuclei_grayscale)
    listing3 = os.listdir(path3)    
    
    for item in listing3:
        print("Checking file for star counting " + item)
        filename3, extension3 = os.path.splitext(path3+item) # separate the filename from the extension
        if extension3 == ".tif":
            print("Processing file for star counting " + item)
            if os.path.isfile(path3 + item):
                full_image = plt.imread(path3 + item)
                # TODO: make sure it's greyscale first 
                # look up how to remove alpha channel from image?
                image = full_image[:,:].astype(float)

                #subtract median from image
                #image -= np.median(image)
                
                mean, median, std = sigma_clipped_stats(image, sigma=3.0, iters=5)
                backgroundmedian_inverted.append(std)

                #get median absolute deviation
                bkg_sigma = mad_std(image)

                #use function to dectect stars. output (sources) is a table like object
                sources = daofind(image, fwhm=8, threshold=6*bkg_sigma)

                #get x-y coordinates of each "star" or light soruce
                positions = (sources['xcentroid'], sources['ycentroid'])

                #print number of detected sources
                num_stars_inverted = len(sources['xcentroid'])
                numberofnuclei_inverted.append(num_stars_inverted)
                #print("Number of detected sources: %s" %num_stars)

                #define a circule around each source
                #apertures = CircularAperture(positions, r=4.)

                #find brightest source
                #phot_table = aperture_photometry(image, apertures)
                #brightest_source_id = phot_table['aperture_sum'].argmax()

                #plot image with blue circles around each detected "star"
                #plt.figure()
                #plt.imshow(image, cmap='gray_r', origin='lower')
                #apertures.plot(color='blue', lw=1.5, alpha=0.5)
        else:
            print("Skipping file for star counting: " + item)

    print("Number of detected sources from inverted: %s" %numberofnuclei_inverted)    
    
    sumofboth = [x + y for x, y in zip(numberofnuclei_grayscale, numberofnuclei_inverted)]
    print("Number of detected sources sum: %s" %sumofboth)
    
    sumofboth_corrected = [0.7245*(x + y) - 52.452 for x, y in zip(numberofnuclei_grayscale, numberofnuclei_inverted)]
    
    sample = '22.6_27_linearfit'
    omega = 22.6        
    
    x = 1

    y = 697 #Number of images 
    
    time = 116.167 * 60 #Input time in minutes 

    timeofexp = np.array(np.linspace(x, time, y))    #start from 1, go to amount of time in s, broken into an equal number of pieces
    numberofnuclei_grayscale = np.array(numberofnuclei_grayscale)  #making an array - tried to normalize for size of microscope image,but can't multiply by non-integers here
    numberofnuclei_inverted = np.array(numberofnuclei_inverted)
    sumofboth = np.array(sumofboth)
    sumofboth_corrected = np.array(sumofboth_corrected)
    
    
    #print("Number of images: %s" %numimages)
    
    plt.figure()    #Plotting figure with all three (grayscale, inverted, and sum) nuclei counts vs time
    #fig.scatter(numimages, numimages2, color = 'blue', edgecolor = 'none' )
    #fig.set_aspect(1./ax1.get_data_ratio())
    plt.plot(timeofexp, numberofnuclei_grayscale, 'ro', timeofexp, numberofnuclei_inverted, 'g--', timeofexp, sumofboth, 'bs', timeofexp, sumofboth_corrected, 'ms') 
    plt.xlabel('Time (s)')
    plt.ylabel('Nuclei')    
    plt.show  
    plt.savefig('/Users/Marta/Documents/Python/Plots/22.6_27_fig_corrected_clippedstats.png')
    
    z = 30          #Input time in minutes where you want to slice for the linearfit
    w = 15        #Input time in minutes where you want the slice for the linearfit to start    
    
    Wheretoslice = (y/time)*(z*60) #knowing where to slice the nucleation list - input the time in minutes where you want to slice, and divide by the step unit of time
    Wheretostart = (y/time)*(w*60)    
    linearfitx = timeofexp[Wheretostart:Wheretoslice]    
    lineary = sumofboth_corrected[Wheretostart:Wheretoslice]
    lineary = np.array(lineary)  
    linearfitx = np.array(linearfitx)
    
    plt.figure()
    m, b = np.polyfit(linearfitx, lineary, deg=1)
    plt.plot(linearfitx, lineary,'s',color='m')
    plt.plot(linearfitx, m*linearfitx + b, 'r-')
    plt.plot(timeofexp, sumofboth_corrected, 'ms')
    plt.xlabel('Time (s)')
    plt.ylabel('Number of Nuclei') 
    plt.show
    plt.savefig('/Users/Marta/Documents/Python/Plots/22.6_27_linearfit_corrected_clippedstats.png')
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(linearfitx,lineary)
    print("r-squared:", r_value**2)
    print("slope:", slope)
    print("intercept:", intercept)

    style0 = xlwt.easyxf('font:name Times New Roman, color-index red, bold on')
    #style1 = xlwt.easyxf(num_format_str='#,##0.00')

    #wb = xlwt.Workbook('/Users/Marta/Documents/Python/PhotutilsOutput_Copy.xls')
    w = copy(open_workbook('/Users/Marta/Documents/Python/PhotutilsOutput_Copy.xls'))
    ws = w.add_sheet('Output22.6_27')
    
    ws.write(0,0, 'Sample' , style0)
    ws.write(1,0, sample)
    ws.write(0,1, 'Omega', style0)
    ws.write(1,1, omega) 
    ws.write(0,2, 'RSquared', style0)
    ws.write(1,2, r_value**2)
    ws.write(0,3, 'Slope', style0)
    ws.write(1,3, slope)
    ws.write(0,4, 'Intercept', style0)
    ws.write(1,4, intercept)
    
    
    #wb.save('/Users/Marta/Documents/Python/PhotutilsOutput_Copy.xls')
    w.save('/Users/Marta/Documents/Python/PhotutilsOutput_Copy.xls')        
    
    print("Background Median: %s" %backgroundmedian_grayscale)
    
process_fpath()    
    

#p = Pool(4)
#p.map(process_fpath, listing)
#p.close()
#p.join()
