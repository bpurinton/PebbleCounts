# PebbleCounts - a tool for measuring gravel-bed river grain-size
#       Developed by Ben Purinton (purinton[at]uni-potsdam.de)
#       26 November 2018
#       See the help manual for instructions on use   

# =============================================================================
# Load the modules
# =============================================================================

import os, csv, time, argparse, sys
import cv2
import numpy as np
import PCfunctions as func
from osgeo import gdal
from scipy import ndimage as ndi
from skimage import color
from sklearn import cluster as clust
from skimage import measure as meas
from skimage import segmentation as segm
from skimage import feature as feat
from skimage import morphology as morph
from skimage import filters as filt
import matplotlib.pyplot as plt
# ignore some warnings thrown by sci-kit
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# =============================================================================
# Configure the run with command line arguments
# =============================================================================

# Get the arguments
parser = argparse.ArgumentParser()
parser.add_argument("-im", type=str, 
                    help="The image to use including the full path and extension.")
parser.add_argument("-ortho", type=str, 
                    help="'y' if geo-referenced ortho-image, 'n' if not. Supply input resolution if 'n'.")
parser.add_argument("-input_resolution", type=float, 
                    help="If image is not ortho-image, input the calculated resolution from calculate_camera_resolution.py")
parser.add_argument("-lithologies", type=int, 
                    help="What is the expected number of lithologies with distinct colors? DEFAULT=1", default=1)
parser.add_argument("-maxGS", type=float, 
                    help="Maximum expected longest axis grain size in meters. DEFAULT=0.3", default=0.3)
parser.add_argument("-cutoff", type=int, 
                    help="Cutoff factor in pixels for inclusion of pebble in final count. DEFAULT=9", default=9)
parser.add_argument("-min_sz_factors", nargs='+', type=float, 
                    help="Factors to multiply cutoff value by at each scale. Used to clean-up the masks for easier clicking. The default values are good for ~1 mm/pixel imagery but should be doubled for sub-millimeter or halved for centimeter resolution imagery. DEFAULT=[100, 10, 2]", default=[100, 10, 2])
parser.add_argument("-win_sz_factors", nargs='+', type=float, 
                    help="Factors to multiply maximum grain-size (in pixels) by at each scale. The default values are good for millimeter and sub-millimeter imagery, but should be doubled for coarser centimeter imagery. DEFAULT=[10, 2, 0.5]", default=[10, 2, 0.5])
parser.add_argument("-improvement_ths", nargs='+', type=float, 
                    help="Improvement threshold values for each window scale that tells k-means when to halt. DEFAULT=[0.01, 0.1, 0.1]", default=[0.01, 0.1, 0.1])
parser.add_argument("-coordinate_scales", nargs='+', type=float, 
                    help="Fraction to scale X/Y coordinates by in k-means. DEFAULT=[0.5, 0.5, 0.5]", default=[0.5, 0.5, 0.5])
parser.add_argument("-overlaps", nargs='+', type=float, 
                    help="Fraction of overlap between windows at the different scales. DEFAULT=[0.5, 0.3, 0.1]", default=[0.5, 0.3, 0.1])
parser.add_argument("-nl_means_chroma_filts", nargs='+', type=int, 
                    help="Nonlocal means chromaticity filtering strength for the different scales. DEFAULT=[3, 2, 1]", default=[3, 2, 1])
parser.add_argument("-bilat_filt_szs", nargs='+', type=int, 
                    help="Size of bilateral filtering windows for the different scales. DEFAULT=[9, 5, 3]", default=[9, 5, 3])
parser.add_argument("-tophat_th", type=int, 
                    help="Top percentile threshold to take from tophat filter for edge detection. DEFAULT=90", default=90)
parser.add_argument("-sobel_th", type=int, 
                    help="Top percentile threshold to take from sobel filter for edge detection. DEFAULT=90", default=90)
parser.add_argument("-canny_sig", type=int, 
                    help="Canny filtering sigma value for edge detection. DEFAULT=2", default=2)
parser.add_argument("-resize", type=float, 
                    help="Value to resize windows by should be between 0 and 1. DEFAULT=0.8", default=0.8)
args = parser.parse_args()

# assign the arguments
resize = args.resize
im = args.im
# exit if there was no image supplied or the file doesn't exist
if im == None:
    print("\nSupply and image with the -im command line argument")
    sys.exit()
elif not os.path.exists(im):
    print("\nThe image doesn't exist, check the path and name")
    sys.exit()
else:
    pass
# also exit if there was no ortho option supplied or if there was no input resolution
ortho = args.ortho
if ortho == 'n':
    ortho=False
elif ortho == 'y':
    ortho=True
else:
    print("\nIs the input image an ortho? Use '-ortho y' or '-ortho n'")
    sys.exit()
input_resolution = args.input_resolution
if not ortho and input_resolution == None:
        print("\nSupply an input image resolution with '-input_resolution' calculated with 'calculate_camera_resolution.py'")
        sys.exit()

lithologies = args.lithologies
min_sz_factors = args.min_sz_factors
win_sz_factors = args.win_sz_factors
improvement_ths = args.improvement_ths
coordinate_scales = args.coordinate_scales
cutoff = args.cutoff
maxGS = args.maxGS
overlaps = args.overlaps
nl_means_chroma_hs = args.nl_means_chroma_filts
bilat_filt_szs = args.bilat_filt_szs
tophat_th = args.tophat_th
canny_sig = args.canny_sig
sobel_th = args.sobel_th
# there needs to be three scales
if not len(improvement_ths)==3 or not len(min_sz_factors)==3 or not len(win_sz_factors)==3 or not len(coordinate_scales)==3 or not len(overlaps)==3 or not len(nl_means_chroma_hs)==3 or not len(bilat_filt_szs)==3:
    print("\nPlease supply three values to '-min_sz_factors', '-win_sz_factors', '-improvement_ths', '-coordinate_scales', '-overlaps', '-nl_means_chroma_hs', and '-bilat_filt_szs'.")
    print("Or just let PebbleCounts use their default values.")
    sys.exit()
# output files
base_name = os.path.splitext(os.path.basename(im))[0]
base_dir = os.path.dirname(im)
csv_out = os.path.join(base_dir, base_name + "_PebbleCounts_CSV.csv")
im_out = os.path.join(base_dir, base_name + "_PebbleCounts_LABELS.tif")
fig_out = os.path.join(base_dir, base_name + "_PebbleCounts_FIGURE.png")
# check existence of outputs and warn the user
overwrite=None
if os.path.exists(csv_out) or os.path.exists(im_out) or os.path.exists(fig_out):
    print("\nLooks like some or all of the output files exist already:\n{}\n\n{}\n\n{}\n".format(csv_out, im_out, fig_out))
    while True:
        overwrite = input("Do you want to overwrite these? ('y' or 'n'): ")
        if overwrite=='y' or overwrite=='n':
            break
        else:
            print("incorrect input, should be 'y' or 'n'")
if overwrite=='y':
    print("\nFiles will be overwritten\n")
elif overwrite=='n':
    print("\nEnding the script, move the files to a different directory or rename them to proceed\n")
    sys.exit()
else:
    print("No output files exist for the image yet, proceeding\n")
# open the datset and get the step size if ortho image
if ortho:
    ds = gdal.Open(im)
    gt = ds.GetGeoTransform()
    ds = None
    step = gt[1]
# otherwise take the given input_resolution
if not ortho:
    step = input_resolution
    # convert to meters
    step /= 1000
# set the window sizes and minimum sizes based on the maximum expected GS and resolution
windowSizes = [int(np.round((maxGS/step)*win_sz_factors[0])), int(np.round((maxGS/step)*win_sz_factors[1])), 
               int(np.round((maxGS/step))*win_sz_factors[2])]
min_sizes = [int(cutoff*min_sz_factors[0]), int(cutoff*min_sz_factors[1]), 
                 int(cutoff*min_sz_factors[2])]

# get start time
start = time.time()

# =============================================================================
# Pre-process the entire image (subsetting, shadow masking, color masking)
# =============================================================================

# subset the image?
while True:
    subset = input("subset the image? (y/n): ")
    if subset=='y' or subset=='n':
        break
    else:
        print("incorrect input, should be 'y' or 'n'")
while subset=='y':
    img = cv2.imread(im)
    win_name = "ROI Selector ('spacebar' to end)"
    cv2.startWindowThread()
    cv2.namedWindow(win_name, cv2.WINDOW_NORMAL)
    cv2.moveWindow(win_name, 0, 0)
    cv2.resizeWindow(win_name, func.resizeWin(img, resize)[0], 
                     func.resizeWin(img, resize)[1])
    r = cv2.selectROI(win_name, img, False, False)
    break
cv2.destroyAllWindows()
try:
    if r[2] < 10 or r[3] < 10:
        print("\nBad ROI: {}\nThe image will not be subset".format(str(r)))
        subset='n'
except:
    pass

# read the image
if subset=='y':
    bgr = cv2.imread(im)[int(r[1]):int(r[1]+r[3]), int(r[0]):int(r[0]+r[2])]
else:
    bgr = cv2.imread(im)

# get the X/Y UTM grid if it's an ortho image
if ortho:
    if subset=='y':
        xgrid = func.getXYgrid(im)[0][int(r[1]):int(r[1]+r[3]), int(r[0]):int(r[0]+r[2])]
        ygrid = func.getXYgrid(im)[1][int(r[1]):int(r[1]+r[3]), int(r[0]):int(r[0]+r[2])]
    else:
        xgrid, ygrid = func.getXYgrid(im)

# get the x/y grid (indices)
rgrid = np.arange(0, bgr.shape[0], 1)
cgrid = np.arange(0, bgr.shape[1], 1)
cgrid, rgrid = np.meshgrid(cgrid, rgrid)

# create gray and do otsu thresholding to remove shadow regions
gray = cv2.cvtColor(bgr, cv2.COLOR_BGR2GRAY)
otsu_th, _ = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY+cv2.THRESH_OTSU)
otsu_thresholding = func.otsu_threshold()
while True:
    otsu_thresholding.percent_of_otsu(otsu_th)
    otsu_thresholding.apply_threshold(gray, bgr, otsu_th, resize)
    if otsu_thresholding.thresh != None:
        ignore_mask = gray > otsu_th*(otsu_thresholding.thresh/100)
        break

# instantiate percentage sand for colormask
perc_sand = 0

# interactively select the color for masking
while True:
    do_masking = input("\ncreate a color mask (for sand, vegetation, etc.)? (y/n): ")
    if do_masking=='y':
        # copy the image to prevent modification when adding more masks
        img = bgr.copy()
        # create hsv for colorpicking
        hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV_FULL)
        break
    elif do_masking=='n':
        break
    else:
        print("incorrect input, should be 'y' or 'n'")
while do_masking == 'y':    
    # instantiate coordinate class for storing the clicks
    coords = func.pick_colors()
    # also create a copy of the current shadow mask to pass to the function
    shadow_mask = np.invert(ignore_mask.copy().astype(np.uint8))
    shadow_mask[shadow_mask==254] = 0    
    shadow_mask = np.dstack((shadow_mask, shadow_mask, shadow_mask))
    img = cv2.addWeighted(shadow_mask, 0.6, img, 1, 0)
    # create a window and set the callback function
    win_name = "Color Selector ('q' to close, 'r' to flash image)"
    cv2.namedWindow(win_name, cv2.WINDOW_NORMAL)
    cv2.setMouseCallback(win_name, coords.clicker, param=[img, hsv, resize])
    while cv2.getWindowProperty(win_name, 0) >= 0:
        cv2.imshow(win_name, img)
        cv2.moveWindow(win_name, 0, 0)
        cv2.resizeWindow(win_name, func.resizeWin(img, resize)[0], 
                         func.resizeWin(img, resize)[1])
        k = cv2.waitKey(1)
        if k == ord('q') & 0xFF:
            cv2.destroyAllWindows()
            break
        elif k == ord('r') & 0xFF:
            timeout = time.time() + 0.5
            while time.time() < timeout:
                cv2.namedWindow("Image Overlay", cv2.WINDOW_NORMAL)
                cv2.imshow("Image Overlay", bgr)
                cv2.moveWindow("Image Overlay", 0, 0)
                cv2.resizeWindow("Image Overlay", func.resizeWin(bgr, resize)[0], 
                             func.resizeWin(bgr, resize)[1])
                cv2.waitKey(1)
            cv2.destroyWindow("Image Overlay")
        elif coords.closeWin == True:
            cv2.destroyAllWindows()
            break
    try:
        # apply the clicked range
        lower = coords.lower
        upper = coords.upper
        color_mask = cv2.inRange(hsv, lower, upper)
        color_mask = np.invert(color_mask).astype(bool)
        # clean it up
        color_mask = morph.remove_small_holes(color_mask, area_threshold=10, connectivity=2)
        color_mask = morph.opening(color_mask, selem=morph.selem.disk(1))
        color_mask = morph.closing(color_mask, selem=morph.selem.disk(1))
        # add the hue mask to the full ignore mask
        ignore_mask = np.logical_and(ignore_mask, color_mask)
        # the "color mask" represents percentage of pixels that are sand
        perc_sand = np.sum(np.invert(color_mask.astype(bool)))/color_mask.size
        # make the color mask three channel image for stacking
        color_mask = color_mask.astype(np.uint8)
        color_mask[color_mask==0] = 255
        color_mask[color_mask==1] = 0
        color_mask = np.dstack((color_mask, color_mask, color_mask))
        # stack the color mask for subsequent runs
        img = cv2.addWeighted(color_mask, 0.6, img, 1, 0)
        # do we repeat the color masking?
        while True:
            do_masking = input("add another color mask? (y/n): ")
            if do_masking=='y' or do_masking=='n':
                break
            else:
                print("incorrect input, should be 'y' or 'n'")
    except:
        while True:
            do_masking = input("add another color mask? (y/n): ")
            if do_masking=='y' or do_masking=='n':
                break
            else:
                print("incorrect input, should be 'y' or 'n'")

# instantiate the empty grains to fill with all the region props below
grains = []
all_labels = np.ones(bgr.shape[0:2])*0

# =============================================================================
# Run through the windows at the different scales
# =============================================================================

print("\nBeginning k-means segmentation\n")

# loop over all window sizes
for index in range(len(windowSizes)):
    improvement_th = improvement_ths[index]
    windowSize = int(windowSizes[index])
    scale = coordinate_scales[index]
    min_size = int(min_sizes[index])
    color_filt = nl_means_chroma_hs[index]
    filt_sz = bilat_filt_szs[index]
    overlap = overlaps[index]
    
    print("\nScale {:d} of {:d}\n".format(index+1, len(windowSizes)))
    
    # choose the windows based based on chosen window size, step size is taken as % of window size
    windows = []
    for (x, y, sz) in func.sliding_window(bgr[:,:,0], np.int16(windowSize*(1-overlap)), windowSize):
        windows.append((x, y, sz[0], sz[1]))
    
    # now loop through the windows
    for winNumb, window in enumerate(windows):
        print("\nWindow {:d} of {:d}\n".format(winNumb+1, len(windows)))
        
        ulx, uly, lrx, lry = window[0], window[1], window[0]+window[2], window[1]+window[3]
        BGR = bgr[ulx:lrx, uly:lry]
        # pass if empty
        if BGR.size == 0:
            print("\nEmpty window, skipping\n")
            continue
        RGRID = rgrid[ulx:lrx, uly:lry]
        CGRID = cgrid[ulx:lrx, uly:lry]
        MASK = ignore_mask[ulx:lrx, uly:lry]
        GRAY = gray[ulx:lrx, uly:lry]
        
        # do non-local means for denoising
        print("\nNon-local means filtering")
        BGR = cv2.fastNlMeansDenoisingColored(BGR, None, color_filt, 1, 3, 9)
        
        # get the hsv from the de-noised image for lithology estimation
        HSV = cv2.cvtColor(BGR, cv2.COLOR_BGR2HSV_FULL)
        hue = HSV[:,:,0]
        sat = HSV[:,:,1]
        
        print("Bilateral filtering")
        # bilateral filter (preserved edges) in opencv on CIELab
        LAB = cv2.cvtColor(BGR, cv2.COLOR_BGR2Lab)
        a_blur = cv2.bilateralFilter(LAB[:,:,1], filt_sz, 75, 75)
        b_blur = cv2.bilateralFilter(LAB[:,:,2], filt_sz, 75, 75)
        
        # tophat edges
        print("Black tophat edge detection")
        tophat = morph.black_tophat(GRAY, selem=morph.selem.disk(1))
        tophat = tophat < np.percentile(tophat, tophat_th)
        tophat = morph.remove_small_holes(tophat, area_threshold=cutoff, connectivity=2)
        foo = func.featAND_fast(MASK, tophat)
        MASK = np.logical_and(foo, MASK)
        # canny edges
        print("Canny edge detection")
        canny = feat.canny(GRAY, sigma=canny_sig)
        canny = np.invert(canny)
        foo = func.featAND_fast(MASK, canny)
        MASK = np.logical_and(foo, MASK)
        # sobel edges
        print("Sobel edge detection")
        sobel = filt.sobel(GRAY)
        sobel = sobel < np.percentile(sobel, sobel_th)
        sobel = morph.remove_small_holes(sobel, area_threshold=cutoff, connectivity=2)
        sobel = morph.thin(np.invert(sobel))
        sobel = np.invert(sobel)
        foo = func.featAND_fast(MASK, sobel)
        MASK = np.logical_and(foo, MASK)
        
        # find the remaining pixels in the mask
        idx = np.where(MASK == True)
        
        # skip if there's only a small number of pixels left
        # as this will lead to errors if the number of k-means clusters
        # becomes greater than the number of pixels
        if len(idx[0]) < 100:
            print("\nEmpty window, skipping\n")
            continue
        
        # get X/Y vectors
        rgrid_ = RGRID[idx]
        cgrid_ = CGRID[idx]        
        
        # get the color vectors
        a_blur_ = a_blur[idx]
        b_blur_ = b_blur[idx]        
        
        # rescale color between 0 and 1 and X/Y between 0 and scaling factor
        a_blur_ = (((a_blur_ - a_blur_.min()) / (a_blur_.max() - a_blur_.min()))*1)
        b_blur_ = (((b_blur_ - b_blur_.min()) / (b_blur_.max() - b_blur_.min()))*1)
        cgrid_scaled = (((cgrid_ - cgrid_.min()) / (cgrid_.max() - cgrid_.min()))*scale)
        rgrid_scaled = (((rgrid_ - rgrid_.min()) / (rgrid_.max() - rgrid_.min()))*scale)
        
        # create kmeans vector
        X = np.column_stack((a_blur_.reshape(-1, 1), b_blur_.reshape(-1, 1),
                             rgrid_scaled.reshape(-1, 1), cgrid_scaled.reshape(-1, 1)))
        
        # run kmeans
        print("Running k-means")
        # dummy variables for looping
        inertias = [0]
        iteration = 1
        n_clusters = 1
        # run it once
        kmean=clust.MiniBatchKMeans(n_clusters=n_clusters, batch_size=1000).fit(X)
        inertias.append(kmean.inertia_)
        # stop when new inertia is < improvement threshold
        while abs(inertias[iteration]-inertias[iteration-1]) > inertias[iteration-1]*improvement_th:
            n_clusters += 1
            kmean=clust.MiniBatchKMeans(n_clusters=n_clusters, batch_size=1000).fit(X)
            inertia = kmean.inertia_
            centroids = kmean.cluster_centers_
            print("Current number of clusters: {:d}, total inertia: {:0.3f}".format(n_clusters, inertia))
            inertias.append(inertia)
            iteration += 1

        # get the labels
        labels = kmean.labels_
        # reshape into image
        rc_reduced = list(zip(rgrid_ - ulx, cgrid_ - uly, labels))
        im_cluster = np.ones(GRAY.shape)*np.nan
        for r, c, l in rc_reduced:
            im_cluster[r, c] = l
        
        # add one for indexing
        im_cluster += 1

        # loop over masks and clean them up, creating a master list
        print("Cleaning up k-means mask")
        master_mask = np.zeros((im_cluster.shape[0], im_cluster.shape[1], 3)).astype(np.uint8)
        for i in range(1, n_clusters+1):
            # pull out each cluster as mask
            mask = np.invert(im_cluster.copy().astype(bool))
            mask[im_cluster == i] = True
            mask[mask != True] = False
            # binary operations to clean the mask a bit
            mask = morph.remove_small_objects(mask, min_size=min_size, connectivity=2)
            mask = morph.erosion(mask, selem=morph.selem.square(3))
            mask = morph.dilation(mask, selem=morph.selem.square(2))
            mask = segm.clear_border(mask)
            mask = morph.remove_small_objects(mask, min_size=min_size, connectivity=2)
            # make sure we didn't add any pixels back in from the original mask
            mask[MASK==False] = False
            # label the mask and give each region a random rgb color
            color_labels, _ = ndi.label(mask)
            color_choice = [list(np.random.choice(range(255), size=3)) for x in range(40)]
            mask_color = color.label2rgb(color_labels, colors=color_choice, bg_label=0, bg_color=[0, 0, 0])
            # add to master mask
            master_mask = master_mask + mask_color.astype(np.uint8)
        
        # skip the master mask if it is empty after cleaning
        idx = np.where(master_mask != [0, 0, 0])
        if len(idx[0]) == 0:
            print("\nEmpty window, skipping\n")
            continue

        # click seed points with OpenCV
        img = master_mask.copy()    
        # overlay the color grain mask and original image
        img = cv2.addWeighted(img, 1, BGR, 0.5, 0)
        # instantiate coordinate class for storing the clicks
        coords = func.select_grains()
        # create a window, call it open and click through it
        win_name = "KMeans ('r' see image, 'q' close)"
        cv2.namedWindow(win_name, cv2.WINDOW_NORMAL)
        cv2.setMouseCallback(win_name, coords.clicker, param=[img])
        while cv2.getWindowProperty(win_name, 0) >= 0:
            cv2.imshow(win_name, img)
            cv2.moveWindow(win_name, 0, 0)
            cv2.resizeWindow(win_name, func.resizeWin(img, resize)[0], 
                             func.resizeWin(img, resize)[1])
            k = cv2.waitKey(1)
            # create a call to overlapping window of RGB image if 'r' is pressed
            if k == ord('r') & 0xFF:
                while cv2.getWindowProperty(win_name, 0) >= 0:
                    cv2.namedWindow("Image Overlay, 'r' close", cv2.WINDOW_NORMAL)
                    cv2.imshow("Image Overlay, 'r' close", BGR)
                    cv2.moveWindow("Image Overlay, 'r' close", 0, 0)
                    cv2.resizeWindow("Image Overlay, 'r' close", func.resizeWin(BGR, resize)[0], 
                                     func.resizeWin(BGR, resize)[1])
                    l = cv2.waitKey(1)
                    if l == ord('r') & 0xFF:
                        cv2.destroyWindow("Image Overlay, 'r' close")
                        break
            # close all windows and go to next mask
            if k == ord('q') & 0xFF:
                cv2.destroyAllWindows()
                break

        # get the region properties of clicked areas
        print("Getting properties of clicked grains")
        master_mask = master_mask[:,:,0]
        if len(coords.clicks) != 0:
            master_mask, _ = ndi.label(master_mask)
            labels = np.zeros(master_mask.shape).astype(np.uint8) 
            for click in coords.clicks:
                if not master_mask[click] == 0:
                    labels[master_mask==master_mask[click]] = 1
            labels, _ = ndi.label(labels)      
            for grain in meas.regionprops(labels, coordinates='xy'):
                # only append if the grain is greater than the cutoff size
                if grain.minor_axis_length < float(cutoff) or grain.major_axis_length < float(cutoff):
                    continue
                else:
                    # get the average color from HSV color space in the grain mask
                    Havg = np.mean(np.asarray([hue[tuple(idx)] for idx in grain.coords]))
                    Savg = np.mean(np.asarray([sat[tuple(idx)] for idx in grain.coords]))
                    # append the grain and color parameters
                    grains.append((grain.centroid[0]+ulx, grain.centroid[1]+uly, 
                               grain.minor_axis_length, grain.major_axis_length, 
                               grain.orientation, grain.area, Havg, Savg))

            # add the chosen grains to the ignore mask
            labels[labels != 0] = 1
            labels = labels.astype(bool)
            ignore_mask[ulx:lrx, uly:lry] = np.logical_and(ignore_mask[ulx:lrx, uly:lry], 
                                                           np.invert(labels))
            
            # also add it to an "all labels" mask if there are any labels in it
            if np.sum(labels) != 0:
                all_labels[ulx:lrx, uly:lry] = all_labels[ulx:lrx, uly:lry] + labels
        
        else:
            pass

# =============================================================================
# Get the rock colors based on expected number of lithologies
# =============================================================================

# get the lithology labels using kmeans on the colors
hues, sats = [], []
for grain in grains:    
    hue, sat = grain[-2], grain[-1]
    hues.append(hue)
    sats.append(sat)
hues = np.asarray(hues)
sats = np.asarray(sats)
# re-run kmeans separation on the color channels
X = np.column_stack((hues.reshape(-1, 1), sats.reshape(-1, 1)))
kmean=clust.KMeans(n_clusters=lithologies).fit(X)
color_labels = kmean.labels_
# order colors from low to high values with lookup table
# so the numbering is consistent for different images from the
# same field site
idx = np.argsort(kmean.cluster_centers_.sum(axis=1))
lut = np.zeros_like(idx)
lut[idx] = np.arange(lithologies)
color_labels = lut[color_labels]

# =============================================================================
# Output a final plot, .csv of grain-size data, and label mask
# =============================================================================

# make a plot of the fit grains 
plt.figure(figsize=(10,10))
plt.imshow(cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB))
for kmeans_label, grain in enumerate(grains):    
    litho = color_labels[kmeans_label]
    y0, x0 = grain[0], grain[1]
    a, b = grain[3], grain[2]
    orientation = grain[4]
    x1 = x0 + np.cos(orientation) * .5 * a
    y1 = y0 - np.sin(orientation) * .5 * a
    x2 = x0 - np.sin(orientation) * .5 * b
    y2 = y0 - np.cos(orientation) * .5 * b
    # also plot the ellipse
    phi = np.linspace(0,2*np.pi,50)
    x = x0 + a/2 * np.cos(phi) * np.cos(-orientation) - b/2 * np.sin(phi) * np.sin(-orientation)
    y = y0 + a/2 * np.cos(phi) * np.sin(-orientation) + b/2 * np.sin(phi) * np.cos(-orientation)
    plt.plot((x0, x1), (y0, y1), '-r', linewidth=1)
    plt.plot((x0, x2), (y0, y2), '-r', linewidth=1)
    plt.plot(x0, y0, '.g', markersize=2)
    plt.plot(x, y, 'r--', linewidth=0.7)
    plt.text(x0, y0, "{:d}".format(litho), color='m', fontsize=7)
    plt.axis('off')
plt.savefig(fig_out, dpi=300)
plt.close()

# what is the percent of the image not measured (fines or unfound rocks)
perc_nongrain = np.sum(np.invert(ignore_mask.astype(bool)))/ignore_mask.size
# output the measured grains as a csv
with open(csv_out, "w") as csv_file:
    writer=csv.writer(csv_file, delimiter=",",lineterminator="\n",)
    if ortho:
        writer.writerow(["perc. not meas.", "perc. background color", 
                         "UTM X (m)", "UTM Y (m)", "a (px)", "b (px)", 
                         "a (m)", "b (m)", "area (px)", "area (m2)", 
                         "orientation", "lithology"])
    if not ortho:
        writer.writerow(["perc. not meas.", "perc. background color", 
                         "a (px)", "b (px)", "a (m)", "b (m)", 
                         "area (px)", "area (m2)", 
                         "orientation", "lithology"])
    
    for kmeans_label, grain in enumerate(grains):
        litho = color_labels[kmeans_label]
        y0, x0 = grain[0], grain[1]
        a, b = grain[3], grain[2]
        orientation = grain[4]
        area = grain[5]
        
        if ortho:
            x_coord = xgrid[np.round(y0).astype(int), np.round(x0).astype(int)]
            y_coord = ygrid[np.round(y0).astype(int), np.round(x0).astype(int)]
            writer.writerow([perc_nongrain, perc_sand, x_coord, y_coord, a, b, 
                             a*step, b*step, area, area*step**2, orientation, litho])
    
        if not ortho:
            writer.writerow([perc_nongrain, perc_sand, a, b, 
                             a*step, b*step, area, area*step**2, orientation, litho])
                
    csv_file.close()

# save out as raster or image
labels, _ = ndi.label(all_labels)
if ortho:
    func.array2rast(labels, im, im_out, xgrid, ygrid)
if not ortho:
    labels = (color.label2rgb(labels, bg_label=0, bg_color=[1, 1, 1])*255).astype(np.uint8)
    cv2.imwrite(im_out, labels)
    
# get end time    
end = time.time()
print("\nThat took about {:.0f} minutes, you counted {:d} pebbles!\n".format(end/60-start/60, len(grains)))