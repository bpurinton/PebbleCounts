# Calculate the approximate image resolution in mm/pixel for plugging into PebbleCounts
#       Developed by Ben Purinton (purinton[at]uni-potsdam.de)
#       26 November 2018
# 
# Usage:
#       python calculate_camera_resolution.py -focal [focal length in mm] \
#                                               -height [caputre height in m] \
#                                               -sensorHW [sensor height in mm] [sensor width in mm]
#                                               -imageHW [image height in pixels] [image width in pixels]

import numpy as np
import argparse

# Get the arguments
parser = argparse.ArgumentParser()
parser.add_argument("-focal", type=float, 
                    help="Camera focal length in millimeters", default=35)
parser.add_argument("-height", type=float, 
                    help="Photo capture height in meters", default=1.6)
parser.add_argument("-sensorHW", nargs='+', type=float,
                    help="The height and width of the internal camera sensor in millimeters", default=[15.6, 23.5])
parser.add_argument("-imageHW", nargs='+', type=int,
                    help="The height and width of the photography in pixels", default=[4000, 6000])
args = parser.parse_args()
    

def calculate_camera_res(focal_length_mm, height_m, sensorH_mm=15.6, sensorW_mm=23.5, 
                         pixelsH=4000, pixelsW=6000):
    fovH_m = (sensorH_mm/focal_length_mm)*height_m
    fovW_m = (sensorW_mm/focal_length_mm)*height_m
    print("\nThe field of view is {:0.2f} by {:0.2f} m\n".format(fovH_m, fovW_m))
    return np.round(fovH_m*1000/pixelsH, 4), np.round(fovW_m*1000/pixelsW, 4) 

try:
    print("\nFocal length {:0.2f} mm; Shot from {:0.2f} m; Sensor size ({:0.2f}, {:0.2f}) mm; Image size ({:d}, {:d}) pixels:".format(args.focal, args.height, args.sensorHW[0], 
                                        args.sensorHW[1], args.imageHW[0], args.imageHW[1]))
    X_res, Y_res = calculate_camera_res(args.focal, args.height, args.sensorHW[0], 
                                        args.sensorHW[1], args.imageHW[0], args.imageHW[1])
    average_res = np.round(np.mean([X_res, Y_res]), 4)
    print("approximate (x,y) resolution in mm/pixel = ({:0.4f}, {:0.4f})".format(X_res, Y_res))
    print("average resolution in mm/pixel = {:0.4f}\n".format(average_res))
except:
    print("\nOops that didn't work try 'python calculate_camera_resolution.py --help' and check your inputs!")