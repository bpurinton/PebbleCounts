#!/usr/bin/zsh
###############################################################################
# buildOpenCV.zsh -                                                           #
#                                                                             #
# Script to build OpenCV 3.4.2 from source into an Anaconda virtual           #
# environment				                                      #
###############################################################################

###############################################################################
# SOURCE: https://gist.github.com/clamytoe/3424d0201bba0073ff1313e80ffc6328   #
# Modified by Ben Purinton, December 2018				      #
###############################################################################

# Install requirements
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install build-essential cmake unzip pkg-config
sudo apt-get install libjpeg-dev libpng-dev libtiff-dev
sudo apt-get install libavcodec-dev libavformat-dev libswscale-dev libv4l-dev
sudo apt-get install libxvidcore-dev libx264-dev
sudo apt-get install libgtk-3-dev
sudo apt-get install libatlas-base-dev gfortran
sudo apt-get install python3-dev

# Move to Downloads folder
cd ~/Downloads
wget -O opencv.zip -c https://github.com/opencv/opencv/archive/3.4.2.zip
wget -O opencv_contrib.zip -c https://github.com/opencv/opencv_contrib/archive/3.4.2.zip
unzip opencv.zip
unzip opencv_contrib.zip

# Create virtual environment
conda create -y --name pebblecounts numpy
conda activate pebblecounts

# Generate make files and install
cd opencv-3.4.2
mkdir build
cd build
cmake -D CMAKE_BUILD_TYPE=RELEASE -D CMAKE_INSTALL_PREFIX=~/miniconda3/envs/pebblecounts -D INSTALL_PYTHON_EXAMPLES=ON -D INSTALL_C_EXAMPLES=OFF -D OPENCV_EXTRA_MODULES_PATH=~/opencv_contrib-3.4.2/modules -D PYTHON_EXECUTABLE=~/miniconda3/envs/pebblecounts/bin/python -D WITH_GTK=ON -D WITH_OPENGL=ON -D WITH_QT=ON -D WITH_TBB=ON -D WITH_V4L=ON -D BUILD_EXAMPLES=ON ..
make
make install
ldconfig -n ~/miniconda3/envs/pebblecounts/lib

# Verify that it worked
python -c "import cv2; print(cv2.__version__)"

# Cleanup
cd ~/Downloads
rm opencv.zip opencv_contrib.zip
rm -rf opencv-3.4.2 opencv_contrib-3.4.2
cd

# Deactivate virtual environment
conda deactivate
