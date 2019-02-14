# Install requirements
sudo apt-get -y update
sudo apt-get -y upgrade
sudo apt-get -y install build-essential cmake unzip pkg-config
sudo apt-get -y install libjpeg-dev libpng-dev libtiff-dev
sudo apt-get -y install libavcodec-dev libavformat-dev libswscale-dev libv4l-dev
sudo apt-get -y install libxvidcore-dev libx264-dev
sudo apt-get -y install libgtk-3-dev
sudo apt-get -y install libatlas-base-dev gfortran
sudo apt-get -y install python3-dev

# Move to Downloads folder and get OpenCV source code from GITHUB
cd ~/Downloads
wget -O opencv.zip -c https://github.com/opencv/opencv/archive/3.4.2.zip
wget -O opencv_contrib.zip -c https://github.com/opencv/opencv_contrib/archive/3.4.2.zip
unzip opencv.zip
unzip opencv_contrib.zip

# Create virtual environment called "pebblecounts" with packages (besides OpenCV) installed and activate the environment
conda create -y --name pebblecounts python=3.6.6 numpy scikit-image scikit-learn gdal scipy matplotlib tk
conda activate pebblecounts

# configure opencv with cmake
cd opencv-3.4.2
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=RELEASE \
    -DCMAKE_INSTALL_PREFIX=~/miniconda3/envs/pebblecounts \
    -DINSTALL_PYTHON_EXAMPLES=ON \
    -DINSTALL_C_EXAMPLES=OFF \
    -DOPENCV_EXTRA_MODULES_PATH=~/Downloads/opencv_contrib-3.4.2/modules \
    -DPYTHON_EXECUTABLE=~/miniconda3/envs/pebblecounts/bin/python \
    -DBUILD_EXAMPLES=ON ..

# make it! this command might take 15-20 minutes
make
make install

# Verify that it worked
python -c "import cv2; print(cv2.__version__)"

# if this doesn't work, a clean install is done with "conda env remove -n pebblecounts" and "rm -rf build"

# Cleanup
cd ~/Downloads
rm -rf opencv*
cd

# Deactivate virtual environment
source deactivate
