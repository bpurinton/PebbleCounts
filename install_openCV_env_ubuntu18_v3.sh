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
wget -O opencv.zip -c https://github.com/opencv/opencv/archive/3.4.4.zip
wget -O opencv_contrib.zip -c https://github.com/opencv/opencv_contrib/archive/3.4.4.zip
unzip opencv.zip
unzip opencv_contrib.zip
mv opencv-3.4.4 opencv
mv opencv_contrib-3.4.4 opencv_contrib

# install pip
wget https://bootstrap.pypa.io/get-pip.py
sudo python3 get-pip.py

# install virtualenv and virtualenvwrapper
sudo pip install virtualenv virtualenvwrapper
sudo rm -rf ~/Downloads/get-pip.py ~/Downloads/.cache/pip

# add lines to bashrc
echo -e "\n# virtualenv and virtualenvwrapper" >> ~/.bashrc
echo "export WORKON_HOME=$HOME/.virtualenvs" >> ~/.bashrc
echo "export VIRTUALENVWRAPPER_PYTHON=/usr/bin/python3" >> ~/.bashrc
echo "source /usr/local/bin/virtualenvwrapper.sh" >> ~/.bashrc
source ~/.bashrc

# create virtual environment for pebblecounts
mkvirtualenv pebblecounts -p python3
deactivate
workon pebblecounts
pip install numpy scikit-image scikit-learn scipy matplotlib
# GDAL!!!
# 

# configure opencv with cmake
cd opencv
mkdir build
cd build
cmake -D CMAKE_BUILD_TYPE=RELEASE \
	-D CMAKE_INSTALL_PREFIX=/usr/local \
	-D INSTALL_PYTHON_EXAMPLES=ON \
	-D INSTALL_C_EXAMPLES=OFF \
	-D OPENCV_ENABLE_NONFREE=ON \
	-D OPENCV_EXTRA_MODULES_PATH=~/Downloads/opencv_contrib/modules \
	-D PYTHON_EXECUTABLE=~/.virtualenvs/pebblecounts/bin/python3 \
	-D BUILD_EXAMPLES=ON ..
make
sudo make install
sudo ldconfig

# finishing
cd /usr/lib/python3/dist-packages
sudo mv cv2.cpython-36m-x86_64-linux-gnu.so cv2.so
cd ~/.virtualenvs/pebblecounts/lib/python3.6/site-packages/
ln -s /usr/lib/python3/dist-packages/cv2.so cv2.so

# Verify that it worked
python -c "import cv2; print(cv2.__version__)"

# Cleanup
rm -rf ~/Downloads/opencv*
