# PebbleCounts GUI - Installation Guide

This guide helps you install PebbleCounts GUI with different configurations depending on your needs.

## Quick Decision Tree

**Do you need to process georeferenced ortho-images (GeoTIFFs)?**
- **NO** → Use [Basic Installation](#basic-installation-without-gdal) (easier, no GDAL needed)
- **YES** → Use [Full Installation with GDAL](#full-installation-with-gdal) (requires conda)

---

## Basic Installation (Without GDAL)

**Use this if you only process non-orthorectified images (JPG/PNG).**

### Step 1: Create Virtual Environment
```bash
cd PebbleCounts
python3 -m venv pebblecounts-gui
source pebblecounts-gui/bin/activate  # Mac/Linux
# OR
pebblecounts-gui\Scripts\activate  # Windows
```

### Step 2: Install Basic Requirements
```bash
pip install -r requirements-basic.txt
```

### Step 3: Launch GUI
```bash
python pebblecounts_gui.py
```

**That's it!** The GUI will work for all non-georeferenced images.

---

## Full Installation with GDAL

**Use this if you need to process georeferenced ortho-images (GeoTIFFs).**

GDAL is notoriously difficult to install via pip. The recommended approach is to use **conda**, which handles the C library dependencies automatically.

### Option A: Conda Installation (RECOMMENDED)

```bash
# Create conda environment with Python and GDAL
conda create -n pebblecounts-gui python=3.9 gdal shapely -c conda-forge

# Activate environment
conda activate pebblecounts-gui

# Install remaining dependencies
pip install numpy scipy opencv-python scikit-image scikit-learn matplotlib Pillow gradio

# Launch GUI
python pebblecounts_gui.py
```

### Option B: Homebrew (macOS Only)

```bash
# Install GDAL via Homebrew
brew install gdal

# Create virtual environment
python3 -m venv pebblecounts-gui
source pebblecounts-gui/bin/activate

# Install Python bindings matching your GDAL version
GDAL_VERSION=$(gdal-config --version)
pip install gdal==$GDAL_VERSION

# Install other dependencies
pip install -r requirements-basic.txt

# Launch GUI
python pebblecounts_gui.py
```

### Option C: Ubuntu/Debian

```bash
# Install GDAL system libraries
sudo apt-get update
sudo apt-get install gdal-bin libgdal-dev

# Create virtual environment
python3 -m venv pebblecounts-gui
source pebblecounts-gui/bin/activate

# Install Python bindings matching your GDAL version
GDAL_VERSION=$(gdal-config --version)
pip install gdal==$GDAL_VERSION

# Install other dependencies
pip install -r requirements-basic.txt

# Launch GUI
python pebblecounts_gui.py
```

---

## Troubleshooting

### "GDAL version mismatch" Error

This happens when the Python GDAL bindings don't match your system GDAL library.

**Solution:** Install the matching version:
```bash
# Find your system GDAL version
gdal-config --version  # Example output: 3.7.1

# Install matching Python bindings
pip install gdal==3.7.1
```

### "Can't find gdal-config" Error

This means GDAL isn't installed on your system.

**Solution:** Use conda to install GDAL first:
```bash
conda install -c conda-forge gdal
```

### "SSL module not available" (Conda)

This is a conda/OpenSSL issue on macOS.

**Solution:** Use pip + venv instead of conda (see [Basic Installation](#basic-installation-without-gdal))

### "ImportError: No module named 'osgeo'"

This means GDAL Python bindings aren't installed.

**Solution:**
```bash
# If using conda:
conda install -c conda-forge gdal

# If using pip (after installing system GDAL):
pip install gdal==$(gdal-config --version)
```

---

## Testing Your Installation

### Test Basic GUI (No GDAL)
```bash
python -c "import gradio; import cv2; import numpy; import sklearn; print('✓ Basic installation OK')"
```

### Test Full Installation (With GDAL)
```bash
python -c "from osgeo import gdal; print('✓ GDAL available')"
```

### Launch Test
```bash
python pebblecounts_gui.py
```

If the web interface opens at `http://127.0.0.1:7860`, you're good to go!

---

## Which Dependencies Do I Actually Need?

| Dependency | Required For | Can Skip If... |
|------------|--------------|----------------|
| numpy, scipy, scikit-learn | All processing | Never |
| opencv-python | Image I/O, GUI | Never |
| scikit-image | Segmentation, filtering | Never |
| matplotlib | Output figures | Never |
| Pillow | Image handling | Never |
| shapely | Geometry operations | Never |
| gradio | Web GUI | Using CLI only |
| **gdal** | **Georeferenced images** | **Only processing non-ortho images** |

---

## Recommended Setup for Different Users

### Student/Beginner
```bash
# Just use basic installation - works for 90% of use cases
pip install -r requirements-basic.txt
python pebblecounts_gui.py
```

### Researcher (Need GeoTIFF Support)
```bash
# Use conda for hassle-free GDAL
conda create -n pebblecounts-gui python=3.9 gdal -c conda-forge
conda activate pebblecounts-gui
pip install numpy scipy opencv-python scikit-image scikit-learn matplotlib Pillow gradio shapely
```

### Developer/Advanced User
```bash
# Use existing environment, install selectively
pip install -r requirements-basic.txt
# Install GDAL only if needed via your preferred method
```

---

## Still Having Issues?

1. **Check your Python version**: `python --version` (Need 3.7+)
2. **Try in a fresh environment**: Delete old venv and recreate
3. **Use conda for GDAL**: It's much more reliable than pip
4. **Skip GDAL if possible**: Most users don't need georeferenced image support

For more help, see:
- [PebbleCounts Manual](docs/PebbleCounts_Manual.pdf)
- [GUI User Guide](GUI_GUIDE.md)
- [GitHub Issues](https://github.com/bpurinton/PebbleCounts/issues)
