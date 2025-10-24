# PebbleCounts GUI - Conda Setup Guide

## Quick Start (If Your Conda Works)

```bash
cd PebbleCounts
conda env create -f environment.yml
conda activate pebblecounts-gui
python pebblecounts_gui.py
```

That's it! If the above works, you're done. If you get an SSL error, continue below.

---

## Debugging the SSL Error

You're seeing this error:
```
ImportError: Can't connect to HTTPS URL because the SSL module is not available.
```

This means your conda installation's Python was built without SSL support. This is NOT a package problem—it's a conda installation issue.

### Check Your Conda's SSL Module

```bash
# Test if SSL module works
/Users/ben/miniconda3/bin/python -c "import ssl; print(ssl.OPENSSL_VERSION)"
```

**If this fails**, your conda needs to be reinstalled or repaired.

### Check OpenSSL Libraries

```bash
# Check if OpenSSL is available on your system
which openssl
openssl version

# Check if conda's Python can find SSL libraries
otool -L /Users/ben/miniconda3/bin/python | grep ssl
```

---

## Solutions (In Order of Preference)

### Solution 1: Reinstall Miniconda (RECOMMENDED)

Your conda installation appears corrupted. A fresh install is the cleanest fix.

```bash
# 1. Backup your conda environments list (optional)
conda env list > ~/conda_environments_backup.txt

# 2. Remove old miniconda
rm -rf ~/miniconda3

# 3. Download latest Miniconda for macOS
# For Intel Mac:
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh

# For Apple Silicon (M1/M2/M3):
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh

# 4. Install (follow prompts)
bash Miniconda3-latest-MacOSX-*.sh

# 5. Restart terminal or source
source ~/.bash_profile  # or ~/.zshrc

# 6. Verify SSL works
conda --version
python -c "import ssl; print('SSL OK')"

# 7. Create PebbleCounts environment
cd PebbleCounts
conda env create -f environment.yml
conda activate pebblecounts-gui
```

### Solution 2: Repair Conda's SSL

If you don't want to reinstall, try repairing:

```bash
# Update conda base environment
conda update -n base conda --yes

# Reinstall OpenSSL in base
conda install -n base openssl --yes

# Test again
python -c "import ssl; print('SSL OK')"
```

### Solution 3: Use Mamba Instead of Conda

Mamba is a faster, more reliable conda replacement:

```bash
# Install mamba in base environment
conda install mamba -n base -c conda-forge

# Use mamba instead of conda
mamba env create -f environment.yml
mamba activate pebblecounts-gui
```

### Solution 4: Bypass Conda Notices (Temporary Workaround)

This won't fix the underlying issue but lets you create environments:

```bash
# Disable the notice system that's failing
export CONDA_NO_PLUGINS=true

# Or use offline mode
conda env create -f environment.yml --offline

# Or use --no-plugins flag
conda env create -f environment.yml --no-plugins
```

---

## What Caused This?

Common causes on macOS:

1. **macOS System Update**: You're on macOS 15.5 (Sequoia). System updates can break conda's library links.

2. **Outdated Conda**: Your conda is version 23.7.2 (from August 2023). There have been several SSL-related fixes since then.

3. **Conflicting Python Installations**: You have multiple Python-related paths in your PATH:
   - `/Users/ben/miniconda3`
   - `/Users/ben/.rbenv/shims` (Ruby environment manager)
   - `/Users/ben/.nvm` (Node version manager)
   - System Python

   These can sometimes interfere with each other.

---

## After Reinstalling Conda

Once you have a fresh conda installation with working SSL:

```bash
# Navigate to PebbleCounts
cd /Users/ben/Library/CloudStorage/Dropbox/GITHUB/PebbleCounts

# Create environment from yml file
conda env create -f environment.yml

# Activate
conda activate pebblecounts-gui

# Verify everything installed
python -c "from osgeo import gdal; import gradio; print('✓ All packages OK')"

# Launch GUI
python pebblecounts_gui.py
```

The GUI will open at http://127.0.0.1:7860

---

## Environment Details

The `environment.yml` file installs:

- **Python 3.9** (stable, well-tested)
- **Scientific stack**: numpy, scipy, scikit-learn, scikit-image
- **Image processing**: OpenCV, Pillow
- **Visualization**: matplotlib
- **Geospatial**: GDAL, shapely (conda handles C dependencies)
- **GUI**: Gradio (installed via pip within conda env)

All dependencies are pinned to major versions for stability.

---

## Verifying Your Installation

### Step-by-step verification:

```bash
# 1. Check conda works
conda --version

# 2. Check SSL module
python -c "import ssl; print(ssl.OPENSSL_VERSION)"

# 3. Create environment
conda env create -f environment.yml

# 4. Activate
conda activate pebblecounts-gui

# 5. Check GDAL (most problematic package)
python -c "from osgeo import gdal; print(f'GDAL {gdal.__version__}')"

# 6. Check Gradio (GUI framework)
python -c "import gradio; print(f'Gradio {gradio.__version__}')"

# 7. Check OpenCV
python -c "import cv2; print(f'OpenCV {cv2.__version__}')"

# 8. Launch GUI
python pebblecounts_gui.py
```

If all steps pass, you're ready to process images!

---

## Alternative: Use Homebrew Python + pip

If conda continues to give you trouble, you can use Homebrew:

```bash
# Install Python and GDAL via Homebrew
brew install python gdal

# Create virtual environment
python3 -m venv pebblecounts-gui
source pebblecounts-gui/bin/activate

# Install Python packages
pip install numpy scipy opencv-python scikit-image scikit-learn matplotlib Pillow shapely gradio

# Install GDAL matching your system version
GDAL_VERSION=$(gdal-config --version)
pip install gdal==$GDAL_VERSION

# Launch
python pebblecounts_gui.py
```

---

## Still Having Issues?

### Collect diagnostic information:

```bash
# Save this output to share for debugging
conda info --all > conda_debug.txt
python -c "import sys; print(sys.version)" >> conda_debug.txt
which python >> conda_debug.txt
echo $PATH >> conda_debug.txt
ls -la $(which python) >> conda_debug.txt
otool -L $(which python) | grep ssl >> conda_debug.txt
```

### Common Questions

**Q: Why not just use pip?**
A: GDAL is a C library with complex dependencies. Conda handles the entire stack (C libraries + Python bindings). Pip only handles Python packages and expects system libraries to exist.

**Q: Can I skip GDAL?**
A: Only if you never process georeferenced ortho-images (GeoTIFFs). For regular photos (JPG/PNG), GDAL isn't needed, but it's better to have it installed.

**Q: Why Python 3.9 instead of 3.11?**
A: Better compatibility with older GDAL versions in conda-forge. Python 3.9 is mature and stable.

**Q: What if I need a newer Python version?**
A: Edit `environment.yml` and change `python=3.9` to `python=3.10` or `python=3.11`. Test carefully.

---

## Next Steps

Once your environment is working:

1. Read [GUI_GUIDE.md](GUI_GUIDE.md) for usage instructions
2. Try the example data in `example_data/`
3. Explore the Manual and Automated processing modes
4. Check out the Resolution Calculator for non-ortho images

---

## Contact

If you continue having SSL/conda issues after trying the above:
- Check if it's a known conda bug: https://github.com/conda/conda/issues
- Consider using Homebrew + pip instead (more manual but more reliable on macOS)
- Email: purinton@uni-potsdam.de
