# PebbleCounts GUI User Guide

## Overview

The PebbleCounts GUI is a modern, web-based graphical interface that makes grain-size analysis accessible to users without command-line experience. Built with Gradio, it provides an intuitive interface while maintaining all the power of the original PebbleCounts algorithms.

## Installation

### Prerequisites
- Python 3.7 or higher
- pip package manager

### Steps

1. **Clone or download the repository**
   ```bash
   git clone https://github.com/bpurinton/PebbleCounts.git
   cd PebbleCounts
   ```

2. **Install dependencies**
   ```bash
   pip install -r requirements.txt
   ```

   This will install all required packages including:
   - Core scientific libraries (numpy, scipy, scikit-image, scikit-learn)
   - Image processing (opencv-python, matplotlib, Pillow)
   - Geospatial tools (gdal, shapely)
   - GUI framework (gradio)

3. **Launch the GUI**
   ```bash
   python run_gui.py
   ```

   The launcher will:
   - Check for missing dependencies
   - Offer to install them automatically if needed
   - Launch the web interface at http://127.0.0.1:7860

## Using the GUI

### Interface Overview

The GUI consists of four main tabs:

#### 1. Manual Mode (Interactive)
Interactive k-means segmentation with manual grain selection via mouse clicks.

**When to use:**
- When you need highest accuracy
- For validation datasets
- When working with complex grain arrangements
- For research-grade measurements

**Key Features:**
- Multi-scale k-means clustering
- Interactive Otsu thresholding
- Manual grain selection for quality control
- Color masking for sand/vegetation removal

**Important Note:** This mode will open OpenCV windows outside the web interface for interactive selection. Monitor your terminal for instructions during processing.

#### 2. Automated Mode
Fully automatic segmentation without user interaction.

**When to use:**
- Batch processing multiple images
- When speed is more important than perfect accuracy
- For preliminary analysis
- When consistent (non-interactive) processing is required

**Key Features:**
- Fully automated workflow
- Automatic overlap and misfit filtering
- No manual intervention required
- Suitable for scripting and batch jobs

#### 3. Resolution Calculator
Calculate pixel resolution for non-orthorectified images.

**Inputs needed:**
- **Focal Length** (mm): Your camera's focal length
- **Flight Height** (meters): Height above the water surface when photo was taken
- **Sensor Width** (mm): Camera sensor width (check camera specs)
- **Sensor Height** (mm): Camera sensor height (check camera specs)
- **Image Width** (pixels): Width of your image in pixels
- **Image Height** (pixels): Height of your image in pixels

**Output:**
Ground Sample Distance (GSD) in mm/pixel to use as the `input_resolution` parameter.

#### 4. About
Complete documentation, feature descriptions, and citation information.

---

## Processing Workflow

### For Manual Mode:

1. **Upload Your Image**
   - Click "Upload Image" in the Manual Mode tab
   - Select a GeoTIFF (ortho) or JPG/PNG (non-ortho) file

2. **Set Basic Parameters**
   - **Georeferenced Ortho-Image?**: Check if your image is a georeferenced GeoTIFF
   - **Input Resolution**: If not ortho, provide the mm/pixel value (use Resolution Calculator)
   - **Interactively Subset Image?**: Check to select a region of interest
   - **Maximum Grain Size**: Expected largest grain in meters (default: 0.3m)
   - **Minimum B-axis Length**: Smallest grain to detect in pixels (default: 20)
   - **Otsu Threshold %**: Leave empty for interactive selection, or provide value (50-100)

3. **Adjust Advanced Parameters** (optional)
   - **Initial Denoising Strength**: Non-local means filter strength (1-10, default: 5)
   - **Edge Detection Thresholds**: Tophat, Sobel percentiles (80-100, default: 90)
   - **Canny Sigma**: Edge detection sensitivity (1-5, default: 2)

4. **Run Processing**
   - Click "Run Manual Processing"
   - Monitor the progress bar and status messages
   - **IMPORTANT**: Watch your terminal/console for interactive prompts:
     - Otsu threshold selection window (if not pre-specified)
     - Color masking selection (y/n prompt)
     - Grain selection windows (click grains, press 'q' to continue)

5. **Review Results**
   - Result visualization appears in the right panel
   - CSV file with measurements is available for download
   - Check the status box for file locations and processing details

### For Automated Mode:

1. **Upload Your Image**
   - Click "Upload Image" in the Automated Mode tab

2. **Set Parameters**
   - Basic parameters (same as Manual Mode)
   - **Max Ellipse Overlap %**: Maximum allowed grain overlap (default: 15)
   - **Max Area Misfit %**: Maximum ellipse fit error (default: 30)
   - **Min Grain Area**: Minimum grain area in pixels (default: 10)

3. **Run Processing**
   - Click "Run Automated Processing"
   - Processing runs without interaction
   - Monitor progress bar

4. **Review Results**
   - Result image and CSV download appear automatically

---

## Parameter Guide

### Critical Parameters

| Parameter | Range | Default | Description |
|-----------|-------|---------|-------------|
| **Input Resolution** | 0.5-2.0 mm/px | - | Pixel size in mm (critical for accurate measurements) |
| **Max Grain Size** | 0.1-1.0 m | 0.3 m | Expected largest grain long-axis |
| **Cutoff** | 5-50 px | 20 px | Minimum b-axis length to detect |
| **Otsu Threshold** | 50-100 % | 85% | Shadow masking sensitivity |

### Advanced Parameters

| Parameter | Range | Default | Effect |
|-----------|-------|---------|--------|
| **Initial Denoising** | 1-10 | 5 | Noise reduction strength (higher = smoother) |
| **Tophat Threshold** | 80-100 | 90 | Edge detection sensitivity for tophat filter |
| **Sobel Threshold** | 80-100 | 90 | Edge detection sensitivity for Sobel filter |
| **Canny Sigma** | 1-5 | 2 | Canny edge detector smoothing |

### Automated Mode Parameters

| Parameter | Range | Default | Purpose |
|-----------|-------|---------|---------|
| **Percent Overlap** | 5-30 % | 15% | Reject grains with >X% overlap with neighbors |
| **Misfit Threshold** | 10-50 % | 30% | Reject grains with poor ellipse fit |
| **Min Size Threshold** | 5-20 px | 10 px | Minimum grain area |

---

## Understanding Outputs

### CSV File Structure

The output CSV contains three sections:

1. **Parameters Section**
   - All processing parameters used
   - Allows reproduction of results

2. **Image Statistics**
   - Percentage of image not measured (fines/unfound rocks)
   - Percentage of background color (sand mask)

3. **Grain Measurements**

For **Orthorectified images**:
- UTM X, Y coordinates (meters)
- a-axis, b-axis (pixels and meters)
- Area (pixels and square meters)
- Orientation (radians)
- Ellipse area and percent area difference (misfit)

For **Non-ortho images**:
- a-axis, b-axis (pixels and meters)
- Area (pixels and square meters)
- Orientation (radians)
- Ellipse area and percent area difference (misfit)

### Visualization Figure

The PNG figure shows:
- Original image as background
- Semi-transparent mask of identified grains
- Red ellipse outlines for each grain
- Red lines showing major/minor axes
- Green dots at grain centroids

### Label Raster

For ortho images: GeoTIFF with unique labels for each grain
For non-ortho images: PNG with colored grain labels

---

## Tips and Best Practices

### Image Quality
✓ Use well-lit images with minimal shadows
✓ Capture images perpendicular to water surface
✓ Aim for 0.8-1.2 mm/pixel resolution
✓ Avoid images with excessive glare or water surface reflection

### Processing Strategy
✓ Start with small test regions before processing large areas
✓ Break large images into 2m x 2m tiles
✓ Use Manual mode on representative samples to validate Automated mode results
✓ Save parameter values that work well for similar images

### Interactive Mode Tips
✓ When clicking grains, click near the center
✓ Press 'r' in selection windows to toggle to original image
✓ Press 'q' when done selecting (not spacebar or Enter)
✓ Take breaks! You can't save mid-session, so plan accordingly

### Common Issues

**"No module named 'gradio'"**
- Solution: Run `pip install gradio` or `pip install -r requirements.txt`

**"Error: GDAL not found"**
- Solution: GDAL can be tricky to install. Try:
  - Windows: `conda install -c conda-forge gdal`
  - Mac: `brew install gdal` then `pip install gdal`
  - Linux: `sudo apt-get install gdal-bin libgdal-dev` then `pip install gdal`

**OpenCV windows not appearing (Manual mode)**
- Solution: Make sure you're running the GUI in an environment with display support
- For remote servers, you may need X11 forwarding or to use Automated mode instead

**Processing seems stuck**
- Check your terminal for interactive prompts (y/n questions, window instructions)
- The GUI can't display these prompts, so always monitor your console

---

## Keyboard Shortcuts (During Interactive Selection)

- **'q'**: Quit/close current window and proceed
- **'r'**: Toggle between mask and original image view
- **'spacebar'**: Accept ROI selection (subsetting)
- **Left Click**: Select grain or color for masking
- **Ctrl+C** (in terminal): Emergency stop

---

## Comparison: GUI vs Command Line

| Feature | GUI | Command Line |
|---------|-----|--------------|
| **Ease of Use** | ★★★★★ Very intuitive | ★★☆☆☆ Requires CLI knowledge |
| **Parameter Control** | ★★★★☆ Visual sliders and inputs | ★★★★★ Full control via flags |
| **Batch Processing** | ★★☆☆☆ One at a time | ★★★★★ Easy scripting |
| **Learning Curve** | ★★★★★ Minimal | ★★★☆☆ Moderate |
| **Customization** | ★★★☆☆ Pre-built interface | ★★★★★ Complete flexibility |
| **Progress Feedback** | ★★★★★ Real-time progress bar | ★★★☆☆ Terminal output |
| **Remote Use** | ★★★★☆ Web-based | ★★★★★ SSH-friendly |

**Bottom line**: Use the GUI for learning and single-image analysis. Use the command line for production workflows and batch processing.

---

## Getting Help

- **Documentation**: See the [full manual](docs/PebbleCounts_Manual.pdf)
- **Issues**: Report bugs at https://github.com/bpurinton/PebbleCounts/issues
- **Email**: purinton@uni-potsdam.de
- **Citation**: Purinton & Bookhagen (2019) https://doi.org/10.5194/esurf-7-859-2019

---

## Technical Notes

### Architecture
The GUI uses Gradio to provide a web interface that calls the original PebbleCounts Python scripts via subprocess. This means:
- All processing uses the same validated algorithms
- Results are identical between GUI and CLI
- The CLI scripts remain available for advanced users
- No processing code was modified (only wrapped)

### Security
- The GUI runs locally on your machine (127.0.0.1)
- No data is sent to external servers
- Files are processed on your computer
- Gradio's `share=False` prevents remote access

### Performance
- Processing speed is identical to command-line version
- Large images (>2m x 2m at high resolution) may take 10-30 minutes
- Memory usage scales with image size
- GPU acceleration is not used (CPU only)

---

## Future Enhancements

Potential improvements for future versions:
- [ ] Embedded OpenCV windows within the web interface
- [ ] Real-time parameter preview
- [ ] Batch processing interface for multiple images
- [ ] Progress logging with detailed status updates
- [ ] Result comparison tools
- [ ] Parameter presets for common scenarios
- [ ] Integrated grain-size distribution plotting

Contributions welcome at https://github.com/bpurinton/PebbleCounts
