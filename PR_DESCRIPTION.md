# Pull Request: Add Modern GUI Interface to PebbleCounts

## Summary

This PR introduces a comprehensive graphical user interface (GUI) for PebbleCounts, making grain-size analysis accessible to users without command-line expertise. The GUI is built with Gradio, a modern Python framework for creating web-based interfaces.

## Motivation

While PebbleCounts is a powerful tool, the command-line interface can be intimidating for new users or those unfamiliar with terminal operations. This GUI addresses that barrier while maintaining 100% backward compatibility with the existing CLI scripts.

## What's New

### Core Files Added

1. **`pebblecounts_gui.py`** (550 lines)
   - Full-featured Gradio web interface
   - Four main tabs: Manual Mode, Automated Mode, Resolution Calculator, About
   - Real-time progress tracking
   - Visual parameter controls (sliders, inputs, checkboxes)
   - File upload/download functionality

2. **`run_gui.py`** (70 lines)
   - Smart launcher script with dependency checking
   - Automatic installation prompting for missing packages
   - User-friendly error messages and guidance

3. **`requirements.txt`**
   - Complete dependency list for pip installation
   - Includes core scientific libraries, geospatial tools, and Gradio

4. **`GUI_GUIDE.md`** (400+ lines)
   - Comprehensive user documentation
   - Installation and setup instructions
   - Detailed workflow guides for both modes
   - Parameter reference tables
   - Troubleshooting section
   - Comparison of GUI vs CLI approaches

### Files Modified

1. **`README.md`**
   - Added prominent GUI section at the top
   - Quick start guide for new users
   - Feature overview and comparison
   - Updated installation instructions with GUI option

## Key Features

### User Interface
- **Web-based**: Runs locally in browser (no external servers)
- **Intuitive**: Visual controls for all parameters
- **Responsive**: Real-time progress updates
- **Documented**: Built-in help and tooltips

### Processing Modes
- **Manual Mode Tab**: Interactive k-means segmentation with full parameter control
- **Automated Mode Tab**: Fully automatic processing for batch workflows
- **Resolution Calculator Tab**: Built-in tool for camera resolution calculations
- **About Tab**: Complete documentation and citation information

### Technical Highlights
- Uses subprocess to call original CLI scripts (no algorithm changes)
- Same processing code = identical results to CLI
- Backward compatible - all existing scripts remain unchanged
- Runs on localhost only (127.0.0.1) for security
- Progress bars and status messages for user feedback

## Usage

### Quick Start
```bash
# Install dependencies
pip install -r requirements.txt

# Launch GUI
python run_gui.py
```

The interface opens automatically at `http://127.0.0.1:7860`

### GUI vs CLI

**Use the GUI for:**
- Learning PebbleCounts
- Single image analysis
- Parameter experimentation
- Visual feedback preference

**Use the CLI for:**
- Batch processing
- Automated workflows
- Scripting and integration
- Maximum parameter control

## Testing

The GUI has been designed to:
- ✅ Call the exact same processing code as CLI
- ✅ Accept all parameters available in CLI
- ✅ Produce identical output files (CSV, figures, labels)
- ✅ Handle both ortho and non-ortho images
- ✅ Support interactive features (mouse selection, color picking)

## Documentation

Comprehensive documentation is provided in:
- `GUI_GUIDE.md` - Complete user guide with examples
- `README.md` - Updated with GUI quick start
- In-app documentation via the "About" tab

## Backward Compatibility

**Important**: This PR does NOT modify any existing functionality:
- ✅ `PebbleCounts.py` - unchanged
- ✅ `PebbleCountsAuto.py` - unchanged
- ✅ `PCfunctions.py` - unchanged
- ✅ All CLI arguments work exactly as before
- ✅ Output formats remain identical

The GUI is purely additive - existing users can continue using the CLI without any changes to their workflows.

## Dependencies Added

The only new dependency is **Gradio** (and its sub-dependencies):
```
gradio>=4.0.0
```

All other dependencies were already required by PebbleCounts.

## Future Enhancements

Potential improvements for future versions:
- Embedded OpenCV windows within web interface
- Batch processing UI for multiple images
- Real-time parameter preview
- Result comparison tools
- Parameter presets

## Screenshots

The GUI includes:
- Clean, modern interface with tabbed navigation
- Visual parameter controls with sensible defaults
- Real-time progress tracking
- Integrated file management
- Comprehensive documentation

## Notes for Reviewers

1. **No algorithm changes**: The PR only adds interface code
2. **Fully optional**: Users can ignore the GUI and use CLI as before
3. **Well documented**: Extensive user guide and inline help
4. **Tested locally**: GUI successfully launches and processes test images
5. **Accessible**: Lowers barrier to entry for new users

## Closes

This PR addresses the need for a more user-friendly interface while maintaining the scientific rigor and flexibility of the original command-line tool.

---

**Ready for review and testing!** 🎉
