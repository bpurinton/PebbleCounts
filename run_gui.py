#!/usr/bin/env python3
"""
Simple launcher script for PebbleCounts GUI
This script checks dependencies and launches the GUI application
"""

import sys
import subprocess
import importlib.util


def check_dependency(package_name, import_name=None):
    """Check if a package is installed"""
    if import_name is None:
        import_name = package_name

    spec = importlib.util.find_spec(import_name)
    return spec is not None


def check_dependencies():
    """Check all required dependencies"""
    dependencies = {
        "numpy": "numpy",
        "scipy": "scipy",
        "cv2": "opencv-python",
        "skimage": "scikit-image",
        "sklearn": "scikit-learn",
        "matplotlib": "matplotlib",
        "PIL": "Pillow",
        "osgeo": "gdal",
        "shapely": "shapely",
        "gradio": "gradio",
    }

    missing = []
    for import_name, package_name in dependencies.items():
        if not check_dependency(package_name, import_name):
            missing.append(package_name)

    return missing


def main():
    print("=" * 70)
    print("PebbleCounts GUI Launcher")
    print("=" * 70)
    print("\nChecking dependencies...")

    missing = check_dependencies()

    if missing:
        print("\n⚠️  Missing dependencies detected:")
        for pkg in missing:
            print(f"   - {pkg}")

        print("\n📦 To install missing dependencies, run:")
        print(f"   pip install {' '.join(missing)}")
        print("\nOr install all requirements:")
        print("   pip install -r requirements.txt")

        response = input("\nWould you like to attempt automatic installation? (y/n): ")
        if response.lower() == 'y':
            print("\nInstalling dependencies...")
            try:
                subprocess.check_call([sys.executable, "-m", "pip", "install"] + missing)
                print("\n✓ Dependencies installed successfully!")
            except subprocess.CalledProcessError:
                print("\n❌ Failed to install dependencies automatically.")
                print("   Please install them manually using the commands above.")
                sys.exit(1)
        else:
            print("\nExiting. Please install dependencies and try again.")
            sys.exit(1)

    print("✓ All dependencies are installed!\n")
    print("Launching PebbleCounts GUI...")
    print("=" * 70 + "\n")

    # Import and run the GUI
    try:
        from pebblecounts_gui import main as gui_main
        gui_main()
    except Exception as e:
        print(f"\n❌ Error launching GUI: {e}")
        print("\nPlease check that all files are present:")
        print("   - pebblecounts_gui.py")
        print("   - PebbleCounts.py")
        print("   - PebbleCountsAuto.py")
        print("   - PCfunctions.py")
        sys.exit(1)


if __name__ == "__main__":
    main()
