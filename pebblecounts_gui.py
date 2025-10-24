#!/usr/bin/env python3
"""
PebbleCounts GUI - A graphical interface for gravel grain-size analysis
Developed using Gradio framework
"""

import gradio as gr
import cv2
import numpy as np
import os
import sys
import subprocess
import tempfile
from pathlib import Path
import matplotlib.pyplot as plt
from PIL import Image
import io

# Try to import processing modules
try:
    import PCfunctions as func
    from osgeo import gdal
except ImportError as e:
    print(f"Warning: Some dependencies may not be installed: {e}")


class PebbleCountsGUI:
    """Main GUI application for PebbleCounts"""

    def __init__(self):
        self.temp_dir = tempfile.mkdtemp()

    def run_pebblecounts_manual(
        self,
        image_file,
        is_ortho,
        input_resolution,
        subset_image,
        max_grain_size,
        cutoff,
        otsu_threshold,
        first_nl_denoise,
        tophat_th,
        sobel_th,
        canny_sig,
        progress=gr.Progress()
    ):
        """
        Run PebbleCounts in manual mode
        """
        if image_file is None:
            return None, None, "Please upload an image file."

        try:
            progress(0, desc="Initializing...")

            # Build command
            cmd = ["python", "PebbleCounts.py"]
            cmd.extend(["-im", image_file.name])
            cmd.extend(["-ortho", "y" if is_ortho else "n"])

            if not is_ortho and input_resolution:
                cmd.extend(["-input_resolution", str(input_resolution)])
            elif not is_ortho and not input_resolution:
                return None, None, "Please provide input resolution for non-ortho images."

            cmd.extend(["-subset", "y" if subset_image else "n"])
            cmd.extend(["-maxGS", str(max_grain_size)])
            cmd.extend(["-cutoff", str(cutoff)])

            if otsu_threshold:
                cmd.extend(["-otsu_threshold", str(otsu_threshold)])

            cmd.extend(["-first_nl_denoise", str(first_nl_denoise)])
            cmd.extend(["-tophat_th", str(tophat_th)])
            cmd.extend(["-sobel_th", str(sobel_th)])
            cmd.extend(["-canny_sig", str(canny_sig)])
            cmd.extend(["-resize", "0.6"])  # Smaller for GUI

            progress(0.1, desc="Starting processing...")

            # Run the command
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                stdin=subprocess.PIPE,
                text=True,
                bufsize=1
            )

            # Handle interactive prompts
            progress(0.3, desc="Processing... (check terminal for interactive prompts)")

            stdout, stderr = process.communicate(input="n\n")

            progress(0.9, desc="Finalizing...")

            if process.returncode != 0:
                return None, None, f"Error running PebbleCounts:\n{stderr}\n{stdout}"

            # Get output files
            base_name = os.path.splitext(os.path.basename(image_file.name))[0]
            base_dir = os.path.dirname(image_file.name)
            csv_out = os.path.join(base_dir, base_name + "_PebbleCounts_CSV.csv")
            fig_out = os.path.join(base_dir, base_name + "_PebbleCounts_FIGURE.png")

            # Check if outputs exist
            result_image = None
            csv_file = None

            if os.path.exists(fig_out):
                result_image = fig_out

            if os.path.exists(csv_out):
                csv_file = csv_out

            progress(1.0, desc="Complete!")

            status = f"Processing complete!\n\nResults saved to:\n- {csv_out}\n- {fig_out}\n\n{stdout}"

            return result_image, csv_file, status

        except Exception as e:
            return None, None, f"Error: {str(e)}"

    def run_pebblecounts_auto(
        self,
        image_file,
        is_ortho,
        input_resolution,
        subset_image,
        cutoff,
        percent_overlap,
        misfit_threshold,
        min_size_threshold,
        otsu_threshold,
        first_nl_denoise,
        tophat_th,
        sobel_th,
        canny_sig,
        progress=gr.Progress()
    ):
        """
        Run PebbleCounts in automated mode
        """
        if image_file is None:
            return None, None, "Please upload an image file."

        try:
            progress(0, desc="Initializing...")

            # Build command
            cmd = ["python", "PebbleCountsAuto.py"]
            cmd.extend(["-im", image_file.name])
            cmd.extend(["-ortho", "y" if is_ortho else "n"])

            if not is_ortho and input_resolution:
                cmd.extend(["-input_resolution", str(input_resolution)])
            elif not is_ortho and not input_resolution:
                return None, None, "Please provide input resolution for non-ortho images."

            cmd.extend(["-subset", "y" if subset_image else "n"])
            cmd.extend(["-cutoff", str(cutoff)])
            cmd.extend(["-percent_overlap", str(percent_overlap)])
            cmd.extend(["-misfit_threshold", str(misfit_threshold)])
            cmd.extend(["-min_size_threshold", str(min_size_threshold)])

            if otsu_threshold:
                cmd.extend(["-otsu_threshold", str(otsu_threshold)])

            cmd.extend(["-first_nl_denoise", str(first_nl_denoise)])
            cmd.extend(["-tophat_th", str(tophat_th)])
            cmd.extend(["-sobel_th", str(sobel_th)])
            cmd.extend(["-canny_sig", str(canny_sig)])
            cmd.extend(["-resize", "0.6"])

            progress(0.1, desc="Starting processing...")

            # Run the command
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                stdin=subprocess.PIPE,
                text=True,
                bufsize=1
            )

            # Provide default answers to interactive prompts
            progress(0.3, desc="Processing... (answering prompts automatically)")
            stdout, stderr = process.communicate(input="n\nn\n")

            progress(0.9, desc="Finalizing...")

            if process.returncode != 0:
                return None, None, f"Error running PebbleCounts:\n{stderr}\n{stdout}"

            # Get output files
            base_name = os.path.splitext(os.path.basename(image_file.name))[0]
            base_dir = os.path.dirname(image_file.name)
            csv_out = os.path.join(base_dir, base_name + "_PebbleCountsAuto_CSV.csv")
            fig_out = os.path.join(base_dir, base_name + "_PebbleCountsAuto_FIGURE.png")

            result_image = None
            csv_file = None

            if os.path.exists(fig_out):
                result_image = fig_out

            if os.path.exists(csv_out):
                csv_file = csv_out

            progress(1.0, desc="Complete!")

            status = f"Processing complete!\n\nResults saved to:\n- {csv_out}\n- {fig_out}\n\n{stdout}"

            return result_image, csv_file, status

        except Exception as e:
            return None, None, f"Error: {str(e)}"

    def calculate_resolution(
        self,
        focal_length,
        flight_height,
        sensor_width,
        sensor_height,
        image_width,
        image_height
    ):
        """Calculate camera resolution from parameters"""
        try:
            # Ground sample distance calculation
            gsd_x = (sensor_width * flight_height * 1000) / (focal_length * image_width)
            gsd_y = (sensor_height * flight_height * 1000) / (focal_length * image_height)
            gsd_avg = (gsd_x + gsd_y) / 2

            result = f"""
Camera Resolution Calculation:

Focal Length: {focal_length} mm
Flight Height: {flight_height} m
Sensor Size: {sensor_width} x {sensor_height} mm
Image Size: {image_width} x {image_height} pixels

Ground Sample Distance (GSD):
- X direction: {gsd_x:.3f} mm/pixel
- Y direction: {gsd_y:.3f} mm/pixel
- Average: {gsd_avg:.3f} mm/pixel

Use the average value ({gsd_avg:.3f}) as the input_resolution parameter.
"""
            return result, gsd_avg

        except Exception as e:
            return f"Error calculating resolution: {str(e)}", None

    def create_interface(self):
        """Create the Gradio interface"""

        with gr.Blocks(title="PebbleCounts - Grain Size Analysis", theme=gr.themes.Soft()) as app:
            gr.Markdown("""
            # PebbleCounts - Gravel Grain-Size Analysis Tool

            A tool for identifying and measuring gravel grain sizes from river photographs.
            Supports both **Manual** (interactive) and **Automated** processing modes.

            **Note:** Manual mode requires interactive mouse input. Some features may need to be run in a terminal.
            """)

            with gr.Tabs() as tabs:
                # Manual Mode Tab
                with gr.Tab("Manual Mode (Interactive)"):
                    gr.Markdown("""
                    ### Manual Mode - Interactive K-means Segmentation
                    This mode uses multi-scale k-means clustering with manual grain selection via mouse clicks.
                    **Best for higher accuracy** but requires user interaction.

                    **Note:** This mode will open OpenCV windows for interactive selection.
                    Monitor your terminal/command prompt for instructions.
                    """)

                    with gr.Row():
                        with gr.Column(scale=1):
                            manual_image = gr.File(label="Upload Image", file_types=["image"])
                            manual_ortho = gr.Checkbox(label="Georeferenced Ortho-Image?", value=True)
                            manual_resolution = gr.Number(label="Input Resolution (mm/pixel, if not ortho)", value=0.8)
                            manual_subset = gr.Checkbox(label="Interactively Subset Image?", value=False)

                            with gr.Accordion("Basic Parameters", open=True):
                                manual_maxgs = gr.Slider(0.1, 1.0, value=0.3, step=0.05,
                                                        label="Maximum Grain Size (meters)")
                                manual_cutoff = gr.Slider(5, 50, value=20, step=1,
                                                         label="Minimum B-axis Length (pixels)")
                                manual_otsu = gr.Number(label="Otsu Threshold % (0 for interactive)", value=0, minimum=0, maximum=100)

                            with gr.Accordion("Advanced Parameters", open=False):
                                manual_denoise = gr.Slider(1, 10, value=5, step=1,
                                                          label="Initial Denoising Strength")
                                manual_tophat = gr.Slider(80, 100, value=90, step=1,
                                                         label="Tophat Threshold %")
                                manual_sobel = gr.Slider(80, 100, value=90, step=1,
                                                        label="Sobel Threshold %")
                                manual_canny = gr.Slider(1, 5, value=2, step=1,
                                                        label="Canny Sigma")

                            manual_btn = gr.Button("Run Manual Processing", variant="primary")

                        with gr.Column(scale=1):
                            manual_result_img = gr.Image(label="Result Visualization", type="filepath")
                            manual_result_csv = gr.File(label="Download CSV Results")
                            manual_status = gr.Textbox(label="Status", lines=10, max_lines=20)

                    manual_btn.click(
                        fn=self.run_pebblecounts_manual,
                        inputs=[
                            manual_image, manual_ortho, manual_resolution, manual_subset,
                            manual_maxgs, manual_cutoff, manual_otsu, manual_denoise,
                            manual_tophat, manual_sobel, manual_canny
                        ],
                        outputs=[manual_result_img, manual_result_csv, manual_status]
                    )

                # Automated Mode Tab
                with gr.Tab("Automated Mode"):
                    gr.Markdown("""
                    ### Automated Mode - Fully Automatic Segmentation
                    This mode performs fully automated grain segmentation without user interaction.
                    **Best for batch processing** but may have higher uncertainty in measurements.
                    """)

                    with gr.Row():
                        with gr.Column(scale=1):
                            auto_image = gr.File(label="Upload Image", file_types=["image"])
                            auto_ortho = gr.Checkbox(label="Georeferenced Ortho-Image?", value=True)
                            auto_resolution = gr.Number(label="Input Resolution (mm/pixel, if not ortho)", value=0.8)
                            auto_subset = gr.Checkbox(label="Interactively Subset Image?", value=False)

                            with gr.Accordion("Basic Parameters", open=True):
                                auto_cutoff = gr.Slider(5, 50, value=20, step=1,
                                                       label="Minimum B-axis Length (pixels)")
                                auto_overlap = gr.Slider(5, 30, value=15, step=1,
                                                        label="Max Ellipse Overlap %")
                                auto_misfit = gr.Slider(10, 50, value=30, step=1,
                                                       label="Max Area Misfit %")
                                auto_minsize = gr.Slider(5, 20, value=10, step=1,
                                                        label="Min Grain Area (pixels)")
                                auto_otsu = gr.Number(label="Otsu Threshold % (0 for interactive)", value=0, minimum=0, maximum=100)

                            with gr.Accordion("Advanced Parameters", open=False):
                                auto_denoise = gr.Slider(1, 10, value=5, step=1,
                                                        label="Initial Denoising Strength")
                                auto_tophat = gr.Slider(80, 100, value=90, step=1,
                                                       label="Tophat Threshold %")
                                auto_sobel = gr.Slider(80, 100, value=90, step=1,
                                                      label="Sobel Threshold %")
                                auto_canny = gr.Slider(1, 5, value=2, step=1,
                                                      label="Canny Sigma")

                            auto_btn = gr.Button("Run Automated Processing", variant="primary")

                        with gr.Column(scale=1):
                            auto_result_img = gr.Image(label="Result Visualization", type="filepath")
                            auto_result_csv = gr.File(label="Download CSV Results")
                            auto_status = gr.Textbox(label="Status", lines=10, max_lines=20)

                    auto_btn.click(
                        fn=self.run_pebblecounts_auto,
                        inputs=[
                            auto_image, auto_ortho, auto_resolution, auto_subset,
                            auto_cutoff, auto_overlap, auto_misfit, auto_minsize, auto_otsu,
                            auto_denoise, auto_tophat, auto_sobel, auto_canny
                        ],
                        outputs=[auto_result_img, auto_result_csv, auto_status]
                    )

                # Resolution Calculator Tab
                with gr.Tab("Resolution Calculator"):
                    gr.Markdown("""
                    ### Calculate Camera Resolution
                    Use this tool to calculate the pixel resolution (mm/pixel) for non-orthorectified images.
                    You'll need your camera specifications and flight parameters.
                    """)

                    with gr.Row():
                        with gr.Column():
                            calc_focal = gr.Number(label="Focal Length (mm)", value=35)
                            calc_height = gr.Number(label="Flight Height (meters)", value=2.0)
                            calc_sensor_w = gr.Number(label="Sensor Width (mm)", value=23.6)
                            calc_sensor_h = gr.Number(label="Sensor Height (mm)", value=15.7)
                            calc_img_w = gr.Number(label="Image Width (pixels)", value=4000)
                            calc_img_h = gr.Number(label="Image Height (pixels)", value=3000)

                            calc_btn = gr.Button("Calculate Resolution", variant="primary")

                        with gr.Column():
                            calc_result = gr.Textbox(label="Calculation Results", lines=15)
                            calc_value = gr.Number(label="Resolution to Use (mm/pixel)", interactive=False)

                    calc_btn.click(
                        fn=self.calculate_resolution,
                        inputs=[calc_focal, calc_height, calc_sensor_w, calc_sensor_h,
                               calc_img_w, calc_img_h],
                        outputs=[calc_result, calc_value]
                    )

                # About Tab
                with gr.Tab("About"):
                    gr.Markdown("""
                    ## About PebbleCounts

                    **PebbleCounts** is a Python-based image analysis application for identifying and measuring
                    gravel grain sizes from river photographs.

                    ### Features
                    - **Two Processing Modes:**
                      - Manual: Interactive k-means segmentation with manual grain selection
                      - Automated: Fully automatic segmentation for batch processing

                    - **Supports Multiple Image Types:**
                      - Georeferenced orthorectified images (GeoTIFF with UTM projection)
                      - Non-orthorectified overhead photos with calculated resolution

                    - **Advanced Image Processing:**
                      - Multi-scale edge detection (Tophat, Canny, Sobel)
                      - Shadow masking via Otsu thresholding
                      - Optional color masking for sand/vegetation removal
                      - K-means clustering for grain segmentation

                    - **Comprehensive Outputs:**
                      - CSV files with grain measurements (size, position, orientation, area)
                      - Visualization figures with ellipse overlays
                      - Labeled raster outputs (GeoTIFF or PNG)

                    ### Optimal Image Characteristics
                    - Resolution: 0.8-1.2 mm/pixel
                    - UTM projected for georeferenced images
                    - Clear lighting without excessive shadows
                    - Overhead view perpendicular to water surface

                    ### Citation
                    Developed by Ben Purinton (purinton[at]uni-potsdam.de), 2019

                    See the GitHub repository for more details:
                    [https://github.com/bpurinton/PebbleCounts](https://github.com/bpurinton/PebbleCounts)

                    ### License
                    GNU General Public License v3.0

                    ### GUI Version
                    This graphical interface was created to make PebbleCounts more accessible.
                    The original CLI scripts remain available for advanced users and scripting.
                    """)

        return app


def main():
    """Main entry point for the GUI"""
    gui = PebbleCountsGUI()
    app = gui.create_interface()

    print("\n" + "="*70)
    print("PebbleCounts GUI - Gravel Grain-Size Analysis Tool")
    print("="*70)
    print("\nStarting the GUI application...")
    print("Once launched, open the URL shown below in your web browser.")
    print("\nNote: Manual mode will open OpenCV windows for interactive selection.")
    print("Monitor this terminal for instructions during processing.")
    print("="*70 + "\n")

    app.launch(
        server_name="127.0.0.1",
        server_port=7860,
        share=False,
        show_error=True
    )


if __name__ == "__main__":
    main()
