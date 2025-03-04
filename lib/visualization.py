"""
Visualization utilities for particle tracking.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import os


def plot_frame(frame_data, title, colorbar_label="Value", cmap='gray'):
    """
    Plot a single frame with colorbar.
    
    Args:
        frame_data (numpy.ndarray): 2D array to display
        title (str): Plot title
        colorbar_label (str): Label for the colorbar
        cmap (str): Matplotlib colormap name
    """
    plt.figure(figsize=(6, 6))
    plt.imshow(frame_data, cmap=cmap)
    plt.colorbar(label=colorbar_label)
    plt.title(title)
    plt.show()


def plot_tracking_results(stack_matrix, subtracted_matrix, particle_region_roi, 
                        particle_centroid_roi, start_row, start_col, roi_size, 
                        save_path=None):
    """
    Create a composite figure showing tracking results.
    
    Args:
        stack_matrix (numpy.ndarray): Original image stack
        subtracted_matrix (numpy.ndarray): Background-subtracted stack
        particle_region_roi (numpy.ndarray): Boolean ROI tracking mask
        particle_centroid_roi (numpy.ndarray): Particle centroids
        start_row (int): Y coordinate of ROI center
        start_col (int): X coordinate of ROI center
        roi_size (int): Half-