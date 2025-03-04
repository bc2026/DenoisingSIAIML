"""
Image processing utilities for denoising and background subtraction.
"""

import numpy as np
from skimage.morphology import remove_small_objects
from scipy.ndimage import gaussian_filter, label
import tifffile as tiff


def load_tiff_stack(file_path, start_frame=1, end_frame=None):
    """
    Load a TIFF stack from file.
    
    Args:
        file_path (str): Path to the TIFF file
        start_frame (int): First frame to load (1-indexed)
        end_frame (int): Last frame to load (1-indexed), or None for all frames
        
    Returns:
        numpy.ndarray: Stack of images with shape (frames, height, width)
    """
    stack = tiff.imread(file_path)
    
    frames_total = stack.shape[0]
    if end_frame is None or end_frame > frames_total:
        end_frame = frames_total
    
    # Convert to 0-indexed for array slicing
    start_idx = start_frame - 1
    end_idx = end_frame
    
    return stack[start_idx:end_idx]


def moving_average_subtraction(stack, frame_range):
    """
    Perform background subtraction using moving average.
    
    Args:
        stack (numpy.ndarray): Input image stack with shape (frames, height, width)
        frame_range (int): Number of frames to average for background estimation
        
    Returns:
        tuple: (subtracted_matrix, subtracted_matrix_b, subtracted_matrix_d, subtracted_matrix_double)
            - subtracted_matrix: Absolute difference between frame and background (uint16)
            - subtracted_matrix_b: Positive difference (frame - background) (uint16)
            - subtracted_matrix_d: Positive difference (background - frame) (uint16)
            - subtracted_matrix_double: Raw difference (frame - background) (float64)
    """
    stack_frames, height, width = stack.shape
    
    # Pre-allocate matrices
    subtracted_matrix_double = np.zeros_like(stack, dtype=np.float64)
    subtracted_matrix = np.zeros_like(stack, dtype=np.int16)
    subtracted_matrix_b = np.zeros_like(stack, dtype=np.int16)
    subtracted_matrix_d = np.zeros_like(stack, dtype=np.int16)
    
    half_range = frame_range // 2
    
    for i in range(stack_frames):
        # Calculate moving average
        start_frame = max(i - half_range, 0)
        end_frame = min(start_frame + frame_range, stack_frames)
        if (end_frame - start_frame) < (frame_range - 1):
            start_frame = max(end_frame - frame_range, 0)
            
        average_image = np.mean(stack[start_frame:end_frame], axis=0).astype(np.uint16)
        
        # Calculate subtracted image as a double to preserve negative values
        image_matrix_double = stack[i].astype(np.float64)
        subtracted_matrix_double[i] = image_matrix_double - average_image.astype(np.float64)
        
        # Create different representations for analysis
        subtracted_matrix[i] = np.abs(subtracted_matrix_double[i]).astype(np.uint16)
        subtracted_matrix_b[i] = np.maximum(stack[i] - average_image, 0).astype(np.uint16)
        subtracted_matrix_d[i] = np.maximum(average_image - stack[i], 0).astype(np.uint16)
        
    return subtracted_matrix, subtracted_matrix_b, subtracted_matrix_d, subtracted_matrix_double


def identify_particles(subtracted_matrix, subtracted_matrix_b, subtracted_matrix_d, 
                       thresh_stack=10, particle_size_thresh=20):
    """
    Identify potential particles in subtracted images.
    
    Args:
        subtracted_matrix (numpy.ndarray): Absolute difference matrix
        subtracted_matrix_b (numpy.ndarray): Positive difference matrix (frame - background)
        subtracted_matrix_d (numpy.ndarray): Positive difference matrix (background - frame)
        thresh_stack (float): Threshold multiplier for standard deviation
        particle_size_thresh (int): Minimum size of particles in pixels
        
    Returns:
        tuple: (particle_region, particle_region_b, particle_region_d)
            - particle_region: Combined boolean matrix of identified particles
            - particle_region_b: Boolean matrix of bright particles
            - particle_region_d: Boolean matrix of dark particles
    """
    stack_frames, height, width = subtracted_matrix.shape
    
    # Pre-allocate matrices
    particle_region = np.zeros_like(subtracted_matrix, dtype=bool)
    particle_region_b = np.zeros_like(subtracted_matrix, dtype=bool)
    particle_region_d = np.zeros_like(subtracted_matrix, dtype=bool)
    
    for j in range(stack_frames):
        # Determine non-background pixels under the assumption that most pixels are background
        average_uint = np.mean(subtracted_matrix[j])
        stdev_uint = np.std(subtracted_matrix[j].astype(np.float64))
        
        # Create logical mask by filtering by intensity and size
        particle_region_d[j] = subtracted_matrix_d[j] > (average_uint + (thresh_stack * stdev_uint))
        particle_region_d[j] = remove_small_objects(particle_region_d[j], min_size=particle_size_thresh)
        
        particle_region_b[j] = subtracted_matrix_b[j] > (average_uint + (thresh_stack * stdev_uint))
        particle_region_b[j] = remove_small_objects(particle_region_b[j], min_size=particle_size_thresh)
        
        particle_region[j] = particle_region_d[j] | particle_region_b[j]
    
    return particle_region, particle_region_b, particle_region_d