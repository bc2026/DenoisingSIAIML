"""
Particle tracking functionality.
"""

import numpy as np
from scipy.ndimage import label, center_of_mass, find_objects
from skimage.measure import regionprops
from scipy.spatial.distance import pdist


def extract_roi(particle_region, start_row, start_col, roi_size):
    """
    Extract region of interest around a particle.
    
    Args:
        particle_region (numpy.ndarray): Boolean matrix of identified particles
        start_row (int): Y coordinate of ROI center
        start_col (int): X coordinate of ROI center
        roi_size (int): Half-width of ROI
        
    Returns:
        numpy.ndarray: Boolean ROI matrix
    """
    stack_frames = particle_region.shape[0]
    particle_region_roi = np.zeros((stack_frames, 2*roi_size + 1, 2*roi_size + 1), dtype=bool)
    
    for j in range(stack_frames):
        # Handle edge cases at image boundaries
        row_start = max(0, start_row - roi_size)
        row_end = min(particle_region.shape[1], start_row + roi_size + 1)
        col_start = max(0, start_col - roi_size)
        col_end = min(particle_region.shape[2], start_col + roi_size + 1)
        
        # Extract ROI
        roi_height = row_end - row_start
        roi_width = col_end - col_start
        
        # Create ROI and center it within the ROI space
        roi_row_start = max(0, roi_size - (start_row - row_start))
        roi_row_end = roi_row_start + roi_height
        roi_col_start = max(0, roi_size - (start_col - col_start))
        roi_col_end = roi_col_start + roi_width
        
        particle_region_roi[j, roi_row_start:roi_row_end, roi_col_start:roi_col_end] = \
            particle_region[j, row_start:row_end, col_start:col_end]
    
    return particle_region_roi


def track_particles(particle_region_roi, particle_default_diameter=10, 
                   jump_thresh=20, min_area=78.5):
    """
    Track particles across frames within an ROI.
    
    Args:
        particle_region_roi (numpy.ndarray): Boolean ROI matrix
        particle_default_diameter (float): Expected particle diameter in pixels
        jump_thresh (float): Maximum allowed distance for particle movement between frames
        min_area (float): Minimum particle area in pixels
        
    Returns:
        numpy.ndarray: Particle centroids in each frame with shape (frames, 2)
    """
    stack_frames, roi_height, roi_width = particle_region_roi.shape
    particle_centroid_roi = np.zeros((stack_frames, 2))
    jump_frames = 0
    
    # Initialize with center of ROI
    previous_coord = np.array([roi_height // 2, roi_width // 2])
    
    for j in range(stack_frames):
        # Label regions in the logical mask
        labeled_image, num_features = label(particle_region_roi[j])
        
        if num_features > 0:
            # Get properties of the possible particle regions
            slices = find_objects(labeled_image)
            centroids = [center_of_mass(particle_region_roi[j], labeled_image, i + 1) for i in range(num_features)]
            particle_props = [{'Centroid': centroids[i], 'BoundingBox': slices[i]} for i in range(num_features)]
            
            # If more than one region identified, choose the one closest to previous centroid
            if num_features >= 2:
                coord_list = np.zeros((num_features, 2))
                distances = np.zeros(num_features)
                
                for h in range(num_features):
                    coord_list[h] = particle_props[h]['Centroid']
                    distances[h] = pdist([previous_coord, coord_list[h]])[0]
                
                chosen_particle = np.argmin(distances)
                particle_props = [particle_props[chosen_particle]]
                particle_region_roi[j] = (labeled_image == (chosen_particle + 1))
            
            # Set the centroid
            centroid = particle_props[0]['Centroid']
            particle_centroid_roi[j] = centroid
            
            # If the particle jumped too far, assume it stayed in place
            distance_from_prev = pdist([previous_coord, particle_centroid_roi[j]])[0]
            if distance_from_prev > jump_thresh:
                jump_frames += 1
                particle_centroid_roi[j] = previous_coord
                
                # Create a circular mask centered at previous location
                yy, xx = np.mgrid[:roi_height, :roi_width]
                circle_mask = ((xx - previous_coord[1])**2 + (yy - previous_coord[0])**2) < (particle_default_diameter/2)**2
                particle_region_roi[j] = circle_mask
        else:
            # No particles detected, keep previous position
            jump_frames += 1
            particle_centroid_roi[j] = previous_coord
            
            # Create a circular mask centered at previous location
            yy, xx = np.mgrid[:roi_height, :roi_width]
            circle_mask = ((xx - previous_coord[1])**2 + (yy - previous_coord[0])**2) < (particle_default_diameter/2)**2
            particle_region_roi[j] = circle_mask
        
        # Ensure minimum particle size
        particle_props_final = regionprops(particle_region_roi[j].astype(int))
        if particle_props_final and particle_props_final[0].area < min_area:
            yy, xx = np.mgrid[:roi_height, :roi_width]
            circle_mask = ((xx - particle_centroid_roi[j][1])**2 + 
                          (yy - particle_centroid_roi[j][0])**2) < (particle_default_diameter/2)**2
            particle_region_roi[j] = circle_mask
        
        # Update previous coordinates for next frame
        previous_coord = particle_centroid_roi[j].copy()
    
    return particle_centroid_roi, particle_region_roi, jump_frames