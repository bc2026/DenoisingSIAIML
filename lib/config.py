"""
Configuration parameters for the denoising and particle tracking system.
"""

# Processing parameters
FRAME_RANGE = 20  # Used in moving average
THRESH_STACK = 10  # Threshold for identifying non-background pixels
THRESH_FRAME = 0.6  # Threshold for identifying non-background pixels within ROI
PARTICLE_SIZE_THRESH = 20  # Minimum number of pixels in particle
JUMP_THRESH = 20  # Longest distance allowed between frames for particle tracking
ROI = 50  # ROI size around particle (pixels)
ROI_SMALL = 18  # Smaller region around particle
PARTICLE_DEFAULT_DIAMETER = 10
MIN_AREA = 3.14159 * (PARTICLE_DEFAULT_DIAMETER/2)**2  # Minimum particle area in pixels

# Analysis mode
ANALYZE_MODE = "track"  # Options: "track", "analyze", etc.

# Camera parameters
EXPOSURE_TIME = 0.100  # seconds
PIXEL_SIZE = 0.065  # microns