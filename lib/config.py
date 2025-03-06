import os

"""
Configuration parameters for the denoising and particle tracking system.
"""

# Processing parameters
FRAME_RANGE = 20
THRESH_STACK = 10
THRESH_FRAME = 0.6
PARTICLE_SIZE_THRESH = 20
JUMP_THRESH = 20
ROI = 50
ROI_SMALL = 18
PARTICLE_DEFAULT_DIAMETER = 10
MIN_AREA = 3.14159 * (PARTICLE_DEFAULT_DIAMETER / 2) ** 2

# Analysis mode
ANALYZE_MODE = "track"

# Camera parameters
EXPOSURE_TIME = 0.100
PIXEL_SIZE = 0.065


def directory_setup(stackPath, stackFile, particleID, startColX, startRowY, X, Y):

    saveID = os.path.splitext(stackFile)[0]
    savingPath = os.path.join(stackPath, "Process", f"Particle {particleID}")

    # Ensure directory exists
    os.makedirs(savingPath, exist_ok=True)

    saveBaseName = os.path.join(
        savingPath, f"{saveID}_Particle{particleID}_FRAME_RANGE{FRAME_RANGE}_X{startColX}_Y{startRowY}"
    )
    saveBaseNameMovie = os.path.join(
        savingPath, f"{saveID}_FRAME_RANGE{FRAME_RANGE}"
    )
    save_base_name = os.path.join(
        savingPath, f"{saveID}_Particle{particleID}_frameRange{FRAME_RANGE}_X{X}_Y{Y}"
    )

    return saveBaseName, saveBaseNameMovie, save_base_name
