"""
Visualization tools for particle tracking and denoising.
"""

import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.animation import FuncAnimation
from matplotlib.colors import LinearSegmentedColormap

class ParticleVisualizer:
    def __init__(self, save_path=None, figsize=(10, 8), dpi=100):
        """
        Initialize visualizer with optional save path.
        
        Parameters:
        -----------
        save_path : str
            Directory to save visualizations
        figsize : tuple
            Figure size (width, height) in inches
        dpi : int
            DPI for saved figures
        """
        self.save_path = save_path
        self.figsize = figsize
        self.dpi = dpi
        
        # Create custom colormaps
        self._create_custom_colormaps()

    def _create_custom_colormaps(self):
        """Create custom colormaps for different visualization types"""
        # Create a heatmap colormap (red-orange-yellow)
        self.heatmap_cmap = LinearSegmentedColormap.from_list(
            "particle_heatmap", ["darkred", "red", "orange", "yellow", "white"]
        )
        
        # Create a diverging colormap for difference images
        self.diff_cmap = plt.cm.RdBu_r
        
    def plot_frame(self, image, title=None, colormap='gray', show_colorbar=True, 
                  roi=None, centroids=None, frame_num=None, save_name=None):
        """
        Plot a single frame with optional ROI overlay and centroids.
        
        Parameters:
        -----------
        image : numpy.ndarray
            2D array with image data
        title : str
            Title for the plot
        colormap : str or matplotlib.colors.Colormap
            Colormap to use
        show_colorbar : bool
            Whether to show a colorbar
        roi : tuple
            Region of interest to highlight (y, x, height, width)
        centroids : list
            List of (y, x) centroids to mark
        frame_num : int
            Frame number to include in title
        save_name : str
            File name for saving the figure
        """
        plt.figure(figsize=self.figsize)
        
        # Plot main image
        plt.imshow(image, cmap=colormap)
        
        # Add colorbar if requested
        if show_colorbar:
            cbar = plt.colorbar(label='Value')
            cbar.ax.tick_params(labelsize=8)
        
        # Build title with frame number if provided
        plot_title = title if title else "Frame"
        if frame_num is not None:
            plot_title = f"{plot_title} (Frame {frame_num})"
        plt.title(plot_title)
        
        # Draw ROI rectangle if provided
        if roi is not None:
            y, x, h, w = roi
            rect = plt.Rectangle((x, y), w, h, linewidth=2, 
                               edgecolor='yellow', facecolor='none')
            plt.gca().add_patch(rect)
        
        # Plot centroids if provided
        if centroids is not None:
            if isinstance(centroids, tuple) or (isinstance(centroids, np.ndarray) and centroids.size == 2):
                # Single centroid
                y, x = centroids
                plt.plot(x, y, 'r+', markersize=10, markeredgewidth=2)
            else:
                # Multiple centroids
                for centroid in centroids:
                    y, x = centroid
                    plt.plot(x, y, 'r+', markersize=10, markeredgewidth=2)
        
        plt.axis('on')
        plt.tight_layout()
        
        # Save if requested and path provided
        if save_name and self.save_path:
            save_path = os.path.join(self.save_path, save_name)
            plt.savefig(save_path, dpi=self.dpi)
        
        plt.show()
    
    def plot_comparison(self, images, titles, colormap='gray', show_colorbar=True, 
                      frame_num=None, save_name=None, layout=None):
        """
        Plot multiple images for comparison.
        
        Parameters:
        -----------
        images : list
            List of image arrays to plot
        titles : list
            List of titles for each subplot
        colormap : str, matplotlib.colors.Colormap, or list
            Colormap(s) to use for each image
        show_colorbar : bool
            Whether to show colorbars
        frame_num : int
            Frame number to include in suptitle
        save_name : str
            File name for saving the figure
        layout : tuple
            Grid layout (rows, cols), if None will be determined automatically
        """
        n_images = len(images)
        
        # Determine layout if not specified
        if layout is None:
            if n_images <= 2:
                layout = (1, n_images)
            elif n_images <= 4:
                layout = (2, 2)
            elif n_images <= 6:
                layout = (2, 3)
            else:
                layout = (3, 3)
        
        rows, cols = layout
        fig, axes = plt.subplots(rows, cols, figsize=self.figsize)
        
        # Convert to array for consistent indexing
        if n_images == 1:
            axes = np.array([axes])
        elif rows == 1 and cols > 1:
            axes = axes.reshape(1, -1)
        
        # Convert colormap to list if it's a single value
        if not isinstance(colormap, list):
            colormap = [colormap] * n_images
        
        # Plot each image
        for i, (img, title, cmap) in enumerate(zip(images, titles, colormap)):
            if i >= rows * cols:
                break
                
            r, c = i // cols, i % cols
            im = axes[r, c].imshow(img, cmap=cmap)
            axes[r, c].set_title(title)
            axes[r, c].axis('on')
            
            if show_colorbar:
                fig.colorbar(im, ax=axes[r, c], fraction=0.046, pad=0.04)
        
        # Turn off unused subplots
        for i in range(n_images, rows * cols):
            r, c = i // cols, i % cols
            axes[r, c].axis('off')
        
        # Add overall title if frame number provided
        if frame_num is not None:
            fig.suptitle(f'Frame {frame_num}', fontsize=16)
            
        plt.tight_layout()
        
        # Save if requested and path provided
        if save_name and self.save_path:
            save_path = os.path.join(self.save_path, save_name)
            plt.savefig(save_path, dpi=self.dpi)
            
        plt.show()
    
    def plot_tracking_results(self, original_frames, processed_frames, centroids, 
                            frame_indices=None, save_name_prefix=None):
        """
        Plot tracking results for selected frames.
        
        Parameters:
        -----------
        original_frames : numpy.ndarray
            Original image stack
        processed_frames : numpy.ndarray
            Processed image stack
        centroids : numpy.ndarray
            Array of centroids (frame, y, x)
        frame_indices : list
            List of frame indices to plot, if None will select evenly spaced frames
        save_name_prefix : str
            Prefix for saved files
        """
        n_frames = original_frames.shape[0]
        
        # Select frames to display if not specified
        if frame_indices is None:
            if n_frames <= 9:
                frame_indices = list(range(n_frames))
            else:
                # Select evenly spaced frames
                frame_indices = np.linspace(0, n_frames-1, 9, dtype=int)
        
        # Plot each selected frame
        for i, frame_idx in enumerate(frame_indices):
            # Create a 2x1 subplot
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize)
            
            # Plot original frame
            ax1.imshow(original_frames[frame_idx], cmap='gray')
            ax1.set_title(f"Original (Frame {frame_idx})")
            ax1.axis('on')
            
            # Plot processed frame with centroid
            ax2.imshow(processed_frames[frame_idx], cmap='viridis')
            ax2.set_title(f"Processed (Frame {frame_idx})")
            
            # Add centroid marker
            if frame_idx < len(centroids):
                y, x = centroids[frame_idx]
                ax2.plot(x, y, 'r+', markersize=12, markeredgewidth=2)
            
            ax2.axis('on')
            
            plt.tight_layout()
            
            # Save if requested and path provided
            if save_name_prefix and self.save_path:
                save_path = os.path.join(
                    self.save_path, f"{save_name_prefix}_frame{frame_idx}.png")
                plt.savefig(save_path, dpi=self.dpi)
                
            plt.show()
    
    def plot_trajectory(self, centroids, start_col, start_row, roi_size, 
                      frame_skip=1, save_name=None):
        """
        Plot the trajectory of a tracked particle.
        
        Parameters:
        -----------
        centroids : numpy.ndarray
            Array of centroids (frame, y, x)
        start_col : int
            Starting column (X) coordinate
        start_row : int
            Starting row (Y) coordinate
        roi_size : int
            Size of ROI
        frame_skip : int
            Plot every n-th frame
        save_name : str
            File name for saving the figure
        """
        plt.figure(figsize=self.figsize)
        
        # Convert ROI-relative coordinates to global coordinates
        global_x = centroids[::frame_skip, 1] + (start_col - roi_size)
        global_y = centroids[::frame_skip, 0] + (start_row - roi_size)
        
        # Create color gradient for time progression
        n_points = len(global_x)
        colors = plt.cm.viridis(np.linspace(0, 1, n_points))
        
        # Plot trajectory with gradient color
        for i in range(n_points-1):
            plt.plot(global_x[i:i+2], global_y[i:i+2], color=colors[i], linewidth=1.5)
        
        # Add start and end markers
        plt.plot(global_x[0], global_y[0], 'go', markersize=8, label='Start')
        plt.plot(global_x[-1], global_y[-1], 'ro', markersize=8, label='End')
        
        # Add colorbar to indicate time progression
        sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis)
        sm.set_array([])
        cbar = plt.colorbar(sm)
        cbar.set_label('Time progression')
        
        plt.title('Particle Trajectory')
        plt.xlabel('X Position (pixels)')
        plt.ylabel('Y Position (pixels)')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Make equal aspect ratio so circles look like circles
        plt.axis('equal')
        plt.tight_layout()
        
        # Save if requested and path provided
        if save_name and self.save_path:
            save_path = os.path.join(self.save_path, save_name)
            plt.savefig(save_path, dpi=self.dpi)
            
        plt.show()
    
    def create_animation(self, frames, centroids=None, roi_size=None, 
                       fps=10, save_name=None, title=None):
        """
        Create an animation of particle tracking.
        
        Parameters:
        -----------
        frames : numpy.ndarray
            3D array with image frames
        centroids : numpy.ndarray
            Array of centroids (frame, y, x)
        roi_size : int
            Size of ROI to zoom in on particle
        fps : int
            Frames per second in animation
        save_name : str
            File name for saving the animation
        title : str
            Title for the animation
        
        Returns:
        --------
        matplotlib.animation.FuncAnimation
            Animation object
        """
        fig, ax = plt.subplots(figsize=self.figsize)
        
        # Function to update the plot for each frame
        def update(frame_idx):
            ax.clear()
            
            # Plot the frame
            ax.imshow(frames[frame_idx], cmap='viridis')
            
            # Add centroid marker if available
            if centroids is not None and frame_idx < len(centroids):
                y, x = centroids[frame_idx]
                ax.plot(x, y, 'r+', markersize=12, markeredgewidth=2)
            
            # Set the title with the frame number
            frame_title = f"Frame {frame_idx}"
            if title:
                frame_title = f"{title} - {frame_title}"
            ax.set_title(frame_title)
            
            # If ROI size provided, zoom in on the particle
            if centroids is not None and roi_size is not None and frame_idx < len(centroids):
                y, x = centroids[frame_idx]
                ax.set_xlim(x - roi_size, x + roi_size)
                ax.set_ylim(y + roi_size, y - roi_size)  # Reverse y-axis
                
            return ax,
        
        # Create the animation
        anim = FuncAnimation(fig, update, frames=len(frames), interval=1000/fps, blit=False)
        
        # Save animation if requested
        if save_name and self.save_path:
            save_path = os.path.join(self.save_path, save_name)
            anim.save(save_path, writer='pillow', fps=fps)
        
        plt.tight_layout()
        plt.show()
        
        return anim
    
    def plot_intensity_profile(self, frames, centroids, roi_small, 
                             frame_indices=None, save_name=None):
        """
        Plot intensity profiles for a particle across frames.
        
        Parameters:
        -----------
        frames : numpy.ndarray
            3D array with image frames
        centroids : numpy.ndarray
            Array of centroids (frame, y, x)
        roi_small : int
            Size of small ROI for intensity profile
        frame_indices : list
            List of frame indices to plot, if None will select evenly spaced frames
        save_name : str
            File name for saving the figure
        """
        n_frames = frames.shape[0]
        
        # Select frames to display if not specified
        if frame_indices is None:
            if n_frames <= 4:
                frame_indices = list(range(n_frames))
            else:
                # Select evenly spaced frames
                frame_indices = np.linspace(0, n_frames-1, 4, dtype=int)
        
        # Create figure with subplots
        fig, axes = plt.subplots(len(frame_indices), 3, figsize=(15, 4*len(frame_indices)))
        
        for i, frame_idx in enumerate(frame_indices):
            # Get frame and centroid
            frame = frames[frame_idx]
            y, x = centroids[frame_idx].astype(int)
            
            # Extract small ROI around particle
            roi_y_min = max(0, y - roi_small)
            roi_y