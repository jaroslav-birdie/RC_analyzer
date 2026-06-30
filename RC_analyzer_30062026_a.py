# =============================================================================
# <RC_analyzer to automatically evaluate RC curve>
#     Copyright (C) <2021>  <Jaroslav Ptacek (WG PET/CT/MRI QC EFOMP)>
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.
# =============================================================================

# =============================================================================
# 07.09.2022    - newer Python version compatibility issues solved
#               - by Matteo Bagnalasta (Nuclear Medicine, Foundation IRCCS Istituto Nazionale Tumori, Milan, Italy)
#
# 26.07.2024    - RC_analyzer calculation of COV corrected by Evon Smyth (University College Dublin, Dublin, Ireland)
# 
# 24.09.2024    - all scripts corrected. A deprecated Pydicom function read_file was replaced by dcmread function.
#               - This prevents any issues with newer Pydicom versions (thanks to Christian Bracco, Mauriziano Hospital, Turin, Italy)
#
# 10.06.2025    - Added RC_peak calculation using 3D convolution with a spherical kernel (according to PERCIST methodology and Siemens White Paper)
#               - Extended GUI and results export to include RC_peak values
#               - by Anna Kufova (Nemocnice AGEL Novy Jicin, Czech Republic); anna.kufova@nnj.agel.cz
# =============================================================================

import os

import pydicom as dcm

import tkinter as tk
import tkinter.filedialog as filedialog

from matplotlib import cm
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator, AutoMinorLocator
from matplotlib.ticker import FormatStrFormatter
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from scipy.signal import convolve 

from PIL import Image, ImageTk

import numpy as np
import pandas as pd

# =============================================================================
# Data
# =============================================================================
class phantomData:
    def __init__(self, directory, data, reference_slice, size, all_setup):
# =============================================================================
# object properties
# =============================================================================
        self.directory = directory
        self.data = data
        self.mask = [] # mask set up
        self.mask_labelled = []
        self.reference_slice = reference_slice # reference slice containing DICOM info
        self.size = size
        self.primary_zoom = []
        self.primary_mip_zoom = []
        self.peak_map = None # Attribute to store the result of the 3D convolution.

        self.all_setup = all_setup

        self.zoom = 1
        self.mask_zoom = 1

        self.position = [0,0] # phantom image position row x column
        self.mip_position = [0,0]

        self.slice = int(size[2]/2) # initial slice to open

        self.seg_thresh = 50 # initial segmentation threshold settings
        self.roi_nr = 6 # initial number of lesions
        # self.roi_nr = [] # initial number of lesions
        self.roi_obj_size = dict() # dictionary containing lesion diameters
        self.roi_bg = [0] # background ROI identification
        self.roi_bg_info = dict()  # dictionary with background info - [0] real bckg, [1] measured bckg, [2] diff(%), [3] COV(%)
        self.roi_bg_value = []
        self.roi_bg_diff = []
        self.roi_bg_COV = []
        # self.roi_bg_nr = 6
        self.roi_position = dict()
        self.roi_RC_max = dict()
        self.roi_RC_max_lbl = dict()
        self.roi_RC_A = dict()
        self.roi_RC_A_lbl = dict()
        self.roi_RC_peak = dict()
        self.roi_RC_peak_lbl = dict()
        # CRC (contrast recovery coefficient) dictionaries - Eq.2 (Sunderland et al., J Nucl Med)
        self.roi_CRC_max = dict()
        self.roi_CRC_max_lbl = dict()
        self.roi_CRC_A = dict()
        self.roi_CRC_A_lbl = dict()
        self.roi_CRC_peak = dict()
        self.roi_CRC_peak_lbl = dict()

        self.image_to_show = []
        self.mip_to_show = []
        self.thresh = 1

# =============================================================================
# object functions
# =============================================================================
    def show_me_slices(self, of_what, which, where, position, zoom, thresh):
        # temp = Image.fromarray(cm.Greys(of_what[:,:,which]/(np.max(of_what[:,:,:])/thresh), bytes=True))
        temp = Image.fromarray(cm.Greys(of_what[:,:,which]/(np.max(of_what[:,:,:])/thresh), bytes=True))     
        temp = temp.resize((int(where.winfo_width()*zoom), int(where.winfo_height()*zoom))) # zoom according to window size col x row
        self.primary_zoom = [self.size[0]/where.winfo_height(), self.size[1]/where.winfo_width()] # size row x col
        self.image_to_show = ImageTk.PhotoImage(temp)
        where.create_image(position[0],position[1], anchor='nw', image=self.image_to_show) # this is col x row

    def show_me_mip(self, of_what, where, position, zoom):
        self.mip_to_show = []
        temp = [np.sum(of_what, axis=2), np.rot90(np.sum(of_what, axis=0),1)]
        for i in range(len(temp)):
            # division by 0 exemption
            if np.max(temp[i]) != 0:
                mip = Image.fromarray(cm.Greys(temp[i]/np.max(temp[i]), bytes = True))
            else:
                mip = Image.fromarray(cm.Greys(temp[i], bytes = True))

            mip = mip.resize((int(zoom*where[i].winfo_width()),int(zoom*where[i].winfo_height()))) # zoom according to window size col x row
            self.primary_mip_zoom = [self.size[0]/where[i].winfo_height(), self.size[1]/where[i].winfo_width()] # size row x col
            self.mip_to_show.append(ImageTk.PhotoImage(mip))
            where[i].create_image(position[0],position[1], anchor='nw', image=self.mip_to_show[i])  # this is col x row

    def slice_zoom(self, event, what, zoom_point, where):
    # zoom[0] mouse scroll button rotation, zoom[1] position row, zoom[2] position col, where = where to zoom
        zoom_step = 0.5
        if (zoom_point[0] > 0):
            if self.zoom >= 5: # zoom setup = object property
                ...
            else:
                self.zoom += (120/zoom_point[0])*zoom_step
                self.position = [self.position[0]-zoom_point[1]*zoom_step, self.position[1]-zoom_point[2]*zoom_step] # position setup = object property
                self.show_me_slices(what, self.slice, where, self.position, self.zoom, self.thresh)
        elif (zoom_point[0] < 0):
            if self.zoom <= 1:
                ...
            else:
                self.position = [self.position[0]+zoom_point[1]*zoom_step, self.position[1]+zoom_point[2]*zoom_step]
                self.zoom += (120/zoom_point[0])*zoom_step
                self.show_me_slices(what, self.slice, where, self.position, self.zoom, self.thresh)

    def mip_zoom(self, event, what, zoom_point, where):
        # zoom[0] mouse scroll button rotation, zoom[1] position row, zoom[2] position col, where = where to zoom
        zoom_step = 0.5
        if (zoom_point[0] > 0):
            if self.mask_zoom >= 5: # zoom setup = object property
                ...
            else:
                self.mask_zoom += (120/zoom_point[0])*zoom_step # MIP projection created using the regiongrowed mask
                self.mip_position = [self.mip_position[0]-zoom_point[1]*zoom_step, self.mip_position[1]-zoom_point[2]*zoom_step] # position setup = object property
                self.show_me_mip(what, where, self.mip_position, self.mask_zoom)
        elif (zoom_point[0] < 0):
            if self.mask_zoom <= 1:
                ...
            else:
                self.mip_position = [self.mip_position[0]+zoom_point[1]*zoom_step, self.mip_position[1]+zoom_point[2]*zoom_step]
                self.mask_zoom += (120/zoom_point[0])*zoom_step
                self.show_me_mip(what, where, self.mip_position, self.mask_zoom)

    def slice_list(self, event, what, index, where):
    # index = scroll button rotation, where = where to zoom
        self.slice += int(120/index) # slice setup = property

        if self.slice <= 0:
            self.slice = 0

        elif self.slice >= self.size[2]-1:
            self.slice = self.size[2]-1

        self.show_me_slices(what, self.slice, where, self.position, self.zoom, self.thresh)

    def slice_click(self, event, position): # coordinates conversion - from iameged to real
        position_r = (position[0]-self.position[1])/self.zoom*self.primary_zoom[0]
        position_c = (position[1]-self.position[0])/self.zoom*self.primary_zoom[1]
        position_s = self.slice
        return(position_r, position_c, position_s)

    def slice_click_reversed(self, event, position): # coordinates conversion - from real to imaged
        position_r = position[0]*self.zoom/self.primary_zoom[0]+self.position[1]
        position_c = position[1]*self.zoom/self.primary_zoom[1]+self.position[0]
        position_s = self.slice
        return(position_r, position_c, position_s)

    def get_bg_for_segmentation(self):  
        background = self.roi_position.get('bg')
        if (self.roi_bg == [0]) or (background[0][0] == background[1][0]):
            ...
        else:
            r = np.arange(0, self.size[0])
            c = np.arange(0, self.size[1])
            mask_array = np.zeros((self.size[0], self.size[1]))

            cr = (background[1][0]+background[0][0])/2 # centre coordinates - rows
            cc = (background[1][1]+background[0][1])/2 # centre coordinates - cols

            a = (background[0][0]-background[1][0])/2 # a axis lenght - rows
            b = (background[0][1]-background[1][1])/2 # b axis lenght - cols

            mask = (r[:,np.newaxis]-cr)**2/a**2 + (c[np.newaxis,:]-cc)**2/b**2 < 1 # ellipse equation
            mask_array[mask] = 1

            temp = []

            for i in [-2, -1, 0, 1, 2]: # +-2 slices from the one where ellipse was drawn will be evaluated
                data = self.data[:,:,background[0][2]+i][mask_array == 1]
                temp.append(data)
            self.roi_bg_value = np.mean(sum(temp)/len(temp))
            self.roi_bg_info[1].set(str(np.round(self.roi_bg_value,1)))

        # calculate the volume activity at the time of measurement 
        self.ph_vol_act_src, ph_vol_act_bg = self.get_calculations()

        # convert to Bq/ml
        self.ph_vol_act_bg = 1000000*ph_vol_act_bg/float(self.all_setup['background_activity_vol'])

        self.roi_bg_diff = 100*self.roi_bg_value/self.ph_vol_act_bg-100
        self.roi_bg_diff = np.mean(100*self.roi_bg_value/self.ph_vol_act_bg-100)
        self.roi_bg_COV = 100*np.mean([(np.std(i)/np.mean(i)) for i in temp])
        
        self.roi_bg_info[0].set(str(np.round(self.ph_vol_act_bg, 1)))
        self.roi_bg_info[2].set(str(np.round(self.roi_bg_diff, 2)))
        self.roi_bg_info[3].set(str(np.round(self.roi_bg_COV, 2)))


    def segment_me_with_A(self, seed_position, background, thresh, label):
    # =============================================================================
    # - seed_position = vector containing x, y, s coordinates 
    # - background = background mean value
    # - thresh = selected threshold
    # - label = which lesion?
    # output:
    # - mask
    # - maximum position
    # =============================================================================
        im_row=self.size[0] # number of rows determination
        im_col=self.size[1] # number of dols determination
        im_slice=self.size[2] # number of slices determination

        # seed_position already exists? ... to avoid situation when regiongrow cannot be started
        try:
            seed_value=self.data[int(seed_position[0]),int(seed_position[1]),int(seed_position[2])] # seed_value setup ... will be used for the threshold calculation
        except:
            return
   
        # regiongrow neigbourhood setup - from 12 o'clock position, CW direction + up + down
        neighb=[[-1, 0, 0], [0, +1, 0], [+1, 0, 0], [0, -1 , 0], [0, 0, -1], [0, 0, 1]]

        cnt_sx=np.array([int(seed_position[0])]) # creation of the starting container for x position (rows)
        cnt_sy=np.array([int(seed_position[1])]) # creation of the starting container for y position (cols)
        cnt_sz=np.array([int(seed_position[2])]) # creation of the starting container for z position (slices)

        self.max_position = [int(seed_position[0]),int(seed_position[1]),int(seed_position[2])] # output value for the real lesion maximum search
        # the real max can be in a different positon compared to the initial seed selection

        growth=1 # initial growth setup - the regiongrowing method will finish when growth == 0

        while (growth !=0):
            temp=[]
            tempx=np.array(temp, dtype=int) # empty the tested pixel list - x position (rows)
            tempy=np.array(temp, dtype=int) # empty the tested pixel list - y position (cols)
            tempz=np.array(temp, dtype=int) # empty the tested pixel list - z position (slices)

            for i in range(0,len(cnt_sx)):
                x_t=[0, 0, 0, 0, 0, 0] # x_t variable definition - working coordinate for the pixel surrounding
                y_t=[0, 0, 0, 0, 0, 0] # y_t variable definition - working coordinate for the pixel surrounding 
                z_t=[0, 0, 0, 0, 0, 0] # z_t variable definition - working coordinate for the pixel surrounding 

                for j in range(0,6): # range 0-6 for the connectivity = 6
                    x_t[j]=cnt_sx[i]+neighb[j][0] # new x pixel coordinates calculation
                    y_t[j]=cnt_sy[i]+neighb[j][1] # new y pixel coordinates calculation
                    z_t[j]=cnt_sz[i]+neighb[j][2] # new z pixel coordinates calculation
                    # a condition for the pixel position in the image, preceding pixel impulses (volume activity) testing
                    # let's try if the comparison can be done - check the correctness of the threshold, seed point - if not correct = no regiongrowing
                    # region_grow will hit the image border, which it shouldn't - cannot continue further
                    try:
                        # tested pixel inside image
                        # tested pixel equal or larger than 0
                        # pixel was not tested before and it's value is equal or larger than the threshold condition
                        if x_t[j] <= im_row and y_t[j] <= im_col and z_t[j] <= im_slice \
                        and x_t[j] >= 0 and y_t[j] >= 0 and z_t[j] >= 0 \
                        and self.mask[x_t[j], y_t[j] , z_t[j]] == 0 and self.data[x_t[j], y_t[j], z_t[j]] >= \
                        thresh*(seed_value+self.roi_bg_value):
                            self.mask[x_t[j], y_t[j], z_t[j]]=1
                            self.mask_labelled[x_t[j], y_t[j], z_t[j]]= label # if condition is met a number corresponding to the lesion label is recorded in the mask_labelled
                            tempx=np.append(tempx,[x_t[j]]) # bin for pixels fulfilling the condition prolonged - x coordinate
                            tempy=np.append(tempy,[y_t[j]]) # bin for pixels fulfilling the condition prolonged - y coordinate
                            tempz=np.append(tempz,[z_t[j]]) # bin for pixels fulfilling the condition prolonged - z coordinate
                            if self.data[x_t[j], y_t[j], z_t[j]] > seed_value: # check if the initial seed was the max the lesion
                                self.max_position = [x_t[j], y_t[j], z_t[j]]
                                seed_value = self.data[x_t[j], y_t[j], z_t[j]] # if not, the new max is the new seed
                    except:
                        tk.messagebox.showwarning("Warning", "Something went wrong - reset segmentation and continue!")
                        return
                    
        

            growth=len(tempx) # growth determination - how many new pixels were found?

            cnt_sx=tempx.copy() # new starting container now contains the last recorded pixels from the mask
            cnt_sy=tempy.copy()
            cnt_sz=tempz.copy()

        return self.max_position
    
    
    
#--------------------------------------------------------------------------
# Creates and normalizes 3D convolution kernel representing a 1cm^3 sphere.
#
# This method uses a dynamic sub-sampling approach to accurately account for
# partial volume effects at the sphere's boundary. The resulting kernel is
# designed to be convolved with an entire PET image volume to efficiently 
# generate a complete peak map.
#
# The algorithm is based on the methodology described in the Siemens 
# Healthineers whitepaper on SUVpeak calculation.
#--------------------------------------------------------------------------

    def create_spherical_kernel(self, dx, dy, dz):
    
    # 1. Preparing a zero 3D kernel 
        
        radius_mm = (3 * 1000.0 / (4 * np.pi))**(1/3.0) # radius of 1 cm^3 sphere
        radius_sq_mm = radius_mm**2 # optimization for further distance calculations

        # Define the kernel's voxel-bounding box, initialize it, and find its center (+1 for exactly one central voxel)
        # dx, dy, dy - voxel size
        k_size_x = int(np.ceil(radius_mm / dx)) * 2 + 1
        k_size_y = int(np.ceil(radius_mm / dy)) * 2 + 1
        k_size_z = int(np.ceil(radius_mm / dz)) * 2 + 1

        # Convolution 3D-kernel creation, filled with zeros
        kernel = np.zeros((k_size_y, k_size_x, k_size_z))

        # Get the indices for the geometric center of the kernel 
        center_y, center_x, center_z = (k_size_y - 1) / 2, (k_size_x - 1) / 2, (k_size_z - 1) / 2

    # 2. Supersampling for edge voxels
       
        # Determine subdivision count by repeatedly halving voxel dimensions until the target sub-voxel size (<= 0.5 mm) is met
        N_sub_x, N_sub_y, N_sub_z = 1, 1, 1
        while (dx / N_sub_x) > 0.5:
            N_sub_x *= 2 
        while (dy / N_sub_y) > 0.5:
            N_sub_y *= 2
        while (dz / N_sub_z) > 0.5:
            N_sub_z *= 2
            
        N_sub_total = N_sub_x * N_sub_y * N_sub_z # Total number of subvoxels per voxel
        
        # Subvoxel size
        sub_dx = dx / N_sub_x
        sub_dy = dy / N_sub_y
        sub_dz = dz / N_sub_z


    # 3. Generate coordinates for subvoxel centers along each axis

        x_coords = np.linspace(-0.5*dx + 0.5*sub_dx, 0.5*dx - 0.5*sub_dx, N_sub_x)
        y_coords = np.linspace(-0.5*dy + 0.5*sub_dy, 0.5*dy - 0.5*sub_dy, N_sub_y)
        z_coords = np.linspace(-0.5*dz + 0.5*sub_dz, 0.5*dz - 0.5*sub_dz, N_sub_z)
        
        
        # Create 3D coordinate grids from the 1D axes, then flatten them for efficient vectorization.
        grid_x, grid_y, grid_z = np.meshgrid(x_coords, y_coords, z_coords, indexing='ij')

        sub_x_offsets_mm = grid_x.flatten()
        sub_y_offsets_mm = grid_y.flatten()
        sub_z_offsets_mm = grid_z.flatten()

    # 4. Kernel Weight Calculation 

        # Iterate over each voxel in the kernel's grid.
        for r in range(k_size_y):
            for c in range(k_size_x):
                for s in range(k_size_z):

                    # Get the physical center of the current kernel voxel, relative to the origin.
                    voxel_center_x_mm = (c - center_x) * dx
                    voxel_center_y_mm = (r - center_y) * dy
                    voxel_center_z_mm = (s - center_z) * dz

                    # Calculate the absolute coordinates for all corresponding sub-voxel centers.
                    points_x = voxel_center_x_mm + sub_x_offsets_mm
                    points_y = voxel_center_y_mm + sub_y_offsets_mm
                    points_z = voxel_center_z_mm + sub_z_offsets_mm

                    # Calculate the squared distance of each subvoxel from the kernel's origin.
                    dist_sq_points = points_x**2 + points_y**2 + points_z**2
                    
                    # The weight of the voxel is the fraction of subvoxels whose centers are inside the sphere.
                    n_inside = np.sum(dist_sq_points <= radius_sq_mm)
                    weight = n_inside / N_sub_total
                    
                    # Overwrite the zero kernel
                    if weight > 0:
                        kernel[r, c, s] = weight

        # Kernel normalization so that the sum of the weights is 1
        kernel_sum = np.sum(kernel)
        if kernel_sum > 0:
            kernel /= kernel_sum
        
        return kernel
    
    
#-------------------------------------------------------------------------------------------
# Pre-computes the peak map for the entire image volume using a single 3D convolution.
# The result is cached in self.peak_map for fast subsequent access.
#-------------------------------------------------------------------------------------------
    
    def precompute_peak_map(self):

        # 1. Getting pixel dimensions from DICOM header
        try:
            if not (hasattr(self, 'reference_slice') and self.reference_slice and
                    hasattr(self.reference_slice, 'PixelSpacing') and len(self.reference_slice.PixelSpacing) == 2 and
                    hasattr(self.reference_slice, 'SliceThickness')):
                tk.messagebox.showerror("DICOM Error", "Missing voxel dimensions. Cannot pre-compute SUVpeak map.")
                return

            dy_vox_dim = float(self.reference_slice.PixelSpacing[0])
            dx_vox_dim = float(self.reference_slice.PixelSpacing[1])
            dz_vox_dim = float(self.reference_slice.SliceThickness)

        # 2. Create the spherical convolution 3D-kernel.
            kernel = self.create_spherical_kernel(dx_vox_dim, dy_vox_dim, dz_vox_dim)
            
        # 3. Perform the 3D convolution and cache the resulting map.
            print("Performing 3D convolution... (This may take a moment)")
            self.peak_map = convolve(self.data, kernel, mode='same', method='auto')
            print("Peak map has been successfully pre-computed.")

        except Exception as e:
            self.peak_map = None 
            import traceback
            error_msg = f"Failed to pre-compute peak map: {e}\n{traceback.format_exc()}"
            print(error_msg)
            tk.messagebox.showerror("Computation Error", error_msg)

        
#-------------------------------------------------------------------------------------------
# This method calculates all recovery coefficients (Max, A, Peak) for a single
# segmented lesion and updates the GUI with the results.
# 'where': The full 3D image data (self.data).
# 'mask_labelled': The 3D mask with unique integer labels for each lesion.
# 'which_label': The integer label of the lesion currently being processed.
#-------------------------------------------------------------------------------------------

    def get_results(self, where, mask_labelled, which_label):

    # 1. Initial Data Extraction and Validation 
        lesion_mask = (mask_labelled == which_label)
        lesion_values = where[lesion_mask]

        if lesion_values.size == 0:
            tk.messagebox.showwarning("Warning", f"Segmented lesion {which_label} is empty!")
            self.roi_RC_max[which_label] = np.nan
            self.roi_RC_A[which_label] = np.nan
            self.roi_RC_peak[which_label] = np.nan
            self.roi_CRC_max[which_label] = np.nan
            self.roi_CRC_A[which_label] = np.nan
            self.roi_CRC_peak[which_label] = np.nan
            if which_label in self.roi_RC_max_lbl: self.roi_RC_max_lbl[which_label].set("NaN")
            if which_label in self.roi_RC_A_lbl: self.roi_RC_A_lbl[which_label].set("NaN")
            if which_label in self.roi_RC_peak_lbl: self.roi_RC_peak_lbl[which_label].set("NaN")
            if which_label in self.roi_CRC_max_lbl: self.roi_CRC_max_lbl[which_label].set("NaN")
            if which_label in self.roi_CRC_A_lbl: self.roi_CRC_A_lbl[which_label].set("NaN")
            if which_label in self.roi_CRC_peak_lbl: self.roi_CRC_peak_lbl[which_label].set("NaN")
            return

    # 2. Basic Metric Calculation 
        av_max_val = np.max(lesion_values) # maximal volume activity in the segmented lesion
        av_A = np.mean(lesion_values) # mean
        av_bg = self.roi_bg_value # mean background
        av_peak = np.nan # peak

    # 3. SUVpeak Calculation from Pre-computed Map
        if self.peak_map is None:
            tk.messagebox.showerror("Error", "Peak map is not available. Cannot calculate RC_peak.")
        else:    
            qualifying_voxels_mask = lesion_mask 
            av_peak = np.max(self.peak_map[qualifying_voxels_mask])

    # 4. Calculation of real volume activity ratio between lesions and background (ph_ratio)
        ph_ratio = np.nan 

        if not (self.all_setup.get('sources_time') and len(str(self.all_setup['sources_time'])) == 6 and
                self.all_setup.get('sources_residual_time') and len(str(self.all_setup['sources_residual_time'])) == 6 and
                self.all_setup.get('background_time') and len(str(self.all_setup['background_time'])) == 6 and
                self.all_setup.get('background_residual_time') and len(str(self.all_setup['background_residual_time'])) == 6):
            tk.messagebox.showwarning("Warning", "Wrong time format in setup - RC evaluation might be incorrect or impossible!")
        else:
            try:
                #  volume activities of lesions / background at the time of the measurement calculation
                ph_vol_act_src, ph_vol_act_bg = self.get_calculations()
                
                src_activity_vol = float(self.all_setup.get('sources_activity_vol', 0))
                bg_activity_vol = float(self.all_setup.get('background_activity_vol', 0))

                if src_activity_vol == 0 or bg_activity_vol == 0 or ph_vol_act_bg == 0:
                    print(f"Warning for lesion {which_label}: Source/BG activity volume or calculated BG activity is zero. ph_ratio will be NaN.")
                
                else:
                    # lesion to background ratio as filled
                    ph_ratio = (ph_vol_act_src / src_activity_vol) / (ph_vol_act_bg / bg_activity_vol)
            except Exception as e_ratio:
                print(f"Error calculating ph_ratio for lesion {which_label}: {e_ratio}")

    # 5. Final RC and CRC Calculation
    # RC (recovery coefficient)  = (C_sphere / C_Bkg) / (A_sphere / A_Bkg)
    # CRC (contrast recovery)    = ((C_sphere - C_Bkg) / C_Bkg) / ((A_sphere - A_Bkg) / A_Bkg)   ... Eq.2
    #                            = ((C_sphere / C_Bkg) - 1) / (ph_ratio - 1)
    # with C_sphere = av_max_val / av_A / av_peak, C_Bkg = av_bg and ph_ratio = A_sphere / A_Bkg
        valid_av_bg = isinstance(av_bg, (int, float)) and av_bg != 0
        # CRC denominator - the known (true) contrast of the phantom filling; ph_ratio == 1 means no contrast
        crc_denom = ph_ratio - 1 if not np.isnan(ph_ratio) else np.nan
        if np.isnan(ph_ratio) or ph_ratio == 0 or not valid_av_bg:
            self.roi_RC_max[which_label] = np.nan
            self.roi_RC_A[which_label] = np.nan
            self.roi_RC_peak[which_label] = np.nan
        else:
            self.roi_RC_max[which_label] = (av_max_val / av_bg / ph_ratio)
            self.roi_RC_A[which_label] = (av_A / av_bg / ph_ratio)
            self.roi_RC_peak[which_label] = (av_peak / av_bg / ph_ratio) if not np.isnan(av_peak) else np.nan

        if np.isnan(crc_denom) or crc_denom == 0 or not valid_av_bg:
            self.roi_CRC_max[which_label] = np.nan
            self.roi_CRC_A[which_label] = np.nan
            self.roi_CRC_peak[which_label] = np.nan
        else:
            self.roi_CRC_max[which_label] = (av_max_val / av_bg - 1) / crc_denom
            self.roi_CRC_A[which_label] = (av_A / av_bg - 1) / crc_denom
            self.roi_CRC_peak[which_label] = ((av_peak / av_bg - 1) / crc_denom) if not np.isnan(av_peak) else np.nan

    # 6. GUI Update
        rc_max_to_display = self.roi_RC_max.get(which_label, np.nan)
        rc_a_to_display = self.roi_RC_A.get(which_label, np.nan)
        rc_peak_to_display = self.roi_RC_peak.get(which_label, np.nan)

        if which_label in self.roi_RC_max_lbl: self.roi_RC_max_lbl[which_label].set(f"{rc_max_to_display:.4f}" if not np.isnan(rc_max_to_display) else "NaN")
        if which_label in self.roi_RC_A_lbl: self.roi_RC_A_lbl[which_label].set(f"{rc_a_to_display:.4f}" if not np.isnan(rc_a_to_display) else "NaN")
        if which_label in self.roi_RC_peak_lbl: self.roi_RC_peak_lbl[which_label].set(f"{rc_peak_to_display:.4f}" if not np.isnan(rc_peak_to_display) else "NaN")

        crc_max_to_display = self.roi_CRC_max.get(which_label, np.nan)
        crc_a_to_display = self.roi_CRC_A.get(which_label, np.nan)
        crc_peak_to_display = self.roi_CRC_peak.get(which_label, np.nan)

        if which_label in self.roi_CRC_max_lbl: self.roi_CRC_max_lbl[which_label].set(f"{crc_max_to_display:.4f}" if not np.isnan(crc_max_to_display) else "NaN")
        if which_label in self.roi_CRC_A_lbl: self.roi_CRC_A_lbl[which_label].set(f"{crc_a_to_display:.4f}" if not np.isnan(crc_a_to_display) else "NaN")
        if which_label in self.roi_CRC_peak_lbl: self.roi_CRC_peak_lbl[which_label].set(f"{crc_peak_to_display:.4f}" if not np.isnan(crc_peak_to_display) else "NaN")

    def get_calculations(self):
        # calculate sources activity

        time_delta = (int(self.all_setup['series_time'][0:2])*3600+int(self.all_setup['series_time'][2:4])*60+int(self.all_setup['series_time'][4:6])) - \
        (int(self.all_setup['sources_time'][0:2])*3600+int(self.all_setup['sources_time'][2:4])*60+int(self.all_setup['sources_time'][4:6]))

        time_delta_res = (int(self.all_setup['series_time'][0:2])*3600+int(self.all_setup['series_time'][2:4])*60+int(self.all_setup['series_time'][4:6])) - \
        (int(self.all_setup['sources_residual_time'][0:2])*3600+int(self.all_setup['sources_residual_time'][2:4])*60+int(self.all_setup['sources_residual_time'][4:6]))

        res_activity = float(self.all_setup['sources_residual_activity'])*np.exp(-np.log(2)*time_delta_res/float(self.all_setup['halflife']))

        ph_vol_act_src = float(self.all_setup['sources_activity'])*np.exp(-np.log(2)*time_delta/float(self.all_setup['halflife'])) - res_activity

        # calculate background activity
        time_delta = (int(self.all_setup['series_time'][0:2])*3600+int(self.all_setup['series_time'][2:4])*60+int(self.all_setup['series_time'][4:6])) - \
        (int(self.all_setup['background_time'][0:2])*3600+int(self.all_setup['background_time'][2:4])*60+int(self.all_setup['background_time'][4:6]))

        time_delta_res = (int(self.all_setup['series_time'][0:2])*3600+int(self.all_setup['series_time'][2:4])*60+int(self.all_setup['series_time'][4:6])) - \
        (int(self.all_setup['background_residual_time'][0:2])*3600+int(self.all_setup['background_residual_time'][2:4])*60+int(self.all_setup['background_residual_time'][4:6]))

        res_activity = float(self.all_setup['background_residual_activity'])*np.exp(-np.log(2)*time_delta_res/float(self.all_setup['halflife']))

        ph_vol_act_bg = float(self.all_setup['background_activity'])*np.exp(-np.log(2)*time_delta/float(self.all_setup['halflife'])) - res_activity

        return(ph_vol_act_src, ph_vol_act_bg)

    def show_me_results(self, where):
        # data preparation for graph creation - only already segmented ROI
        obj_diam = np.zeros(self.roi_nr)
        RC_max = obj_diam.copy()
        RC_A = obj_diam.copy()
        RC_peak = obj_diam.copy()
        CRC_max = obj_diam.copy()
        CRC_A = obj_diam.copy()
        CRC_peak = obj_diam.copy()
        # list prints all items in the corresponding dictionary
        for i in list(self.roi_RC_max):

            obj_diam[i-1] = self.roi_obj_size[i].get()
            RC_max[i-1] = self.roi_RC_max[i]
            RC_A[i-1] = self.roi_RC_A[i]
            RC_peak[i-1] = self.roi_RC_peak[i]
            CRC_max[i-1] = self.roi_CRC_max.get(i, 0)
            CRC_A[i-1] = self.roi_CRC_A.get(i, 0)
            CRC_peak[i-1] = self.roi_CRC_peak.get(i, 0)

        # x-axis (diameters) shared by RC and CRC - keep only segmented lesions
        seg = obj_diam != 0
        obj_diam = obj_diam[seg]
        RC_max = RC_max[seg]
        RC_A = RC_A[seg]
        RC_peak = RC_peak[seg]
        CRC_max = CRC_max[seg]
        CRC_A = CRC_A[seg]
        CRC_peak = CRC_peak[seg]

        GraphRC(where, (RC_analyzer.Xsize, RC_analyzer.Ysize), [0,0], self.seg_thresh, [obj_diam,RC_max],[obj_diam,RC_A], [obj_diam,RC_peak], \
                [obj_diam,CRC_max], [obj_diam,CRC_A], [obj_diam,CRC_peak], RC_analyzer.RC[0], RC_analyzer.RC[1], RC_analyzer.RC[2], \
                RC_analyzer.config, RC_analyzer.RC_limit_max, RC_analyzer.RC_limit_A, RC_analyzer.RC_limit_peak)

    def get_out_results(self, file_path):
        # export of results
        temp = [['RC_A' + str(self.seg_thresh), 'RC_max', 'RC_peak', \
                 'CRC_A' + str(self.seg_thresh), 'CRC_max', 'CRC_peak', \
                 'diameter', 'bg_real', 'bg_measured', 'bg_diff(%)', 'bg_COV(%)']]
        for i in range(1, len(self.roi_RC_A)+1):
            try:
                temp.append([self.roi_RC_A[i], self.roi_RC_max[i], self.roi_RC_peak[i], \
                             self.roi_CRC_A[i], self.roi_CRC_max[i], self.roi_CRC_peak[i], \
                             self.roi_obj_size[i].get()])
            except:
                tk.messagebox.showwarning("Warning", 'Segmentation of ROI ' + str(i) + ' is missing!')
                return

        temp.append(['', '', '', '', '', '', self.ph_vol_act_bg, self.roi_bg_value, self.roi_bg_diff, self.roi_bg_COV])

        file_name = self.all_setup['device'] + '_' + self.all_setup['study_date'][0:2] + self.all_setup['study_date'][3:5] + \
            self.all_setup['study_date'][6:10] + 'result.csv'
        file_name = file_name.replace(" ", "") # delete spaces
        file_name = file_name.replace("/", "_") # delete /

        file_dir = filedialog.asksaveasfilename(parent=root, initialfile=file_name, initialdir="/",title='Please select a directory to save')

        if file_dir: # check if filename was stated
            pd.DataFrame(temp).to_csv(file_dir, sep=';')
            tk.messagebox.showinfo("Export", "Export finished!")
        else:
            ...


# =============================================================================
# Graph creation
# =============================================================================
class GraphRC:
    def __init__(self, which_frame, graph_size, graph_position, seg_thresh, graph_RCmax, graph_RCA, graph_RCpeak, graph_CRCmax, graph_CRCA, graph_CRCpeak, RC_max, RC_A, RC_peak, config, RC1, RC2, RC3):

        fig = Figure(figsize=(graph_size[0]/200, graph_size[1]/(200)), dpi = 200) 
        axs = fig.subplots(3) 
        fig.subplots_adjust(hspace=0.3, wspace = 0, top=0.95, bottom=0.15, left = 0.15, right=0.98)
        idx_1, idx_2, idx_3 = 0, 0, 0 

        # initial values of graph setup 
        RC1_obj_diam, RC1_low_max, RC1_up_max, RC2_obj_diam, RC2_low_max, RC2_up_max, RC3_obj_diam, RC3_low_max, RC3_up_max  = [0],[0],[0],[0],[0],[0],[0],[0],[0]

        # definiton of the selected RC curve region - from the config.cfg file
        if RC1.get():  
            for i in range(0, len(config)):
                if config[i][0] == RC1.get():
                    idx_1 = i
                if config[i][0] == RC2.get():
                    idx_2 = i
                if config[i][0] == RC3.get():
                    idx_3 = i

            # conver to float
            try: 
                RC1_obj_diam = list(map(float,config[idx_1][1].split(";")))
                RC1_low_max = list(map(float,config[idx_1][2].split(";")))
                RC1_up_max = list(map(float,config[idx_1][3].split(";")))

                RC2_obj_diam = list(map(float,config[idx_2][1].split(";")))
                RC2_low_max = list(map(float,config[idx_2][2].split(";")))
                RC2_up_max = list(map(float,config[idx_2][3].split(";")))

                RC3_obj_diam = list(map(float,config[idx_3][1].split(";")))
                RC3_low_max = list(map(float,config[idx_3][2].split(";")))
                RC3_up_max = list(map(float,config[idx_3][3].split(";")))

            except:
                tk.messagebox.showwarning("Warning", "Wrong format of RC boundary values! Correct configuration file and load PET data again")
                return

        if len(RC1_obj_diam) == len(RC1_low_max) and len(RC1_low_max) == len(RC1_up_max) and \
            len(RC2_obj_diam) == len(RC2_low_max) and len(RC2_low_max) == len(RC2_up_max) and \
            len(RC3_obj_diam) == len(RC3_low_max) and len(RC3_low_max) == len(RC3_up_max):  

            title_font = {'fontname':'Arial', 'size':'5'}  
            axis_font = {'fontname':'Arial', 'size':'5'} 
            axis_tick = {'labelsize':'4'} 
            data_line = {'color':'black', 'marker':'o', 'markersize':'2', 'linestyle':'-', 'linewidth':'0.5'}
            crc_line = {'color':'blue', 'marker':'s', 'markersize':'2', 'linestyle':'-', 'linewidth':'0.5'}
            RC_line = { 'marker':'', 'linestyle':'--', 'linewidth':'1'}

            # RC (black circles) and CRC (blue squares) plotted together on the same axes
            axs[0].plot(graph_RCmax[0],graph_RCmax[1], label='RC', **data_line)
            axs[1].plot(graph_RCA[0],graph_RCA[1], label='RC', **data_line)
            axs[2].plot(graph_RCpeak[0],graph_RCpeak[1], label='RC', **data_line)

            axs[0].plot(graph_CRCmax[0],graph_CRCmax[1], label='CRC', **crc_line)
            axs[1].plot(graph_CRCA[0],graph_CRCA[1], label='CRC', **crc_line)
            axs[2].plot(graph_CRCpeak[0],graph_CRCpeak[1], label='CRC', **crc_line)

            axs[0].set_title('RC_max / CRC_max', y=0.95,  **title_font)
            axs[1].set_title('RC_A' + str(seg_thresh) + ' / CRC_A' + str(seg_thresh), y=0.95, **title_font)
            axs[2].set_title('RC_peak / CRC_peak', y=0.95, **title_font)
            axs[2].set_xlabel('lesion diameter (mm)', **axis_font)

            if RC_max.get() == 1: # is RC_max limits selected?
                axs[0].plot(RC1_obj_diam,RC1_low_max, color="red", **RC_line)
                axs[0].plot(RC1_obj_diam,RC1_up_max, color="red", **RC_line)

            if RC_A.get() == 1: # is RC_A limits selected?
                axs[1].plot(RC2_obj_diam,RC2_low_max, color="red", **RC_line)
                axs[1].plot(RC2_obj_diam,RC2_up_max, color="red", **RC_line)

            if RC_peak.get() == 1: # is RC_peak limits selected?
                axs[2].plot(RC3_obj_diam,RC3_low_max, color="red", **RC_line)
                axs[2].plot(RC3_obj_diam,RC3_up_max, color="red", **RC_line) 
    

            for i in range(0, len(axs)):  
                axs[i].autoscale(enable=True) 
                axs[i].set_ylabel('RC / CRC',**axis_font)
                axs[i].tick_params(axis = 'both', **axis_tick)
                axs[i].legend(loc='lower right', fontsize=4)
                axs[i].xaxis.set_major_locator(MaxNLocator(integer=True)) 
                axs[i].xaxis.set_major_formatter(FormatStrFormatter('%.0f')) 
                axs[i].xaxis.set_minor_locator(AutoMinorLocator()) 
                axs[i].yaxis.set_minor_locator(AutoMinorLocator())

                axs[i].set_xlim(left = 0) 
                axs[i].set_ylim(bottom = 0) 

                axs[i].set_ylim(bottom = 0) 
                axs[i].grid(which = 'major', color = 'grey', linestyle = '-', linewidth = '0.25') 
                axs[i].grid(which = 'minor', color = 'grey', linestyle = '--', linewidth = '0.20') 

            axs[0].axes.xaxis.set_ticklabels([]) 
            axs[1].axes.xaxis.set_ticklabels([]) 
            axs[2].tick_params(axis='both', which='major', pad=1) 

            self.canvas = FigureCanvasTkAgg(fig, master = which_frame) 
            # navigation toolbar
            frame = tk.Frame(which_frame) 
            frame.grid(row = graph_position[0]+1, column=0) 
            toolbar = NavigationToolbar2Tk(self.canvas, frame).update() 
            self.canvas.get_tk_widget().grid(row = graph_position[0], column = graph_position[1]) 
        else:
             tk.messagebox.showwarning("Warning", 'Wrong lenght of RC boundary values! Correct configuration file and load PET data again')
             return 

# =============================================================================
# MainWindow
# =============================================================================
             
class mainWindow:
    def __init__(self, master):
          
        self.master = master

        self.master.geometry('+1915+0') # this can be use to place the window on a particular coordinates of the computer workspace (dual monitors, etc.)        
        self.master.title("RC analyzer") 
        self.master.resizable(0,0) # window is not resizable 

        self.Xpix = master.winfo_screenwidth() # get the resolution of the whole desktop x (šířka) (be aware - if two same sized monitors are used this means double the size of one)
        self.Ypix = master.winfo_screenheight() # get the resolution of the whole desktop y (výška)

        self.Xsize = (self.Xpix-350)/3 # set the window size according to the desktop resolution x - šířka okna
        self.Ysize = self.Xsize # set the window size according to the desktop resolution y - délka okna

        self.RC = [tk.IntVar(),tk.IntVar(), tk.IntVar()] # setup the RC curve boundaries lists 
        self.config = [] # preprare for opening the config.cfg file 
        self.RC_limit_max = tk.StringVar()
        self.RC_limit_A = tk.StringVar()
        self.RC_limit_peak = tk.StringVar()

        self.plus = tk.PhotoImage(file = r"plus.png") # open icon files - plusko a minusko pro změnu prahu zobrazení
        self.minus = tk.PhotoImage(file = r"minus.png")

        # inititate Tkinter frames
      
        self.frame_a = tk.LabelFrame(master, text='Data') 
        self.frame_b = tk.LabelFrame(master, text='MIPs')
        self.frame_c = tk.LabelFrame(master, text='RC curves')
        self.frame_d = tk.LabelFrame(master, text='Load data, DICOM info')
        self.frame_e = tk.LabelFrame(master, text='Phantom, segmentation and RC boundaries configuration')
        self.frame_f = tk.LabelFrame(master, text='Experiment information')
        self.frame_g = tk.LabelFrame(master, text='ROI data')
        self.frame_h = tk.LabelFrame(master, text='Image threshold')

        # =============================================================================
        # frame_a 
        # =============================================================================
        self.phantom_slice_view = tk.Canvas(self.frame_a, width=self.Xsize, height=self.Ysize, borderwidth=2, bg='white') 
        self.phantom_slice_view.config(highlightbackground='black', highlightthickness='2') 
        self.phantom_slice_view.grid(row=0, column=0)  

        self.phantom_slice_lbl = tk.Label(self.frame_a, text='slice: ??', anchor='nw', justify='left')
        self.phantom_slice_lbl.grid(row=0, column=0, sticky='nw') 

        # =============================================================================
        # frame_b 
        # =============================================================================
           
        self.mip_UP_view = tk.Canvas(self.frame_b, width=self.Xsize/2.5, height=self.Ysize/2.5, borderwidth=0, bg='white')
        self.mip_UP_view.config(highlightbackground='black', highlightthickness='1') 
        self.mip_UP_view.grid(row=0, column=1, sticky='nsew') 

        self.mip_LEFT_view = tk.Canvas(self.frame_b, width=self.Xsize/2.5, height=self.Ysize/2.5, borderwidth=0, bg='white') 
        self.mip_LEFT_view.config(highlightbackground='black', highlightthickness='1') 
        self.mip_LEFT_view.grid(row=0, column=2, sticky='nsew') 
        

        # =============================================================================
        # frame_c 
        # =============================================================================
        # startovni zobrazeni grafu 
        GraphRC(self.frame_c, (self.Xsize, self.Ysize), [0,0], ' ',  [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], self.RC[0], self.RC[1], self.RC[2], self.config, self.RC_limit_max, self.RC_limit_A, self.RC_limit_peak)


        # =============================================================================
        # frame_d 
        # =============================================================================
        button = tk.Button(self.frame_d, text="Open data", command=self.open_data)
        button.grid(row=1, column=0, sticky='w') 
        self.data_size = tk.Label(self.frame_d, text='\n\n\n\n\n\n\n\n\n\n\n', anchor='sw', justify='left') # místo pro DICOM tagy
        self.data_size.grid(row=2, column=0)

        # =============================================================================
        # frame_e
        # =============================================================================
        self.roi_nr_lbl = tk.Label(self.frame_e, width = 1)
        self.roi_nr_lbl.grid(row=0, column=0, sticky='nsew')
        # self.roi_nr_lbl['text'] = '???'
        self.roi_nr_lbl['text'] = '6'
        label = tk.Label(self.frame_e, text=' lesions (1-12)', anchor='sw', justify='left').grid(row=0, column=2, sticky='nsew')

        button = tk.Button(self.frame_e, image = self.minus, width="20", height="20", command=lambda where = self.frame_g: self.phantom_rois_nr(-1, where)).grid(row=0, column=3)
        button = tk.Button(self.frame_e, image = self.plus, width="20", height="20", command=lambda where = self.frame_g: self.phantom_rois_nr(1, where)).grid(row=0, column=4)

        self.segmentation_thresh_lbl = tk.Label(self.frame_e, width = 1)
        self.segmentation_thresh_lbl.grid(row=1, column=0, sticky='nsew')
        self.segmentation_thresh_lbl['text'] = '???'
        label = tk.Label(self.frame_e, text=' regiongrow thresh (%)', anchor='sw', justify='left').grid(row=1, column=2, sticky='nsew')

        button =  tk.Button(self.frame_e, image = self.minus, width="20", height="20", command=lambda where = self.frame_g: self.seg_thresh(-1, where)).grid(row=1, column=3)
        button =  tk.Button(self.frame_e, image = self.plus, width="20", height="20", command=lambda where = self.frame_g: self.seg_thresh(1, where)).grid(row=1, column=4)

        label = tk.Label(self.frame_e, text='').grid(row=2, column=0, columnspan=4)
        label = tk.Label(self.frame_e, text='Show RC curve boundaries').grid(row=3, column=0, columnspan=4)

        ch_button = tk.Checkbutton(self.frame_e, text = "RC_max", variable = self.RC[0], onvalue = 1, offvalue = 0, command = self.RCs_redraw).grid(row=4, column = 0, columnspan=2)
        ch_button = tk.Checkbutton(self.frame_e, text = "RC_A50", variable = self.RC[1], onvalue = 1, offvalue = 0, command = self.RCs_redraw).grid(row=4, column = 2, columnspan=2)
        ch_button = tk.Checkbutton(self.frame_e, text = "RC_peak", variable = self.RC[2], onvalue = 1, offvalue = 0, command = self.RCs_redraw).grid(row=4, column = 4, columnspan=2)

        self.RC_limit_spin_max = tk.Spinbox(master=self.frame_e, textvariable=self.RC_limit_max, wrap=True, command = self.RCs_redraw)
        self.RC_limit_spin_max.grid(row=5, column = 0, columnspan = 2, sticky='w')
        self.RC_limit_spin_A = tk.Spinbox(master=self.frame_e, textvariable=self.RC_limit_A, wrap=True, command = self.RCs_redraw)
        self.RC_limit_spin_A.grid(row=5, column = 2, columnspan = 2, sticky='w')
        self.RC_limit_spin_peak = tk.Spinbox(master=self.frame_e, textvariable=self.RC_limit_peak, wrap=True, command = self.RCs_redraw)
        self.RC_limit_spin_peak.grid(row=5, column = 4, columnspan = 2, sticky='w')



        label = tk.Label(self.frame_e, text='').grid(row=6, column=0, columnspan=4)

        button = tk.Button(self.frame_e, text='IEC Body Phantom', command = lambda: self.roi_character_autofill(self.phantom, 6))
        button.grid(row=7, column=0, columnspan=5, sticky='nsew')
        
        button = tk.Button(self.frame_e, text='IEC_12 Body Phantom', command = lambda: self.roi_character_autofill(self.phantom, 12))
        button.grid(row=8, column=0, columnspan=5, sticky='nsew')

        # =============================================================================
        # frame_f
        # =============================================================================
        width = 10

        label = tk.Label(self.frame_f, text='sources', anchor='center', justify='center').grid(row=0, column=0, columnspan=2)
        label = tk.Label(self.frame_f, text='background', anchor='center', justify='center').grid(row=0, column=2, columnspan=2)

        label = tk.Label(self.frame_f, text="activity (MBq)", anchor='w', justify = 'left').grid(row=1, column=0, sticky="nsew")
        self.data_srcs_act = tk.Entry(self.frame_f, width = width, exportselection=0)
        self.data_srcs_act.insert(0, "???")
        # self.data_srcs_act.insert(0, "10")

        label = tk.Label(self.frame_f, text="time [hhmmss]", anchor='w', justify = 'left').grid(row=2, column=0, sticky="nsew")
        self.data_srcs_time = tk.Entry(self.frame_f, width = width, exportselection=0)
        self.data_srcs_time.insert(0, "hhmmss")
        # self.data_srcs_time.insert(0, "100000")

        label = tk.Label(self.frame_f, text="residual activity (MBq)", anchor='w', justify = 'left').grid(row=3, column=0, sticky="nsew")
        self.data_srcs_res_act = tk.Entry(self.frame_f, width = width, exportselection=0)
        self.data_srcs_res_act.insert(0, "0")

        label = tk.Label(self.frame_f, text="residual time [hhmmss]", anchor='w', justify = 'left').grid(row=4, column=0, sticky="nsew")
        self.data_srcs_res_time = tk.Entry(self.frame_f, width = width, exportselection=0)
        self.data_srcs_res_time.insert(0, "hhmmss")
        # self.data_srcs_res_time.insert(0, "100000")

        label = tk.Label(self.frame_f, text="total dilution volume (ml)", anchor='w', justify = 'left').grid(row=5, column=0, sticky="nsew")
        self.data_srcs_volume_tot = tk.Entry(self.frame_f, width = width, exportselection=0)
        self.data_srcs_volume_tot.insert(0, "1000")

        self.data_srcs_act.grid(row=1, column=1)
        self.data_srcs_time.grid(row=2, column=1)
        self.data_srcs_res_act.grid(row=3, column=1)
        self.data_srcs_res_time.grid(row=4, column=1)
        self.data_srcs_volume_tot.grid(row=5, column=1)

        label = tk.Label(self.frame_f, text="   activity (MBq)", anchor='w', justify = 'left').grid(row=1, column=2, sticky="nsew")
        self.data_bcg_act = tk.Entry(self.frame_f, width = width, exportselection=0)
        self.data_bcg_act.insert(0, "???")
        # self.data_bcg_act.insert(0, "10")

        label = tk.Label(self.frame_f, text="   time [hhmmss]", anchor='w', justify = 'left').grid(row=2, column=2, sticky="nsew")
        self.data_bcg_time = tk.Entry(self.frame_f, width = width, exportselection=0)
        self.data_bcg_time.insert(0, "hhmmss")
        # self.data_bcg_time.insert(0, "100000")

        label = tk.Label(self.frame_f, text="   residual acitivity (MBq)", anchor='w', justify = 'left').grid(row=3, column=2, sticky="nsew")
        self.data_bcg_res_act = tk.Entry(self.frame_f, width = width, exportselection=0)        # self.data_bcg_act.insert(0, "0")
        self.data_bcg_res_act.insert(0, "0")

        label = tk.Label(self.frame_f, text="   residual time [hhmmss]", anchor='w', justify = 'left').grid(row=4, column=2, sticky="nsew")
        self.data_bcg_res_time = tk.Entry(self.frame_f, width = width, exportselection=0)
        self.data_bcg_res_time.insert(0, "hhmmss")
        # self.data_bcg_res_time.insert(0, "100000")

        label = tk.Label(self.frame_f, text="   total volume (ml)", anchor='w', justify = 'left').grid(row=5, column=2, sticky="nsew")
        self.data_bcg_volume_tot = tk.Entry(self.frame_f, width = width, exportselection=0)
        self.data_bcg_volume_tot.insert(0, "9800")

        self.data_bcg_act.grid(row=1, column=3)
        self.data_bcg_time.grid(row=2, column=3)
        self.data_bcg_res_act.grid(row=3, column=3)
        self.data_bcg_res_time.grid(row=4, column=3)
        self.data_bcg_volume_tot.grid(row=5, column=3)

        button = tk.Button(self.frame_f, text='Reset segmentation', command = self.reset_segmentation).grid(row = 6, column=0, columnspan=4)
        label = tk.Label(self.frame_f, text=' ', anchor='w', justify = 'left').grid(row=7, column=0, sticky="nsew")
        button = tk.Button(self.frame_f, text='Export results', command = self.get_out_results).grid(row = 8, column=0, columnspan=4)
        
        # =============================================================================
        # frame_g - empty when initiated - no data opened yet
        # =============================================================================

        # =============================================================================
        # frame_h
        # =============================================================================
        
        self.thresh_set = tk.Label(self.frame_h, width = 10)
        self.thresh_set.grid(row=0, column=0, sticky='nsew')
        self.thresh_set['text'] = '1'
        label = tk.Label(self.frame_h, text=' image threshold', anchor='sw', justify='left').grid(row=0, column=1, sticky='nsew')

        button = tk.Button(self.frame_h, image = self.minus, width="20", height="20", command=lambda: self.threshold_img(-1)).grid(row=0, column=3)
        button = tk.Button(self.frame_h, image = self.plus, width="20", height="20", command=lambda: self.threshold_img(1)).grid(row=0, column=4)

        # =============================================================================
        # master.grid
        # =============================================================================
        self.frame_a.grid(row=0,rowspan=2, column=0, columnspan=2)
        self.frame_b.grid(row=0,rowspan=1, column=2, columnspan=1)
        self.frame_c.grid(row=0,rowspan=2, column=3)
        self.frame_d.grid(row=2, column=0, sticky='nsew')
        self.frame_e.grid(row=2, column=1, sticky='nsew')
        frame_switch(self.frame_e, 0) # inactive until data opened - frame_switch = function to activate/deactivate frames
        self.frame_f.grid(row=2, column=2, sticky='nsew')
        frame_switch(self.frame_f, 0) # inactive until data opened - frame_switch = function to activate/deactivate frames
        self.frame_g.grid(row=2, column=3, sticky='nsew')
        self.frame_h.grid(row=1, column=2, columnspan=1)
        frame_switch(self.frame_h, 0) # inactive until data opened - frame_switch = function to activate/deactivate frames

    def open_data(self):
        # dirname = filedialog.askdirectory(parent=root,initialdir="/",title='Please select a directory') # select the directory with data
        dirname = filedialog.askdirectory(parent=root,initialdir="d:/Python_RC_curves_analysis_tool/c_vision",title='Please select a directory') # select the directory with data
        # what happens if no directory selected - nothing
        if not dirname:
            return

        else:
            PT_lst = []  # create an empty list for PT data

            for item in os.listdir(dirname): # lists the selected directory
                if os.path.isfile(os.path.join(dirname,item)):
                    PT_lst.append(os.path.join(dirname,item))

        temp_data = {} # using dictionary to store data about slices

        # opening the files one after another and storing them into the dictionary 
        # no need for preallocation of the memory and sorting according to the slice position
        # sorting the dictionary is easy
        for i in range(0, len(PT_lst)):
            # check whether it's a DICOM file - if not nothing happens and next file is opened
            try:
                temp = dcm.dcmread(PT_lst[i]) # slice opened
            except:  # if any error in dcm.read accurs, the slice is just skipped
                ...
                      
            modality = str(temp[0x0008, 0x0060].value) # modality information
            number = int(temp[0x0020, 0x0013].value) # slice number information
            position = int(temp[0x0020, 0x0032].value[2]) # recorded in the dictionary based on the slice position

            # reference slice with image DICOM informaiton
            if modality == 'PT' and number == 1: # must be modality PT and first slice in the series
                ref_PT = dcm.dcmread(PT_lst[i], stop_before_pixels=True)
    
            # if calculated like this, with DecayCorrection = START, volume activity is always scaled by the PET itself to SeriesTime
            # SeriesTime is used to correct the injected activity for SUV calculations
            # volume activity at the time of measurement (which can differ from SeriesTime in dynamic studies) can be obtained by
            # dividing temp_data by a corresponding DecayFactor
            temp_data[position] = (temp.pixel_array*temp.RescaleSlope + temp.RescaleIntercept)

        data_info = {} # data_info is a dictionary containing PT series information extracted from the rederrence slice DICOM data

         # data_info as a dictionary
        try:
             data_info = {'rows':ref_PT[0x0028, 0x0010].value, 'cols':ref_PT[0x0028, 0x0011].value, 'slices':len(temp_data), \
                          'slice_thickness':str(ref_PT.SliceThickness), \
                          'pixel_size':str(ref_PT.PixelSpacing[0]), \
                          'study_date': ref_PT.StudyDate[6:8]+'.'+ref_PT.StudyDate[4:6]+'.'+ref_PT.StudyDate[0:4], \
                          'series_date': ref_PT.SeriesDate[6:8]+'.'+ref_PT.SeriesDate[4:6]+'.'+ref_PT.SeriesDate[0:4], \
                          'series_time': ref_PT.SeriesTime[0:2] + ref_PT.SeriesTime[2:4] + ref_PT.SeriesTime[4:6], \
                          'radionuclide': ref_PT[0x0054,0x0016][0][0x0054, 0x0300][0][0x0008,0x0104].value[1::], \
                          'halflife': str(ref_PT[0x0054,0x0016][0][0x0018, 0x1075].value), \
                          'injected_activity': str(ref_PT[0x0054,0x0016][0][0x0018, 0x1074].value/1000000), \
                          'injection_date': str(ref_PT[0x0054,0x0016][0][0x0018, 0x1078].value[6:8])+ \
                                                str(ref_PT[0x0054,0x0016][0][0x0018, 0x1078].value[4:6])+ \
                                                    str(ref_PT[0x0054,0x0016][0][0x0018, 0x1078].value[0:4]), \
                          'injection_time': str(ref_PT[0x0054,0x0016][0][0x0018, 0x1072].value[0:6]), \
                          'device': ref_PT[0x0008,0x1090].value}
         # to be used if [0x0054,0x0016][0][0x0018, 0x1078] is missing
        except:
             data_info = {'rows':ref_PT[0x0028, 0x0010].value, 'cols':ref_PT[0x0028, 0x0011].value, 'slices':len(temp_data), \
                          'slice_thickness':str(ref_PT.SliceThickness), \
                          'pixel_size':str(ref_PT.PixelSpacing[0]), \
                          'study_date': ref_PT.StudyDate[6:8]+'.'+ref_PT.StudyDate[4:6]+'.'+ref_PT.StudyDate[0:4], \
                          'series_date': ref_PT.SeriesDate[6:8]+'.'+ref_PT.SeriesDate[4:6]+'.'+ref_PT.SeriesDate[0:4], \
                          'series_time': ref_PT.SeriesTime[0:2] + ref_PT.SeriesTime[2:4] + ref_PT.SeriesTime[4:6], \
                          'radionuclide': ref_PT[0x0054,0x0016][0][0x0054, 0x0300][0][0x0008,0x0104].value[1::], \
                          'halflife': str(ref_PT[0x0054,0x0016][0][0x0018, 0x1075].value), \
                          'injected_activity': str(ref_PT[0x0054,0x0016][0][0x0018, 0x1074].value/1000000), \
                          'injection_date': str(ref_PT[0x0008, 0x0022].value[6:8])+ \
                                                str(ref_PT[0x0008, 0x0022].value[4:6])+ \
                                                    str(ref_PT[0x0008, 0x0022].value[0:4]), \
                          'injection_time': str(ref_PT[0x0054,0x0016][0][0x0018, 0x1072].value[0:6]), \
                          'device': ref_PT[0x0008,0x1090].value}            

        self.data_size['text'] = 'rows: ' + str(data_info['rows']) + '\n' \
                            'columns: ' + str(data_info['cols']) + '\n' \
                            'slices: ' + str(data_info['slices']) + '\n' \
                            'slices thickness: ' + data_info['slice_thickness'] + ' mm' +'\n' \
                            'pixel size: ' + data_info['pixel_size'] + ' mm' +'\n' \
                            'study date: ' + data_info['study_date'] + '\n' \
                            'series date: ' + data_info['series_date'] + '\n' \
                            'series time: ' + data_info['series_time'] + '\n' \
                            'injection date:' + data_info['injection_date'] + '\n' \
                            'radionuclide: ' + data_info['radionuclide'] + '\n' \
                            'halflife: ' + data_info['halflife'] + ' seconds' + '\n' \
                            'device: ' + data_info['device']

        # dictionary sorting
        temp_items = temp_data.items() # list all dictionary records according to the key for sorting (slice position in this case)
        temp_data = sorted(temp_items)

        # creation of 3D matrix as dictionary is not possible to be visualized in 2D from different views
        PT_temp = np.zeros((data_info['rows'], data_info['cols'], data_info['slices']), dtype=int)
        i = 0
        for x in temp_data:
            PT_temp[:,:,i] = x[1]
            i += 1

        # exception for Philips DICOM heads - additional filtration not present
        if ref_PT.get([0x0018,0x1210]) == None:
            ...
        else:
            data_info['filtration'] = ref_PT[0x0018,0x1210].value
            self.data_size['text'] = self.data_size['text'] + '\n' + \
                'Post-recon filter FWHM: ' + data_info['filtration'] + ' mm'

        # insert data in phantomData object
        self.phantom = phantomData(dirname, PT_temp, ref_PT, [data_info['rows'], data_info['cols'], data_info['slices']], data_info)

        
        # =============================================================================
        # CALLING THE CONVOLUTION CALCULATION
        # =============================================================================
        self.phantom.precompute_peak_map() 
        # =============================================================================
        
        # phantom visualization - self.phantom.show_me_slices((self, of_what, which, where, position, zoom))
        self.phantom.show_me_slices(self.phantom.data, self.phantom.slice, self.phantom_slice_view, self.phantom.position, self.phantom.zoom, self.phantom.thresh)

        # segmentation masks creation - zeroes only
        self.phantom.mask = np.zeros([data_info['rows'], data_info['cols'], data_info['slices']])
        self.phantom.mask_labelled = np.zeros([data_info['rows'], data_info['cols'], data_info['slices']])

        # number of displayed slice setup
        self.phantom_slice_lbl['text']= 'slice: '+ str(self.phantom.slice)

        # initial graph display
        GraphRC(self.frame_c, (self.Xsize, self.Ysize), [0,0], self.phantom.seg_thresh, [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], self.RC[0], self.RC[1], self.RC[2], self.config, self.RC_limit_max, self.RC_limit_A, self.RC_limit_peak)

        # =============================================================================
        # switching additional frames on
        # =============================================================================
        frame_switch(self.frame_e, 1)
        frame_switch(self.frame_f, 1)
        frame_switch(self.frame_h, 1)
        
        # how many roi after reopening a file
        if int(self.roi_nr_lbl['text']) == 6:
            ...
        else:    
            ...
            
        self.segmentation_thresh_lbl['text'] = '50'

        # config.cfg file read
        self.RC_analyzer_config_read()

        # =============================================================================
        # bind mouse events
        # =============================================================================
        self.phantom_slice_view.bind("<Shift-MouseWheel>", self.mouse_zoom)
        self.phantom_slice_view.bind("<MouseWheel>", self.mouse_list)
        self.phantom_slice_view.bind("<ButtonPress-1>", lambda event, on_off=1, roi_list=self.phantom.roi_bg: self.mouse_draw_start(event, on_off, roi_list))
        self.phantom_slice_view.bind("<ButtonPress-2>", self.mouse_click_reset)
        self.phantom_slice_view.bind("<B1-Motion>", lambda event, canvas=self.phantom_slice_view: self.mouse_draw(event, canvas))
        self.phantom_slice_view.bind("<ButtonRelease-1>",lambda event, on_off=0, roi_list=self.phantom.roi_bg: self.mouse_draw_start(event, on_off, roi_list))
        self.mip_UP_view.bind("<Shift-MouseWheel>", self.mask_zoom)
        self.mip_LEFT_view.bind("<Shift-MouseWheel>", self.mask_zoom)

        # =============================================================================
        # number of lesions context menu
        # =============================================================================
        self.rois_popup_menu = tk.Menu(master = self.frame_a, tearoff=0)

        # new roi generation
        
        self.phantom.roi_nr = int(self.roi_nr_lbl['text'])
        
        rois =[str(x) for x in range(1,self.phantom.roi_nr+1)] # roi_nr must be +1, background 1 and 6 rois for IEC Body - labeling of rois

        for i in range(0, len(rois)):
            # if in cycle, lambda i=i must be called to work with actual value
            self.rois_popup_menu.add_command(label=str(rois[i]), command= lambda i=i: self.rois_pos_save(i))

        self.phantom_slice_view.bind("<Button-3>", self.rois)

        # =============================================================================
        # create a list of rois
        # =============================================================================
        self.roi_character_create(self.phantom.roi_nr, self.frame_g, self.phantom.roi_obj_size, \
                                  self.phantom.roi_RC_max_lbl, self.phantom.roi_RC_A_lbl, self.phantom.roi_RC_peak_lbl, \
                                  self.phantom.roi_CRC_max_lbl, self.phantom.roi_CRC_A_lbl, self.phantom.roi_CRC_peak_lbl, \
                                  self.phantom.roi_bg_info, str(self.phantom.seg_thresh))

    # =============================================================================
    #  main windows functions section
    # =============================================================================

    def RC_analyzer_config_read(self):
        with open('RC_config.cfg') as f:
            temp = f.read()

            # divide config file into particular main_block(s) 
            main_block = temp.split("##")

            name = []
            for i in range(1, len(main_block)):
                # divide main_block(s) into rows + delete empty inputs using (filter(None,x))
                config = (list(filter(None, main_block[i].splitlines())))
                self.config.append(config)
                name.append(config[0])

        self.RC_limit_spin_max.configure(values=name[::-1])
        self.RC_limit_spin_A.configure(values=name[::-1])
        self.RC_limit_spin_peak.configure(values=name[::-1])
        self.RC_limit_max.set(name[2])
        self.RC_limit_A.set(name[3])
        self.RC_limit_peak.set(name[4])

    # popup menu behaviour - where it occurs and what it does
    def rois(self, event):
        # check - roi sizes defined?
        test = []
        # position 0 belongs to the background - always
        for i in range(1, self.phantom.roi_nr):
            test.append(self.phantom.roi_obj_size[i].get())

        if ((self.phantom.roi_bg[0] == 0) or (self.data_srcs_act.get() == '0') \
            or (self.data_srcs_time.get() == 'hhmmss') or (self.data_bcg_act.get() == '0') \
                or (self.data_bcg_time.get() == 'hhmmss') or (0 in test)):
            tk.messagebox.showwarning("Warning", "Delineate background ROI and/or fill Experiment information and/or specify ROI diameters!")
        else:
            self.all_setup_fill(self.phantom.all_setup) # add experiment setup values
            self.phantom.position_click = self.phantom.slice_click(self, (event.y, event.x)) # check the position of mouse click row x col!
            try:
                self.rois_popup_menu.tk_popup(event.x_root+15, event.y_root, 0) # show the menu
            finally:
                self.rois_popup_menu.grab_release()

    def rois_pos_save(self, i):
        self.phantom.roi_position[i] = self.phantom.position_click
        # pre-segmentation to find the maximum value
        new_seed = self.phantom.segment_me_with_A(self.phantom.roi_position.get(i), self.phantom.get_bg_for_segmentation(), float(self.phantom.seg_thresh)/100, i+1) # proved pre-segmentaci k nalezeni pozice maxima
        # main segmentation using the result from a previous step
        self.phantom.segment_me_with_A(new_seed, self.phantom.get_bg_for_segmentation(), float(self.phantom.seg_thresh)/100, i+1) # proved finalni segmentaci
        # display MIPs (self, of_what, which, where, position, zoom)
        self.phantom.show_me_mip(self.phantom.mask_labelled, [self.mip_UP_view, self.mip_LEFT_view], self.phantom.mip_position, self.phantom.mask_zoom)
        # get the result
        self.phantom.get_results(self.phantom.data, self.phantom.mask_labelled, i+1)
        # display the result
        self.phantom.show_me_results(self.frame_c)

    # changing the number of rois
    def phantom_rois_nr(self, i, where):
        self.phantom.roi_nr += i      

        # =============================================================================
        # create a list of rois
        # =============================================================================
        if self.phantom.roi_nr <= 1: self.phantom.roi_nr = 1
        elif self.phantom.roi_nr >= 12: self.phantom.roi_nr = 12

        self.roi_nr_lbl['text'] = str(self.phantom.roi_nr)

        for roi_character_lbl in where.winfo_children():
                roi_character_lbl.destroy() # detele all segROI checkbuttons
                rois = []

        # new roi list generation
        self.roi_character_create(self.phantom.roi_nr, where, self.phantom.roi_obj_size, \
                                      self.phantom.roi_RC_max_lbl, self.phantom.roi_RC_A_lbl, self.phantom.roi_RC_peak_lbl, \
                                  self.phantom.roi_CRC_max_lbl, self.phantom.roi_CRC_A_lbl, self.phantom.roi_CRC_peak_lbl, \
                                  self.phantom.roi_bg_info, str(self.phantom.seg_thresh))

        # new roi popup menu generation
        self.rois_popup_menu = tk.Menu(master = self.frame_a, tearoff=0)

        # new roi generation
        rois =[str(x) for x in range(1,self.phantom.roi_nr+1)] # roi_nr must be +1, to have 1 backgroung and 6 rois for IEC Body - labelling - label nr used in subsequent functions

        for i in range(0, len(rois)):
            # if in cycle, lambda i=i must be called to work with actual value
            self.rois_popup_menu.add_command(label=str(rois[i]), command= lambda i=i: self.rois_pos_save(i))
            
    # changing the image threshold
    def threshold_img(self, i):
        self.phantom.thresh += i 
        
        if self.phantom.thresh <= 1: self.thresh = 1
        elif self.phantom.thresh >= 10: self.phantom.thresh = 10
        
        self.thresh_set['text'] = str(self.phantom.thresh)
        
        # phantom visualization - self.phantom.show_me_slices((self, of_what, which, where, position, zoom))
        self.phantom.show_me_slices(self.phantom.data, self.phantom.slice, self.phantom_slice_view, self.phantom.position, self.phantom.zoom, self.phantom.thresh)
        
    # changing the segmentation threshold
    def seg_thresh(self, i, where):
        self.phantom.seg_thresh += i

        if self.phantom.seg_thresh <= 1: self.phantom.seg_thresh = 1
        elif self.phantom.seg_thresh >= 100: self.phantom.seg_thresh = 100

        self.segmentation_thresh_lbl['text'] = str(self.phantom.seg_thresh)

        # zeroing of already existing masks
        self.phantom.mask = np.zeros([self.phantom.size[0], self.phantom.size[1], self.phantom.size[2]])
        self.phantom.mask_labelled = np.zeros([self.phantom.size[0], self.phantom.size[1], self.phantom.size[2]])

        # repeated generation of rois list
        self.roi_character_create(self.phantom.roi_nr, where, self.phantom.roi_obj_size, \
                                      self.phantom.roi_RC_max_lbl, self.phantom.roi_RC_A_lbl, self.phantom.roi_RC_peak_lbl, \
                                  self.phantom.roi_CRC_max_lbl, self.phantom.roi_CRC_A_lbl, self.phantom.roi_CRC_peak_lbl, \
                                  self.phantom.roi_bg_info, str(self.phantom.seg_thresh))

        GraphRC(self.frame_c, (self.Xsize, self.Ysize), [0,0], self.phantom.seg_thresh, [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], self.RC[0], self.RC[1], self.RC[2], self.config, self.RC_limit_max, self.RC_limit_A, self.RC_limit_peak)

    def all_setup_fill(self, dict_to_fill):
        dict_to_fill['sources_activity'] = self.data_srcs_act.get().replace(',','.')
        dict_to_fill['sources_residual_activity'] = self.data_srcs_res_act.get().replace(',','.')
        dict_to_fill['background_activity'] = self.data_bcg_act.get().replace(',','.')
        dict_to_fill['background_residual_activity'] = self.data_bcg_res_act.get().replace(',','.')
        dict_to_fill['sources_time'] = self.data_srcs_time.get()
        dict_to_fill['sources_residual_time'] = self.data_srcs_res_time.get()
        dict_to_fill['background_time'] = self.data_bcg_time.get()
        dict_to_fill['background_residual_time'] = self.data_bcg_res_time.get()
        dict_to_fill['sources_activity_vol'] = self.data_srcs_volume_tot.get()
        dict_to_fill['background_activity_vol'] = self.data_bcg_volume_tot.get()

    def mouse_zoom(self, event, *kwargs):
        zoom = [event.delta, event.x, event.y]
        self.phantom.slice_zoom(self, self.phantom.data, zoom, self.phantom_slice_view)
        self.bg_roi_redraw(self.phantom_slice_view)

    def mouse_list(self, event, *kwargs):
        index = event.delta
        self.phantom.slice_list(self, self.phantom.data, index, self.phantom_slice_view)
        self.phantom_slice_lbl['text']= 'slice: '+ str(self.phantom.slice)
        self.bg_roi_redraw(self.phantom_slice_view)

    def mask_zoom(self, event):
        zoom = [event.delta, event.x, event.y]
        self.phantom.mip_zoom(self, self.phantom.mask, zoom, [self.mip_UP_view, self.mip_LEFT_view])

    def roi_character_create(self, how_many, where, roi_size, RC_max_lbl, RC_A_lbl, RC_peak_lbl, CRC_max_lbl, CRC_A_lbl, CRC_peak_lbl, bcg_info, thresh):
        # new rois generation
        end = how_many+1
        rois = ['empty'] + [str(x) for x in range(1,end)]
        width = 8
        # roi list generation
        for r in range(0,end):

            if r == 0: # this is only for the background
                label = tk.Label(where, text = ' ')
                label.grid(row=r, column=0)

            else:
                label = tk.Label(where, text = 'roi ' + rois[r])
                label.grid(row=r, column=0)
                roi_size[r] = tk.IntVar()
                roi_character_size = tk.Entry(where, textvariable = roi_size[r], width = width, exportselection = 0)
                roi_character_size.grid(row=r, column=1)

                label = tk.Label(where, text = ' mm')
                label.grid(row=r, column=2)

                label = tk.Label(where, text = 'RC_max = ')
                label.grid(row=r, column=3)
                RC_max_lbl[r] = tk.StringVar()
                roi_RC_max = tk.Label(where, width = width,textvariable = RC_max_lbl[r])
                roi_RC_max.grid(row=r, column=4)

                label = tk.Label(where, text = 'RC_A' + thresh + ' = ')
                label.grid(row=r, column=5)
                RC_A_lbl[r] = tk.StringVar()
                roi_RC_A = tk.Label(where, width = width,textvariable = RC_A_lbl[r], anchor='w')
                roi_RC_A.grid(row=r, column=6)

                label = tk.Label(where, text = 'RC_peak = ')
                label.grid(row=r, column=7)
                RC_peak_lbl[r] = tk.StringVar()
                roi_RC_peak = tk.Label(where, width = width,textvariable = RC_peak_lbl[r], anchor='w')
                roi_RC_peak.grid(row=r, column=8)

                label = tk.Label(where, text = 'CRC_max = ')
                label.grid(row=r, column=9)
                CRC_max_lbl[r] = tk.StringVar()
                roi_CRC_max = tk.Label(where, width = width,textvariable = CRC_max_lbl[r], anchor='w')
                roi_CRC_max.grid(row=r, column=10)

                label = tk.Label(where, text = 'CRC_A' + thresh + ' = ')
                label.grid(row=r, column=11)
                CRC_A_lbl[r] = tk.StringVar()
                roi_CRC_A = tk.Label(where, width = width,textvariable = CRC_A_lbl[r], anchor='w')
                roi_CRC_A.grid(row=r, column=12)

                label = tk.Label(where, text = 'CRC_peak = ')
                label.grid(row=r, column=13)
                CRC_peak_lbl[r] = tk.StringVar()
                roi_CRC_peak = tk.Label(where, width = width,textvariable = CRC_peak_lbl[r], anchor='w')
                roi_CRC_peak.grid(row=r, column=14)



        background_s = ['background real volume activity: '] + ['background measured volume activity: '] + ['background real / measured difference: '] + ['backgroung COV: ']
        background_e = [' Bq/ml'] + [' Bq/ml'] + [' %'] + [' %']

        start = how_many+2
        for r in range(0,4):

            label = tk.Label(where, text = background_s[r])
            label.grid(row=start+r, column=0, columnspan=4, sticky = 'w')

            bcg_info[r] = tk.StringVar()
            roi_RC_max = tk.Label(where, width = width, textvariable = bcg_info[r])
            roi_RC_max.grid(row=start+r, column=4)
            label = tk.Label(where, text = background_e[r])
            label.grid(row=start+r, column=5)

    # automatically fill NEMA IEC Body phantom data -  what = fantom obj
    def roi_character_autofill(self, what, spec):
        if spec == 6:
            # IEC Body Phantom parameters
            IEC = ['37', '28', '22', '17', '13', '10']
            if what.roi_nr == 6:
                for i in range(0,6):
                    what.roi_obj_size[i+1].set(IEC[i])
            else:
                ...
        elif spec == 12:
                # IEC_12 Body Phantom parameters
            IEC = ['37', '28', '22', '17', '13', '10', '15.4', '12.4', '7.9', '6.2', '5', '4']
            if what.roi_nr == 12:
                for i in range(0,12):
                    what.roi_obj_size[i+1].set(IEC[i])
            else:
                ...

    # reset button
    def mouse_click_reset(self, event):
        # phantom display -  self.phantom.show_me_slices(which, where, zoom)
        self.phantom.show_me_slices(self.phantom.data, self.phantom.slice, self.phantom_slice_view, [0,0], 1, self.phantom.thresh)
        self.phantom.position=[0,0]
        self.phantom.zoom = 1
        self.bg_roi_redraw(self.phantom_slice_view)

    # get mouse click position to start drawing the background 
    def mouse_draw_start(self, event, on_off, roi_list):

        # check whether there is already a background roi existing
        if self.phantom.roi_bg[0] == 0: # it is not
            ...
        else:
            self.phantom_slice_view.delete(self.phantom.roi_bg[0]) # existing background roi is deleted and new one will be created

        if on_off == 1: # start drawing the roi
            # check whether all data are prefilled - after the roi is created, the automatic evaluation of its data starts
            # but data about the experiment must be available
            if (self.data_srcs_act.get() == '0') or (self.data_srcs_time.get() == 'hhmmss') or \
                (self.data_bcg_act.get() == '0') or (self.data_bcg_time.get() == 'hhmmss'):
                    tk.messagebox.showwarning("Warning", "Fill Experiment information!")
            else:
                self.x_click = event.x
                self.y_click = event.y
                # record the position into self.phantom.roi_position - position in real image coordinates system
                self.phantom.roi_position['bg'] = self.phantom.slice_click(self, (self.y_click, self.x_click))
                self.drawn = None

        else: # drawing finish
            if (self.data_srcs_act.get() == '0') or (self.data_srcs_time.get() == 'hhmmss') or \
                (self.data_bcg_act.get() == '0') or (self.data_bcg_time.get() == 'hhmmss'):
                    ...
            else:
                roi_list[0] = self.drawn
                self.x_click = event.x
                self.y_click = event.y
                # record the position into self.phantom.roi_position
                self.phantom.roi_position['bg'] = self.phantom.roi_position.get('bg') , self.phantom.slice_click(self, (self.y_click, self.x_click))
                self.all_setup_fill(self.phantom.all_setup) # add experiment values data
                # signal roi bg created - only after mouse button released
                if self.phantom.roi_position['bg'][0][0] == self.phantom.roi_position['bg'][1][0]: # single point background roi check - unintentional single click
                    ...
                else:
                    self.phantom.get_bg_for_segmentation()

    # background roi drawn
    def mouse_draw(self, event, canvas):
        if (self.data_srcs_act.get() == '0') or (self.data_srcs_time.get() == 'hhmmss') or \
            (self.data_bcg_act.get() == '0') or (self.data_bcg_time.get() == 'hhmmss'):
                ...
        else: # delete all so far exisiting 'roi'
            canvas.delete('roi')
            canvas.create_oval(self.x_click, self.y_click, event.x, event.y, outline = 'red', width = 2, tags = 'roi')

    def bg_roi_redraw(self, canvas):
        # # delete all so far exisiting 'roi'
        canvas.delete('roi')
        if self.phantom.roi_bg[0] == 0:
            ...
        else:
            if abs(self.phantom.slice - self.phantom.roi_position.get('bg')[0][2]) <=2: # background roi drawn on central and +-2 neighbouring slices
                # self.phantom.roi_position.get('bg') # this coordinates muse be converted into the actual image coordinates
                start = self.phantom.slice_click_reversed(self, (self.phantom.roi_position.get('bg')[0][0], self.phantom.roi_position.get('bg')[0][1]))
                end = self.phantom.slice_click_reversed(self, (self.phantom.roi_position.get('bg')[1][0], self.phantom.roi_position.get('bg')[1][1]))
                # re-creating of the oval object in correct coordinates - event.x vs event.y vs rows vs cols 
                self.phantom_slice_view.create_oval(start[1],start[0],end[1],end[0], outline = 'red', width = 2, tags='roi')
            else:
                ...

    def RCs_redraw(self):
        self.phantom.show_me_results(self.frame_c)

    def get_out_results(self):
        self.phantom.get_out_results(self.phantom.directory)

    def reset_segmentation(self):
        self.phantom.mask = np.zeros([self.phantom.size[0], self.phantom.size[1], self.phantom.size[2]])
        self.phantom.mask_labelled = np.zeros([self.phantom.size[0], self.phantom.size[1], self.phantom.size[2]])
        GraphRC(self.frame_c, (self.Xsize, self.Ysize), [0,0], self.phantom.seg_thresh, [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], self.RC[0], self.RC[1], self.RC[2], self.config, self.RC_limit_max, self.RC_limit_A, self.RC_limit_peak)
        self.phantom.show_me_mip(self.phantom.mask_labelled, [self.mip_UP_view, self.mip_LEFT_view], self.phantom.mip_position, self.phantom.mask_zoom)

#%% global functions
# function for switching frames on and off
def frame_switch(frame_name, on_off):
    if on_off == 0:
        for child in frame_name.winfo_children():
            child.configure(state='disabled')
    else:
        for child in frame_name.winfo_children():
            child.configure(state='normal')

root = tk.Tk()
RC_analyzer = mainWindow(root)
root.mainloop()

