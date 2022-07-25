"""""
 @Purpose: Stack images to improve SNR on RAW images taken by Falcon 16
 @Input: Directory where RAW images are
 @Output: Stacked images/filter (5 images)
 @Written by: Hector Erives
 @Date: Summer 2022
"""""

import os
import numpy as np
from astropy.io import fits
import PySimpleGUI as sg
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder
from photutils.aperture import aperture_photometry, CircularAperture, CircularAnnulus


def main():

    # Define header information.
    frame_dict = {
        'Light': 'Light Frame',
        'Dark': 'Dark Frame',
        'Bias': 'Bias Frame',
        'Flat': 'Flat Field'
    }
    filter_dict = {
        'B': 'B',
        'V': 'V',
        'R': 'R',
        'H': '100 lines/mm'
    }
    header_dict = {
        'ImgType': 'IMAGETYP',
        'Filter': 'FILTER',
        'ExpTime': 'EXPOSURE',
        'ObjName': 'OBJECT',
        'NumXpx': 'NAXIS1',
        'NumYpx': 'NAXIS2',
        'SetTemp': 'SET - TEMP',
        'ActTemp': 'CCD - TEMP',
        'Hist': 'HISTORY'
    }
    ext = '.fit'
    fwhm = 8
    r_aper = 1.50
    r_inn = 2.50
    r_out = (r_inn**2 + r_aper**2)**0.5
    snr = 10
    filt = 'B'
    filt = 'V'
    filt = 'R'
    #filt = 'H'  # Use SNR=10 for 'H' filter

    # read images and stack them (only 5 because satellite drift)
    folder_path = sg.popup_get_folder("Please enter the parent folder for reduction.")
    stacked_image = stack_images(
        folder_path,
        frame_dict['Light'],
        ext,
        filter_dict[filt],  # Define the filter to analyze
        header_dict,
        filt
    )

    # plot stacked image
    plt.figure(1)
    plt.imshow(stacked_image, cmap='gray_r', norm=LogNorm(), origin='lower')
    plt.colorbar()

    # extract sources and plot them
    source_extraction(
        stacked_image,
        fwhm,
        r_aper,
        r_inn,
        r_out,
        snr
    )


def stack_images(
    directory,
    exp_type,
    ext,
    filter_type,
    header_dict,
    filt
):
    head_ExpTime = 'EXPOSURE'
    dir_separ = "/"
    files = []
    total_exp_time = 0
    header_array = []

    # scan directory and record file names
    for item in os.scandir(directory):
        if ext in os.path.split(item)[1] and item.is_file():
            with fits.open(item) as img:
                if img[0].header[header_dict['ImgType']] == exp_type and \
                        img[0].header[header_dict['Filter']] == filter_type:
                    files.append(item)

    # stack only 5 images per filter
    indx = np.array([149, 150, 151, 152, 153])
    img_data = fits.getdata(files[indx[0]]).byteswap().newbyteorder()
    stacked_image = np.zeros(np.shape(img_data))
    for idx in range(len(indx)):
        img_data = fits.getdata(files[indx[idx]]).byteswap().newbyteorder()
        stacked_image += img_data
        with fits.open(files[indx[idx]]) as hdulist:
            exp_time = hdulist[0].header[head_ExpTime]
            total_exp_time += exp_time

    # generate header file for the stacked image
    header_array.append(hdulist[0].header)
    temp_header = header_array[0]
    temp_header['T'] = (str(total_exp_time), "seconds")
    temp_header['NSI'] = (str(len(files)), "images")
    header_array[0] = temp_header

    # generate name of output file
    stack_exp_time = str(total_exp_time)
    Indv_exp_time = str(exp_time)
    output_name = filt + Indv_exp_time + "-TotalEXPtime-" + stack_exp_time + ext

    # Check if stacked folder exists , if not make one
    sub_path = directory + dir_separ + "Stacked"
    if not os.path.exists(sub_path):
        os.makedirs(sub_path)

    # Write file to path
    file_path = sub_path + dir_separ + output_name
    if not os.path.exists(file_path):
        fits.writeto(file_path, stacked_image, header_array[0])
    print(file_path)

    # return stacked file name
    return stacked_image


def source_extraction(
    stacked_image,
    fwhm,
    r_aper,
    r_inn,
    r_out,
    snr
):
    bg_mean, bg_median, bg_std = sigma_clipped_stats(
        stacked_image,
        sigma=3.0,
        maxiters=5,
        std_ddof=0,
    )
    print("\nMean: {0:.4f}, Median: {1:.4f}, STD: {2:.4f}".format(bg_mean, bg_median, bg_std))
    dao_find = DAOStarFinder(
        fwhm=fwhm,
        threshold=snr*bg_std,
        peakmax=None
    )

    # Print and plot out source
    sources = dao_find(stacked_image - bg_median)
    print("\nDaofind Source Table ")
    for col in sources.colnames:
        sources[col].info.format = '%.4g'
    sources.pprint_all()

    # Plot circles on the identified sources
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    apertures = CircularAperture(positions, r=r_aper*fwhm)
    annulus_apertures = CircularAnnulus(positions, r_in=r_inn*fwhm, r_out=r_out*fwhm)
    apertures.plot(color='green', lw=1.5, alpha=0.75)
    annulus_apertures.plot(color='red', lw=0.8, alpha=0.38)
    plt.xlabel("X Axis")
    plt.ylabel("Y Axis")
    plt.show(block=True)

    #       id xcentroid ycentroid sharpness roundness1 roundness2 npix sky peak       flux   mag
    #       --- --------- --------- --------- ---------- ---------- ---- --- ----      ----- -------
    # (B)  2       703     165.4    0.3478     0.3486     0.1077  121   0 2731 14.94  -2.936
    # (V)  2     702.9     165.4    0.4064     0.3731     0.1445  121   0 1.04e+04 54.23 -4.336
    # (R)  2     701.6     165.9    0.3909     0.3938   -0.01214  121   0 2.173e+04   125 -5.242


# Define main method
if __name__ == '__main__':
    main()
