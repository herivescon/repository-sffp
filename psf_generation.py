"""""
 @Purpose: Generate PSF and photometric spectrum from images
 @Input: Directory where stacked images are created
 @Output: PSF and spectrum of hyperspectral
 @Written by: Hector Erives
 @Date: Summer 2022
"""""

import os
import numpy as np
from astropy.io import fits
from astropy.modeling import models, fitting
import PySimpleGUI as sg
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
import pandas as pd
import math


header_dict = {'ImgType': 'IMAGETYP ', 'Filter': 'FILTER', 'Exposure': 'EXPOSURE'}
frame_dict = {'Dark': 'Dark Frame', 'Light': 'Light Frame', 'Flat': 'Flat Field', 'Bias': 'Bias Frame'}
filter_dict = {'B': 'B', 'V': 'V', 'R': 'R', 'H': '100 lines/mm'}
dir_separ = '/'


def main():

    fwhm = 8
    # Information to be used to generate PSFs for different filters
    #       id xcentroid ycentroid sharpness roundness1 roundness2 npix sky peak       flux   mag
    #       --- --------- --------- --------- ---------- ---------- ---- --- ----      ----- -------
    # (B)  2     703       165.4    0.3478     0.3486     0.1077   121   0 2731         14.94 -2.936
    # (V)  2     702.9     165.4    0.4064     0.3731     0.1445   121   0 1.04e+04     54.23 -4.336
    # (R)  2     701.6     165.9    0.3909     0.3938    -0.01214  121   0 2.173e+04    125   -5.242
    xcent_B = 703
    ycent_B = 166
    npix_B = 121
    xcent_V = 703
    ycent_V = 166
    npix_V = 121
    xcent_R = 702
    ycent_R = 166
    npix_R = 121

    # read stacked images
    folder_path = sg.popup_get_folder("Please enter the parent folder where stacked images are.")
    if folder_path is None or folder_path == '':
        sg.PopupCancel('Cancelled - No valid folder entered')
        return
    for folder_item in os.listdir(folder_path):
        if folder_item == 'Photometric':
            stacked_dir = folder_path + dir_separ + folder_item + dir_separ + 'Stacked'
            for file_item in os.listdir(stacked_dir):
                file = stacked_dir + dir_separ + file_item
                with fits.open(file) as img:
                    if img[0].header[header_dict['Filter']] == 'B':
                        psf_B = generate_psf("Band B", file, xcent_B, ycent_B, npix_B)
                    if img[0].header[header_dict['Filter']] == 'V':
                        psf_V = generate_psf("Band V", file, xcent_V, ycent_V, npix_V)
                    if img[0].header[header_dict['Filter']] == 'R':
                        psf_R = generate_psf("Band R", file, xcent_R, ycent_R, npix_R)
        if folder_item == 'Spectral':
            stacked_dir = folder_path + dir_separ + folder_item + dir_separ + 'Stacked'
            for file_item in os.listdir(stacked_dir):
                file = stacked_dir + dir_separ + file_item
                img_H = fits.getdata(file).byteswap().newbyteorder()
                wv, spectrum = extract_spectrum(img_H)

    # define the directory where data will be saved
    filename = folder_path + dir_separ + 'BVR_psf.csv'
    data = np.column_stack((psf_B/np.sum(psf_B), psf_V/np.sum(psf_V), psf_R/np.sum(psf_R)))
    header = []
    dx, dy = psf_B.shape
    for i in range(dx):
        header.append('B'+str(i))
    for i in range(dx):
        header.append('V'+str(i))
    for i in range(dx):
        header.append('R'+str(i))
    data = pd.DataFrame(data)
    data.to_csv(filename, index=False, header=False)

    filename = folder_path + dir_separ + 'measured_diffraction.csv'
    header = ['Wavelength', 'Transmission']
    data = np.column_stack((wv, spectrum))
    data = pd.DataFrame(data, columns=header)
    data.to_csv(filename, index=False)


def generate_psf(
    band,
    file,
    xcent, ycent,
    npix
):
    img = fits.getdata(file).byteswap().newbyteorder()
    #crop_x = np.arange(xcent - math.ceil((npix ** 0.5)/2), xcent + math.ceil((npix ** 0.5)/2)-1, 1)
    #crop_y = np.arange(ycent - math.ceil((npix ** 0.5)/2), ycent + math.ceil((npix ** 0.5)/2)-1, 1)
    crop_x = np.arange(xcent - math.ceil((npix ** 0.5)), xcent + math.ceil((npix ** 0.5))-1, 1)
    crop_y = np.arange(ycent - math.ceil((npix ** 0.5)), ycent + math.ceil((npix ** 0.5))-1, 1)
    box = img[np.ix_(crop_y, crop_x)]; box = box - np.min(box); box = box/np.max(box)
    yp, xp = box.shape

    y, x, = np.mgrid[:yp, :xp]
    amp = 1
    x0 = xp/2; y0 = yp/2
    sigx = xp/4; sigy = yp/4
    f_init = models.Gaussian2D(amp, x0, y0, sigx, sigy)
    fit_f = fitting.LevMarLSQFitter()
    fun = fit_f(f_init, x, y, box)
    psf = fun(x, y)

    # plot point source, psf, and fitting error
    plt.figure()
    plt.subplot(1, 3, 1)
    plt.imshow(box)
    plt.title(band)
    plt.subplot(1, 3, 2)
    plt.imshow(psf)
    plt.title("Best Gaussian fit of band "+band)
    plt.subplot(1, 3, 3)
    plt.imshow(box - psf)
    plt.title("Residual Error")
    plt.show()
    plt.show(block=True)

    return psf


def extract_spectrum(img):
    # define wavelengths
    m1 = (760-400)/(629-400)
    b1 = 400 - m1 * 400
    pix = np.linspace(0, 1024, num=1024)
    wv = pix * m1 + b1

    # compute slope and offset to get diffraction grating from image
    x1 = 700; y1 = 1024 - 165
    x2 = 690; y2 = 1024 - 560
    m = (y2 - y1) / (x2 - x1); b = y1 - m * x1

    # extract counts from image
    spectrum = np.empty(1024, dtype=object)
    for y in range(1024):
        xx = round((y - b)/m)
        yy = 1023 - y
        spectrum[1023-y] = img[yy, xx]
    spectrum = savgol_filter(spectrum, 35, 5)
    wv = wv[350:800]
    spectrum = spectrum[350:800]
    wv_sel = np.arange(400, 1000, 0.1)
    interp = interp1d(wv, spectrum, kind="linear")
    spectrum_sel = interp(wv_sel)
    plt.figure()
    plt.plot(wv_sel, spectrum_sel, color='red')
    plt.title('Average spectrum')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Counts')
    plt.show()
    plt.show(block=True)

    return wv_sel, spectrum_sel


# Define main method
if __name__ == '__main__':
    main()
    