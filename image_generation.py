"""""
 @Purpose: Run DIRSIG
 @Input: none
 @Output: none
 @Written by: Hector Erives
 @Date: Summer 2022
"""""

import os
import numpy as np
import pandas as pd
from spectral import *
import matplotlib.pyplot as plt
from astropy.convolution import convolve
import math
from array import *


def main():
    # Define constants as described by ref. in Journal of the Astronautial Sciences
    # X: airmass, k: extinction coefficient, m_zp: zero point
    X = 1.4
    B_k = -0.25; B_mzp = 21.47
    V_k = -0.19; V_mzp = 21.62
    R_k = -0.12; R_mzp = 21.63

    # read PSF for conditions when images were acquired
    psf_file = 'BVR_psf.csv'
    PSF = pd.read_csv(psf_file)
    #PSF_B = PSF.iloc[:, 0:31]; PSF_V = PSF.iloc[:, 32:53]; PSF_R = PSF.iloc[:, 54:75]
    # PSF_B = PSF.iloc[:, 0:11]; PSF_V = PSF.iloc[:, 11:22]; PSF_R = PSF.iloc[:, 22:33]
    PSF_B = PSF.iloc[:, 0:21]; PSF_V = PSF.iloc[:, 21:42]; PSF_R = PSF.iloc[:, 42:63]

    # define file to open and initialize values
    DIRSIG_filename = "photo_watts-t0000-c0000.img.hdr"
    hours = list(range(5, 9, 1)); minutes = list(range(0, 60, 10))
    #         5:00 5:10 5:20 5:30 5:40 5:50 6.00 6:10 6.20 6:30 6:40 6.50 7:00 7.10 7:20 7:30 7:40 7:50 8:00 8.10 8:20 8:30 8:40 8:50
    #angles = list([[0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0],
    #            [-45, -45, -45, -45, -45, -45, -45, -45, -45, -30, -20, -10,   0,  10,  20,  30,  45,  45,  45,  45,  45,  45,  45,  45]])
    #hours = list(range(5, 9, 1)); minutes = list(range(0, 60, 30))
    #               5:00  6.00    7:00   8:00
    angles = list([[   0,    0,      0,     0],
                   [ -45,  -45,      0,    45]])

    times = np.zeros(len(hours)*len(minutes))
    cal_val_B = np.zeros([len(hours)*len(minutes), 2])
    cal_val_V = np.zeros([len(hours)*len(minutes), 2])
    cal_val_R = np.zeros([len(hours)*len(minutes), 2])

    path_DIRSIG = "C:/Users/herivescon/Downloads/DIRSIG/DIRSIG_SSA1/Ssa1/"; path_geometry = path_DIRSIG + "geometry"
    os.chdir(path_geometry)
    # input_file = " --input_filename=directv_10_0.obj";      output_file = " --output_filename=directv_10.obj"
    # input_panels = " --input_filename=simple_panels_scaled.obj"; output_panels = " --output_filename=simple_panels.obj"
    # input_bus = " --input_filename=simple_bus_scaled.obj";       output_bus = " --output_filename=simple_bus.obj"
    input_panels = " --input_filename=panels_DIRSIG_0.obj"
    output_panels = " --output_filename=panels_DIRSIG.obj"
    input_bus = " --input_filename=bus_DIRSIG_0.obj"
    output_bus = " --output_filename=bus_DIRSIG.obj"
    rotatex = " --rotatex=-5"; scale = " --scale=0.3"

    for idx in range(2):
        t = 0
        for h, hour in enumerate(hours):
            for m, minute in enumerate(minutes):
                # rotate object
                os.chdir(path_geometry)
                rotatey = " --rotatey=" + str(angles[idx][t])
                # os.system("object_tool" + input_file + output_file + rotatey)
                os.system("object_tool" + input_bus + output_bus + rotatex + scale)
                os.system("object_tool" + input_panels + output_panels + rotatex + rotatey + scale)

                # define time of acquisition, and angle and run DIRSIG
                generate_tasks(path_DIRSIG, hour, minute)

                # compute instrument magnitudes
                Photo_B, Photo_V, Photo_R, Photo_B_psf, Photo_V_psf, Photo_R_psf = convolve_images(DIRSIG_filename, PSF_B, PSF_V, PSF_R)
                x, y = np.shape(Photo_B); x0 = y0 = math.floor(x / 2); roix = roiy = math.floor(x0 / 4)
                B_minst, V_minst, R_minst = compute_maginst(x0, y0, roix, roiy, Photo_B_psf, Photo_V_psf, Photo_R_psf)

                # compute calibrated magnitude
                B_mcal = B_k * X + B_mzp + B_minst
                V_mcal = V_k * X + V_mzp + V_minst
                R_mcal = R_k * X + R_mzp + R_minst

                # save data to arrays
                times[t] = hour + minute / 60
                cal_val_B[t, idx] = B_mcal[roix, roiy]
                cal_val_V[t, idx] = V_mcal[roix, roiy]
                cal_val_R[t, idx] = R_mcal[roix, roiy]
                t += 1
                plot_images(Photo_B, Photo_V, Photo_R, Photo_B_psf, Photo_V_psf, Photo_R_psf, B_mcal, V_mcal, R_mcal)

    plot_calvals(times, cal_val_B, cal_val_V, cal_val_R)


def generate_tasks(path_DIRSIG, hour, minute):
    os.chdir(path_DIRSIG)
    file = open("demo.tasks", "w")
    with open('demo.tasks.template') as f:
        for i, line in enumerate(f):
            if i == 2:
                time = list(line); time[42] = str(hour)
                if minute < 10:
                    time[45] = str(minute)
                else:
                    parts = list(str(minute)); time[44] = parts[0]; time[45] = parts[1]
                line = ''.join(time)
            file.write(line)
    file.close()
    os.system("dirsig4 demo.sim")


def convolve_images(DIRSIG_filename, PSF_B, PSF_V, PSF_R):
    # read DIRSIG file name
    hdr_file = DIRSIG_filename
    Photo_B = open_image(hdr_file).read_band(0)
    Photo_B_psf = convolve(Photo_B, PSF_B)
    Photo_V = open_image(hdr_file).read_band(1)
    Photo_V_psf = convolve(Photo_V, PSF_V)
    Photo_R = open_image(hdr_file).read_band(2)
    Photo_R_psf = convolve(Photo_R, PSF_R)
    return Photo_B, Photo_V, Photo_R, Photo_B_psf, Photo_V_psf, Photo_R_psf


def compute_maginst(x0, y0, roix, roiy, Photo_B_psf, Photo_V_psf, Photo_R_psf):
    # subselect psf images
    Photo_B_psf_sel = Photo_B_psf[x0 - roix: x0 + roix, y0 - roiy: y0 + roiy]
    Photo_V_psf_sel = Photo_V_psf[x0 - roix: x0 + roix, y0 - roiy: y0 + roiy]
    Photo_R_psf_sel = Photo_R_psf[x0 - roix: x0 + roix, y0 - roiy: y0 + roiy]

    # compute to magnitude
    B_minst = -2.5 * np.log10(Photo_B_psf_sel)
    V_minst = -2.5 * np.log10(Photo_V_psf_sel)
    R_minst = -2.5 * np.log10(Photo_R_psf_sel)
    return B_minst, V_minst, R_minst


# plot results
def plot_images(Photo_B, Photo_V, Photo_R, Photo_B_psf, Photo_V_psf, Photo_R_psf, B_mcal, V_mcal, R_mcal):
    plt.figure()
    plt.subplot(3, 4, 1)
    plt.imshow(Photo_B); plt.colorbar(); plt.title('Band B ($W/m^{2}$)')
    plt.subplot(3, 4, 2)
    plt.imshow(Photo_V); plt.colorbar(); plt.title('Band V ($W/m^{2}$)')
    plt.subplot(3, 4, 3)
    plt.imshow(Photo_R); plt.colorbar(); plt.title('Band R ($W/m^{2}$)')
    plt.subplot(3, 4, 5)
    plt.imshow(Photo_B_psf); plt.colorbar(); plt.title('Band B ($W/m^{2}$)')
    plt.subplot(3, 4, 6)
    plt.imshow(Photo_V_psf); plt.colorbar(); plt.title('Band V ($W/m^{2}$)')
    plt.subplot(3, 4, 7)
    plt.imshow(Photo_R_psf); plt.colorbar(); plt.title('Band R ($W/m^{2}$)')
    plt.subplot(3, 4, 9)
    plt.imshow(B_mcal); plt.colorbar(); plt.title('Cal Mag B')
    plt.subplot(3, 4, 10)
    plt.imshow(V_mcal); plt.colorbar(); plt.title('Cal Mag V')
    plt.subplot(3, 4, 11)
    plt.imshow(R_mcal); plt.colorbar(); plt.title('Cal Mag R')
    plt.show()
    plt.show(block=True)


def plot_calvals(times, cal_val_B, cal_val_V, cal_val_R):
    color = ['bo-', 'go-', 'ro-', 'co-', 'mo-', 'yo-', 'ko-']
    plt.figure()
    plt.ylim(max(cal_val_B[:, 1]) + 5, 0); plt.title('B filter')
    for a in range(2):
        plt.plot(times, cal_val_B[:, a], color[a])
    plt.xlabel('Time (UTC)'); plt.ylabel('Cal Mag')
    plt.figure()
    plt.ylim(max(cal_val_V[:, 1]) + 5, 0); plt.title('V filter')
    for a in range(2):
        plt.plot(times, cal_val_V[:, a], color[a])
    plt.xlabel('Time (UTC)'); plt.ylabel('Cal Mag')
    plt.figure()
    plt.ylim(max(cal_val_R[:, 1]) + 5, 0); plt.title('R filter')
    for a in range(2):
        plt.plot(times, cal_val_R[:, a], color[a])
    plt.xlabel('Time (UTC)'); plt.ylabel('Cal Mag')
    plt.show()
    plt.show(block=True)


# define main
if __name__ == '__main__':
    main()
