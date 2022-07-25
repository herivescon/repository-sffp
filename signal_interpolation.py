"""""
 @Purpose: Input JC Filter transmission and QE of U47 CCD and interpolate to wavelengths desired
 @Input: Files in *.csv format
 @Output: Interpolated spectra
 @Written by: Hector Erives
 @Date: Summer 2022
"""""

import numpy as np
import PySimpleGUI as sg
import pandas as pd
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

dir_separ = '/'


def main():

    # resample to a wv_sel set of wevelengths
    wv_sel = np.arange(400, 1000, 1)

    # read solar panel reflectances and resample
    filename = popup_get_file("Select solar panel SiK7 reflectance file")
    if filename is None or filename == '':
        sg.PopupCancel('Cancelled - No valid folder entered')
        return
    sik7_cell = pd.read_csv(filename)
    sik7_sel = wv_sel
    for col in range(3):
        wv = sik7_cell.iloc[:, 0]
        ref = sik7_cell.iloc[:, col + 1]
        interp = interp1d(wv, ref, kind="linear")
        ref_int = interp(wv_sel)
        sik7_sel = np.column_stack((sik7_sel, ref_int))
    sik7_avg = sik7_sel[:, [1, 2, 3]].mean(axis=1)

    filename = popup_get_file("Select solar panel GaAsGe reflectance file")
    if filename is None or filename == '':
        sg.PopupCancel('Cancelled - No valid folder entered')
        return
    gaasge_cell = pd.read_csv(filename)
    gaasge_sel = wv_sel
    for col in range(3):
        wv = gaasge_cell.iloc[:, 0]
        ref = gaasge_cell.iloc[:, col + 1]
        interp = interp1d(wv, ref, kind="linear")
        ref_int = interp(wv_sel)
        gaasge_sel = np.column_stack((gaasge_sel, ref_int))
    gaasge_avg = gaasge_sel[:, [1, 2, 3]].mean(axis=1)

    # read filter transmission and resample
    # resample filter transmission
    filename = popup_get_file("Select ORIGINAL Filter transmission file")
    if filename is None or filename == '':
        sg.PopupCancel('Cancelled - No valid folder entered')
        return
    Filters = pd.read_csv(filename)
    wv = Filters.Wavelength; B = Filters.B
    V = Filters.V; R = Filters.R
    interp = interp1d(wv, B, kind="linear")
    B_sel = interp(wv_sel)
    interp = interp1d(wv, V, kind="linear")
    V_sel = interp(wv_sel)
    interp = interp1d(wv, R, kind="linear")
    R_sel = interp(wv_sel)

    # resample QE and resample
    filename = popup_get_file("Select ORIGINAL Quantum Efficiency file")
    if filename is None or filename == '':
        sg.PopupCancel('Cancelled - No valid folder entered')
        return
    CCD = pd.read_csv(filename)
    wv = CCD.Wavelength; QE = CCD.QE
    interp = interp1d(wv, QE, kind="linear")
    QE_sel = interp(wv_sel)
    B_sel = B_sel*QE_sel
    V_sel = V_sel*QE_sel
    R_sel = R_sel*QE_sel

    # request place where transmission will be saved for DIRSIG
    folder_path = sg.popup_get_folder("Please enter folder where filter transmission will be saved.")
    if folder_path is None or folder_path == '':
        sg.PopupCancel('Cancelled - No valid folder entered')
        return

    # plot and save resampled filter transmission
    plt.figure()
    plt.plot(wv_sel, B_sel, color='blue')
    plt.plot(wv_sel, V_sel, color='orange')
    plt.plot(wv_sel, R_sel, color='red')
    plt.title('Filter transmission with QE')
    plt.xlabel('Wavelength (nm)'); plt.ylabel('Filter transmission (%)')
    plt.show()
    plt.show(block=True)
    # save resampled filter transmission
    # B transmission
    filename = folder_path + dir_separ + 'B_filter_transmission.csv'
    header = ['Wavelength', 'B']
    data = np.column_stack((wv_sel, B_sel))
    data = pd.DataFrame(data, columns=header)
    data.to_csv(filename, index=False)

    # V filter transmission
    filename = folder_path + dir_separ + 'V_filter_transmission.csv'
    header = ['Wavelength', 'V']
    data = np.column_stack((wv_sel, V_sel))
    data = pd.DataFrame(data, columns=header)
    data.to_csv(filename, index=False)

    # R filter transmission
    filename = folder_path + dir_separ + 'R_filter_transmission.csv'
    header = ['Wavelength', 'R']
    data = np.column_stack((wv_sel, R_sel))
    data = pd.DataFrame(data, columns=header)
    data.to_csv(filename, index=False)

    # write a xml file with the new wavelengths
    filename2 = folder_path + dir_separ + 'xml_transmission.txt'
    with open(filename2, 'w') as f:
        f.write('              <spectralresponse>\n')
        f.write('                <bandpass spectralunits="nanometers">\n')
        f.write('                  <minimum>400</minimum>\n')
        f.write('                  <maximum>1000</maximum>\n')
        f.write('                  <delta>1</delta>\n')
        f.write('                </bandpass>\n')
        f.write('                <channellist>\n')
        f.write('                  <channel shape="tabulated" bias="0" name="B Channel" gain="1" normalize="true">\n')
        for i, cent in enumerate(wv_sel):
            f.write('                    <entry>\n')
            f.write('                      <spectralpoint>%0.1f</spectralpoint>\n' % cent)
            f.write('                      <value>%0.6f</value>\n' % B_sel[i])
            f.write('                    </entry>\n')
        f.write('                    <polarizer type="none"/>\n')
        f.write('                    <dataoutput type="total"/>\n')
        f.write('                  </channel>\n')
        f.write('                  <channel shape="tabulated" bias="0" name="V Channel" gain="1" normalize="true">\n')
        for i, cent in enumerate(wv_sel):
            f.write('                    <entry>\n')
            f.write('                      <spectralpoint>%0.1f</spectralpoint>\n' % cent)
            f.write('                      <value>%0.6f</value>\n' % V_sel[i])
            f.write('                    </entry>\n')
        f.write('                    <polarizer type="none"/>\n')
        f.write('                    <dataoutput type="total"/>\n')
        f.write('                  </channel>\n')
        f.write('                  <channel shape="tabulated" bias="0" name="R Channel" gain="1" normalize="true">\n')
        for i, cent in enumerate(wv_sel):
            f.write('                    <entry>\n')
            f.write('                      <spectralpoint>%0.1f</spectralpoint>\n' % cent)
            f.write('                      <value>%0.6f</value>\n' % R_sel[i])
            f.write('                    </entry>\n')
        f.write('                    <polarizer type="none"/>\n')
        f.write('                    <dataoutput type="total"/>\n')
        f.write('                  </channel>\n')
        f.write('                </channellist>\n')
        f.write('               </spectralresponse>\n')

    # generate resampled diffraction gratings and plot them
    sig = 50; wv = np.arange(400, 1050, 50)
    plt.figure()
    for i, cent in enumerate(wv):
        func_sel = np.exp(-(wv_sel - cent) ** 2 / (2 * sig ** 2))*QE_sel
        filename = folder_path + dir_separ + str(cent)+'nm_hyperspec_transmission.csv'
        header = ['Wavelength', 'Transmission']
        data = np.column_stack((wv_sel, func_sel))
        data = pd.DataFrame(data, columns=header)
        data.to_csv(filename, index=False)
        plt.plot(wv_sel, func_sel, color='red')
    plt.title('Diffraction grating transmission with QE')
    plt.xlabel('Wavelength (nm)'); plt.ylabel('Transmission (%)')
    plt.show()
    plt.show(block=True)

    # plot and generate resampled solar cell reflectance
    plt.figure()
    for i in range(3):
        plt.plot(wv_sel, sik7_sel[:, i+1], color='blue')
    plt.plot(wv_sel, sik7_avg, color='red')
    for i in range(3):
        plt.plot(wv_sel, gaasge_sel[:, i+1], color='green')
    plt.plot(wv_sel, gaasge_avg, color='red')
    # save SiK7 reflectance
    filename = folder_path + dir_separ + 'SiK7_reflectance.csv'
    header = ['Wavelength', 'Reflectance']
    data = np.column_stack((wv_sel, sik7_avg))
    data = pd.DataFrame(data, columns=header)
    data.to_csv(filename, index=False)
    # save GaAsGe reflectance
    filename = folder_path + dir_separ + 'GaAsGe_reflectance.csv'
    data = np.column_stack((wv_sel, gaasge_avg))
    data = pd.DataFrame(data, columns=header)
    data.to_csv(filename, index=False)
    # plot reflectance
    plt.title('Solar cell reflectance')
    plt.xlabel('Wavelength (nm)'); plt.ylabel('Reflectance (%)')
    plt.ylim(0, 0.25)
    plt.xlim(400, 1000)
    plt.grid()
    plt.show()
    plt.show(block=True)


def popup_get_file(message, title=None):

    layout = [
        [sg.Text(message)],
        [sg.Input(key='-INPUT-'), sg.FilesBrowse('Browse')],
        [sg.Button('Ok'), sg.Button('Cancel')],
    ]
    window = sg.Window(title if title else message, layout)
    event, values = window.read(close=True)
    return values['-INPUT-'] if event == 'Ok' else None


# Define main method
if __name__ == '__main__':
    main()