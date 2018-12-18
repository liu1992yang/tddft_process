import sys
import os
import subprocess
import pandas as pd
import math
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    print("usage: python extract-tddft-process.py td-dft.log")
    sys.exit()

HWHM = 18
S2_FILTER = 2.6
FILE = sys.argv[1]
SCALER = 1
WAVE_START = 200
WAVE_END = 1100


def grep_excited():
  """
  input str
  output list of string
  """
  output=subprocess.check_output(["grep", "^ Excited State ", FILE]).decode('utf-8')
  states=output.splitlines()
  return states


def clean_data(data_strlist):
  """
  input str
  return pd.df
  """
  states=[i.strip().split() for i in data_strlist ]
  df=pd.DataFrame(states)
  df=df.iloc[:,[2,4,6,8,9]]
  df.columns=['Excite_state','Energy(eV)', 'wavelength(nm)', 'Oscilation_Strength', '<S**2>']
  df['Excite_state']=df['Excite_state'].str[:-1].astype(int)
  df['wavelength(nm)']=pd.to_numeric(df['wavelength(nm)'], downcast='float')
  df['Energy(eV)']=pd.to_numeric(df['Energy(eV)'], downcast='float')
  df['Oscilation_Strength']=df['Oscilation_Strength'].str[2:].astype(float)
  df['<S**2>']=df['<S**2>'].str[7:].astype(float)
  return df

def filter_real_excited(df):
  return df[df['<S**2>'] <= S2_FILTER]
  

  
def lorentzian_pdf_intensity(wl, peak, intensity):
  """
  input int wl, float peak, float instensity 
  output float
  """
  return intensity*SCALER/(1+math.pow(((wl-peak)/HWHM),2.0))
  

  
def map_reduce_peaks(wl, peaks, intensities):
  """
  input int wl,
  input list of float peaks
  input list of float intensities
  return list of float
  """
  densities = list(map(lambda a, b: lorentzian_pdf_intensity(wl, a, b), peaks, intensities))
  return sum(densities)

def gen_spec(wls, peaks, intensities):
  """
  input list of float peaks
  input list of float intensities
  """
  specs = [map_reduce_peaks(wl, peaks, intensities) for wl in wls]
  return specs

def get_T(specs):
  return [spec/max(specs) for spec in specs]
  
def get_abs(specs):
  return [-math.log10(1-x) for x in specs]
  

def gen_whole_set(df):
  """
  df pd.dataFrame
  rtype pd.dataFrame
  """
  peaks = list(df['wavelength(nm)'])
  intensities = list(df['Oscilation_Strength'])
  wls = list(range(WAVE_START,WAVE_END+1))
  specs = gen_spec(wls,peaks, intensities)
  whole=pd.concat(map(lambda x: pd.Series(x), [wls, specs, get_T(specs),get_abs(specs)]),axis=1)
  whole.columns=['wavelength_number(nm)','spectrum','T','Absorbance']
  return whole

def plot_vline(xvalues, yvalues):
  plt.stem(xvalues,yvalues,basefmt=' ', markerfmt=' ')
  
def plot_specs(wavelength_number, spectrum):
  plt.plot(wavelength_number,spectrum, color='black')

if __name__ == "__main__":  
  states = grep_excited()
  cleaned_df = clean_data(states)
  cleaned_df_fn = FILE.replace('.log', '_td_excited_states.csv')
  print(cleaned_df_fn)
  cleaned_df.to_csv(cleaned_df_fn,index=False)
  real_df=filter_real_excited(cleaned_df)
  whole_set= gen_whole_set(real_df)
  #plot 
  plt.figure()
  plot_vline(real_df['wavelength(nm)'],real_df['Oscilation_Strength'])
  plot_specs(whole_set['wavelength_number(nm)'],whole_set['spectrum'])
  plt.title(FILE.replace('.log','spectrum')+'HWHM'+str(HWHM))
  plt.axis([WAVE_START,800, -0.02,whole_set['spectrum'].max()*1.1]) #up to 800 since experiment only run up to 700nm
  plt.savefig(FILE.replace('.log','_td_spectrum.png'))
  final_data=pd.concat([whole_set, real_df['Excite_state'],real_df['wavelength(nm)'],real_df['Oscilation_Strength']], axis=1)
  final_data.to_csv(FILE.replace('.log','_td_spectrum_')+'HWHM'+str(HWHM)+'.csv',index=False)
  
