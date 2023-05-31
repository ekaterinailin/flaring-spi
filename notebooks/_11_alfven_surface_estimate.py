"""
UTF-8, Python 3

------------------
Flaring SPI
------------------

Judy Chebly, 2023, MIT License

Estimation of the alfven surface size in au using the results table. 

The input parameters are the following:

- The magnetic field strength derived from Rossby number(column name: B_G)
- The maximum magnetic field strength (column name: B_G_high)
- The minimum magnetic field strength (column name: B_G_low)
- Star radius (Column name: st_rad)
- Star radius minimum error (Column: st_rad_err1)
- Star radius maximum error (Column: st_rad_err2)
- Star system (Column name: ID) 

Output:

- Alfven surface (AS) estimate in au with errors. In the errors we take in 
consideration the error on the magetic field 
strength, on the radius, and the slope and intercept (in the equation).
- The code will generate  a table with 4 columns corresponding respectively 
to the star's;
name, AS in au, AS min error, AS max error, 15% of the magnetic field strength 
derived from Rossby number

"""

import numpy as np
import pandas as pd

import os

### Put the columns name as input parameter

def Alfven_surface_size_au (B_G, B_G_high, B_G_low, st_rad, ID, st_rad_err1, st_rad_err2):
    
    #### Reading the table an defining the different parameters
    
    path = '../results/'
    input_file = 'results'
    input_file = pd.read_csv(f"{path}{input_file}.csv")

    
    ########### Remove rows corresponding to systems with unknown radius and 
    # magnetic field
    B_G_              = input_file[B_G]
    radius_           = input_file[st_rad]
    input_file_cleaned  = input_file [radius_ .notna()]
    input_file_cleaned  = input_file [B_G_ .notna()]
    pd.set_option('display.max_columns', None)
    input_file.head()
    print(input_file_cleaned)
    
    ##########  Work with the cleaned table (After removing the rows with nan 
    # values in column radius and magnetic field)
    B_G              = input_file_cleaned[B_G]
    B_G_high         = input_file_cleaned[B_G_high]
    B_G_low          = input_file_cleaned[B_G_low]
    radius           = input_file_cleaned[st_rad]
    ID               = input_file_cleaned[ID]
    radius_min_error = input_file_cleaned[st_rad_err1]
    radius_max_error = input_file_cleaned[st_rad_err2]
    
    print(len(B_G))
    print(len(radius))
    
    ######################## Calculate B_G errors###############################
    B_G_min_error = [a-b for a,b in zip(B_G, B_G_low)]
    print(f'Minimum B_G error  = {B_G_min_error}\n')
    
    B_G_max_error = [a-b for a,b in zip(B_G_high, B_G) ]
    print(f'Maximum B_G error  = {B_G_max_error}\n')
    
    
    ##################### Calculate 15% of B_G, B_G_min_error, B_G_max_error###
    #Command to get the percentage of a number
    percent = lambda part, whole:float(whole) / 100 * float(part) 
    
    B_G_15perc = [percent(15,a) for a in B_G] 
    print(f'15% of B_G  = {B_G_15perc}\n')
    
    B_G_min_error_15perc = [percent(15,a) for a in B_G_min_error] 
    print(f'15% of B_G_min_error  = {B_G_min_error_15perc}\n')
    
    B_G_max_error_15perc = [percent(15,a) for a in B_G_max_error] 
    print(f'15% of B_G_max_error  = {B_G_max_error_15perc}\n')
    
    ####### Calculating the Alfven surface in stellar radius with the errors####
    a = 0.44
    a_min_error = -0.05
    a_max_error = 0.05
    
    b = 0.54
    b_min_error = -0.08
    b_max_error = 0.08
    
    # Calculate AS in stellar radius
    log_AS_Rstar = [a*np.log10(c) + b for c in B_G_15perc]
    AS_Rstar     = [10**(d) for d in log_AS_Rstar]
    print(f'AS_Rstar = {AS_Rstar}\n')
    
    # Calculate AS min in stellar radius
    log_AS_min = [(a+a_min_error)*np.log10(c-f) + (b+b_min_error) for c,f in zip(B_G_15perc, B_G_min_error_15perc)]
    AS_min     = [10**(d) for d in log_AS_min]
    print(f'AS_min = {AS_min}\n')
    
    # Calculate AS max in stellar radius
    log_AS_max = [(a+a_max_error)*np.log10(c+f) + (b+b_max_error) for c,f in zip(B_G_15perc, B_G_max_error_15perc)]
    AS_max     = [10**(d) for d in log_AS_max]
    print(f'AS_max = {AS_max}\n')


  ########## Calculating the Alfven surface in au with the errors ##############


    ###### Convert stellar radius to astronomical unit (au)

    def Rstar_AU_list(distance_Rstar, star_radius):
        Rsun = 215
        Rstar_to_Rsun   = [a*b for a,b in zip(distance_Rstar, star_radius)]
        Rsun_to_AU_calc = [d/Rsun for d in Rstar_to_Rsun]
        return  Rsun_to_AU_calc

    
    # Calculate AS in au        
    AS_au = Rstar_AU_list(AS_Rstar , radius)
    print(f'AS in au = {AS_au}\n')
    print(len(AS_au))
    
    #Convert AS min from stellar radius to au
    radius_min = [a-b for a,b in zip(radius, radius_min_error)]
    AS_min_au = Rstar_AU_list(AS_min, radius_min)
    print(f'AS min in au = {AS_min_au}\n')
    
    #Convert AS max from stellar radius to au
    radius_max = [a+b for a,b in zip(radius, radius_max_error)]
    AS_max_au = Rstar_AU_list(AS_max, radius_max)
    print(f'AS max in au = {AS_max_au}\n')

    #Calculate AS min error in au
    AS_min_error_au = [a- b for a,b in zip(AS_au , AS_min_au)]
    print(f'AS_min_error_au = {AS_min_error_au}\n')
    
    #Calculate AS max error in au
    AS_max_error_au = [a - b for a,b in zip(AS_max_au , AS_au)]
    print(f'AS_max_error_au = {AS_max_error_au}\n')
    

    ## Create csv table to store the estimated Alfven surface size with the errors######
    
    output_table = "AS_estimation.txt"
    OutputData_values = pd.read_csv(output_table) if os.path.exists(output_table) else pd.DataFrame([])
    OutputData_values = [ID, AS_au,  AS_min_error_au ,   AS_max_error_au , B_G_15perc ]
    data = {'ID': OutputData_values[0],'AS_au':OutputData_values[1], 
            'AS_au_min_error':OutputData_values[2] ,'AS_au_max_error':OutputData_values[3],
            'B_G_15perc':OutputData_values[4]}
    
    Dataframe = pd.DataFrame(data)

    Dataframe.to_csv(f'{path}AS_estimation.csv', index=False)
    path_to_paper = "/home/ekaterina/Documents/002_writing/flaring-spi-paper/src/data/"
    Dataframe.to_csv(f'{path_to_paper}AS_estimation.csv', index=False)

if __name__ == '__main__':
    
    Alfven_surface_size_au('B_G', 'B_G_high', 'B_G_low','st_rad',
                           'ID', 'st_rad_err1', 'st_rad_err2')