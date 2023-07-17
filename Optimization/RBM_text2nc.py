import pandas as pd
import subprocess
import os

#outlets_df = pd.read_csv('/glade/scratch/dblaskey/RBM/outlet_list.csv')
comids = [81000402, 81000433, 81027459, 81027546, 81030789, 81030912, 81035777, 81035794] #81020021, 81014324

alpha_it = list(range(7, 16, 2))
beta_it = list(range(5, 13, 2))
gamma_it = [0.28, 0.38, 0.89]
mu_it = [0.001, 0.117, 5.001]

previous_yr = 2013 # initialize previous string with what is currently in file
previous_id = 81020021
previous_var = '7_9_0.001_0.89'

for x in range(len(alpha_it)):
    alpha = alpha_it[x]
    for j in range(len(beta_it)):
        beta = beta_it[j]
        for k in range(len(mu_it)):
            mu = mu_it[k]
            for l in range(len(gamma_it)):
                gamma = gamma_it[l]
                
                folder = '%s_%s_%s_%s'%(alpha, beta, mu, gamma)
        
                for outlet_id in comids:

                    # Convert output to NetCDF file
                    years = range(2013,2018)

                    
                    #### Change COMID #######
                    
                    # creating a variable and storing the text
                    # that we want to search
                    search_text = str(previous_id)

                    previous_id = outlet_id

                    # creating a variable and storing the text
                    # that we want to add
                    replace_text = str(outlet_id)

                    # Opening our text file in read only
                    # mode using the open() function
                    with open(r'/glade/u/home/dblaskey/RBM/Calibration/postprocess_rbm2.cfg', 'r') as file:

                        # Reading the content of the file
                        # using the read() function and storing
                        # them in a new variable
                        data = file.read()

                        # Searching and replacing the text
                        # using the replace() function
                        data = data.replace(search_text, replace_text)

                    # Opening our text file in write only
                    # mode to write the replaced content
                    with open(r'/glade/u/home/dblaskey/RBM/Calibration/postprocess_rbm2.cfg', 'w') as file:

                        # Writing the replaced data in our
                        # text file
                        file.write(data)

                    #### Change Variable #######
                    
                    # creating a variable and storing the text
                    # that we want to search
                    search_text = str(previous_var)

                    previous_var = folder

                    # creating a variable and storing the text
                    # that we want to add
                    replace_text = str(folder)

                    # Opening our text file in read only
                    # mode using the open() function
                    with open(r'/glade/u/home/dblaskey/RBM/Calibration/postprocess_rbm2.cfg', 'r') as file:

                        # Reading the content of the file
                        # using the read() function and storing
                        # them in a new variable
                        data = file.read()

                        # Searching and replacing the text
                        # using the replace() function
                        data = data.replace(search_text, replace_text)

                    # Opening our text file in write only
                    # mode to write the replaced content
                    with open(r'/glade/u/home/dblaskey/RBM/Calibration/postprocess_rbm2.cfg', 'w') as file:

                        # Writing the replaced data in our
                        # text file
                        file.write(data)
                    
                    for year in years:

                        # creating a variable and storing the text
                        # that we want to search
                        search_text = str(previous_yr)

                        previous_yr = year

                        # creating a variable and storing the text
                        # that we want to add
                        replace_text = str(year)

                        # Opening our text file in read only
                        # mode using the open() function
                        with open(r'/glade/u/home/dblaskey/RBM/Calibration/postprocess_rbm2.cfg', 'r') as file:

                            # Reading the content of the file
                            # using the read() function and storing
                            # them in a new variable
                            data = file.read()

                            # Searching and replacing the text
                            # using the replace() function
                            data = data.replace(search_text, replace_text)

                        # Opening our text file in write only
                        # mode to write the replaced content
                        with open(r'/glade/u/home/dblaskey/RBM/Calibration/postprocess_rbm2.cfg', 'w') as file:

                            # Writing the replaced data in our
                            # text file
                            file.write(data)

                        # Specify path
                        file_path = '/glade/scratch/dblaskey/RBM/Output/Calibration/%s/%s_%s.nc'%(folder, outlet_id, year)

                        if os.path.isfile(file_path) == True:
                            print("Skip ", folder)
                        
                        else:    
                            # Call function
                            args ='python /glade/u/home/dblaskey/RBM/Functions/convert_rbm_to_nc.py /glade/u/home/dblaskey/RBM/Calibration/postprocess_rbm2.cfg' 
                            subprocess.call(args, shell=True)
