import pandas as pd
import subprocess
import os

#outlets_df = pd.read_csv('/glade/scratch/dblaskey/RBM/outlet_list.csv')
comids =  [81014324, 81000005, 81000402, 81000069, 81000433, 81005554, 81027459, 81030789, 81027546, 81030912, 81035777, 81035794]

df_runs = pd.read_csv("/glade/u/home/dblaskey/RBM/Optimization/Opt_runs.csv")

years = range(2013,2018)

previous_yr = 2017 # initialize previous string with what is currently in file
previous_id = 81035777
previous_var = 'pe_basin_0_0188'

for i in range(200):  
    folder = df_runs.iloc[i,0]
                
    # Specify path
    path = '/glade/scratch/dblaskey/RBM/Output/Optimize/%s/'%folder

    # Check whether the specified
    # path exists or not
    if os.path.exists(path) == True:
        print("Skip Creating ", folder)            
    else:
        print("Creating ", folder)
        os.mkdir(path)

    for outlet_id in comids:

        # Specify path
        file_path = '/glade/scratch/dblaskey/RBM/Output/Optimize/%s/%s_2017.temp'%(folder, outlet_id)

        if os.path.isfile(file_path) == True:
            print("Skip ", outlet_id, " in ", folder)
                        
        else:
            print("Running ", outlet_id, " in ", folder)
            # Call function to run RBM
            args ='/glade/work/dblaskey/RBM/src/rbm10_mizu /glade/scratch/dblaskey/RBM/RBM_Input/Optimize/%s/%s /glade/scratch/dblaskey/RBM/Output/Optimize/%s/%s'%(folder, outlet_id, folder, outlet_id) 
            subprocess.call(args, shell=True)

        # Specify path
        file_path = '/glade/scratch/dblaskey/RBM/Output/Optimize/%s/%s_2017.nc'%(folder, outlet_id)
        # Convert output to NetCDF file                  
        if os.path.isfile(file_path) == True:
            print("Skip making nc for ", outlet_id, " in ", folder)
                        
        else:
            print("Running text2nc for ", outlet_id, " in ", folder)

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
            with open(r'/glade/u/home/dblaskey/RBM/Optimization/postprocess_rbm.cfg', 'r') as file:

                # Reading the content of the file
                # using the read() function and storing
                # them in a new variable
                data = file.read()

                            # Searching and replacing the text
                            # using the replace() function
                data = data.replace(search_text, replace_text)

                        # Opening our text file in write only
                        # mode to write the replaced content
            with open(r'/glade/u/home/dblaskey/RBM/Optimization/postprocess_rbm.cfg', 'w') as file:

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
            with open(r'/glade/u/home/dblaskey/RBM/Optimization/postprocess_rbm.cfg', 'r') as file:

                            # Reading the content of the file
                            # using the read() function and storing
                            # them in a new variable
                data = file.read()

                            # Searching and replacing the text
                            # using the replace() function
                data = data.replace(search_text, replace_text)

                        # Opening our text file in write only
                        # mode to write the replaced content
            with open(r'/glade/u/home/dblaskey/RBM/Optimization/postprocess_rbm.cfg', 'w') as file:

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
                with open(r'/glade/u/home/dblaskey/RBM/Optimization/postprocess_rbm.cfg', 'r') as file:

                                # Reading the content of the file
                                # using the read() function and storing
                                # them in a new variable
                    data = file.read()

                                # Searching and replacing the text
                                # using the replace() function
                    data = data.replace(search_text, replace_text)

                            # Opening our text file in write only
                            # mode to write the replaced content
                with open(r'/glade/u/home/dblaskey/RBM/Optimization/postprocess_rbm.cfg', 'w') as file:

                                # Writing the replaced data in our
                                # text file
                    file.write(data)

                # Call function
                args ='python /glade/u/home/dblaskey/RBM/Functions/convert_rbm_to_nc.py /glade/u/home/dblaskey/RBM/Optimization/postprocess_rbm.cfg' 
                subprocess.call(args, shell=True)

