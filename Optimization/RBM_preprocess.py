import pandas as pd
import subprocess
import os

# outlets_df = pd.read_csv('/glade/scratch/dblaskey/RBM/outlet_list.csv')

outlets = [81014324, 81000005, 81000402, 81000069, 81000433, 81005554, 81027459, 81030789, 81027546, 81030912, 81035777, 81035794]

previous_id = 81030789 # initialize previous string with what is currently in file

for outlet_id in outlets:  #outlets_df.Outlet.values 
    # creating a variable and storing the text
    # that we want to search
    search_text = str(previous_id)
    
    previous_id = outlet_id

    # creating a variable and storing the text
    # that we want to add
    replace_text = str(outlet_id)

    # Opening our text file in read only
    # mode using the open() function
    with open(r'/glade/u/home/dblaskey/RBM/Optimization/prepare_rbm.cfg', 'r') as file:

        # Reading the content of the file
        # using the read() function and storing
        # them in a new variable
        data = file.read()

        # Searching and replacing the text
        # using the replace() function
        data = data.replace(search_text, replace_text)

    # Opening our text file in write only
    # mode to write the replaced content
    with open(r'/glade/u/home/dblaskey/RBM/Optimization/prepare_rbm.cfg', 'w') as file:

        # Writing the replaced data in our
        # text file
        file.write(data)
    
    # Call function
    args ='python /glade/u/home/dblaskey/RBM/Optimization/generate_RBM_input_files_hru_annual.py /glade/u/home/dblaskey/RBM/Optimization/prepare_rbm.cfg' 
    subprocess.call(args, shell=True)

