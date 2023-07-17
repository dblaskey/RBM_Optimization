import pandas as pd
import subprocess

outlets_df = pd.read_csv('/glade/scratch/dblaskey/RBM/outlet_list.csv')

previous_id = 81014496 # initialize previous string with what is currently in file

for outlet_id in outlets_df.Outlet.values[3:]: # remove brackets if you want all of them
    
    # creating a variable and storing the text
    # that we want to search
    search_text = str(previous_id)
    
    previous_id = outlet_id

    # creating a variable and storing the text
    # that we want to add
    replace_text = str(outlet_id)

    # Opening our text file in read only
    # mode using the open() function
    with open(r'/glade/u/home/dblaskey/RBM/Control/mizuRoute.control', 'r') as file:

        # Reading the content of the file
        # using the read() function and storing
        # them in a new variable
        data = file.read()

        # Searching and replacing the text
        # using the replace() function
        data = data.replace(search_text, replace_text)

    # Opening our text file in write only
    # mode to write the replaced content
    with open(r'/glade/u/home/dblaskey/RBM/Control/mizuRoute.control', 'w') as file:

        # Writing the replaced data in our
        # text file
        file.write(data)
    
    # Call function
    args ='/glade/work/dblaskey/mizuRoute/route/bin/route_runoff /glade/u/home/dblaskey/RBM/Control/mizuRoute.control' 
    subprocess.call(args, shell=True)
    