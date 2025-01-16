import os, datetime

def make_directory(maindir = None, verbose = None):
    """
    Function to create directory in you current working directory.
    The function will have time stamp assigned
    
    Parameters: 
    dirname : name of main directory to hold newly created directories
    
    """
    
#     os.chdir('..') # go one step back to the current dir
    
    if maindir is False or  maindir is True:
        raise ValueError("dirname can't be True or False")
    
    if maindir is None:
        mydir = os.path.join(os.getcwd(),
                     datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
        
    elif maindir is not None:
        mydir = os.path.join(os.getcwd(), maindir,
                     datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
        
    os.makedirs(mydir)
        
    if verbose:
        print("Directory created:", mydir)
        
    return(mydir)
