'''
For several ranges along coast, save times drifters enter and exit from other
ranges.
'''




    # Loop through along-coast boxes to find which other boxes they are connected to
    mat = np.zeros((len(paths),len(paths)))
    years = np.arange(2014,2015)
    months = [1,2,7,8]
    for year in years:
        for month in months:
            matfile = 'calcs/alongcoastconn/BAY-' + str(year) + '-' + str(month).zfill(2) + '.npz'
