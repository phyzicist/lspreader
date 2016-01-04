# Call the actual stitching functions from a specific, local location

#from h5stitch2D import h5fields2D

#folder = r'/media/sf_temp/lsp laser test fields'
#h5path = h5fields2D(folder)


import h5stitch2D as h52D

folder = r'/media/sf_temp/lsp laser test fields'

h5path = h52D.h5fields2Da(folder)
h5path = h52D.h5fields2Db(folder)
