'''
Reader for LSP output xdr files (.p4's)

'''
import xdrlib as xdr
import numpy as np
import gzip # For reading .p4.gz files
import os # For path splitting
import re # For regular expression searches of the history.p4 text file

# Define a basic function for checking (and get rid of misc.py dependency)
test = lambda d,k: k in d and d[k];

#get basic dtypes
def get_int(file,N=1,forcearray=False,lowlev=False):
    # If forcearray is True, the output will always be a numpy array.
    
    dtype = np.dtype('>i4')
    if lowlev:
        ret=np.fromfile(file,dtype=dtype,count=N)    
    else:
        s = file.read(dtype.itemsize * N)  
        ret=np.fromstring(s, dtype=dtype,count=N)
        
    if N==1 and not forcearray:
        return ret[0];
    return ret;


def read_hist(p4dir):
    """ Read probes in the 'history.p4' or 'history.p4.gz' file contained in the given directory 
    Inputs:
        p4dir: string, path to a directory containing a "history.p4" or "history.p4.gz"
    Outputs:
        values: 2D NumPy array: (N x tsteps) If there are N probes, and tsteps time steps in the simulation.
        labels: list, length N, of two-element tuples with each tuple containing ('probe name', 'units')
    
    Note: Reads the entire file twice into memory. Not efficient if for some reason history.p4 is a huge file.
    Note: Not thoroughly tested. Not a super-duper rigorous function (assumes a certain format, doesn't do as many checks as it should), so don't be surprised if it breaks.
    """
    
    ## Decide if it is a .gz or .p4
    if os.path.isfile(os.path.join(p4dir, 'history.p4')):
        histpath = os.path.join(p4dir, 'history.p4')    
        myopen = open
    elif os.path.isfile(os.path.join(p4dir, 'history.p4.gz')):
        histpath = os.path.join(p4dir, 'history.p4.gz')
        myopen = gzip.open
    else:
        raise Exception("history.p4(.gz) file does not exist in directory")
    
    ## Read in the entire file as a string to look at the header (inefficient)
    with myopen(histpath, 'r') as f:
        s = f.read()
    
    ## Look at the header and extract the number of probes and their labels (including their units)
    m1 = re.search(r"^#Number of data items: ([0-9]{1,})$", s, re.MULTILINE) # e.g. reading "Number of data items: 9" and extracting 9
    if not m1:
        raise Exception("Something's wrong with the .p4 file or reader. Couldn't find a 'Number of data items' entry.")
    else:
        nitems = int(m1.group(1)) # This is the number of probes. Assume there are this many header lines, plus four more.
    
    labels = re.findall(r"^#[0-9]{1,}?: (.*?): (.*?)$", s, re.MULTILINE) # e.g. reading "#0: time: ns", extracting the tuple ('time', 'ns')
    del s

    ## Read in the array of probe data
    values = np.genfromtxt(histpath, skip_header = nitems + 4).swapaxes(0,1)[1:] # Fortunately, genfromtext works natively on .gz files as well
    
    ## Make sure we didn't screw something up royally, that we have the same number of probe labels as probe dimensions
    if np.abs(len(labels) - values.shape[0]) > 0.5:
        raise Exception("Error reading in the history.p4 file and its header.")
    
    return values, labels
    
def get_float(file,N=1,forcearray=False,lowlev=False):
    # If forcearray is True, the output will always be a numpy array.
    dtype = np.dtype('>f4')
    if lowlev:
        ret=np.fromfile(file,dtype=dtype,count=N)    
    else:
        s = file.read(dtype.itemsize * N)  
        ret=np.fromstring(s, dtype=dtype,count=N)

    if N==1 and not forcearray:
        return ret[0];
    return ret;

def get_str(file):
    l1 = get_int(file);
    l2 = get_int(file);
    if l1 != l2:
        print("warning, string prefixes are not equal...");
        print("{}!={}".format(l1,l2));
    size=l1;
    if l1%4:
        size+=4-l1%4;
    return xdr.Unpacker(file.read(size)).unpack_fstring(l1);

def get_list(file,fmt):
    '''makes a list out of the fmt from the LspOutput f using the format
       i for int
       f for float
       d for double
       s for string'''
    out=[]
    for i in fmt:
        if i == 'i':
            out.append(get_int(file));
        elif i == 'f' or i == 'd':
            out.append(get_float(file));
        elif i == 's':
            out.append(get_str(file));
        else:
            raise ValueError("Unexpected flag '{}'".format(i));
    return out;
    
def get_dict(file,fmt,keys):
    return dict(
        zip(keys, get_list(file,fmt))
    );

def get_header(file,**kw):
    '''gets the header for the .p4 file, note that this
       Advanced the file position to the end of the header.
       
       Returns the size of the header, and the size of the header, if the header keyword is true.
    '''
    if type(file) == str:
        #if called with a filename, recall with opened file.
        with open(file,'r') as f:
            return get_header(f,**kw);
    if test(kw, "size"):
        size = file.tell();
    header = get_dict(
        file,'iiss',
        ['dump_type','dversion', 'title','revision']);
    if header['dump_type'] == 2 or header['dump_type'] == 3:
        #this is a fields file or a scalars file.
        d = get_dict(file,'fii',['timestamp','geometry','domains']);
        header.update(d);
        #reading quantities
        n = get_int(file);
        names=[get_str(file) for i in range(n)];
        units=[get_str(file) for i in range(n)];
        #print "Now I'm done reading units!"
        header['quantities'] = zip(names,units);
    elif header['dump_type'] == 6:
        #this is a particle movie file
        d = get_dict(file,
            'iiii',
            ['geometry','sflagsx','sflagsy','sflagsz']);
        header.update(d);
        #reading params
        n = get_int(file);
        flags=[bool(get_int(file)) for i in range(n)];
        units=[get_str(file) for i in range(n)];
        labels=['q','x','y','z','ux','uy','uz','E'];
        if n == 8:
            pass;
        elif n == 7:
            labels = labels[:-1];
        elif n == 11:
            labels+=['xi','yi','zi'];
        else:
            raise NotImplementedError('Not implemented for these number of parameters:{}.'.format(n));
        header['params'] = [
            (label,unit) for (label,unit,flag) in zip(labels,units,flags) if flag
        ];
    elif header['dump_type'] == 10:
        #this is a particle extraction file:
        header['geometry'] = get_int(file);
        #reading quantities
        n = get_int(file);
        header['quantities'] = [get_str(file) for i in range(n)];
    else:
        raise ValueError('Unknown dump_type: {}'.format(header['dump_type']));
    if test(kw,'size'):
        return header, file.tell()-size;
    return header;

def read_flds(file, header, var=None, vector=True):
    if vector:
        size=3;
        readin = set();
        for i in var:#we assume none of the fields end in x
            if i[-1] == 'x' or i[-1] == 'y' or i[-1] == 'z':
                readin.add(i[:-1]);
            else:
                readin.add(i);
    else:
        size=1;
        readin = set(var);
    if not var:
        var = header;
        header = get_header(file);
        if header['dump_type'] != 2 or header['dump_type'] != 3:
            print("warning: reading something that doesn't look like a flds or sclr as such!");
    doms = [];
    qs = [i[0] for i in header['quantities']];
    for i in range(header['domains']):
        iR, jR, kR = get_int(file, N=3);
        #getting grid parameters (real coordinates)
        nI = get_int(file); Ip = get_float(file,N=nI);
        nJ = get_int(file); Jp = get_float(file,N=nJ);
        nK = get_int(file); Kp = get_float(file,N=nK);
        nAll = nI*nJ*nK;
        #self.logprint('Dimensions are {}x{}x{}={}.'.format(nI,nJ,nK,nAll));
        d={}
        #self.logprint('Making points.');
        #the way meshgrid works, it has to be in this order.
        d['y'], d['z'], d['x'] = np.vstack(np.meshgrid(Jp,Kp,Ip)).reshape(3,-1);
        for quantity in qs:
            if quantity not in readin:
                file.seek(nAll*4*size,1);
            else:
                #self.logprint('Reading in {}'.format(quantity));
                d[quantity] = get_float(file,N=nAll*size);
                if size==3:
                    data=d[quantity].reshape(nAll,size).T;
                    d[quantity+'x'],d[quantity+'y'],d[quantity+'z']= data;
                    del data, d[quantity];
        doms.append(d);
    #self.logprint('Done! Stringing together.');
    out = { k : np.concatenate([i[k] for i in doms]) for k in doms[0] };
    #converting to little endian
    for k in out:
        out[k] = out[k].astype('<f4');
    return out;

def read_flds2(fname, flds=None):
    # Re-written by Scott to be more python-syntax oriented. Reads fields files and scalar files.
    # Inputs:
    #   fname: a file path (string) for a fields or scalars .p4 file
    #   flds (optional): a list of vector field variables to read in. e.g. flds = ["E","B"] (do NOT put "Ex", etc.), or scalar variables (for a sclr file). If omitted, all fields available will be read in.
    # Outputs:
    #   doms: the data of the fields, as a list of domains with sub-content of the field arrays
    #   header: the header dict, as generated by "get_header()" subroutine
    # Changelog:
    #   2015-12-28  Fixed the header read-in problem, where it needs to seek forward.
    #   2015-02-29  Wrapped up the function such that it takes in the filename
    #               Changed file open argument to 'rb' rather than 'r' for compatibility with linux and windows binary file reads
    #   2015-12-30  Spun off this function, read_flds3, from read_flds2. Got rid of legacy read-in of data X, Y, Z
    #               Renamed this function read_flds2 to avoid extraneous version jumping.
    #   2016-01-15  Added scalar read functionality back in.
    # Issues:
    #   * Unknown if this function will work for 3D field dumps. Specifically, the order of the reshaping of fields may be wrong.
    _, ext = os.path.splitext(fname)
    
    if ext == '.gz':
        myopen = gzip.open # If the file is .p4.gz, open using gunzip
    else:
        myopen = open
	
    with myopen(fname, 'rb') as file:
        # Read in the header, and check that this dump is a vector fields file or scalars file
        header = get_header(file) # Gets the header, and advances the file cursor to after the header
        if header['dump_type'] == 2:
            size = 3 # Fields file type. Each vector has three components.
        elif header['dump_type'] == 3:
            size = 1 # Scalars file type. Each scalar has one component.
        else:
            size = None
            raise TypeError("Invalid fields file: Dump_type of header is not 2 or 3, indicating this is something other than a fields or scalars file.")
            
        # Check that the requested field variables are available in this file (if any is not, abort.)
        qs = [i[0] for i in header['quantities']];
        if not flds: # If no fields specified at input, read in all available fields.
            flds = qs
        else: # If fields were specified at input, check their validity.
            for fld in flds:
                if qs.count(fld) < 1:
                    raise ValueError("Invalid field request: Field '" + fld + "' is not a stored quantity in this file.")
        #print "All requested fields are available! Reading domains."
        
        flds_set = set(flds) # A python set can be an easier type than a list when comparing elements
    
        doms = []
        
        for i in range(header['domains']): # Iterate over domains
            #print "Domain: " + str(i + 1) + " of " + str(header['domains'])
            iR, jR, kR = get_int(file, N=3)
            #getting grid parameters (real coordinates)
            nI = get_int(file)
            Ip = get_float(file,N=nI,forcearray=True)
            nJ = get_int(file)
            Jp = get_float(file,N=nJ,forcearray=True)
            nK = get_int(file)
            Kp = get_float(file,N=nK,forcearray=True)
            nAll = nI*nJ*nK
            #self.logprint('Dimensions are {}x{}x{}={}.'.format(nI,nJ,nK,nAll));
            d = {} # The domain data storage unit
            #print('Making points.')
            #the way meshgrid works, it has to be in this order.
            d['xgv'] = Ip
            d['ygv'] = Jp
            d['zgv'] = Kp

            for quantity in qs:
                if quantity not in flds_set: # Skip file cursor past unwanted fields/quantities
                    file.seek(nAll*4*size,1)
                else: # Read in the desired fields
                    #self.logprint('Reading in {}'.format(quantity));
                    # COULD PROBABLY DO THIS IN ONE STEP, BUT WHAT'S THE POINT?
                    fld_raw = get_float(file,N=nAll*size)
                    data=fld_raw.reshape(nAll,size).T
                    if size == 1:
                        d[quantity] = data.reshape(nK,nJ,nI)
                    else:
                        fld_x,fld_y,fld_z= data
                        d[quantity+'x'] = fld_x.reshape(nK,nJ,nI) # CHECK THE ORDER FOR 3D!!!
                        d[quantity+'y'] = fld_y.reshape(nK,nJ,nI)
                        d[quantity+'z'] = fld_z.reshape(nK,nJ,nI)
                    del data, fld_raw

            doms.append(d)

    return doms, header
    

def read_sclr(file,header,var):
    return read_flds(file, header, var, False);


def iseof(file):
    c = file.tell();
    file.read(1);
    if file.tell() == c:
        return True;
    file.seek(c);
    
def read_movie(file, header):
    params,_  = zip(*header['params']);
    nparams = len(params);
    pbytes = (nparams+1)*4;
    frames=[];
    pos0 = file.tell(); 
    while not iseof(file):
        d=get_dict(file, 'fii',['t','step','pnum']);
        d['pos']=file.tell();
        file.seek(d['pnum']*pbytes,1);
        frames.append(d);
    for i,d in enumerate(frames):
        N = d['pnum'];
        lt = [('ip','>i4')] + list(zip(params,['>f4']*nparams));
        file.seek(d['pos']);
        arr=np.frombuffer(file.read(N*4*len(lt)),dtype=np.dtype(lt),count=N);
        arr.flags.writeable = True;
        frames[i].update({'data':arr});
        del frames[i]['pos'];
    return frames;

def read_movie2(fname):
    """ Wrapper for read_movie() that takes a filename as input. Gzipped .p4.gz allowed.
    Inputs:
        fname: string, path to a file, e.g. fname = 'pext1.p4'
    Outputs:
        frames: The output from read_movie. A NumPy array with multiple records.
    """
    
    _, ext = os.path.splitext(fname)
    
    if ext == '.gz':
        myopen = gzip.open # If the file is .p4.gz, open using gunzip
    else:
        myopen = open
    
    with myopen(fname, 'rb') as file:
        header = get_header(file)
        frames = read_movie(file, header)
    return frames
    
def read_pext(file, header, lowlev=False):
    nparams = len(header['quantities'])
    params = ['t','q','x','y','z','ux','uy','uz']
    if nparams == 9:
        params+=['E']
    elif nparams == 11:
        params+=['xi','yi','zi']
    elif nparams == 12:
        params+=['E','xi','yi','zi']
    #it's just floats here on out
    dt = list(zip(params, ['>f4']*len(params)))
    if lowlev: # more efficient, but does not work with .p4.gz file handles
        out = np.fromfile(file,dtype=dt,count=-1)
    else: # works with all file handles
        s = file.read()
        out=np.fromstring(s, dtype=dt,count=-1)
    return out

def read_pext2(fname):
    """ Wrapper for read_pext() that takes a filename as input. Gzipped .p4.gz allowed.
    Inputs:
        fname: string, path to a file, e.g. fname = 'pext1.p4'
    Outputs:
        out: The output from read_pext. A NumPy array with multiple records.
    """
    
    _, ext = os.path.splitext(fname)
    
    if ext == '.gz':
        myopen = gzip.open # If the file is .p4.gz, open using gunzip
    else:
        myopen = open
    
    with myopen(fname, 'rb') as file:
        header = get_header(file)
        out = read_pext(file, header)
    return out

def read(fname,**kw):
    '''reads an lsp output file into an h5 file.'''
    with open(fname,'rb') as file:
        header = get_header(file);
        if not test(kw, 'var'):
            var=[i[0] for i in header['quantities']];
        else:
            var=kw['var'];
        readers = {#overkill
            2: lambda: read_flds(file,header,var),
            3: lambda: read_sclr(file,header,var),
            10:lambda: read_pext(file,header)
        };
        try:
            d = readers[header['dump_type']]();
        except KeyError:
            d = None;
            raise NotImplementedError("Other file types not implemented yet!");
    return d;

def stitch2D(doms, fld_id, divsp=1, splitax="z"):
    """ Stitch a simple 2D lsp sim together, where we have N domains all built up along the Z dimension, and a flat Y dimension.
      Such as with Chris' 2D sims of the back-reflected plasma
      Inputs:
        doms: an output of fld_reader3, which is a list of N items, containing each the fields data for that domain
        fld_id: the string identifying the field component, e.g. "Ex" or "Bz"
        divsp: integer, divisor by which to reduce the spatial resolution (e.g. divsp = 2 reduces field dimensions from 300x200 to 150x100)
      Outputs:
        fld: A 2D array which contains the field component data
     xgv, zgv: 1D arrays of the X or Z coordinate along the axis, in centimeters
     TODO: Make this generalized to 3D
    """

    if splitax == "z": # ZSPLIT lsp option
        fld_cat = np.squeeze(doms[0][fld_id])
        xgv = doms[0]['xgv'][::divsp]
        zgv_cat = doms[0]['zgv']
        
        for i in range(1,len(doms)): # Domains are concatenated along the z dimension
            fld_tmp = np.squeeze(doms[i][fld_id])[1:,:]
            fld_cat = np.concatenate((fld_cat,fld_tmp),0)
    
            zgv_tmp = doms[i]['zgv'][1:]
            zgv_cat = np.concatenate((zgv_cat,zgv_tmp),0)
        
        zgv = zgv_cat[::divsp]
        fld = fld_cat[::divsp,::divsp]
    elif splitax == "x": # XSPLIT lsp option
        fld_cat = np.squeeze(doms[0][fld_id])
        xgv_cat = doms[0]['xgv']
        zgv = doms[0]['zgv'][::divsp]
        for i in range(1,len(doms)): # Domains are concatenated along the x dimension
            fld_tmp = np.squeeze(doms[i][fld_id])[:,1:]
            fld_cat = np.concatenate((fld_cat,fld_tmp),1)
    
            xgv_tmp = doms[i]['xgv'][1:]
            xgv_cat = np.concatenate((xgv_cat,xgv_tmp),0)
        
        xgv = xgv_cat[::divsp]
        fld = fld_cat[::divsp,::divsp]
    else:
        raise Exception("Unsupported split axis: " + splitax)
    return fld, xgv, zgv

def pseek(file):
    # Print out the current file seek position ("file.tell")
    print("Current file seek position: " + str(file.tell()))
