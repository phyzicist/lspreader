'''
Reader for LSP output xdr files (.p4's)

'''
import xdrlib as xdr
import numpy as np

# Define a basic function for checking (and get rid of misc.py dependency)
test = lambda d,k: k in d and d[k];

#get basic dtypes
def get_int(file,N=1,forcearray=False):
    # If forcearray is True, the output will always be a numpy array.
    ret=np.fromfile(file,dtype='>i4',count=N);
    if N==1 and not forcearray:
        return ret[0];
    return ret;

def get_float(file,N=1,forcearray=False):
    # If forcearray is True, the output will always be a numpy array.
    ret=np.fromfile(file,dtype='>f4',count=N);
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
        d = xdr.get_dict(
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
    # Re-written by Scott to be more python-syntax oriented. Reads fields files only (not scalar files).
    # Inputs:
    #   fname: a file path (string) for a fields .p4 file
    #   flds (optional): a list of vector field variables to read in. e.g. flds = ["E","B"]. If omitted, all fields available will be read in.
    # Outputs:
    #   doms: the data of the fields, as a list of domains with sub-content of the field arrays
    #   header: the header dict, as generated by "get_header()" subroutine
    # Changelog:
    #   2015-12-28  Fixed the header read-in problem, where it needs to seek forward.
    #   2015-02-29  Wrapped up the function such that it takes in the filename
    #               Changed file open argument to 'rb' rather than 'r' for compatibility with linux and windows binary file reads
    #   2015-12-30  Spun off this function, read_flds3, from read_flds2. Got rid of legacy read-in of data X, Y, Z
    #               Renamed this function read_flds2 to avoid extraneous version jumping.
    # Issues:
    #   * Unknown if this function will work for 3D field dumps. Specifically, the order of the reshaping of fields may be wrong.
    
    with open(fname, 'rb') as file:
        # Read in the header, and check that this dump is a vector fields file
        size = 3;    
        header = get_header(file) # Gets the header, and advances the file cursor to after the header
    
        if header['dump_type'] != 2:
            raise TypeError("Invalid fields file: Dump_type of header is not 2, indicating this is something other than a fields file.")
        else:
            #print "Dump type correct for fields file!"
            pass
            
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
        d=get_dict('fii',['t','step','pnum']);
        d['pos']=file.tell();
        file.seek(d['pnum']*pbytes,1);
        frames.append(d);
    for i,d in enumerate(frames):
        N = d['pnum'];
        lt=[('ip','>i4')]+zip(params,['>f4']*nparams);
        file.seek(d['pos']);
        arr=np.fromfile(file,dtype=np.dtype(lt),count=N);
        frames[i].update({'data':arr});
        del frames[i]['pos'];
    return frames;

def read_pext(file, header):
    nparams = len(header['quantities']);
    params = ['t','q','x','y','z','ux','uy','uz'];
    if nparams == 9:
        params+=['E'];
    elif nparams == 11:
        params+=['xi','yi','zi'];
    elif nparams == 12:
        params+=['E','xi','yi','zi'];
    #it's just floats here on out
    dt = list(zip(params, ['>f4']*len(params)));
    out = np.fromfile(file,dtype=dt,count=-1);
    return out;


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

def stitch2D(doms, fld_id):
    # Stitch a simple 2D lsp sim together, where we have N domains all built up along the Z dimension, and a flat Y dimension.
    # Such as with Chris' 2D sims of the back-reflected plasma
    # Inputs:
    #   doms: an output of fld_reader3, which is a list of N items, containing each the fields data for that domain
    #   fld_id: the string identifying the field component, e.g. "Ex" or "Bz"
    # Outputs:
    #   fld: A 2D array which contains the field component data
    #   xgv, zgv: 1D arrays of the X or Z coordinate along the axis, in centimeters

    fld_cat = np.squeeze(doms[0][fld_id])
    xgv = doms[0]['xgv']
    zgv_cat = doms[0]['zgv']
    
    for i in range(1,len(doms)): # Domains are concatenated along the z dimension
        fld_tmp = np.squeeze(doms[i][fld_id])[1:,:]
        fld_cat = np.concatenate((fld_cat,fld_tmp),0)

        zgv_tmp = doms[i]['zgv'][1:]
        zgv_cat = np.concatenate((zgv_cat,zgv_tmp),0)
    
    zgv = zgv_cat
    fld = fld_cat
    return fld, xgv, zgv

def pseek(file):
    # Print out the current file seek position ("file.tell")
    print "Current file seek position: " + str(file.tell())
