'''
Reader for LSP output xdr files (.p4's)

'''
import xdrlib as xdr;
import numpy as np;
import misc as m;
#get basic dtypes
def get_int(file,N=1):
    return np.fromfile(file,dtype='>i4',counts=N);
def get_float(file,N=1):
    return np.fromfile(file,dtype='>f4',counts=N);
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
            out.append(get_ints(self.file,1));
        elif i == 'f' or i == 'd':
            out.append(get_floats(self.file,1));
        elif i == 's':
            out.append(self.get_str());
        else:
            raise ValueError("Unexpected flag '{}'".format(i));
    return out;
def get_dict(file,fmt,keys):
    return dict(
        zip(keys,l), get_list(file,fmt)
    );

def get_header(file,**kw):
    '''gets the header for the .p4 file, note that this
       Advanced the file position to the end of the header.
    '''
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
        header['quantities'] = zip(names,units);
    elif header['dump_type'] == 6:
        #this is a particle movie file
        d = xdr.get_dict('iiii',
                         ['geometry','sflagsx','sflagsy','sflagsz']);
        header.update(d);
        #reading params
        n = get_int(file);
        flags=[bool(get_int(file)) for i in range(n)];
        units=[get_str(file) for i in range(n)];
        labels=['q','x','y','z','ux','uy','uz','E']
        if n == 8:
            pass;
        elif n == 7:
            labels = labels[:-1];
        elif n == 11:
            labels+=['xi','yi','zi'];
        else:
            raise NotImplementedError('Not implemented for these number of parameters:{}.'.format(n));
        header['params'] = [
            (i[0],i[1]) for i in zip(labels,units,flags) if i[2]
        ];
    elif header['dump_type'] == 10:
        #this is a particle extraction file:
        header['geometry'] = xdr.get_int();
        #reading quantities
        n = xdr.get_int();
        header['quantities'] = [self.get_str() for i in range(n)];
    else:
        raise ValueError('Unknown dump_type: {}'.format(self.header['dump_type']));
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
        nI = get_int(file); Ip = get_float(file,count=nI);
        nJ = get_int(file); Jp = get_float(file,count=nJ);
        nK = get_int(file); Kp = get_float(file,count=nK);
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
                d[quantity] = get_float(file,count=nAll*size);
                if size==3:
                    data=d[quantity].reshape(nAll,size).T;
                    d[quantity+'x'],d[quantity+'y'],d[quantity+'z']= data;
                    del data;
        doms.append(d);
    #self.logprint('Done! Stringing together.');
    out = { k : np.concatenate([i[k] for i in doms]) for k in doms[0] };
    return out;
def read_sclr(file,header,var):
    return read_flds(file, header, var, False);
    
def read(fname,**kw):
    '''reads an lsp output file into an h5 file.'''
    with open(fname,'r') as file:
        header = get_header(file);
        if test(kw, 'var'):
            var=[i[0] for i in self.header['quantities']];
        else:
            var=kw['var'];
        if self.header['dump_type'] == 2:
            d=read_flds(file,header,var);
        elif self.header['dump_type'] == 3:
            d=read_sclr(file,header,var);
        else:
            raise NotImplementedError("Other file types not implemented yet!");
    return d;
