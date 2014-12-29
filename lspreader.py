'''
Reader for LSP output xdr files (.p4's)

'''
import xdrlib as xdr;
import numpy as np;

class LspOutput(object):
    '''represents an lsp output file on call,
       reads the header on open'''
    def __init__(self,filename,verbose=False,prefix='',buffering=-1):
        self.filename = filename;
        self.buffering = buffering;
        self.verbose = verbose;
        self.prefix = prefix;
    
    def __enter__(self):
        self.open();
        return self;
    
    def __exit__(self, type, value, traceback):
        self.close();
        
    def open(self):
        self.file = open(self.filename, 'rb',self.buffering);
        self._get_header();
        
    def close(self):
        self.file.close();
        
    def logprint(self,s):
        if self.verbose:
            if self.prefix == '':
                print('{}'.format(s));
            else:
                print('{}: {}'.format(self.prefix,s));
    def get_chunk(self, n):
        return self.file.read(n);
    def get_int(self):
        return xdr.Unpacker(self.file.read(4)).unpack_int();
    def get_uint(self):
        return xdr.Unpacker(self.file.read(4)).unpack_uint();
    def get_float(self):
        return xdr.Unpacker(self.file.read(4)).unpack_float();
    def get_str(self):
        l1 = self.get_int();
        l2 = self.get_int();
        if l1 != l2:
            print("warning, string prefixes are not equal...");
            print("{}!={}".format(l1,l2));
        size=l1;
        if l1%4:
            size+=4-l1%4;
        return xdr.Unpacker(self.get_chunk(size)).unpack_fstring(l1);
    def get_list(self,fmt):
        '''makes a list out of the fmt from the LspOutput f using the format
        i for int
        f for float
        d for double
        s for string'''
        out=[]
        for i in fmt:
            if i == 'i':
                out.append(self.get_int());
            elif i == 'f' or i == 'd':
                out.append(self.get_float());
            elif i == 's':
                out.append(self.get_str());
            else:
                raise ValueError("Unexpected flag '{}'".format(i));
        return out;
    
    def get_dict(self,fmt,keys):
        l=self.get_list(fmt);
        return dict(zip(keys,l));
    
    def _get_header(self):
        '''gets the header for the .p4 file'''
        self.header = self.get_dict('iiss',
                                    ['dump_type','dversion',
                                     'title','revision']);
        if self.header['dump_type'] == 2 or self.header['dump_type'] == 3:
            #this is a fields file or a scalars file.
            d = self.get_dict('fii',['timestamp','geometry','domains']);
            self.header.update(d);
            #reading quantities
            n = self.get_int();
            names=[self.get_str() for i in range(n)];
            units=[self.get_str() for i in range(n)];
            self.header['quantities'] = zip(names,units);
        elif self.header['dump_type'] == 6:
            #this is a particle movie file
            d = self.get_dict('iiii',
                              ['geometry','sflagsx','sflagsy','sflagsz']);
            self.header.update(d);
            #reading params
            n = self.get_int();
            flags=[bool(self.get_int()) for i in range(n)];
            units=[self.get_str() for i in range(n)];
            labels=['q','x','y','z','ux','uy','uz','E']
            if n == 8:
                pass;
            elif n == 7:
                labels = labels[:-1];
            elif n == 11:
                labels+=['xi','yi','zi'];
            else:
                raise ValueError('Incorrect number of parameters: {}'.format(n));
            self.header['params'] = [(i[0],i[1]) for i in zip(labels,units,flags)
                                     if i[2]];
        elif self.header['dump_type'] == 10:
            #this is a particle extraction file:
            self.header['geometry'] = self.get_int();
            #reading quantities
            n = self.get_int();
            self.header['quantities'] = [self.get_str() for i in range(n)];
        else:
            raise ValueError('Unknown dump_type: {}'.format(self.header['dump_type']));
        return;
    ###################
    #data processing    
    def _getfields(self, var, vector=True):
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
        doms = [];
        qs = [i[0] for i in self.header['quantities']];
        self.logprint('Reading positions and making buffers');
        for i in range(self.header['domains']):
            self.logprint('Reading domain {}.'.format(i));
            iR, jR, kR = self.get_int(),self.get_int(),self.get_int();
            #getting grid parameters (real coordinates)
            nI = self.get_int(); Ip = np.fromfile(self.file,dtype='>f4',count=nI);
            nJ = self.get_int(); Jp = np.fromfile(self.file,dtype='>f4',count=nJ);
            nK = self.get_int(); Kp = np.fromfile(self.file,dtype='>f4',count=nK);
            nAll = nI*nJ*nK;
            self.logprint('Dimensions are {}x{}x{}={}.'.format(nI,nJ,nK,nAll));
            d={}
            self.logprint('Making points.');
            #the way meshgrid works, it has to be in this order.
            d['y'], d['z'], d['x'] = np.vstack(np.meshgrid(Jp,Kp,Ip)).reshape(3,-1);
            for quantity in qs:
                if quantity not in readin:
                    self.file.seek(nAll*4*size,1);
                else:
                    self.logprint('Reading in {}'.format(quantity));
                    d[quantity] = np.fromfile(self.file,dtype='>f4',count=nAll*size);
                    if size==3:
                        data=d[quantity].reshape(nAll,size).T;
                        d[quantity+'x'],d[quantity+'y'],d[quantity+'z']= data;
                        del data;
            doms.append(d);
        self.logprint('Done! Stringing together.');
        out = { k : np.concatenate([i[k] for i in doms]) for k in doms[0] };
        return out;

    def _getmovie(self):
        params,_  = zip(*self.header['params']);
        nparams = len(params);
        pbytes = (nparams+1)*4;
        frames=[];
        while True:
            c=self.tell();
            self.read(1); #python, y u no eof?
            if self.tell() == c:
                break;
            self.seek(c);
            d=self.get_dict('fii',['t','step','pnum']);
            self.logprint('scanning frame at lsp step {}'.format(d['step']));
            d['pos']=self.tell();
            self.seek(d['pnum']*pbytes,1);
            frames.append(d);
        self.logprint('converting frames');
        for i,d in enumerate(frames):
            N = d['pnum'];
            lt=[('ip','>i4')]+zip(params,['>f4']*nparams);
            dt=np.dtype(lt);
            self.seek(d['pos']);
            arr=np.fromfile(self.file,dtype=dt,count=N);
            d['data']=arr;
            del d['pos'];
            self.logprint('done!');
            #checking eof
            c=self.tell();
            self.read(1); #python, y u no eof?
            self.logprint("Do we have eof? {}".format(self.tell() == c));            
            frames[i] = d;
        return frames;

    def _getpext(self):
        nparams = len(self.header['quantities']);
        params = ['t','q','x','y','z','ux','uy','uz'];
        if nparams == 9:
            params+=['E'];
        elif nparams == 11:
            params+=['xi','yi','zi'];
        elif nparams == 12:
            params+=['E','xi','yi','zi'];
        #it's just floats here on out
        dt = list(zip(params, ['>f4']*len(params)));
        out = np.fromfile(self.file,dtype=dt,count=-1);
        return out;
        
    def get_data(self,var=None):
        if not var and (self.header['dump_type']== 2 or self.header['dump_type']== 3):
            var=[i[0] for i in self.header['quantities']];
        if self.header['dump_type'] == 2:
            return self._getfields(var,vector=True);
        elif self.header['dump_type'] == 3:
            return self._getfields(var,vector=False);
        elif self.header['dump_type'] == 6:
            return self._getmovie();
        elif self.header['dump_type'] == 10:
            return self._getpext();
        else:
            return None;
    pass;
