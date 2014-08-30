import xdrlib as xdr;
import itertools as itools;
import multiprocessing;

class LazyIter(object):
    def __init__(self,doms):
        self.doms = doms;
        pass;
    def _next_dom(self):
        self.dom = self.domsi.next();#should throw StopIteration at the last dom
        self.curi = 0;
        self.keys = [i for i in self.dom
                     if  (i != 'xp' and i != 'yp' and i != 'zp')];
        self.l = len(self.dom[self.keys[0]]);
        self.tmp1 = len(self.dom['xp']);
        self.tmp2 = self.tmp1*len(self.dom['yp']);
        pass;
    def __iter__(self):
        self.domsi = self.doms.__iter__();
        self._next_dom();
        return self;
    def next(self):      
        if self.curi == self.l:
            self._next_dom();
        out={};
        for k in self.keys:
            out[k] = self.dom[k][self.curi];
        #lazy evaluate point
        out['z']=self.dom['zp'][self.curi/self.tmp2];
        out['y']=self.dom['yp'][self.curi%self.tmp2/self.tmp1];
        out['x']=self.dom['xp'][self.curi%self.tmp1];
        self.curi+=1;
        return out;

class Callable(object):
    def __init__(self,d,call_func):
        self.d=d;
        self.call=call_func
    def __call__(self,i):
        c=self.call;
        return c(self.d,i);
    pass;
def make_points(d):
    d['x']=[]; d['y']=[]; d['z']=[];
    tmp1 = len(d['xp']);
    tmp2 = tmp1 * len(d['yp']);
    for i in range(d['nAll']):
        d['z'].append(d['zp'][i/tmp2]);
        d['y'].append(d['yp'][i%tmp2/tmp1]);
        d['x'].append(d['xp'][i%tmp1]);
    del d['zp'],d['yp'],d['xp'];
    return d;
    
def read_fields(d):
    for quantity in d['dqs']:
        d[quantity+'x']=[];
        d[quantity+'y']=[];
        d[quantity+'z']=[];
        for k in range(d['nAll']):
            d[quantity+'x'].append(d['x'+quantity].unpack_float());
            d[quantity+'y'].append(d['x'+quantity].unpack_float());
            d[quantity+'z'].append(d['x'+quantity].unpack_float());
        del d['x'+quantity];
    del d['dqs'],d['nAll'];
    return d;

def read_scalars(d):
    for quantity in d['dqs']:
        d[quantity]=[ d['x'+quantity].unpack_float() for k in range(d['nAll'])];
        del d['x'+quantity];
    del d['dqs'],d['nAll'];
    return d;

def convert_frame(d):
    d['ip']=[];
    for param,units in d['params']:
        d[param] = [];
    for i in range(d['pnum']):
        d['ip'].append(d['xdr'].unpack_int())
        for param,units in d['params']:
            d[param].append(d['xdr'].unpack_float());
    del d['xdr'];
    return d;

def convert_particles(d,i):
    d['xdr'].set_position(d['pbytes']*i);
    p={};
    p['ip'] = d['xdr'].unpack_int();
    for param,units in d['params']:
        p[param] = d['xdr'].unpack_float();
    return p;

class LspOutput(file):
    '''represents an lsp output file on call,
       reads the header on open'''
    def __init__(self,filename,verbose=False,prefix='',buffering=-1):
        file.__init__(self,filename,"rb",buffering);
        self._get_header();
        self.verbose = verbose;
        self.prefix = prefix;
    def logprint(self,s):
        if self.verbose:
            if self.prefix == '':
                print('{}'.format(s));
            else:
                print('{}: {}'.format(self.prefix,s));
    def get_chunk(self, n):
        return self.read(n);
    def get_int(self):
        return xdr.Unpacker(self.read(4)).unpack_int();
    def get_uint(self):
        return xdr.Unpacker(self.read(4)).unpack_uint();
    def get_float(self):
        return xdr.Unpacker(self.read(4)).unpack_float();
    def get_str(self):
        l1 = self.get_int();
        l2 = self.get_int();
        if l1 != l2:
            print("warning, string prefixes are not equal...");
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
    def _getfields(self, var, pool_size,lazy,vector=True):
        if vector:
            size=3; call=read_fields;
            readin = set();
            for i in var:#we assume none of the fields end in x
                if i[-1] == 'x' or i[-1] == 'y' or i[-1] == 'z':
                    readin.add(i[:-1]);
                else:
                    readin.add(i);
        else:
            size=1; call=read_scalars;
            readin = set(var);
        doms = [];
        qs = [i[0] for i in self.header['quantities']];
        pool=multiprocessing.Pool(pool_size);
        self.logprint('reading positions and making buffers');
        for i in range(self.header['domains']):
            self.logprint('reading domain {}'.format(i));
            iR, jR, kR = self.get_int(),self.get_int(),self.get_int();
            #getting grid parameters (real coordinates)
            nI = self.get_int(); Ip = [self.get_float() for i in range(nI)];
            nJ = self.get_int(); Jp = [self.get_float() for i in range(nJ)];
            nK = self.get_int(); Kp = [self.get_float() for i in range(nK)];
            self.logprint((nI,nJ,nK));
            nAll = nI*nJ*nK;
            d={}
            dqs=[];
            for quantity in qs:
                if quantity not in readin:
                    self.seek(nAll*4*size,1);
                else:
                    d.update({'x'+quantity:xdr.Unpacker(self.read(nAll*4*size))})
                    dqs.append(quantity);
            d.update({'nAll':nAll,'dqs': dqs,
                      'xp':Ip,'yp':Jp,'zp':Kp});
            doms.append(d);
        self.logprint('making points');
        if not lazy:
            #making points
            doms[:] = pool.map(make_points,doms);
        self.logprint('converting buffers');
        doms[:] = pool.map(call,doms);
        pool.close();
        self.logprint('done! stringing together');
        if lazy:
            return LazyIter(doms);
        else:
            #string together the domains;
            for dom in doms[1:]:
                for k in doms[0]:
                    doms[0][k].extend(dom[k]);
                del dom;
            return doms[0];
        pass;
    
    def _getmovie(self,pool_size,skip=1):
        params  = self.header['params'];
        nparams = len(params);
        pbytes = (nparams+1)*4;
        frames=[];
        cur = 0;
        while True:
            c=self.tell();
            self.read(1); #python, y u no eof?
            if self.tell() == c:
                break;
            self.seek(c);
            d=self.get_dict('fii',['time','step','pnum']);
            if (cur % skip) == 0:
                self.logprint('reading in frame at lsp step {}'.format(d['step']));
                d['xdr'] = xdr.Unpacker(self.read(pbytes*d['pnum']));
                self.logprint('done reading');
                frames.append(d);
            else:
                self.seek(d['pnum']*p_bytes,1);
        pool=multiprocessing.Pool(pool_size);
        self.logprint('converting frames');
        for i,d in enumerate(frames):
            d['params'] = params;
            d['pbytes'] = pbytes;
            f=Callable(d,convert_particles);
            l = range(d['pnum']);
            particles[:] = pool.map(f,l);
            del d['xdr'], d['pnum'], d['pbytes'];
            for k in d['params']:
                d[k] = [];
            for i,p in enumerate(particles):
                self.logprint('doing particle {}'.format(i));
                for k in frame['params']:
                    d[k].append(p[k]);
            del particles;
            frames[i] = d;
        pool.close();
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
        buf = self.read();
        nfloats = len(buf)/4;
        x = xdr.Unpacker(buf);
        del buf;
        data = x.unpack_farray(nfloats,x.unpack_float);
        del x;
        out={};
        for i,key in enumerate(params):
            out[key] = [j for j in data[i::nparams]];
        return out;
        
    def get_data(self,var=None,pool_size=16,lazy=False):
        if not var and (self.header['dump_type']== 2 or self.header['dump_type']== 3):
            var=[i[0] for i in self.header['quantities']];
        if self.header['dump_type'] == 2:
            return self._getfields(var,pool_size,lazy,vector=True);
        elif self.header['dump_type'] == 3:
            return self._getfields(var,pool_size,lazy,vector=False);
        elif self.header['dump_type'] == 6:
            return self._getmovie(pool_size);
        elif self.header['dump_type'] == 10:
            return self._getpext();
        else:
            return None;
    pass;
def burst_pmovie(f,outfmt,skip=1):
    if f.header['dump_type'] != 6:
        raise ValueError("Can't burst non-pmovie type dump");
    nparams=len(f.header['params']);
    p_bytes = (nparams+1)*4;
    for i in itools.count():
        c=f.tell();
        f.read(1); #python, y u no eof?
        if f.tell() == c:
            break;
        f.seek(c);
        d=f.get_dict('fii',['time','step','pnum']);
        if (i % skip) != 0:
            f.seek(d['pnum']*p_bytes,1);
            continue;
        d['xdr']=xdr.Unpacker(f.read(d['pnum']*p_bytes));
        d['params']=f.header['params'];
        print('converting {} for {} bytes'.format(i,d['pnum']*p_bytes));
        d=convert_frame(d);
        outname = outfmt.format(i);
        print('outputting to {}'.format(outname));
        with open(outname,'w') as f:
            cPickle.dump(d,f,2);
    pass;
