import xdrlib as xdr;
import multiprocessing;

class LazyIter(object):
    def __init__(self,doms):
        self.doms = doms;
        pass;
    def _next_dom(self):
        self.dom = self.domsi.next();#should throw StopIteration at the last dom
        self.curi = 0;
        self.keys = [i for i in self.dom if  (i != 'xp' and i != 'yp' and i != 'zp')];
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
    for quantity in d['qs']:
        d[quantity+'x']=[];
        d[quantity+'y']=[];
        d[quantity+'z']=[];
        for k in range(d['nAll']):
            d[quantity+'x'].append(d['x'+quantity].unpack_float());
            d[quantity+'y'].append(d['x'+quantity].unpack_float());
            d[quantity+'z'].append(d['x'+quantity].unpack_float());
        del d['x'+quantity];
    del d['qs'],d['nAll'];
    return d;

def read_scalars(d):
    for quantity in d['dqs']:
        d[quantity]=[ d['x'+quantity].unpack_float() for k in range(d['nAll'])];
        del d['x'+quantity];
    del d['dqs'],d['nAll'];
    return d;

def convert_frame(d):
    particles=[];
    for i in range(d['pnum']):
        p={}
        p['ip']=d['xdr'].unpack_int();
        for i in d['params']:
            p[i[0]] = d['xdr'].unpack_float();
        particles.append(p);
    del d['xdr'];
    d['particles']=particles;
    #print('.');
    return d;

class LspOutput(file):
    '''represents an lsp output file on call,
       reads the header on open'''
    #static variables
    ip=['xi','yi','zi'];
    movie_labels=['q','x','y','z','vx','vy','vz','E']
    movie_labels_w_ip=movie_labels+ip;
    def __init__(self,filename):
        file.__init__(self,filename,"rb");
        self._get_header();
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
            print(units);
            if n == len(LspOutput.movie_labels):
                labels=LspOutput.movie_labels
            elif n == len(LspOutput.movie_labels_w_ip):
                labels=LspOutput.movie_labels_w_ip;
            else:
                raise ValueError('Incorrect number of parameters: {}'.format(n));
            self.header['params'] = [(i[0],i[1]) for i in zip(labels,units,flags) if i[2]];
        elif self.header['dump_type'] == 10:
            #this is a particle extraction file:
            self.header['geometry'] = self.get_int();
            #reading quantities
            n = self.get_int();
            units = [self.get_str() for i in range(n)];
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
        print('reading positions and making buffers');
        for i in range(self.header['domains']):
            iR, jR, kR = self.get_int(),self.get_int(),self.get_int();
            #getting grid parameters (real coordinates)
            nI = self.get_int(); Ip = [self.get_float() for i in range(nI)];
            nJ = self.get_int(); Jp = [self.get_float() for i in range(nJ)];
            nK = self.get_int(); Kp = [self.get_float() for i in range(nK)];
            nAll = nI*nJ*nK;
            d={}
            dqs=[];
            for quantity in qs:
                if quantity not in readin:
                    self.seek(nAll*4*size,1);
                else:
                    d.update({'x'+quantity : xdr.Unpacker(self.read(nAll*4*size))})
                    dqs.append(quantity);
            d.update({'nAll':nAll,'dqs': dqs,
                      'xp':Ip,'yp':Jp,'zp':Kp});
            doms.append(d);
        print('making points');
        if not lazy:
            #making points
            doms[:] = pool.map(make_points,doms);
        print('converting buffers');
        doms[:] = pool.map(call,doms);
        pool.close();
        print('done! stringing together');
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
    
    def _getmovie(self,pool_size):    
        nparams=len(self.header['params']);
        p_bytes = (nparams+1)*4;
        frames=[];
        while True:
            c=self.tell();
            self.read(1); #python, y u no eof?
            if self.tell() == c:
                break;
            self.seek(c);
            d=self.get_dict('fii',['time','step','pnum']);
            d['filepos']=self.tell();
            d['xdr']=xdr.Unpacker(self.read(d['pnum']*p_bytes));
            d['params']=self.header['params'];
            frames.append(d);
        pool=multiprocessing.Pool(pool_size);
        frames[:] = pool.map(convert_frame,frames);
        pool.close();
        return frames;
    
    def _getpext(self):
        nparams = len(self.header['quantities']);
        params = ['t','q','x','y','z','vx','vy','vz'];
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
        if not var:
            var=[i[0] for i in self.header['quantities']];
        if self.header['dump_type'] == 2:
            return self._getfields(var,pool_size,lazy,vector=True);
        elif self.header['dump_type'] == 3:
            return self._getfields(var,pool_size,lazy,vector=False);
        elif self.header['dump_type'] == 6:
            return self._getmovie(pool);
        elif self.header['dump_type'] == 10:
            return self._getpext();
        else:
            return None;
    pass;
