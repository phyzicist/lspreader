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
    #reading points first
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
            d[quantity+'x'].append(d['xdr'].unpack_float());
            d[quantity+'y'].append(d['xdr'].unpack_float());
            d[quantity+'z'].append(d['xdr'].unpack_float());
    del d['xdr'],d['qs'],d['nAll'];
    return d;

def read_scalars(d):
    for quantity in d['qs']:
        d[quantity]=[ d['xdr'].unpack_float() for k in range(d['nAll'])];
    del d['xdr'],d['qs'],d['nAll'];
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
    #globals
    movie_labels=['q','x','y','z','vx','vy','vz','E']
    movie_labels_w_ip=movie_labels+['xi','yi','zi'];
    def __init__(self,filename):
        file.__init__(self,filename,"rb");
        self._get_header();
    def get_chunk(self, n):
        return self.read(n);
        #bad hack, make it better!
        return ''.join([self.get_four() for i in range(n/4)]);
        # r='';
        # if self.i+n > self._mi:
        #     r = self.buf[self.i : self.i-self._mi];
        #     n -= self.i-self._mi;
        #     self._get_more();
        # r += self.buf[self.i:self.i+n];
        # self.i+=n;
        # return r;
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
        else:
            raise ValueError('Unknown dump_type: {}'.format(self.header['dump_type']));
        return;
    ###################
    #data processing    
    def _getfields(self,pool_size,lazy,vector=True):
        if vector:
            size=3; call=read_fields;
        else:
            size=1; call=read_scalars;
        doms = [];
        qs = [i[0] for i in self.header['quantities']];
        pool=multiprocessing.Pool(pool_size);
        print('reading positions and scanning');
        for i in range(self.header['domains']):
            iR, jR, kR = self.get_int(),self.get_int(),self.get_int();
            #getting grid parameters (real coordinates)
            #nI = f.get_int(); f.seek(nI*4,1); #skip grid parameters, they are floats
            nI = self.get_int(); Ip = [self.get_float() for i in range(nI)];
            nJ = self.get_int(); Jp = [self.get_float() for i in range(nJ)];
            nK = self.get_int(); Kp = [self.get_float() for i in range(nK)];
            nAll = nI*nJ*nK;
            doms.append({'filepos':self.tell(),'nAll':nAll,'xp':Ip,'yp':Jp,'zp':Kp});
            
            self.seek(nAll*4*size*len(qs),1);
        print('making points');
        if not lazy:
            #making points
            doms[:] = pool.map(make_points,doms);
        print('making buffers');
        #making buffers
        for dom in doms:
            self.seek(dom['filepos']);
            del dom['filepos'];
            dom.update( {'xdr':xdr.Unpacker(self.read(dom['nAll']*4*size*len(qs))),'qs':qs} );
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

    def get_data(self, pool_size=16,lazy=True):
        if self.header['dump_type'] == 2:
            return self._getfields(pool_size,lazy,vector=True);
        elif self.header['dump_type'] == 3:
            return self._getfields(pool_size,lazy,vector=False);
        elif self.header['dump_type'] == 6:
            return self._getmovie(pool);
    pass;
