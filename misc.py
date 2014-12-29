'''
Miscellaneous definitions.
'''
def conv(arg,default=None,func=None):
    if func:
        return func(arg) if arg else default;
    else:
        return arg if arg else default;

def read(filename,dictlabel='s', dumpfull=True):
    with open(filename,'r') as f:
        d=pickle.load(f);
    if type(d) == np.ndarray or dumpfull:
        return d;
    elif type(d) == dict:
        return d[dictlabel];
    else:
        s = str(type(d));
        errstr='Unknown pickle type "{}" loaded from file "{}".'.format(s,filename);
        raise IOError(errstr);
