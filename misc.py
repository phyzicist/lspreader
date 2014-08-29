'''
Miscellaneous definitions.
'''
def conv(arg,default=None,func=None):
    if func:
        return func(arg) if arg else default;
    else:
        return arg if arg else default;
