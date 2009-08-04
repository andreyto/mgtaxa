### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from types import *
import re
import string, os
from copy import copy, deepcopy

class Struct(object):
    """Class to create 'struct's on the fly.
    Example: o = Struct()
             o.i = 2
             o.x = 'ababab'
             a = Struct({'i':2,'x':'ababab'})
             b = Struct(i=2,x='ababab')
             In all three cases, the result will be the same.
             __str__ method is redefined so that printing the object of this type
             will show all attributes which are not represented as 'instance at 0xXXXXXX'.
             """

    strStyle = "p" # p - print "pretty", s - print as one string

    def __init__(self,*lw,**kw):
        for dict in lw:
            for key in dict.keys():
                setattr(self,key,dict[key])
        for key in kw.keys():
            setattr(self,key,kw[key])

    def __str__(self):
        if self.strStyle == "p":
            return self.strPretty()
        else:
            return self.strDense()

    def __getitem__(self,key):
        try:
            return getattr(self,key)
        except AttributeError:
            raise KeyError(key)

    def __setitem__(self,key,value):
        setattr(self,key,value)
    
    def __len__(self):
        return len(self.__dict__)

    def __delitem__(self,key):
        try:
            delattr(self,key)
        except AttributeError:
            raise KeyError(key)

    def __iter__(self):
        return self.__dict__.iterkeys()

    def iterkeys(self):
        return self.__dict__.iterkeys()

    def setdefault(self,*l):
        try:
            return getattr(self,l[0])
        except AttributeError:
            if len(l) >= 2:
                setattr(self,l[0],l[1])
                return getattr(self,l[0])
            else:
                raise KeyError(l[0])

    def pop(self,key):
        return self.__dict__.pop(key)

    def update(self,other):
        if isinstance(other,Struct):
            o = other.__dict__
        else:
            o = other
        self.__dict__.update(o)

    def updateOtherMissing(self,other):
        """Update in other keys that are not yet present where"""
        if isinstance(other,Struct):
            o = other.__dict__
        else:
            o = other
        s = self.asDict()
        for (key,value) in s.items():
            o.setdefault(key,value)

    def updateFromOtherExisting(self,other):
        """Update in self keys that already present in self"""
        if isinstance(other,Struct):
            o = other.__dict__
        else:
            o = other
        s = self.asDict()
        for key in list(s.keys()):
            if key in o:
                s[key] = o[key]

    def asDict(self):
        return self.__dict__

    def keys(self):
        return self.__dict__.keys()
    
    def has_key(self,key):
        return self.__dict__.has_key(key)

    def get(self,*l):
        try:
            return getattr(self,*l)
        except AttributeError:
            raise KeyError(l[0])
            
    def strDense(self):
        keys = self.keys()
        keys.sort()
        pairs = []
        for key in keys:
            obj = self.__dict__[key]
            s_obj = str(obj)
            if self.isPrintable(s_obj):
                pairs.append((key,s_obj))
        return 'Struct('+`pairs`+')'

    def __repr__(self):
        return self.__str__()

    def isPrintable(self,reprObj):
        return not re.match("^\<.+ instance at 0x[0-9a-z]+\>$",reprObj)

    def strPretty(self):
        keys = self.keys()
        keys.sort()
        s = '\n'
        for key in keys:
            obj = self.__dict__[key]
            s_obj = str(obj)
            if self.isPrintable(s_obj):
                # add to TAB to all rows of attribute's representation
                lines = s_obj.split('\n')
                s_obj = '\n\t'.join(lines)
                s = s + key + '\t=\t' + s_obj + '\n'
        return s

    def scalars(self):
        """Return dictionary mapping names of "scalar" attributes to values.
        "Scalar" attributes are non-sequence primitive types, such as Int, Float, String, None."""
        r = {}
        for key in self.keys():
            val = self.__dict__(key)
            if type(val) in (NoneType,BooleanType,IntType,LongType,FloatType,StringType):
                r[key] = val
        return r

    def copy(self):
        return copy(self)




class Options(Struct):
    
    def copy(self):
        """Deep copy semantics"""
        return deepcopy(self)

    def keys(self):
        """Will ignore all attributes that start with _"""
        return [ k for k in Struct.keys(self) if not k.startswith("_") ]

    def freeze(self):
        """Make this object read-only"""
        Struct.__setattr__(self,"_is_frozen",True)
        for name in self.keys():
            val = getattr(self,name)
            if isinstance(val,Options):
                val.freeze()

    def unfreeze(self):
        """Make this object mutable again after previous call to freeze()"""
        try:
            Struct.__delattr__(self,"_is_frozen")
        except AttributeError:
            pass
        for name in self.keys():
            val = getattr(self,name)
            if isinstance(val,Options):
                val.unfreeze()

    def __setattr__(self,name,value):
        if getattr(self,"_is_frozen",False):
            raise AttributeError(name)
        else:
            Struct.__setattr__(self,name,value)

    def __delattr__(self,name):
        if getattr(self,"_is_frozen",False):
            raise AttributeError(name)
        else:
            Struct.__delattr__(self,name,value)

