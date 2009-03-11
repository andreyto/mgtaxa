"""Factory: Object Oriented Currying"""

# Copyright (c) 2008, Peter Fein
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
# Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
# Neither the name of the Python Factory project nor the names of its
# contributors may be used to endorse or promote products derived from this
# software without specific prior written permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

__author__ = "Peter Fein"
__email__ = "pfein@pobox.com"
__url__ = "http://pypi.python.org/pypi/Factory/"
__revision__ = "$Rev: 9 $"
__version__ = (1, 2)

import types
import inspect
import itertools
import new
import functools

__all__ = ['Factory', 'FactoryMixin', 'returnFactory', 'factoryAttribute',
           'factoryDescriptor', 'Bunch', 'ObjectTemplate']

class Factory(object):
    """Object oriented partial function application"""

    # the interpreter will expand '__args' in __slots__ for us, but explicitly
    # spelling out private names allows __setattr__ to be simpler
    __slots__ = ['_Factory__callable', '_Factory__args', '_Factory__varargs',
                 '_Factory__permit_varargs', '_Factory__permit_varkwargs',
                 '__dict__']

    def __init__(self, callee, *args, **kwargs):
        """
        callee
            the function to bind

        *args
            varargs to pass to ``bind``

        **kwargs
            keyword args to pass to ``bind``
        """
        ## if we're given a Factory as a callee, copy off it and skip the rest of init
        ## XXX this is broken if callee is a subclass of Factory. Oh well.
        if isinstance(callee, Factory):
            for s in self.__slots__:
                setattr(self, s, getattr(callee, s))
            self.__dict__ = callee.__dict__.copy()
            self.bind(*args, **kwargs)
            return

        self.__callable = callee
        self.__args = set()
        self.__varargs = []

        # We do explicit type inspection instead of using
        # callable(), because we need to know where to find
        # the arguments.

        if isinstance(callee, (types.BuiltinFunctionType, types.BuiltinMethodType)):
            # getargspec (below) doesn't support builtin functions
            inspectables = []
        elif isinstance(callee, types.FunctionType):
            inspectables = [callee]
        elif isinstance(callee, type):
            inspectables = [c.__init__.im_func for c in callee.__mro__
                            if hasattr(c.__init__,'im_func')]
        elif isinstance(callee, types.MethodType):
            inspectables = [callee]
        elif hasattr(callee, '__call__'):
            inspectables = [callee.__call__]
        else:
            raise TypeError("must provide known callable type, not %r" % callee)

        if inspectables:
            # We accept variable numbers of arguments if the first inspectable accepts a
            # variable number of arguments; even if later inspectables accept them, we
            # wouldn't get past the first if it doesn't.
            self.__permit_varargs = inspect.getargspec(inspectables[0])[1]

            # We accept variable keyword arguments (i.e., we don't check that a bound keyword
            # arg matches some specification of keyword args) if *all* the inspectables accept
            # **kwargs.
            self.__permit_varkwargs = all(inspect.getargspec(i)[2] for i in inspectables)
        else:
            # no inspectables. This happens for builtin functions, types and
            # subclasses of builtin types. Since we can't examine what
            # parameters are allowed, we accept any args or kwargs. As a
            # result though, we're unable to detect bad args/attrs at
            # assignment time, and will raise a TypeError at call time.
            self.__permit_varargs = self.__permit_varkwargs = True

        for inspectable in inspectables:
            (the_args, _, the_kwargs, _) = inspect.getargspec(inspectable)
            for arg in the_args:
                if arg != 'self':
                    self.__args.add(arg)
            if the_kwargs is None:
                # We can't go any farther than the first class which doesn't
                # cooperatively accept **kwargs (presumably to pass to the super
                # constructor)
                break
        self.bind(*args, **kwargs)

    def __setattr__(self, attr, value):
        if attr not in self.__slots__ and \
               attr not in self.__args and \
               not self.__permit_varkwargs:
            raise AttributeError('No such argument %s' % attr)
        else:
            super(Factory, self).__setattr__(attr, value)

    def generateArgs(self, *args, **kwargs):
        """return the (args, kwargs) that would be used when the Factory is called"""
        args_dict = self.__dict__.copy()
        args_dict.update(kwargs)
        args = list(itertools.chain(self.__varargs, args))
        return (args, args_dict)

    def __call__(self, *args, **kwargs):
        (args, args_dict) = self.generateArgs(*args, **kwargs)
        return self.getCallable()(*args, **args_dict)

    def bind(self, *args, **kwargs):
        """update bound args & kwargs and return self"""
        if args:
            if self.__permit_varargs:
                self.__varargs.extend(args)
            else:
                raise TypeError('Initial callee does not accept *args')
        for k, v in kwargs.iteritems():
            setattr(self,  k, v)
        return self

    def getCallable(self):
        """return the original callable"""
        return self.__callable

    def __repr__(self):
        return "<Factory(%r) at %s>" % (self.__callable, hex(id(self)))

class FactoryMixin(object):
    """a mixin providing a ``factory`` method, which produces Factories"""

    @classmethod
    def factory(cls, *args, **kwargs):
        """return a Factory for the class, binding kwargs"""
        return Factory(cls).bind(*args, **kwargs)

def returnFactory(func):
    """a decorator which *replaces* a function with its Factory-producing equivalent."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return Factory(func).bind(*args, **kwargs)
    return wrapper

def factoryAttribute(func):
    """adds a factory attribute to the decorated func"""
    def factory(*args, **kwargs):
        return Factory(func).bind(*args, **kwargs)
    func.factory = factory
    return func

class factoryDescriptor(object):
    """a descriptor that produces instance methods with a factory attribute.

    Inside classes, use this descriptor instead of factoryAttribute. This class may be used as a decorator.
    """
    __slots__ = ['func']

    def __init__(self, func):
        self.func = func

    def __get__(self, obj, cls):
        return FactoryInstanceMethod(self.func, obj, cls)

class FactoryInstanceMethod(object):
    """an instancemethod-like object that adds a factory attribute. See factoryAttribute"""

    __slots__ = ['im_func', 'im_self', 'im_class', '__call__']

    def __init__(self, func, obj, cls):
        self.im_func = func
        self.im_self = obj
        self.im_class = cls
        self.__call__ = new.instancemethod(func, obj, cls)

    def factory(self, *args, **kwargs):
        return Factory(self.__call__).bind(*args, **kwargs)

class Bunch(object):
    """just a bunch of attributes.

    Calling a Bunch returns a new copy.
    """
    def __init__(self, **kwargs):
        """kwargs are turned into attrs"""
        self.__dict__.update(**kwargs)

    def __call__(self):
        return self.__class__(**self.__dict__)
    
    def get(self, name, default=None):
        return getattr(self, name, default)

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __ne__(self, other):
        return not (self == other)

    def __repr__(self):
        return '%s(%s)' % (self.__class__.__name__,
                           ', '.join('%s=%r' % kv for kv in self.__dict__.iteritems()))

class ObjectTemplate(object):
    """a template for creating objects"""

    # we make bunchClass a slot to make iterating through __dict__ easier
    # need to specify __dict__ so that the interpreter creates one for us
    __slots__ = ['bunchClass', '__dict__']

    def __init__(self, bunchClass=Bunch, **kwargs):
        """
        bunchClass
            the type of object to create.  Defaults to Bunch.
        """
        self.bunchClass = bunchClass
        self.__dict__.update(**kwargs)

    def __call__(self):
        # build a "list" of 2-tuples (key, bunch_value) where bunch_value is
        # the orig_value, unless the orig_value is callable, in which case
        # bunch_value is the result of calling orig_value
        return self.bunchClass(**dict((k, (v() if callable(v) else v))
                                      for k, v in self.__dict__.iteritems()))