"""Contains class that implements in-memory cache of objects using weak references"""
from MGT.UUID import *
import weakref

__all__ = ["ObjectCache"]

class ObjectCache(weakref.WeakValueDictionary):
    """Cache based on WeakValueDictionary.
    Once no hard references any longer exist,
    key-value pairs are up for garbage collection.
    Provides thread safe method that returns existing
    value for key or else creates new value.
    If you want to share cached object among multiple
    calls, use the same value of key.
    If you want to generate new object, use
    object_cache.get_produce(key=object_cache.gen_key(),produce=produce)
    """

    def __init__(self,*l,**kw):
        weakref.WeakValueDictionary.__init__(self,*l,**kw)
        import threading
        self.lock = threading.Lock()
    
    @classmethod
    def gen_key(klass):
        """Generate and return new unique key"""
        return genId()

    def get_produce(self,key,produce=None):
        """Return value for key if it exists, else assign from producer and return.
        This is done thread-safe.
        @param key dictionary key to look up
        @param produce if not None [default is None] and key does not exists,
        produce(key) should generate new value
        @return value corresponding to key in dictionary
        @post if did not exist, new value will be generated and assigned to key in 
        dictionary"""
        try:
            return self[key]
        except KeyError:
            if produce is None:
                raise
            try:
                self.lock.acquire()
                #another thread might have produced and stored
                #the value in the meantime, so try one more time
                #to get it
                try:
                    return self[key]
                except KeyError:
                    #Cannot use self[key] = producer() 
                    #because it might be instantly garbage-collected -
                    #no strong refs to value yet.
                    val = produce(key)
                    self[key] = val
                    return self[key]
            finally:
                self.lock.release()
            
